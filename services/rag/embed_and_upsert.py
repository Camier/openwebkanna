"""Embedding and Qdrant upsert orchestration for materialized evidence."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol

from .materialize_evidence import MaterializedRecord, make_qdrant_point_id
from .qdrant_schema import RagCollectionConfig


class DenseTextEmbedder(Protocol):
    def embed_texts(self, texts: list[str]) -> list[list[float]]:
        """Return one dense vector per input text."""


class SparseTextEmbedder(Protocol):
    def embed_texts(self, texts: list[str]) -> list[dict[str, list[float] | list[int]]]:
        """Return one sparse vector mapping per input text."""


class VisionEmbedder(Protocol):
    def embed_images(self, image_uris: list[str]) -> list[list[list[float]]]:
        """Return one multivector per input image."""


@dataclass(frozen=True)
class EmbedAndUpsertConfig:
    """Execution-time flags for embedding / upsert."""

    collection_name: str = "rag_evidence"
    batch_size: int = 64
    require_vision_for_figures: bool = True


@dataclass
class EmbedAndUpsertResult:
    """Counts emitted by an embed-and-upsert run."""

    points_attempted: int = 0
    points_upserted: int = 0
    points_skipped: int = 0
    page_points: int = 0
    figure_points: int = 0


class EmbedAndUpsertPipeline:
    """Convert materialized records into Qdrant points and upsert them."""

    def __init__(
        self,
        *,
        qdrant_client: Any,
        collection_cfg: RagCollectionConfig,
        runtime_cfg: EmbedAndUpsertConfig | None = None,
        text_dense_embedder: DenseTextEmbedder | None = None,
        text_sparse_embedder: SparseTextEmbedder | None = None,
        vision_embedder: VisionEmbedder | None = None,
    ) -> None:
        self.qdrant_client = qdrant_client
        self.collection_cfg = collection_cfg
        self.runtime_cfg = runtime_cfg or EmbedAndUpsertConfig(collection_name=collection_cfg.collection_name)
        self.text_dense_embedder = text_dense_embedder
        self.text_sparse_embedder = text_sparse_embedder
        self.vision_embedder = vision_embedder

    def run(self, records: list[MaterializedRecord]) -> EmbedAndUpsertResult:
        """Embed all records and upsert them in batches."""

        result = EmbedAndUpsertResult(points_attempted=len(records))
        for batch in _chunked(records, self.runtime_cfg.batch_size):
            points = self._build_points(batch, result)
            if not points:
                continue
            self.qdrant_client.upsert(
                collection_name=self.runtime_cfg.collection_name,
                points=points,
            )
            result.points_upserted += len(points)
        result.points_skipped = result.points_attempted - result.points_upserted
        return result

    def _build_points(
        self,
        records: list[MaterializedRecord],
        result: EmbedAndUpsertResult,
    ) -> list[dict[str, Any]]:
        text_dense_map = self._embed_text_dense(records)
        text_sparse_map = self._embed_text_sparse(records)
        vision_map = self._embed_vision(records)

        points: list[dict[str, Any]] = []
        for record in records:
            object_type = str(record.payload["object_type"])
            vector_payload: dict[str, Any] = {}

            dense_vec = text_dense_map.get(record.point_id)
            sparse_vec = text_sparse_map.get(record.point_id)
            vision_vec = vision_map.get(record.point_id)

            if dense_vec is not None:
                vector_payload["text_dense"] = dense_vec
            if sparse_vec is not None:
                vector_payload["text_sparse"] = sparse_vec
            if vision_vec is not None:
                vector_payload["vision_li"] = vision_vec

            if not self._record_is_eligible(record, vector_payload):
                continue

            points.append(
                {
                    "id": make_qdrant_point_id(record.point_id),
                    "payload": {
                        "point_id": record.point_id,
                        **record.payload,
                    },
                    "vector": vector_payload,
                }
            )
            if object_type == "page":
                result.page_points += 1
            elif object_type == "figure":
                result.figure_points += 1

        return points

    def _embed_text_dense(self, records: list[MaterializedRecord]) -> dict[str, list[float]]:
        if self.text_dense_embedder is None:
            return {}
        eligible = [(record.point_id, record.text_for_embedding) for record in records if record.text_for_embedding]
        if not eligible:
            return {}
        vectors = self.text_dense_embedder.embed_texts([text for _, text in eligible])
        return {point_id: vector for (point_id, _), vector in zip(eligible, vectors)}

    def _embed_text_sparse(self, records: list[MaterializedRecord]) -> dict[str, dict[str, list[float] | list[int]]]:
        if self.text_sparse_embedder is None:
            return {}
        eligible = [
            (record.point_id, record.text_for_embedding)
            for record in records
            if record.text_for_embedding and record.payload["object_type"] in {"page", "figure"}
        ]
        if not eligible:
            return {}
        vectors = self.text_sparse_embedder.embed_texts([text for _, text in eligible])
        return {point_id: vector for (point_id, _), vector in zip(eligible, vectors)}

    def _embed_vision(self, records: list[MaterializedRecord]) -> dict[str, list[list[float]]]:
        if self.vision_embedder is None:
            return {}
        eligible = [
            (record.point_id, record.image_uri)
            for record in records
            if record.payload["object_type"] == "figure" and record.image_uri
        ]
        if not eligible:
            return {}
        vectors = self.vision_embedder.embed_images([image_uri for _, image_uri in eligible])
        return {point_id: vector for (point_id, _), vector in zip(eligible, vectors)}

    def _record_is_eligible(self, record: MaterializedRecord, vector_payload: dict[str, Any]) -> bool:
        object_type = str(record.payload["object_type"])
        if object_type == "page":
            return "text_dense" in vector_payload and "text_sparse" in vector_payload
        if object_type == "figure":
            if "text_dense" not in vector_payload or "text_sparse" not in vector_payload:
                return False
            if self.runtime_cfg.require_vision_for_figures and "vision_li" not in vector_payload:
                return False
            return True
        return False


def _chunked(items: list[MaterializedRecord], batch_size: int) -> list[list[MaterializedRecord]]:
    return [items[i : i + batch_size] for i in range(0, len(items), batch_size)]
