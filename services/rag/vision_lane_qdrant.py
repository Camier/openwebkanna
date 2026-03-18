"""Qdrant-backed visual retrieval lane for ColModernVBERT multivectors.

This follows Qdrant's documented multivector query path:
- store late-interaction embeddings as a named multivector
- query that named vector directly
- constrain scope with payload filters
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .colmodernvbert_profile import ColModernVBERTTuningProfile, get_default_colmodernvbert_profile
from .vision_colmodernvbert import ColModernVBERTLateInteractionRuntime


@dataclass(frozen=True)
class QdrantVisionLaneConfig:
    """Runtime settings for figure-only visual retrieval."""

    collection_name: str = "rag_evidence"
    vector_name: str = "vision_li"
    object_type: str = "figure"
    candidate_pool: int = get_default_colmodernvbert_profile().qdrant_candidate_pool_per_lane
    final_top_k: int = get_default_colmodernvbert_profile().final_top_k
    prioritized_figure_kinds: tuple[str, ...] = field(
        default_factory=lambda: get_default_colmodernvbert_profile().prioritized_figure_kinds
    )
    prioritized_kind_boost: float = 0.05


class QdrantColModernVBERTVisualLane:
    """Figure-scoped visual retrieval lane using ColModernVBERT query embeddings."""

    def __init__(
        self,
        *,
        qdrant_client: Any,
        runtime: ColModernVBERTLateInteractionRuntime,
        config: QdrantVisionLaneConfig | None = None,
        tuning_profile: ColModernVBERTTuningProfile | None = None,
    ) -> None:
        self.qdrant_client = qdrant_client
        self.runtime = runtime
        self.config = config or QdrantVisionLaneConfig()
        self.tuning_profile = tuning_profile or get_default_colmodernvbert_profile()

    def __call__(self, *, query: str, top_k: int = 10) -> list[dict[str, Any]]:
        """Embed a text query and retrieve matching figure points."""

        query = (query or "").strip()
        if not query:
            return []

        query_embeddings = self.runtime.embed_queries([query])
        if not query_embeddings:
            return []

        limit = max(top_k, self.config.final_top_k)
        candidate_pool = max(limit, self.config.candidate_pool)
        raw_hits = self._query_qdrant(query_embeddings[0], candidate_pool)
        normalized_hits = [self._normalize_hit(item) for item in raw_hits]
        rescored_hits = self._boost_prioritized_figure_kinds(normalized_hits)
        rescored_hits.sort(key=lambda item: float(item.get("score", 0.0)), reverse=True)
        return rescored_hits[:limit]

    def _query_qdrant(self, query_embedding: list[list[float]], limit: int) -> list[Any]:
        query_filter = self._build_filter()

        if hasattr(self.qdrant_client, "query_points"):
            result = self.qdrant_client.query_points(
                collection_name=self.config.collection_name,
                query=query_embedding,
                using=self.config.vector_name,
                query_filter=query_filter,
                limit=limit,
            )
            return list(getattr(result, "points", result))

        if hasattr(self.qdrant_client, "search"):
            try:
                return list(
                    self.qdrant_client.search(
                        collection_name=self.config.collection_name,
                        query_vector=(self.config.vector_name, query_embedding),
                        query_filter=query_filter,
                        limit=limit,
                    )
                )
            except TypeError:
                return list(
                    self.qdrant_client.search(
                        collection_name=self.config.collection_name,
                        query_vector=query_embedding,
                        vector_name=self.config.vector_name,
                        query_filter=query_filter,
                        limit=limit,
                    )
                )

        raise RuntimeError("Qdrant client does not expose a supported official query API.")

    def _build_filter(self) -> Any:
        try:
            from qdrant_client.http.models import FieldCondition, Filter, MatchValue
        except ImportError as exc:
            raise RuntimeError("Qdrant filter models are required for visual lane filtering.") from exc

        return Filter(
            must=[
                FieldCondition(
                    key="object_type",
                    match=MatchValue(value=self.config.object_type),
                )
            ]
        )

    def _normalize_hit(self, item: Any) -> dict[str, Any]:
        payload = dict(getattr(item, "payload", {}) or {})
        return {
            "id": str(getattr(item, "id")),
            "object_type": payload.get("object_type", self.config.object_type),
            "score": float(getattr(item, "score", 0.0)),
            "document_id": payload.get("document_id"),
            "page_id": payload.get("page_id"),
            "page_number": payload.get("page_number"),
            "figure_id": payload.get("figure_id"),
            "image_uri": payload.get("image_uri"),
            "figure_kind": payload.get("figure_kind"),
            "payload": payload,
        }

    def _boost_prioritized_figure_kinds(self, hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
        prioritized = set(self.config.prioritized_figure_kinds)
        boosted: list[dict[str, Any]] = []
        for hit in hits:
            score = float(hit.get("score", 0.0))
            if hit.get("figure_kind") in prioritized:
                score += self.config.prioritized_kind_boost
            updated = dict(hit)
            updated["score"] = score
            boosted.append(updated)
        return boosted
