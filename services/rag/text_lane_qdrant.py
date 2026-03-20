"""Qdrant-native text hybrid lane for the one-collection RAG.

This lane uses:
- native dense+sparse query encoding from `services.rag.text_query_runtime`
- Qdrant's official hybrid query API with dense+sparse prefetch and RRF fusion

It targets the canonical `rag_evidence` collection rather than the legacy
page-only collection.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from qdrant_client import models

from services.rag.text_query_runtime import QueryEncodingRuntimeError


@dataclass(frozen=True)
class QdrantTextHybridLaneConfig:
    collection_name: str = "rag_evidence"
    dense_vector_name: str = "text_dense"
    sparse_vector_name: str = "text_sparse"
    allowed_object_types: tuple[str, ...] = ("page", "figure")
    prefetch_k: int = 100
    sparse_prefetch_k: int = 100
    fusion_k: int = 100
    group_by: str = "page_id"
    score_threshold: float | None = None


@dataclass
class QdrantTextHybridLane:
    qdrant_client: Any
    encoder_runtime: Any
    config: QdrantTextHybridLaneConfig = field(default_factory=QdrantTextHybridLaneConfig)
    last_raw_result: dict[str, Any] = field(default_factory=dict)

    def __call__(self, query: str, top_k: int = 10) -> list[dict[str, Any]]:
        raw = self.retrieve(query=query, top_k=top_k)
        self.last_raw_result = raw
        hits = raw.get("reranked_hits") or raw.get("candidate_hits") or []
        normalized: list[dict[str, Any]] = []
        for hit in hits:
            payload = hit.get("payload") or {}
            point_id = hit.get("point_id") or payload.get("point_id") or payload.get("page_id") or hit.get("id")
            if not point_id:
                continue
            normalized.append(
                {
                    "id": point_id,
                    "point_id": point_id,
                    "object_type": payload.get("object_type") or "page",
                    "score": float(hit.get("score") or 0.0),
                    "match_type": "text_hybrid",
                    "payload": payload,
                    "document_id": payload.get("document_id") or payload.get("doc_id"),
                    "page_id": payload.get("page_id"),
                    "page_number": payload.get("page_number"),
                    "figure_id": payload.get("figure_id"),
                    "rank": hit.get("rank"),
                }
            )
        return normalized

    def retrieve(self, *, query: str, top_k: int = 10) -> dict[str, Any]:
        query = (query or "").strip()
        diagnostics: dict[str, Any] = {
            "lane": "qdrant_text_hybrid",
            "collection_name": self.config.collection_name,
            "dense_vector_name": self.config.dense_vector_name,
            "sparse_vector_name": self.config.sparse_vector_name,
            "group_by": self.config.group_by,
            "allowed_object_types": list(self.config.allowed_object_types),
            "warnings": [],
        }
        if not query:
            return {
                "query": query,
                "mode": "hybrid",
                "candidate_hits": [],
                "reranked_hits": [],
                "diagnostics": diagnostics,
            }

        try:
            encoded = self.encoder_runtime.encode_query(query)
        except QueryEncodingRuntimeError as exc:
            diagnostics["warnings"].append(str(exc))
            return {
                "query": query,
                "mode": "hybrid",
                "candidate_hits": [],
                "reranked_hits": [],
                "diagnostics": diagnostics,
            }

        if not encoded.dense:
            diagnostics["warnings"].append("Dense query encoding unavailable for text lane.")
            return {
                "query": query,
                "mode": "hybrid",
                "candidate_hits": [],
                "reranked_hits": [],
                "diagnostics": diagnostics,
            }

        sparse_query = encoded.sparse if encoded.sparse is not None else None
        if sparse_query is None:
            diagnostics["warnings"].append("Sparse query encoding unavailable; text lane ran dense-only.")
        if encoded.diagnostics:
            diagnostics["encoding"] = encoded.diagnostics

        response = self._query_qdrant(
            dense_query=encoded.dense,
            sparse_query=sparse_query,
            limit=max(top_k, self.config.fusion_k),
        )
        candidate_hits = [
            self._normalize_hit(
                point,
                rank=index,
                stage="candidate",
                vector_name="fusion_rrf" if sparse_query is not None else self.config.dense_vector_name,
            )
            for index, point in enumerate(response, start=1)
        ]
        reranked_hits = self._group_hits(candidate_hits, top_k=top_k)
        diagnostics["candidate_count"] = len(candidate_hits)
        diagnostics["reranked_count"] = len(reranked_hits)
        diagnostics["dense_enabled"] = True
        diagnostics["sparse_enabled"] = sparse_query is not None

        return {
            "query": query,
            "mode": "hybrid",
            "candidate_hits": candidate_hits,
            "reranked_hits": reranked_hits,
            "diagnostics": diagnostics,
        }

    def _query_qdrant(
        self,
        *,
        dense_query: list[float],
        sparse_query: models.SparseVector | None,
        limit: int,
    ) -> list[Any]:
        query_filter = self._build_filter()
        if sparse_query is not None:
            result = self.qdrant_client.query_points(
                collection_name=self.config.collection_name,
                prefetch=[
                    models.Prefetch(
                        query=dense_query,
                        using=self.config.dense_vector_name,
                        limit=max(limit, self.config.prefetch_k),
                    ),
                    models.Prefetch(
                        query=sparse_query,
                        using=self.config.sparse_vector_name,
                        limit=max(limit, self.config.sparse_prefetch_k),
                    ),
                ],
                query=models.FusionQuery(fusion=models.Fusion.RRF),
                limit=limit,
                query_filter=query_filter,
                with_payload=True,
                with_vectors=False,
                score_threshold=self.config.score_threshold,
            )
            return list(getattr(result, "points", result))

        result = self.qdrant_client.query_points(
            collection_name=self.config.collection_name,
            query=dense_query,
            using=self.config.dense_vector_name,
            limit=max(limit, self.config.prefetch_k),
            query_filter=query_filter,
            with_payload=True,
            with_vectors=False,
            score_threshold=self.config.score_threshold,
        )
        return list(getattr(result, "points", result))

    def _build_filter(self) -> models.Filter:
        return models.Filter(
            should=[
                models.FieldCondition(
                    key="object_type",
                    match=models.MatchValue(value=object_type),
                )
                for object_type in self.config.allowed_object_types
            ]
        )

    def _normalize_hit(
        self,
        item: Any,
        *,
        rank: int,
        stage: str,
        vector_name: str,
    ) -> dict[str, Any]:
        payload = dict(getattr(item, "payload", {}) or {})
        point_id = str(getattr(item, "id"))
        return {
            "point_id": point_id,
            "score": round(float(getattr(item, "score", 0.0)), 6),
            "rank": rank,
            "stage": stage,
            "doc_id": payload.get("document_id") or payload.get("doc_id"),
            "page_number": payload.get("page_number"),
            "title": payload.get("document_title") or payload.get("title") or payload.get("file_name"),
            "citation_key": payload.get("citation_key"),
            "group_id": payload.get(self.config.group_by) if self.config.group_by.lower() != "none" else None,
            "source_pdf_sha256": payload.get("source_pdf_sha256"),
            "retrieval_mode": "hybrid",
            "vector_name": vector_name,
            "base_score": float(getattr(item, "score", 0.0)),
            "payload": payload,
        }

    def _group_hits(self, hits: list[dict[str, Any]], *, top_k: int) -> list[dict[str, Any]]:
        group_by = (self.config.group_by or "").strip().lower()
        if group_by == "none":
            return hits[:top_k]

        grouped: list[dict[str, Any]] = []
        seen: set[str] = set()
        for hit in hits:
            payload = hit.get("payload") or {}
            group_value = (
                payload.get(self.config.group_by)
                or payload.get("page_id")
                or payload.get("document_id")
                or hit.get("point_id")
            )
            group_key = str(group_value)
            if group_key in seen:
                continue
            seen.add(group_key)
            grouped.append(hit)
            if len(grouped) >= top_k:
                break
        return grouped
