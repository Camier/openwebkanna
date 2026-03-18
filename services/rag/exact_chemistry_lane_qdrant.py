"""Qdrant-backed exact chemistry lane for one-collection molecule points.

This lane uses documented Qdrant payload filtering rather than vector
search. Exact chemical identity is determined by normalized payload
fields:
- `inchikey`
- `canonical_smiles`
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class QdrantExactChemistryLaneConfig:
    """Runtime settings for exact molecule retrieval."""

    collection_name: str = "rag_evidence"
    object_type: str = "molecule"
    review_status: str = "parsed"
    inchikey_score: float = 1.0
    canonical_smiles_score: float = 0.99


def probe_exact_chemistry_runtime() -> list[str]:
    """Return readiness warnings for the exact chemistry lane."""

    warnings: list[str] = []

    try:
        from rdkit import Chem  # noqa: F401
    except ImportError:
        warnings.append(
            "Exact chemistry lane requires RDKit for SMILES normalization in the query path."
        )

    try:
        from qdrant_client.http.models import FieldCondition, Filter, MatchValue  # noqa: F401
    except ImportError:
        warnings.append(
            "Exact chemistry lane requires qdrant-client filter models "
            "(FieldCondition, Filter, MatchValue)."
        )

    return warnings


class QdrantExactChemistryLane:
    """Retrieve exact molecule hits by InChIKey and canonical SMILES."""

    def __init__(self, *, qdrant_client: Any, config: QdrantExactChemistryLaneConfig | None = None) -> None:
        self.qdrant_client = qdrant_client
        self.config = config or QdrantExactChemistryLaneConfig()

    def __call__(
        self,
        *,
        canonical_smiles: str | None,
        inchikey: str | None,
        top_k: int = 10,
    ) -> list[dict[str, Any]]:
        hits: list[dict[str, Any]] = []
        seen: set[str] = set()

        if inchikey:
            for point in self._scroll_exact_match(field_name="inchikey", value=inchikey, limit=top_k):
                hit = self._normalize_hit(point, score=self.config.inchikey_score, match_type="exact_inchikey")
                if hit["id"] in seen:
                    continue
                seen.add(hit["id"])
                hits.append(hit)

        remaining = max(top_k - len(hits), 0)
        if remaining > 0 and canonical_smiles:
            for point in self._scroll_exact_match(field_name="canonical_smiles", value=canonical_smiles, limit=remaining):
                hit = self._normalize_hit(point, score=self.config.canonical_smiles_score, match_type="exact_canonical_smiles")
                if hit["id"] in seen:
                    continue
                seen.add(hit["id"])
                hits.append(hit)

        hits.sort(key=lambda item: float(item.get("score", 0.0)), reverse=True)
        return hits[:top_k]

    def _scroll_exact_match(self, *, field_name: str, value: str, limit: int) -> list[Any]:
        try:
            from qdrant_client.http.models import FieldCondition, Filter, MatchValue
        except ImportError as exc:
            raise RuntimeError("Qdrant filter models are required for exact chemistry retrieval.") from exc

        result = self.qdrant_client.scroll(
            collection_name=self.config.collection_name,
            scroll_filter=Filter(
                must=[
                    FieldCondition(key="object_type", match=MatchValue(value=self.config.object_type)),
                    FieldCondition(key="review_status", match=MatchValue(value=self.config.review_status)),
                    FieldCondition(key=field_name, match=MatchValue(value=value)),
                ]
            ),
            limit=limit,
            with_payload=True,
            with_vectors=False,
        )

        if isinstance(result, tuple):
            return list(result[0])
        return list(result)

    def _normalize_hit(self, point: Any, *, score: float, match_type: str) -> dict[str, Any]:
        payload = dict(getattr(point, "payload", {}) or {})
        return {
            "id": str(getattr(point, "id")),
            "object_type": payload.get("object_type", self.config.object_type),
            "score": score,
            "match_type": match_type,
            "document_id": payload.get("document_id"),
            "page_id": payload.get("page_id"),
            "page_number": payload.get("page_number"),
            "figure_id": payload.get("figure_id"),
            "molecule_id": payload.get("molecule_id"),
            "payload": payload,
        }
