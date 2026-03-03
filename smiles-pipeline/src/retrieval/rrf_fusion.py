"""
Reciprocal Rank Fusion (RRF) utilities for hybrid text + SMILES retrieval.

This module keeps channel fusion deterministic and dependency-light so it can
be reused by API endpoints, evaluation scripts, and tests.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set


def normalize_ranked_docs(docs: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """
    Normalize ranked hits into a stable shape.

    Expected keys per input hit:
    - doc_id (preferred), or id/document_id
    - optional score, metadata
    """
    normalized: List[Dict[str, Any]] = []
    for rank, hit in enumerate(docs, start=1):
        doc_id = hit.get("doc_id") or hit.get("id") or hit.get("document_id")
        if not doc_id:
            continue
        normalized.append(
            {
                "doc_id": str(doc_id),
                "rank": rank,
                "score": hit.get("score"),
                "metadata": hit.get("metadata"),
            }
        )
    return normalized


def fuse_rankings_rrf(
    channel_rankings: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    top_k: int = 20,
    rrf_k: int = 60,
    channel_weights: Optional[Mapping[str, float]] = None,
) -> List[Dict[str, Any]]:
    """
    Fuse multiple ranked lists with Reciprocal Rank Fusion.

    RRF score(doc) = sum_c weight_c / (rrf_k + rank_c(doc))
    """
    if rrf_k <= 0:
        raise ValueError("rrf_k must be > 0")
    if top_k <= 0:
        raise ValueError("top_k must be > 0")

    weights = {name: 1.0 for name in channel_rankings}
    if channel_weights:
        for name, value in channel_weights.items():
            if name in weights:
                weights[name] = float(value)

    fused: MutableMapping[str, Dict[str, Any]] = {}

    for channel_name, raw_hits in channel_rankings.items():
        hits = normalize_ranked_docs(raw_hits)
        weight = weights.get(channel_name, 1.0)
        for hit in hits:
            doc_id = hit["doc_id"]
            rank = int(hit["rank"])
            contribution = weight / float(rrf_k + rank)
            if doc_id not in fused:
                fused[doc_id] = {
                    "doc_id": doc_id,
                    "rrf_score": 0.0,
                    "channels": {},
                }
            fused[doc_id]["rrf_score"] += contribution
            fused[doc_id]["channels"][channel_name] = {
                "rank": rank,
                "score": hit.get("score"),
                "weight": weight,
                "rrf_contribution": contribution,
            }

    ranked = sorted(fused.values(), key=lambda row: row["rrf_score"], reverse=True)
    return ranked[:top_k]


def compute_recall_at_k(
    ranked_doc_ids: Iterable[str],
    relevant_doc_ids: Iterable[str],
    k: int,
) -> float:
    """Compute Recall@k for one query."""
    if k <= 0:
        raise ValueError("k must be > 0")

    relevant: Set[str] = {str(doc_id) for doc_id in relevant_doc_ids}
    if not relevant:
        return 0.0

    top_k = list(ranked_doc_ids)[:k]
    hits = sum(1 for doc_id in top_k if str(doc_id) in relevant)
    return hits / float(len(relevant))
