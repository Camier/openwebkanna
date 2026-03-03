"""Retrieval utilities for multi-channel fusion and evaluation."""

from .rrf_fusion import (
    compute_recall_at_k,
    fuse_rankings_rrf,
    normalize_ranked_docs,
)

__all__ = [
    "compute_recall_at_k",
    "fuse_rankings_rrf",
    "normalize_ranked_docs",
]
