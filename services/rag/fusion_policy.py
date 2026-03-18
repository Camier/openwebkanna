"""Fusion policy helpers for one-collection retrieval.

These helpers keep the exact-chemistry lane hard-prioritized while
allowing soft lanes to be weighted differently by query type.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class FusionPolicy:
    """Weighted RRF configuration for soft retrieval lanes."""

    rrf_k: int = 60
    text_lane_weight: float = 1.0
    visual_lane_weight_text: float = 1.1
    visual_lane_weight_mixed: float = 1.2
    visual_lane_weight_smiles: float = 1.3
    chem_lane_weight: float = 1.05


@dataclass(frozen=True)
class LaneWeights:
    """Resolved lane weights for a specific query."""

    text: float = 1.0
    visual: float = 1.0
    chem: float = 1.0


def weights_for_query_type(query_type: str, policy: FusionPolicy | None = None) -> LaneWeights:
    """Resolve lane weights from query type."""

    policy = policy or FusionPolicy()
    if query_type == "smiles":
        visual_weight = policy.visual_lane_weight_smiles
    elif query_type == "mixed":
        visual_weight = policy.visual_lane_weight_mixed
    else:
        visual_weight = policy.visual_lane_weight_text

    return LaneWeights(
        text=policy.text_lane_weight,
        visual=visual_weight,
        chem=policy.chem_lane_weight,
    )


def weighted_rrf_merge(
    *,
    text_hits: list[dict[str, Any]],
    visual_hits: list[dict[str, Any]],
    chem_hits: list[dict[str, Any]],
    weights: LaneWeights | None = None,
    rrf_k: int = 60,
) -> list[dict[str, Any]]:
    """Reciprocal rank fusion with per-lane weighting."""

    weights = weights or LaneWeights()
    merged: dict[str, dict[str, Any]] = {}
    lane_specs = (
        ("text_hybrid", text_hits, weights.text),
        ("visual", visual_hits, weights.visual),
        ("chem_similarity", chem_hits, weights.chem),
    )

    for lane_name, lane_hits, lane_weight in lane_specs:
        for rank, hit in enumerate(lane_hits, start=1):
            point_id = str(hit["id"])
            weighted_rrf = lane_weight * (1.0 / (rrf_k + rank))
            record = merged.setdefault(point_id, dict(hit))
            existing_lanes = set(record.get("matched_lanes", []))
            existing_lanes.add(lane_name)
            record["matched_lanes"] = sorted(existing_lanes)
            record["rrf_score"] = float(record.get("rrf_score", 0.0)) + weighted_rrf
            record["score"] = max(float(record.get("score", 0.0)), float(hit.get("score", 0.0)))

    return sorted(
        merged.values(),
        key=lambda item: (float(item.get("rrf_score", 0.0)), float(item.get("score", 0.0))),
        reverse=True,
    )
