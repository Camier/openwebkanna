from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from .normalize_smiles import normalize_smiles


def _weighted_rrf(
    lane_hits: list[tuple[str, list[dict[str, Any]], float]],
    *,
    k: int = 60,
) -> list[dict[str, Any]]:
    by_id: dict[str, dict[str, Any]] = {}
    scores: dict[str, float] = {}

    for lane_name, hits, weight in lane_hits:
        for rank, hit in enumerate(hits, start=1):
            point_id = hit.get("point_id") or hit.get("id")
            if not point_id:
                continue
            if point_id not in by_id:
                enriched = dict(hit)
                enriched.setdefault("lanes", [])
                by_id[point_id] = enriched
                scores[point_id] = 0.0
            by_id[point_id]["lanes"].append(lane_name)
            scores[point_id] += weight * (1.0 / (k + rank))

    fused = []
    for point_id, hit in by_id.items():
        merged = dict(hit)
        merged["score"] = scores[point_id]
        fused.append(merged)

    fused.sort(key=lambda item: item.get("score", 0.0), reverse=True)
    return fused


def _merge_exact_and_soft_hits(
    *,
    exact_hits: list[dict[str, Any]],
    soft_hits: list[dict[str, Any]],
    top_k: int,
) -> list[dict[str, Any]]:
    merged: list[dict[str, Any]] = []
    seen: set[str] = set()

    for rank, hit in enumerate(
        sorted(exact_hits, key=lambda item: float(item.get("score", 0.0)), reverse=True),
        start=1,
    ):
        point_id = str(hit.get("point_id") or hit.get("id") or "")
        if not point_id or point_id in seen:
            continue
        enriched = dict(hit)
        enriched.setdefault("point_id", point_id)
        enriched.setdefault("id", point_id)
        enriched["rank"] = rank
        merged.append(enriched)
        seen.add(point_id)

    next_rank = len(merged) + 1
    for hit in soft_hits:
        point_id = str(hit.get("point_id") or hit.get("id") or "")
        if not point_id or point_id in seen:
            continue
        enriched = dict(hit)
        enriched.setdefault("point_id", point_id)
        enriched.setdefault("id", point_id)
        enriched["rank"] = next_rank
        merged.append(enriched)
        seen.add(point_id)
        next_rank += 1
        if len(merged) >= top_k:
            break

    return merged[:top_k]


@dataclass
class OneCollectionRouterResult:
    query_type: str
    normalized_query: dict[str, str | None]
    exact_hits: list[dict[str, Any]]
    text_hits: list[dict[str, Any]]
    visual_hits: list[dict[str, Any]]
    fused_hits: list[dict[str, Any]]
    lanes_used: list[str]
    raw_text_result: dict[str, Any]
    warnings: list[str]


@dataclass
class OneCollectionRouter:
    text_lane: Callable[[str, int], list[dict[str, Any]]]
    visual_lane: Callable[[str, int], list[dict[str, Any]]] | None = None
    exact_lane: Callable[[str | None, str | None, int], list[dict[str, Any]]] | None = None
    rrf_k: int = 60
    text_weight: float = 1.0
    visual_weight_text: float = 1.1
    visual_weight_smiles: float = 1.3

    def retrieve(self, query: str, top_k: int = 10) -> OneCollectionRouterResult:
        normalized = normalize_smiles(query)
        is_smiles = bool(normalized.valid and normalized.canonical_smiles)
        query_type = "smiles" if is_smiles else "text"
        warnings: list[str] = []

        exact_hits: list[dict[str, Any]] = []
        if is_smiles and self.exact_lane is not None:
            try:
                exact_hits = self.exact_lane(
                    canonical_smiles=normalized.canonical_smiles,
                    inchikey=normalized.inchikey,
                    top_k=top_k,
                )
            except Exception as exc:
                warnings.append(f"Exact chemistry lane unavailable: {exc}")

        text_hits = self.text_lane(query, top_k=top_k)

        visual_hits: list[dict[str, Any]] = []
        if self.visual_lane is not None:
            try:
                visual_hits = self.visual_lane(query=query, top_k=top_k)
            except Exception as exc:
                warnings.append(f"ColModernVBERT visual lane unavailable: {exc}")

        soft_hits = _weighted_rrf(
            [
                ("text_hybrid", text_hits, self.text_weight),
                (
                    "visual",
                    visual_hits,
                    self.visual_weight_smiles if is_smiles else self.visual_weight_text,
                ),
            ],
            k=self.rrf_k,
        )
        fused_hits = _merge_exact_and_soft_hits(
            exact_hits=exact_hits,
            soft_hits=soft_hits,
            top_k=top_k,
        )

        lanes_used = ["text_hybrid"]
        if visual_hits:
            lanes_used.append("visual")
        if exact_hits:
            lanes_used.insert(0, "exact_chemistry")

        raw_text_result = getattr(self.text_lane, "last_raw_result", {}) or {}
        text_diagnostics = raw_text_result.get("diagnostics") or {}
        warnings.extend([str(item) for item in text_diagnostics.get("warnings") or []])

        return OneCollectionRouterResult(
            query_type=query_type,
            normalized_query={
                "canonical_smiles": normalized.canonical_smiles,
                "inchikey": normalized.inchikey,
            },
            exact_hits=exact_hits,
            text_hits=text_hits,
            visual_hits=visual_hits,
            fused_hits=fused_hits,
            lanes_used=lanes_used,
            raw_text_result=raw_text_result,
            warnings=warnings,
        )
