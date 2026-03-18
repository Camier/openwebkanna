"""Query parsing and lane orchestration for the one-collection RAG.

This module deliberately separates:
- query parsing
- lane execution
- fusion

The concrete vector search implementations should be injected from the
calling service layer. That keeps the routing contract stable while the
retrieval backend evolves.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from .fusion_policy import FusionPolicy, weighted_rrf_merge, weights_for_query_type
from .normalize_smiles import NormalizedSmiles, normalize_smiles


LaneCallable = Callable[..., list[dict[str, Any]]]


@dataclass(frozen=True)
class ParsedQuery:
    """Normalized representation of query intent."""

    raw_query: str
    query_type: str
    normalized_smiles: NormalizedSmiles | None
    enable_exact_chemistry: bool
    enable_text: bool
    enable_visual: bool
    enable_chem_similarity: bool

    @property
    def canonical_smiles(self) -> str | None:
        return self.normalized_smiles.canonical_smiles if self.normalized_smiles else None

    @property
    def inchikey(self) -> str | None:
        return self.normalized_smiles.inchikey if self.normalized_smiles else None


@dataclass(frozen=True)
class QueryRouterConfig:
    """Feature flags for retrieval lanes."""

    enable_text: bool = True
    enable_visual: bool = True
    enable_chem_similarity: bool = False
    default_group_by: str = "page_id"
    fusion_policy: FusionPolicy = field(default_factory=FusionPolicy)


@dataclass
class LaneBundle:
    """Injected lane callables.

    Each callable should return a list of hit dictionaries. The router
    imposes only a small contract on those hits:
    - `id`: stable point id
    - `object_type`: page | figure | molecule
    - `score`: float
    - grouping field, typically `page_id` or `document_id`
    """

    exact_chemistry: LaneCallable = field(default_factory=lambda: _empty_lane)
    text_hybrid: LaneCallable = field(default_factory=lambda: _empty_lane)
    visual: LaneCallable = field(default_factory=lambda: _empty_lane)
    chem_similarity: LaneCallable = field(default_factory=lambda: _empty_lane)


class QueryRouter:
    """Route a query through exact chemistry, text, visual, and chem lanes."""

    def __init__(self, config: QueryRouterConfig | None = None, lanes: LaneBundle | None = None) -> None:
        self.config = config or QueryRouterConfig()
        self.lanes = lanes or LaneBundle()

    def parse_query(self, raw_query: str) -> ParsedQuery:
        """Classify the query and enable the relevant retrieval lanes."""

        raw_query = (raw_query or "").strip()
        normalized = normalize_smiles(raw_query)
        has_valid_smiles = normalized.is_exact_matchable

        query_type = "smiles" if has_valid_smiles else "text"
        if has_valid_smiles and _has_non_smiles_context(raw_query):
            query_type = "mixed"

        return ParsedQuery(
            raw_query=raw_query,
            query_type=query_type,
            normalized_smiles=normalized if has_valid_smiles else None,
            enable_exact_chemistry=has_valid_smiles,
            enable_text=self.config.enable_text,
            enable_visual=self.config.enable_visual,
            enable_chem_similarity=has_valid_smiles and self.config.enable_chem_similarity,
        )

    def retrieve(self, raw_query: str, top_k: int = 10, group_by: str | None = None) -> dict[str, Any]:
        """Execute all relevant retrieval lanes and fuse the results."""

        parsed = self.parse_query(raw_query)
        group_by = group_by or self.config.default_group_by

        exact_hits: list[dict[str, Any]] = []
        text_hits: list[dict[str, Any]] = []
        visual_hits: list[dict[str, Any]] = []
        chem_hits: list[dict[str, Any]] = []

        if parsed.enable_exact_chemistry:
            exact_hits = self.lanes.exact_chemistry(
                canonical_smiles=parsed.canonical_smiles,
                inchikey=parsed.inchikey,
                top_k=top_k,
            )

        if parsed.enable_text:
            text_hits = self.lanes.text_hybrid(
                query=parsed.raw_query,
                top_k=top_k,
            )

        if parsed.enable_visual:
            visual_hits = self.lanes.visual(
                query=parsed.raw_query,
                top_k=top_k,
            )

        if parsed.enable_chem_similarity:
            chem_hits = self.lanes.chem_similarity(
                canonical_smiles=parsed.canonical_smiles,
                inchikey=parsed.inchikey,
                top_k=top_k,
            )

        fused = fuse_results(
            exact_hits=exact_hits,
            text_hits=text_hits,
            visual_hits=visual_hits,
            chem_hits=chem_hits,
            group_by=group_by,
            query_type=parsed.query_type,
            fusion_policy=self.config.fusion_policy,
        )

        return {
            "query_type": parsed.query_type,
            "normalized_query": {
                "canonical_smiles": parsed.canonical_smiles,
                "inchikey": parsed.inchikey,
            },
            "lanes_used": [
                lane_name
                for lane_name, enabled in (
                    ("exact_chemistry", parsed.enable_exact_chemistry),
                    ("text_hybrid", parsed.enable_text),
                    ("visual", parsed.enable_visual),
                    ("chem_similarity", parsed.enable_chem_similarity),
                )
                if enabled
            ],
            "groups": fused["groups"],
            "ungrouped_hits": fused["hits"],
        }


def fuse_results(
    *,
    exact_hits: list[dict[str, Any]],
    text_hits: list[dict[str, Any]],
    visual_hits: list[dict[str, Any]],
    chem_hits: list[dict[str, Any]],
    group_by: str,
    query_type: str,
    fusion_policy: FusionPolicy | None = None,
) -> dict[str, Any]:
    """Fuse exact and soft hits into grouped evidence bundles.

    Fusion policy:
    - exact chemistry hits always outrank soft matches
    - soft matches are merged conservatively by point id
    - grouping is performed after ranking
    """

    exact_ranked = _dedupe_hits(exact_hits)
    exact_ids = {hit["id"] for hit in exact_ranked}

    resolved_fusion_policy = fusion_policy or FusionPolicy()
    soft_hits = weighted_rrf_merge(
        text_hits=text_hits,
        visual_hits=visual_hits,
        chem_hits=chem_hits,
        weights=weights_for_query_type(query_type, resolved_fusion_policy),
        rrf_k=resolved_fusion_policy.rrf_k,
    )
    soft_ranked = [hit for hit in soft_hits if hit["id"] not in exact_ids]

    ranked = exact_ranked + soft_ranked
    groups = _group_hits(ranked, group_by)
    return {"hits": ranked, "groups": groups}
def _dedupe_hits(hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    deduped: dict[str, dict[str, Any]] = {}
    for hit in hits:
        current = deduped.get(hit["id"])
        if current is None or float(hit.get("score", 0.0)) > float(current.get("score", 0.0)):
            deduped[hit["id"]] = hit
    return sorted(deduped.values(), key=lambda item: float(item.get("score", 0.0)), reverse=True)


def _group_hits(hits: list[dict[str, Any]], group_by: str) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}
    for hit in hits:
        group_key = str(hit.get(group_by) or f"ungrouped::{hit['id']}")
        bucket = grouped.setdefault(
            group_key,
            {
                "group_key": group_key,
                "document_id": hit.get("document_id"),
                "page_id": hit.get("page_id"),
                "page_number": hit.get("page_number"),
                "hits": [],
            },
        )
        bucket["hits"].append(hit)

    for bucket in grouped.values():
        bucket["hits"].sort(key=lambda item: float(item.get("score", 0.0)), reverse=True)

    return list(grouped.values())


def _has_non_smiles_context(raw_query: str) -> bool:
    """Heuristic mixed-query detector.

    We keep this intentionally simple for the skeleton. Replace with a
    proper tokenizer / parser when you wire the real query classifier.
    """

    separators = {" ", "\n", "\t"}
    return any(ch in separators for ch in raw_query.strip())


def _empty_lane(**_: Any) -> list[dict[str, Any]]:
    return []
