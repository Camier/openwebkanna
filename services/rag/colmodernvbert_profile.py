"""Recommended ColModernVBERT operating profile.

This profile keeps the visual lane narrow and figure-first so it
complements, rather than competes with, the existing text hybrid lane.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .fusion_policy import FusionPolicy


@dataclass(frozen=True)
class ColModernVBERTImagePolicy:
    """Image preparation defaults for figure retrieval."""

    convert_to_rgb: bool = True
    max_long_edge: int = 1536
    min_short_edge: int = 448
    upscale_small_images: bool = True


@dataclass(frozen=True)
class ColModernVBERTTuningProfile:
    """End-to-end operating defaults for the visual lane."""

    model_name: str = "ModernVBERT/colmodernvbert"
    index_object_types: tuple[str, ...] = ("figure",)
    enable_page_visual_indexing: bool = False
    require_vision_for_figures: bool = True
    image_policy: ColModernVBERTImagePolicy = field(default_factory=ColModernVBERTImagePolicy)
    embed_batch_size: int = 4
    qdrant_candidate_pool_per_lane: int = 24
    final_top_k: int = 10
    group_by: str = "page_id"
    prioritized_figure_kinds: tuple[str, ...] = ("scheme", "structure_panel")
    fusion_policy: FusionPolicy = field(
        default_factory=lambda: FusionPolicy(
            rrf_k=60,
            text_lane_weight=1.0,
            visual_lane_weight_text=1.1,
            visual_lane_weight_mixed=1.2,
            visual_lane_weight_smiles=1.3,
            chem_lane_weight=1.05,
        )
    )


DEFAULT_COLMODERNVBERT_PROFILE = ColModernVBERTTuningProfile()


def get_default_colmodernvbert_profile() -> ColModernVBERTTuningProfile:
    """Return the default tuned operating profile."""

    return DEFAULT_COLMODERNVBERT_PROFILE
