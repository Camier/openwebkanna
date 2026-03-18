"""Core contracts for the one-collection RAG substrate.

This package is intentionally small. It defines the schema bootstrap,
SMILES normalization, and query routing contracts needed to build the
single-collection multimodal / hybrid / SMILES retrieval layer.
"""

from .normalize_smiles import NormalizedSmiles, normalize_smiles
from .colmodernvbert_profile import (
    ColModernVBERTImagePolicy,
    ColModernVBERTTuningProfile,
    DEFAULT_COLMODERNVBERT_PROFILE,
    get_default_colmodernvbert_profile,
)
from .fusion_policy import FusionPolicy, LaneWeights
from .live_one_collection_router import OneCollectionRouter, OneCollectionRouterResult
from .materialize_evidence import (
    FigureInput,
    MaterializedRecord,
    MoleculeInput,
    PageInput,
    SourceDocument,
    materialize_figures,
    materialize_molecules,
    materialize_pages,
)
from .manifest import EmbeddingModelSpec, IngestionManifest, build_manifest, write_jsonl
from .qdrant_schema import RagCollectionConfig, ensure_payload_indexes, ensure_rag_collection
from .query_router import ParsedQuery, QueryRouter, QueryRouterConfig
from .text_lane_qdrant import QdrantTextHybridLane, QdrantTextHybridLaneConfig
from .text_query_runtime import (
    NativeTextQueryEncoderConfig,
    NativeTextQueryEncoderRuntime,
    QueryEncodingBundle,
    QueryEncodingRuntimeError,
)
from .embed_and_upsert import EmbedAndUpsertConfig, EmbedAndUpsertPipeline, EmbedAndUpsertResult
from .vision_colmodernvbert import (
    ColModernVBERTConfig,
    ColModernVBERTLateInteractionRuntime,
    ColModernVBERTVisionEmbedder,
    probe_colmodernvbert_runtime,
)
from .vision_lane_qdrant import QdrantColModernVBERTVisualLane, QdrantVisionLaneConfig
from .exact_chemistry_lane_qdrant import (
    QdrantExactChemistryLane,
    QdrantExactChemistryLaneConfig,
    probe_exact_chemistry_runtime,
)

__all__ = [
    "ColModernVBERTConfig",
    "ColModernVBERTImagePolicy",
    "ColModernVBERTLateInteractionRuntime",
    "ColModernVBERTTuningProfile",
    "ColModernVBERTVisionEmbedder",
    "DEFAULT_COLMODERNVBERT_PROFILE",
    "EmbeddingModelSpec",
    "EmbedAndUpsertConfig",
    "EmbedAndUpsertPipeline",
    "EmbedAndUpsertResult",
    "FigureInput",
    "FusionPolicy",
    "IngestionManifest",
    "LaneWeights",
    "MaterializedRecord",
    "MoleculeInput",
    "NormalizedSmiles",
    "NativeTextQueryEncoderConfig",
    "NativeTextQueryEncoderRuntime",
    "OneCollectionRouter",
    "OneCollectionRouterResult",
    "PageInput",
    "ParsedQuery",
    "QueryEncodingBundle",
    "QueryEncodingRuntimeError",
    "QdrantExactChemistryLane",
    "QdrantExactChemistryLaneConfig",
    "QdrantColModernVBERTVisualLane",
    "QdrantTextHybridLane",
    "QdrantTextHybridLaneConfig",
    "QdrantVisionLaneConfig",
    "QueryRouter",
    "QueryRouterConfig",
    "RagCollectionConfig",
    "SourceDocument",
    "build_manifest",
    "ensure_payload_indexes",
    "ensure_rag_collection",
    "get_default_colmodernvbert_profile",
    "materialize_figures",
    "materialize_molecules",
    "materialize_pages",
    "normalize_smiles",
    "probe_colmodernvbert_runtime",
    "probe_exact_chemistry_runtime",
    "write_jsonl",
]
