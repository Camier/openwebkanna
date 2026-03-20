"""Core contracts for the lane-based RAG substrate.

This package defines the schema bootstrap, query routing contracts,
and retrieval lane runtime for text and visual evidence.
"""

from .colmodernvbert_profile import (
    ColModernVBERTImagePolicy,
    ColModernVBERTTuningProfile,
    DEFAULT_COLMODERNVBERT_PROFILE,
    get_default_colmodernvbert_profile,
)
from .embed_and_upsert import (
    EmbedAndUpsertConfig,
    EmbedAndUpsertPipeline,
    EmbedAndUpsertResult,
)
from .fusion_policy import FusionPolicy, LaneWeights
from .manifest import EmbeddingModelSpec, IngestionManifest, build_manifest, write_jsonl
from .materialize_evidence import (
    FigureInput,
    make_point_id,
    make_qdrant_point_id,
    materialize_figures,
    materialize_pages,
    MaterializedRecord,
    PageInput,
    SourceDocument,
)
from .qdrant_schema import (
    RagCollectionConfig,
    ensure_payload_indexes,
    ensure_rag_collection,
)
from .text_lane_qdrant import QdrantTextHybridLane, QdrantTextHybridLaneConfig
from .text_query_runtime import (
    NativeTextQueryEncoderConfig,
    NativeTextQueryEncoderRuntime,
    QueryEncodingBundle,
    QueryEncodingRuntimeError,
)
from .vision_colmodernvbert import (
    ColModernVBERTConfig,
    ColModernVBERTLateInteractionRuntime,
    ColModernVBERTVisionEmbedder,
    probe_colmodernvbert_runtime,
)
from .vision_lane_qdrant import QdrantColModernVBERTVisualLane, QdrantVisionLaneConfig

__all__ = [
    "ColModernVBERTConfig",
    "ColModernVBERTImagePolicy",
    "ColModernVBERTLateInteractionRuntime",
    "ColModernVBERTTuningProfile",
    "ColModernVBERTVisionEmbedder",
    "DEFAULT_COLMODERNVBERT_PROFILE",
    "EmbedAndUpsertConfig",
    "EmbedAndUpsertPipeline",
    "EmbedAndUpsertResult",
    "EmbeddingModelSpec",
    "FigureInput",
    "FusionPolicy",
    "IngestionManifest",
    "LaneWeights",
    "make_point_id",
    "make_qdrant_point_id",
    "materialize_figures",
    "materialize_pages",
    "MaterializedRecord",
    "NativeTextQueryEncoderConfig",
    "NativeTextQueryEncoderRuntime",
    "PageInput",
    "QueryEncodingBundle",
    "QueryEncodingRuntimeError",
    "QdrantColModernVBERTVisualLane",
    "QdrantTextHybridLane",
    "QdrantTextHybridLaneConfig",
    "QdrantVisionLaneConfig",
    "RagCollectionConfig",
    "SourceDocument",
    "build_manifest",
    "ensure_payload_indexes",
    "ensure_rag_collection",
    "get_default_colmodernvbert_profile",
    "probe_colmodernvbert_runtime",
    "write_jsonl",
]
