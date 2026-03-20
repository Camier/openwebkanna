from __future__ import annotations

from functools import lru_cache
from typing import Any

from services.multimodal_retrieval_api.contracts import RetrievalBackendInfo
from services.multimodal_retrieval_api.settings import ServiceSettings


def build_one_collection_backend_info(
    *,
    service_settings: ServiceSettings,
) -> RetrievalBackendInfo:
    return RetrievalBackendInfo(
        service="qdrant.text_hybrid",
        qdrant_url=service_settings.qdrant_url,
        qdrant_collection=service_settings.rag_collection_name,
        dense_vector_name=service_settings.text_dense_vector_name,
        sparse_vector_name=service_settings.text_sparse_vector_name,
        late_vector_name=service_settings.vision_vector_name,
    )


@lru_cache(maxsize=1)
def load_qdrant_client() -> Any:
    service_settings = ServiceSettings.load()
    from qdrant_client import QdrantClient

    return QdrantClient(
        url=service_settings.qdrant_url,
        api_key=service_settings.qdrant_api_key,
    )


@lru_cache(maxsize=1)
def load_text_query_runtime() -> object:
    service_settings = ServiceSettings.load()
    from services.rag import NativeTextQueryEncoderConfig, NativeTextQueryEncoderRuntime

    return NativeTextQueryEncoderRuntime(
        NativeTextQueryEncoderConfig(
            model_name_or_path=service_settings.text_query_model_path,
            device=service_settings.text_query_device,
            attn_implementation=service_settings.text_query_attn_implementation,
            sparse_model_name=service_settings.text_sparse_model_name,
        )
    )


@lru_cache(maxsize=1)
def load_text_lane() -> Any:
    service_settings = ServiceSettings.load()
    from services.rag import QdrantTextHybridLane, QdrantTextHybridLaneConfig

    return QdrantTextHybridLane(
        qdrant_client=load_qdrant_client(),
        encoder_runtime=load_text_query_runtime(),
        config=QdrantTextHybridLaneConfig(
            collection_name=service_settings.rag_collection_name,
            dense_vector_name=service_settings.text_dense_vector_name,
            sparse_vector_name=service_settings.text_sparse_vector_name,
        ),
    )


@lru_cache(maxsize=1)
def load_colmodernvbert_visual_lane() -> Any:
    service_settings = ServiceSettings.load()

    from services.rag import (
        ColModernVBERTConfig,
        ColModernVBERTLateInteractionRuntime,
        QdrantColModernVBERTVisualLane,
        QdrantVisionLaneConfig,
    )

    visual_runtime = ColModernVBERTLateInteractionRuntime(
        ColModernVBERTConfig(),
    )
    lane_config = QdrantVisionLaneConfig(
        collection_name=service_settings.rag_collection_name,
        vector_name=service_settings.vision_vector_name,
    )
    return QdrantColModernVBERTVisualLane(
        qdrant_client=load_qdrant_client(),
        runtime=visual_runtime,
        config=lane_config,
    )


@lru_cache(maxsize=1)
def load_retrieval_lanes() -> dict[str, Any]:
    """Return a dict of available retrieval lanes.

    Returns text and visual lanes when configured; exact-chemistry lane removed.
    """
    service_settings = ServiceSettings.load()
    text_lane = load_text_lane()
    visual_lane = None
    try:
        visual_lane = load_colmodernvbert_visual_lane()
    except Exception:
        pass
    return {
        "text": text_lane,
        "visual": visual_lane,
    }
