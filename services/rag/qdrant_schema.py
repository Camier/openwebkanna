"""Qdrant schema bootstrap for the one-collection multimodal RAG."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class RagCollectionConfig:
    """Schema contract for the one-collection RAG index."""

    collection_name: str = "rag_evidence"
    text_dense_dim: int = 1024
    vision_li_dim: int = 128
    chem_dense_dim: int = 0
    enable_chem_dense: bool = False


def ensure_rag_collection(client: Any, cfg: RagCollectionConfig) -> None:
    """Create or validate the canonical RAG collection.

    The collection hosts heterogeneous evidence points:
    - page
    - figure
    """

    models = _models()
    expected_vectors = _expected_vectors(cfg, models)
    expected_sparse_vectors = {"text_sparse": models.SparseVectorParams()}

    if _collection_exists(client, cfg.collection_name):
        _assert_compatible_collection(client, cfg, expected_vectors, expected_sparse_vectors)
        return

    client.create_collection(
        collection_name=cfg.collection_name,
        vectors_config=expected_vectors,
        sparse_vectors_config=expected_sparse_vectors,
    )


def ensure_payload_indexes(client: Any, collection_name: str) -> None:
    """Create the filter-critical payload indexes for retrieval and grouping."""

    models = _models()
    index_fields = [
        ("object_type", models.PayloadSchemaType.KEYWORD),
        ("document_id", models.PayloadSchemaType.KEYWORD),
        ("page_id", models.PayloadSchemaType.KEYWORD),
        ("page_number", models.PayloadSchemaType.INTEGER),
        ("figure_id", models.PayloadSchemaType.KEYWORD),
        ("figure_kind", models.PayloadSchemaType.KEYWORD),
    ]

    for field_name, field_schema in index_fields:
        client.create_payload_index(
            collection_name=collection_name,
            field_name=field_name,
            field_schema=field_schema,
        )


def validate_collection_shape(client: Any, cfg: RagCollectionConfig) -> dict[str, Any]:
    """Return machine-readable schema compatibility details."""

    models = _models()
    collection = client.get_collection(cfg.collection_name)
    config_params = collection.config.params
    actual_vectors = getattr(config_params, "vectors", None) or {}
    actual_sparse_vectors = getattr(config_params, "sparse_vectors", None) or {}

    expected_vector_names = set(_expected_vectors(cfg, models).keys())
    expected_sparse_names = {"text_sparse"}
    actual_vector_names = set(actual_vectors.keys()) if isinstance(actual_vectors, dict) else set()
    actual_sparse_names = set(actual_sparse_vectors.keys()) if isinstance(actual_sparse_vectors, dict) else set()

    return {
        "collection_name": cfg.collection_name,
        "exists": True,
        "expected_vector_names": sorted(expected_vector_names),
        "actual_vector_names": sorted(actual_vector_names),
        "expected_sparse_vector_names": sorted(expected_sparse_names),
        "actual_sparse_vector_names": sorted(actual_sparse_names),
        "vector_names_match": expected_vector_names == actual_vector_names,
        "sparse_vector_names_match": expected_sparse_names == actual_sparse_names,
    }


def _expected_vectors(cfg: RagCollectionConfig, models: Any) -> dict[str, Any]:
    vectors = {
        "text_dense": models.VectorParams(
            size=cfg.text_dense_dim,
            distance=models.Distance.COSINE,
        ),
        "vision_li": models.VectorParams(
            size=cfg.vision_li_dim,
            distance=models.Distance.COSINE,
            multivector_config=models.MultiVectorConfig(
                comparator=models.MultiVectorComparator.MAX_SIM,
            ),
        ),
    }
    return vectors


def _collection_exists(client: Any, collection_name: str) -> bool:
    return any(c.name == collection_name for c in client.get_collections().collections)


def _assert_compatible_collection(
    client: Any,
    cfg: RagCollectionConfig,
    expected_vectors: dict[str, Any],
    expected_sparse_vectors: dict[str, Any],
) -> None:
    collection = client.get_collection(cfg.collection_name)
    params = collection.config.params
    actual_vectors = getattr(params, "vectors", None) or {}
    actual_sparse_vectors = getattr(params, "sparse_vectors", None) or {}

    actual_vector_names = set(actual_vectors.keys()) if isinstance(actual_vectors, dict) else set()
    actual_sparse_names = set(actual_sparse_vectors.keys()) if isinstance(actual_sparse_vectors, dict) else set()
    expected_vector_names = set(expected_vectors.keys())
    expected_sparse_names = set(expected_sparse_vectors.keys())

    if actual_vector_names != expected_vector_names:
        raise ValueError(
            f"incompatible vector names for {cfg.collection_name}: "
            f"expected={sorted(expected_vector_names)} actual={sorted(actual_vector_names)}"
        )
    if actual_sparse_names != expected_sparse_names:
        raise ValueError(
            f"incompatible sparse vector names for {cfg.collection_name}: "
            f"expected={sorted(expected_sparse_names)} actual={sorted(actual_sparse_names)}"
        )


def _models() -> Any:
    from qdrant_client.http import models

    return models
