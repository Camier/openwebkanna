#!/usr/bin/env python3
"""Copy legacy dense/sparse page points into the one-collection RAG index.

This migration avoids a second bespoke text-embedding pipeline:
- read existing page points from the legacy hybrid collection
- rename vector families to the canonical one-collection names
- normalize payloads for `object_type=page`
- upsert into `rag_evidence`

It relies on official Qdrant collection/scroll/upsert behavior and preserves
the existing dense+sparse text vectors already used by the current stack.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Any

from qdrant_client import QdrantClient
from qdrant_client.http import models

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from services.rag.materialize_evidence import make_point_id, make_qdrant_point_id
from services.rag.qdrant_schema import RagCollectionConfig, ensure_payload_indexes, ensure_rag_collection


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qdrant-url", default=os.environ.get("QDRANT_URL", "http://127.0.0.1:6335"))
    parser.add_argument("--qdrant-api-key", default=os.environ.get("QDRANT_API_KEY"))
    parser.add_argument("--source-collection", default=os.environ.get("QDRANT_COLLECTION", "pdf_nemotron_hybrid"))
    parser.add_argument("--target-collection", default="rag_evidence")
    parser.add_argument("--source-dense-vector", default=os.environ.get("QDRANT_VECTOR_NAME") or "dense_prefetch")
    parser.add_argument("--source-sparse-vector", default=os.environ.get("QDRANT_SPARSE_VECTOR_NAME", "sparse_text"))
    parser.add_argument("--target-dense-vector", default="text_dense")
    parser.add_argument("--target-sparse-vector", default="text_sparse")
    parser.add_argument("--vision-li-dim", type=int, default=128)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--limit", type=int, default=0, help="Optional max number of source points to migrate.")
    parser.add_argument("--output-report", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    client = QdrantClient(url=args.qdrant_url, api_key=args.qdrant_api_key)
    dense_dim = resolve_dense_dimension(client, args.source_collection, args.source_dense_vector)

    ensure_rag_collection(
        client,
        RagCollectionConfig(
            collection_name=args.target_collection,
            text_dense_dim=dense_dim,
            vision_li_dim=args.vision_li_dim,
        ),
    )
    ensure_payload_indexes(client, args.target_collection)

    migrated = 0
    cursor: Any = None
    while True:
        if args.limit and migrated >= args.limit:
            break

        points, cursor = client.scroll(
            collection_name=args.source_collection,
            scroll_filter=None,
            limit=min(args.batch_size, args.limit - migrated) if args.limit else args.batch_size,
            with_payload=True,
            with_vectors=True,
            offset=cursor,
        )
        if not points:
            break

        target_points: list[models.PointStruct] = []
        for point in points:
            payload = dict(point.payload or {})
            dense_vector = extract_named_vector(point.vector, args.source_dense_vector)
            sparse_vector = extract_sparse_vector(point.vector, args.source_sparse_vector)
            if dense_vector is None or sparse_vector is None:
                continue
            normalized_payload = normalize_payload(payload)
            page_id = str(normalized_payload.get("page_id") or point.id)
            point_id = str(normalized_payload.get("point_id") or make_point_id("page", page_id))
            normalized_payload["point_id"] = point_id
            target_points.append(
                models.PointStruct(
                    id=make_qdrant_point_id(point_id),
                    payload=normalized_payload,
                    vector={
                        args.target_dense_vector: dense_vector,
                        args.target_sparse_vector: sparse_vector,
                    },
                )
            )

        if target_points:
            client.upsert(collection_name=args.target_collection, points=target_points, wait=True)
            migrated += len(target_points)

        if cursor is None:
            break

    report = {
        "source_collection": args.source_collection,
        "target_collection": args.target_collection,
        "source_dense_vector": args.source_dense_vector,
        "source_sparse_vector": args.source_sparse_vector,
        "target_dense_vector": args.target_dense_vector,
        "target_sparse_vector": args.target_sparse_vector,
        "migrated_points": migrated,
    }
    if args.output_report is not None:
        args.output_report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


def resolve_dense_dimension(client: QdrantClient, collection_name: str, dense_vector_name: str) -> int:
    collection = client.get_collection(collection_name)
    vectors = getattr(collection.config.params, "vectors", None) or {}
    if isinstance(vectors, dict) and dense_vector_name in vectors:
        return int(vectors[dense_vector_name].size)
    raise SystemExit(
        f"Could not resolve dense vector dimension for {dense_vector_name!r} in {collection_name!r}."
    )


def extract_named_vector(vector: Any, vector_name: str) -> list[float] | None:
    if isinstance(vector, dict):
        value = vector.get(vector_name)
        return list(value) if value is not None else None
    return None


def extract_sparse_vector(vector: Any, vector_name: str) -> models.SparseVector | None:
    if not isinstance(vector, dict):
        return None
    value = vector.get(vector_name)
    if value is None:
        return None
    if isinstance(value, dict):
        raw_indices = value.get("indices", [])
        raw_values = value.get("values", [])
    else:
        raw_indices = getattr(value, "indices", [])
        raw_values = getattr(value, "values", [])
    indices = list(raw_indices or [])
    values = list(raw_values or [])
    if not indices or not values:
        return None
    return models.SparseVector(indices=indices, values=values)


def normalize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    normalized = dict(payload)
    normalized["object_type"] = "page"
    document_id = normalized.get("document_id") or normalized.get("doc_id") or normalized.get("source_pdf_sha256")
    normalized["document_id"] = str(document_id) if document_id is not None else None
    if not normalized.get("document_title"):
        normalized["document_title"] = normalized.get("title") or normalized.get("file_name")
    page_number = normalized.get("page_number")
    if normalized.get("page_id") is None and document_id is not None and page_number is not None:
        try:
            normalized["page_id"] = f"{document_id}_p{int(page_number):04d}"
        except (TypeError, ValueError):
            normalized["page_id"] = f"{document_id}_p{page_number}"
    if normalized.get("point_id") is None and normalized.get("page_id") is not None:
        normalized["point_id"] = make_point_id("page", str(normalized["page_id"]))
    return normalized


if __name__ == "__main__":
    raise SystemExit(main())
