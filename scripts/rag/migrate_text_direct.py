#!/usr/bin/env python3
"""Direct migration that bypasses schema validation - for collections with extra vectors."""

import argparse
import json
import os
import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from qdrant_client import QdrantClient
from qdrant_client.http import models


def extract_named_vector(vector, vector_name):
    """Extract named vector from point."""
    if isinstance(vector, dict):
        value = vector.get(vector_name)
        return list(value) if value is not None else None
    return None


def extract_sparse_vector(vector, vector_name):
    """Extract sparse vector from point."""
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


def normalize_payload(payload):
    """Normalize payload for one-collection format."""
    normalized = dict(payload or {})
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
    return normalized


def main():
    parser = argparse.ArgumentParser(description="Direct text migration without schema validation")
    parser.add_argument("--source-collection", default="pdf_nemotron_hybrid")
    parser.add_argument("--target-collection", default="rag_evidence")
    parser.add_argument("--source-dense-vector", default="dense_prefetch")
    parser.add_argument("--source-sparse-vector", default="sparse_text")
    parser.add_argument("--target-dense-vector", default="text_dense")
    parser.add_argument("--target-sparse-vector", default="text_sparse")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    # Setup Qdrant client
    QDRANT_URL = os.environ.get("QDRANT_URL", "http://127.0.0.1:6333")
    with open("/home/miko/.config/qdrant/qdrant.env") as f:
        for line in f:
            if line.startswith("QDRANT"):
                key, val = line.strip().split("=", 1)
                os.environ[key] = val
    QDRANT_API_KEY = os.environ.get("QDRANT__SERVICE__API_KEY")

    client = QdrantClient(url=QDRANT_URL, api_key=QDRANT_API_KEY)

    print(f"=== Migrating text from {args.source_collection} → {args.target_collection} ===")

    # Get source info
    source_info = client.get_collection(args.source_collection)
    print(f"Source: {source_info.points_count} points")

    # Migrate in batches
    cursor = None
    migrated = 0
    batch_num = 0

    while True:
        batch_num += 1
        limit = args.batch_size
        if args.limit:
            remaining = args.limit - migrated
            if remaining <= 0:
                break
            limit = min(limit, remaining)

        results, cursor = client.scroll(
            collection_name=args.source_collection,
            limit=limit,
            offset=cursor,
            with_payload=True,
            with_vectors=True,
        )

        if not results:
            break

        target_points = []
        for point in results:
            # Extract vectors
            dense_vector = extract_named_vector(point.vector, args.source_dense_vector)
            sparse_vector = extract_sparse_vector(point.vector, args.source_sparse_vector)

            if dense_vector is None or sparse_vector is None:
                print(f"  Skipping point {point.id} - missing vectors")
                continue

            # Normalize payload
            normalized_payload = normalize_payload(point.payload)

            target_points.append(
                models.PointStruct(
                    id=point.id,
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
            print(f"Batch {batch_num}: migrated {len(target_points)} points (total: {migrated})")

        if cursor is None:
            break

    print(f"\n=== Migration complete: {migrated} points migrated ===")

    # Verify
    target_info = client.get_collection(args.target_collection)
    print(f"Target collection now has {target_info.points_count} points")

    return 0


if __name__ == "__main__":
    sys.exit(main())
