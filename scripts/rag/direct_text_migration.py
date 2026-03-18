#!/usr/bin/env python3
"""Direct text migration without schema validation - for already-created collections."""

import os
import sys
from pathlib import Path

from qdrant_client import QdrantClient
from qdrant_client.http import models

# Setup
QDRANT_URL = os.environ.get("QDRANT_URL", "http://127.0.0.1:6333")
with open("/home/miko/.config/qdrant/qdrant.env") as f:
    for line in f:
        if line.startswith("QDRANT"):
            key, val = line.strip().split("=", 1)
            os.environ[key] = val

QDRANT_API_KEY = os.environ.get("QDRANT__SERVICE__API_KEY")

SOURCE_COLLECTION = "pdf_nemotron_hybrid"
TARGET_COLLECTION = "rag_evidence"
SOURCE_DENSE = "dense_prefetch"
SOURCE_SPARSE = "sparse_text"
TARGET_DENSE = "text_dense"
TARGET_SPARSE = "text_sparse"
BATCH_SIZE = 64

def main():
    client = QdrantClient(url=QDRANT_URL, api_key=QDRANT_API_KEY)

    print(f"=== Migrating text from {SOURCE_COLLECTION} → {TARGET_COLLECTION} ===")

    # Get source collection info
    source_info = client.get_collection(SOURCE_COLLECTION)
    print(f"Source: {source_info.points_count} points")

    # Migrate in batches
    offset = None
    total_migrated = 0
    batch_num = 0

    while True:
        batch_num += 1
        results, offset = client.scroll(
            collection_name=SOURCE_COLLECTION,
            limit=BATCH_SIZE,
            offset=offset,
            with_vectors=True,
            with_payload=True,
        )

        if not results:
            break

        # Transform points for target collection
        points = []
        for r in results:
            # Build new vector dict with renamed vectors
            new_vectors = {}
            if hasattr(r, 'vector') and r.vector:
                # Handle named vectors
                if isinstance(r.vector, dict):
                    if SOURCE_DENSE in r.vector:
                        new_vectors[TARGET_DENSE] = r.vector[SOURCE_DENSE]
                # Also check for vectors attribute
                elif hasattr(r, 'vectors') and r.vectors:
                    if isinstance(r.vectors, dict):
                        if SOURCE_DENSE in r.vectors:
                            new_vectors[TARGET_DENSE] = r.vectors[SOURCE_DENSE]

            # Build payload - add object_type
            payload = dict(r.payload) if r.payload else {}
            payload["object_type"] = "page"

            points.append(
                models.PointStruct(
                    id=r.id,
                    vector=new_vectors,
                    payload=payload,
                )
            )

        # Upsert to target
        client.upsert(
            collection_name=TARGET_COLLECTION,
            points=points,
        )

        total_migrated += len(points)
        print(f"Batch {batch_num}: migrated {len(points)} points (total: {total_migrated})")

        if offset is None:
            break

    print(f"\n=== Migration complete: {total_migrated} points migrated ===")

    # Verify
    target_info = client.get_collection(TARGET_COLLECTION)
    print(f"Target collection now has {target_info.points_count} points")

if __name__ == "__main__":
    main()
