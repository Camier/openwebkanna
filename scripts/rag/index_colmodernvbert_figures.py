#!/usr/bin/env python3
"""Index figure multivectors with official ColModernVBERT + Qdrant.

This script intentionally stays narrow:
- reads materialized figure JSONL records
- embeds figure images with Hugging Face Transformers ColModernVBERT
- upserts figure points into a named Qdrant multivector

It does not invent a second ingestion framework. It relies on:
- official Transformers ColModernVBERT inference
- official Qdrant named multivectors
- the existing one-collection payload contract in `services.rag`
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, Iterable

from qdrant_client import QdrantClient
from qdrant_client.http.exceptions import UnexpectedResponse
from qdrant_client.http.models import PointStruct

from services.rag.colmodernvbert_profile import get_default_colmodernvbert_profile
from services.rag.materialize_evidence import make_point_id, make_qdrant_point_id
from services.rag.qdrant_schema import RagCollectionConfig, ensure_payload_indexes, ensure_rag_collection
from services.rag.vision_colmodernvbert import ColModernVBERTConfig, ColModernVBERTLateInteractionRuntime


def parse_args() -> argparse.Namespace:
    profile = get_default_colmodernvbert_profile()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, help="Path to a materialization manifest JSON.")
    parser.add_argument("--figures-jsonl", type=Path, help="Path to materialized figure JSONL.")
    parser.add_argument("--qdrant-url", default=os.environ.get("QDRANT_URL", "http://localhost:6333"))
    parser.add_argument("--qdrant-api-key", default=os.environ.get("QDRANT_API_KEY"))
    parser.add_argument("--collection-name", default="rag_evidence")
    parser.add_argument("--vision-vector-name", default="vision_li")
    parser.add_argument("--text-dense-dim", type=int, required=True)
    parser.add_argument("--chem-dense-dim", type=int, default=0)
    parser.add_argument("--model-name", default=profile.model_name)
    parser.add_argument("--batch-size", type=int, default=profile.embed_batch_size)
    parser.add_argument("--upsert-batch-size", type=int, default=8)
    parser.add_argument("--device", default=None)
    parser.add_argument("--torch-dtype", default="float32")
    parser.add_argument("--attn-implementation", default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figures_jsonl = resolve_figures_jsonl(args)
    figure_records = list(read_jsonl(figures_jsonl))
    if not figure_records:
        raise SystemExit(f"No figure records found in {figures_jsonl}")

    runtime = ColModernVBERTLateInteractionRuntime(
        ColModernVBERTConfig(
            model_name=args.model_name,
            batch_size=args.batch_size,
            device=args.device,
            torch_dtype=args.torch_dtype,
            attn_implementation=args.attn_implementation,
        )
    )

    probe_vector = runtime.embed_images([resolve_image_uri(figure_records[0])])[0]
    if not probe_vector or not probe_vector[0]:
        raise SystemExit("ColModernVBERT returned an empty multivector for the probe figure.")

    client = QdrantClient(url=args.qdrant_url, api_key=args.qdrant_api_key)
    ensure_rag_collection(
        client,
        RagCollectionConfig(
            collection_name=args.collection_name,
            text_dense_dim=args.text_dense_dim,
            vision_li_dim=len(probe_vector[0]),
            chem_dense_dim=max(args.chem_dense_dim, 1),
            enable_chem_dense=args.chem_dense_dim > 0,
        ),
    )
    ensure_payload_indexes(client, args.collection_name)

    indexed = 0
    for batch in chunked(figure_records, args.upsert_batch_size):
        image_uris = [resolve_image_uri(record) for record in batch]
        embeddings = runtime.embed_images(image_uris)
        points = [
            PointStruct(
                id=make_qdrant_point_id(resolve_point_id(record)),
                payload=build_payload(record),
                vector={args.vision_vector_name: embedding},
            )
            for record, embedding in zip(batch, embeddings, strict=True)
        ]
        upsert_points_resiliently(
            client,
            collection_name=args.collection_name,
            points=points,
        )
        indexed += len(points)

    print(
        json.dumps(
            {
                "collection_name": args.collection_name,
                "vision_vector_name": args.vision_vector_name,
                "figures_indexed": indexed,
                "figure_source": str(figures_jsonl),
                "model_name": args.model_name,
            },
            indent=2,
        )
    )
    return 0


def resolve_figures_jsonl(args: argparse.Namespace) -> Path:
    if args.figures_jsonl is not None:
        return args.figures_jsonl
    if args.manifest is None:
        raise SystemExit("Pass either --figures-jsonl or --manifest.")

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    artifact_path = (
        manifest.get("artifacts", {}).get("figure_records")
        or manifest.get("artifacts", {}).get("figures")
    )
    if not artifact_path:
        raise SystemExit(f"Manifest does not define figure record artifacts: {args.manifest}")
    return Path(artifact_path)


def read_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            yield json.loads(line)


def resolve_image_uri(record: dict[str, Any]) -> str:
    image_uri = record.get("image_uri") or record.get("payload", {}).get("image_uri")
    if not image_uri:
        raise SystemExit(f"Figure record is missing image_uri: {record.get('point_id')}")
    return str(image_uri)


def resolve_point_id(record: dict[str, Any]) -> str:
    if record.get("point_id"):
        return str(record["point_id"])
    payload = dict(record.get("payload", {}) or {})
    if payload.get("point_id"):
        return str(payload["point_id"])
    figure_id = payload.get("figure_id")
    if figure_id:
        return make_point_id("figure", str(figure_id))
    raise SystemExit(f"Figure record is missing point_id/figure_id: {record}")


def build_payload(record: dict[str, Any]) -> dict[str, Any]:
    payload = dict(record.get("payload", {}) or {})
    payload["point_id"] = resolve_point_id(record)
    return payload


def chunked(items: list[dict[str, Any]], size: int) -> Iterable[list[dict[str, Any]]]:
    for index in range(0, len(items), size):
        yield items[index : index + size]


def upsert_points_resiliently(
    client: QdrantClient,
    *,
    collection_name: str,
    points: list[PointStruct],
) -> None:
    try:
        client.upsert(collection_name=collection_name, points=points, wait=True)
        return
    except UnexpectedResponse as exc:
        if not is_payload_too_large_error(exc) or len(points) == 1:
            raise

    midpoint = max(len(points) // 2, 1)
    upsert_points_resiliently(client, collection_name=collection_name, points=points[:midpoint])
    upsert_points_resiliently(client, collection_name=collection_name, points=points[midpoint:])


def is_payload_too_large_error(exc: UnexpectedResponse) -> bool:
    message = str(exc).lower()
    return "payload error" in message and "larger than allowed" in message


if __name__ == "__main__":
    raise SystemExit(main())
