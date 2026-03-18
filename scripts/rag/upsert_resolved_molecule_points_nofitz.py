#!/usr/bin/env python3
"""Upsert resolved molecule points into the one-collection Qdrant index."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# import fitz
import pandas as pd
from qdrant_client import QdrantClient
from qdrant_client.http import models

from services.multimodal_retrieval_api.settings import ServiceSettings
from services.rag.materialize_evidence import make_point_id, make_qdrant_point_id
from services.rag.normalize_smiles import normalize_smiles


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resolved-records",
        type=Path,
        default=Path("artifacts/rag/smiles_identity/ambiguous_resolved_external_records.parquet"),
        help="Resolved external molecule records parquet.",
    )
    parser.add_argument(
        "--collection-name",
        default=None,
        help="Override target collection name. Defaults to ServiceSettings.rag_collection_name.",
    )
    parser.add_argument(
        "--pdf-dir",
        type=Path,
        default=Path("pdfs"),
        help="Directory of local source PDFs used to resolve file-stem titles to canonical PDF titles.",
    )
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def normalize_title(value: str | None) -> str | None:
    if not isinstance(value, str):
        return None
    normalized = value.replace("—", "-").replace("–", "-").replace("’", "'")
    normalized = " ".join(normalized.strip().split()).casefold()
    return normalized or None


def loose_title_key(value: str | None) -> str | None:
    normalized = normalize_title(value)
    if not normalized:
        return None
    normalized = re.sub(r"[^a-z0-9]+", " ", normalized)
    normalized = " ".join(normalized.split())
    return normalized or None


def iter_title_keys(value: str | None) -> list[str]:
    normalized = normalize_title(value)
    if not normalized:
        return []
    keys: list[str] = [normalized]
    parts = [part.strip() for part in normalized.split(" - ") if part.strip()]
    if len(parts) >= 3:
        for start in range(1, len(parts) - 1):
            candidate = " - ".join(parts[start:])
            if candidate and candidate not in keys:
                keys.append(candidate)
        tail = " - ".join(parts[2:])
        if tail and tail not in keys:
            keys.append(tail)
    for key in list(keys):
        loose_key = loose_title_key(key)
        if loose_key and loose_key not in keys:
            keys.append(loose_key)
    return keys


def build_pdf_title_aliases(pdf_dir: Path) -> dict[str, list[str]]:
    return {}
    for pdf_path in sorted(pdf_dir.glob("*.pdf")):
        stem_key = normalize_title(pdf_path.stem)
        if not stem_key:
            continue
        alias_keys: list[str] = []
        try:
            with fitz.open(pdf_path) as doc:
                metadata = doc.metadata or {}
        except Exception:
            metadata = {}
        for candidate in (pdf_path.stem, metadata.get("title")):
            for key in iter_title_keys(candidate):
                if key and key not in alias_keys:
                    alias_keys.append(key)
        if alias_keys:
            aliases[stem_key] = alias_keys
    return aliases


def parse_block_id(source_record_id: str) -> str | None:
    parts = source_record_id.split("|")
    if len(parts) < 2:
        return None
    block_id = parts[1].strip()
    return block_id or None


def parse_candidate_seq(source_record_id: str) -> int:
    parts = source_record_id.split("|")
    if len(parts) >= 3:
        try:
            return int(parts[2])
        except ValueError:
            pass
    return 1


def parse_page_index(block_id: str | None) -> int | None:
    if not block_id:
        return None
    match = re.search(r"/page/(\d+)/", block_id)
    return int(match.group(1)) if match else None


def parse_block_seq(block_id: str | None) -> int:
    if not block_id:
        return 1
    match = re.search(r"/[^/]+/(\d+)$", block_id)
    return int(match.group(1)) if match else 1


def resolve_dense_dim(client: QdrantClient, collection_name: str, vector_name: str) -> int:
    collection = client.get_collection(collection_name)
    vectors = getattr(collection.config.params, "vectors", None) or {}
    if isinstance(vectors, dict) and vector_name in vectors:
        return int(vectors[vector_name].size)
    raise RuntimeError(f"Missing vector {vector_name!r} in collection {collection_name!r}.")


def build_page_indexes(
    client: QdrantClient,
    collection_name: str,
) -> tuple[dict[tuple[str, int], dict[str, Any]], dict[tuple[str, str], dict[str, Any]]]:
    page_lookup: dict[tuple[str, int], dict[str, Any]] = {}
    figure_lookup: dict[tuple[str, str], dict[str, Any]] = {}
    offset: Any = None

    while True:
        points, offset = client.scroll(
            collection_name=collection_name,
            scroll_filter=models.Filter(
                must=[
                    models.FieldCondition(
                        key="object_type",
                        match=models.MatchValue(value="page"),
                    )
                ]
            ),
            with_payload=True,
            with_vectors=False,
            limit=256,
            offset=offset,
        )
        if not points:
            break
        for point in points:
            payload = dict(getattr(point, "payload", {}) or {})
            title_keys = iter_title_keys(payload.get("document_title") or payload.get("title"))
            page_index = payload.get("page_number")
            if not title_keys or not isinstance(page_index, int):
                continue
            page_index_zero = page_index - 1
            page_entry = {
                "document_id": payload.get("document_id"),
                "document_title": payload.get("document_title"),
                "page_id": payload.get("page_id"),
                "page_number": payload.get("page_number"),
                "source_uri": payload.get("source_uri"),
            }
            for title_key in title_keys:
                page_lookup[(title_key, page_index_zero)] = page_entry
            figure_records = payload.get("figure_records")
            if isinstance(figure_records, list):
                for record in figure_records:
                    if not isinstance(record, dict):
                        continue
                    block_id = record.get("block_id")
                    if not isinstance(block_id, str) or not block_id:
                        continue
                    for title_key in title_keys:
                        figure_lookup[(title_key, block_id)] = {
                            **page_entry,
                            "figure_id": record.get("figure_id"),
                            "block_id": block_id,
                            "bbox": record.get("bbox"),
                            "image_uri": record.get("image_uri") or record.get("image_path"),
                            "figure_kind": record.get("figure_kind"),
                            "caption_text": record.get("caption_text"),
                            "ocr_text": record.get("text"),
                        }
        if offset is None:
            break

    return page_lookup, figure_lookup


def build_molecule_payload(
    row: pd.Series,
    *,
    page_entry: dict[str, Any],
    figure_entry: dict[str, Any] | None,
    pipeline_version: str,
    created_at: str,
) -> dict[str, Any]:
    block_id = parse_block_id(str(row["source_record_id"]))
    block_seq = parse_block_seq(block_id)
    molecule_seq = parse_candidate_seq(str(row["source_record_id"]))
    page_number = int(page_entry["page_number"])
    document_id = str(page_entry["document_id"])
    molecule_id = f"mol_{document_id}_{page_number:04d}_{block_seq:02d}_{molecule_seq:02d}"
    figure_id = figure_entry.get("figure_id") if figure_entry else None
    image_uri = figure_entry.get("image_uri") if figure_entry else None
    bbox = figure_entry.get("bbox") if figure_entry else None
    source_smiles = row.get("raw_smiles") or row.get("canonical_smiles") or ""
    normalized = normalize_smiles(str(source_smiles))

    return {
        "point_id": make_point_id("molecule", molecule_id),
        "object_type": "molecule",
        "document_id": document_id,
        "document_title": page_entry.get("document_title"),
        "page_id": page_entry.get("page_id"),
        "page_number": page_number,
        "parent_id": figure_id or page_entry.get("page_id"),
        "figure_id": figure_id,
        "block_id": block_id,
        "molecule_id": molecule_id,
        "source_uri": page_entry.get("source_uri"),
        "pipeline_version": pipeline_version,
        "created_at": created_at,
        "raw_smiles": normalized.raw_smiles,
        "canonical_smiles": normalized.canonical_smiles,
        "inchikey": normalized.inchikey,
        "formula": None,
        "backend": row["external_resolution_source"],
        "confidence": None,
        "review_status": normalized.review_status,
        "has_smiles": bool(normalized.raw_smiles),
        "normalization_error": normalized.error,
        "source_text": row["candidate_name"],
        "image_uri": image_uri,
        "bbox": bbox,
        "figure_kind": figure_entry.get("figure_kind") if figure_entry else None,
        "caption_text": figure_entry.get("caption_text") if figure_entry else None,
        "ocr_text": figure_entry.get("ocr_text") if figure_entry else None,
        "chembl_ids": row.get("chembl_ids"),
        "pubchem_cids": row.get("pubchem_cids"),
        "external_resolution_source": row["external_resolution_source"],
    }


def ensure_rdkit_normalization_ready() -> None:
    probe = normalize_smiles("C")
    if probe.review_status == "unavailable":
        raise RuntimeError(
            "RDKit normalization is unavailable in this Python runtime. "
            "Re-run molecule upsert with an interpreter that can import rdkit."
        )


def main() -> int:
    args = parse_args()
    ensure_rdkit_normalization_ready()
    settings = ServiceSettings.load()
    client = QdrantClient(url=settings.qdrant_url, api_key=settings.qdrant_api_key)
    collection_name = args.collection_name or settings.rag_collection_name
    dense_dim = resolve_dense_dim(client, collection_name, settings.text_dense_vector_name)
    page_lookup, figure_lookup = build_page_indexes(client, collection_name)
    pdf_title_aliases = build_pdf_title_aliases(args.pdf_dir)

    records = pd.read_parquet(args.resolved_records)
    created_at = utc_now_iso()
    pipeline_version = "rag-v1.0.0-external-molecule-resolution"
    zero_dense = [0.0] * dense_dim
    points: list[models.PointStruct] = []

    matched_pages = 0
    matched_figures = 0
    skipped = 0

    for row in records.to_dict(orient="records"):
        title_keys = iter_title_keys(row.get("source"))
        source_title_key = normalize_title(row.get("source"))
        if source_title_key and source_title_key in pdf_title_aliases:
            for alias_key in pdf_title_aliases[source_title_key]:
                if alias_key not in title_keys:
                    title_keys.append(alias_key)
        block_id = parse_block_id(str(row["source_record_id"]))
        page_index = parse_page_index(block_id)
        if not title_keys or page_index is None:
            skipped += 1
            continue
        page_entry = None
        title_key = None
        for candidate_title_key in title_keys:
            page_entry = page_lookup.get((candidate_title_key, page_index))
            if page_entry is not None:
                title_key = candidate_title_key
                break
        if page_entry is None:
            skipped += 1
            continue
        matched_pages += 1
        figure_entry = figure_lookup.get((title_key, block_id)) if block_id and title_key else None
        if figure_entry is not None:
            matched_figures += 1
        payload = build_molecule_payload(
            row,
            page_entry=page_entry,
            figure_entry=figure_entry,
            pipeline_version=pipeline_version,
            created_at=created_at,
        )
        points.append(
            models.PointStruct(
                id=make_qdrant_point_id(payload["point_id"]),
                payload=payload,
                vector={settings.text_dense_vector_name: zero_dense},
            )
        )

    if args.apply and points:
        for start in range(0, len(points), args.batch_size):
            client.upsert(
                collection_name=collection_name,
                points=points[start : start + args.batch_size],
            )

    summary = {
        "collection_name": collection_name,
        "resolved_records_input": int(len(records)),
        "molecule_points_prepared": int(len(points)),
        "matched_pages": matched_pages,
        "matched_figures": matched_figures,
        "skipped": skipped,
        "apply": bool(args.apply),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
