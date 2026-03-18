#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from functools import lru_cache
from html import unescape
from pathlib import Path
import re
from typing import Any

from qdrant_client import QdrantClient
from qdrant_client.http import models

from services.multimodal_retrieval_api.settings import ServiceSettings
from scripts.rag.local_evidence_store import LocalDocumentEvidence, load_local_evidence_store


def _normalize_text(value: str | None) -> str | None:
    if not isinstance(value, str):
        return None
    normalized = value.lower().replace("—", "-").replace("–", "-")
    normalized = re.sub(r"[^a-z0-9]+", " ", normalized)
    normalized = " ".join(normalized.split())
    return normalized or None


def _normalize_doi(value: str | None) -> str | None:
    if not isinstance(value, str):
        return None
    normalized = value.strip().lower()
    normalized = normalized.removeprefix("https://doi.org/")
    normalized = normalized.removeprefix("http://doi.org/")
    normalized = normalized.removeprefix("doi:")
    return normalized or None


def _coerce_int(value: Any) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _coerce_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _resolve_image_path(repo_root: Path, storage_uri: str | None) -> str | None:
    if not storage_uri:
        return None
    candidate = Path(storage_uri)
    if candidate.is_absolute():
        return str(candidate)
    return str((repo_root / candidate).resolve())


def _html_to_text(value: str | None) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    text = re.sub(r"<[^>]+>", " ", value)
    text = " ".join(unescape(text).split())
    return text or None


def _classify_evidence_type(block_id: str | None) -> str:
    if isinstance(block_id, str) and "/ChemicalBlock/" in block_id:
        return "chemical_block"
    return "figure"


def _best_document_match(documents: list[LocalDocumentEvidence], payload: dict[str, Any]) -> LocalDocumentEvidence | None:
    if not documents:
        return None
    if len(documents) == 1:
        return documents[0]

    target_year = _coerce_int(payload.get("publication_year"))
    if target_year is not None:
        exact_year = [document for document in documents if document.publication_year == target_year]
        if len(exact_year) == 1:
            return exact_year[0]
        if exact_year:
            documents = exact_year

    target_doi = _normalize_doi(
        payload.get("document_doi") if isinstance(payload.get("document_doi"), str) else None
    )
    if target_doi:
        exact_doi = [document for document in documents if document.doi == target_doi]
        if len(exact_doi) == 1:
            return exact_doi[0]
        if exact_doi:
            documents = exact_doi

    target_title = _normalize_text(
        payload.get("document_title") if isinstance(payload.get("document_title"), str) else payload.get("title")
    )
    if target_title:
        exact_title = [document for document in documents if document.normalized_title == target_title]
        if exact_title:
            return exact_title[0]
    return documents[0]


def _resolve_document(local_store: Any, payload: dict[str, Any]) -> LocalDocumentEvidence | None:
    candidate_ids = [
        payload.get("document_id"),
        payload.get("doc_id"),
        payload.get("source_pdf_sha256"),
    ]
    for candidate_id in candidate_ids:
        if isinstance(candidate_id, str):
            document = local_store.get_document(candidate_id)
            if document is not None:
                return document

    candidate_dois = [
        payload.get("document_doi"),
        payload.get("document_source_url"),
    ]
    for candidate_doi in candidate_dois:
        if not isinstance(candidate_doi, str):
            continue
        matches = local_store.get_documents_by_doi(candidate_doi)
        document = _best_document_match(matches, payload)
        if document is not None:
            return document

    candidate_titles = [
        payload.get("title"),
        payload.get("document_title"),
    ]
    for candidate_title in candidate_titles:
        if not isinstance(candidate_title, str):
            continue
        matches = local_store.get_documents_by_normalized_title(candidate_title)
        document = _best_document_match(matches, payload)
        if document is not None:
            return document
    return None


def _resolve_page_index(payload: dict[str, Any]) -> int | None:
    page_index = _coerce_int(payload.get("page_index"))
    if page_index is not None:
        return page_index
    page_number = _coerce_int(payload.get("page_number"))
    if page_number is None:
        page_number = _coerce_int(payload.get("page"))
    if page_number is None:
        return None
    return page_number - 1 if page_number > 0 else page_number


@lru_cache(maxsize=None)
def _load_json(path_str: str) -> dict[str, Any]:
    path = Path(path_str)
    return json.loads(path.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _load_smiles_standardizer():
    try:
        from smiles.core.standardization import standardize_smiles

        return standardize_smiles
    except Exception:
        smiles_repo_root = os.getenv("MULTIMODAL_RETRIEVAL_API_SMILES_REPO_ROOT")
        if smiles_repo_root:
            src_path = str((Path(smiles_repo_root).expanduser().resolve() / "src"))
            if src_path not in sys.path:
                sys.path.insert(0, src_path)
        try:
            from smiles.core.standardization import standardize_smiles

            return standardize_smiles
        except Exception:
            return None


def _standardize_smiles(smiles_value: str) -> dict[str, Any] | None:
    standardize_smiles = _load_smiles_standardizer()
    if standardize_smiles is None:
        return None
    try:
        record = standardize_smiles(smiles_value)
    except Exception:
        return None
    return {
        "canonical_smiles": getattr(record, "canonical_smiles", None),
        "inchikey": getattr(record, "inchikey", None),
        "record_status": getattr(record, "record_status", None),
    }


def build_page_scroll_filter(models_module: Any) -> Any:
    return models_module.Filter(
        must=[
            models_module.FieldCondition(
                key="object_type",
                match=models_module.MatchValue(value="page"),
            )
        ]
    )


@lru_cache(maxsize=None)
def _load_smiles_by_block_id(extraction_dir_str: str) -> dict[str, list[dict[str, Any]]]:
    extraction_dir = Path(extraction_dir_str)
    smiles_path = extraction_dir / "smiles_extracted.json"
    if not smiles_path.exists():
        return {}
    payload = _load_json(str(smiles_path))
    by_block_id: dict[str, list[dict[str, Any]]] = {}
    seen: set[tuple[str, str]] = set()
    for result in payload.get("results", []):
        if not isinstance(result, dict):
            continue
        block_id = result.get("block_id")
        molecule_smiles = result.get("molecule_smiles")
        if not isinstance(block_id, str) or not isinstance(molecule_smiles, str) or not molecule_smiles.strip():
            continue
        standardized = _standardize_smiles(molecule_smiles.strip())
        canonical_smiles = (
            standardized.get("canonical_smiles")
            if isinstance(standardized, dict) and isinstance(standardized.get("canonical_smiles"), str)
            else None
        )
        inchikey = (
            standardized.get("inchikey")
            if isinstance(standardized, dict) and isinstance(standardized.get("inchikey"), str)
            else None
        )
        record_status = (
            standardized.get("record_status")
            if isinstance(standardized, dict) and isinstance(standardized.get("record_status"), str)
            else None
        )
        raw_candidates = result.get("raw_candidates") if isinstance(result.get("raw_candidates"), list) else []
        if canonical_smiles is None:
            for candidate in raw_candidates:
                if not isinstance(candidate, dict):
                    continue
                candidate_smiles = candidate.get("canonical_smiles")
                if isinstance(candidate_smiles, str) and candidate_smiles.strip():
                    canonical_smiles = candidate_smiles.strip()
                    break
        dedupe_key = (block_id, molecule_smiles.strip())
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        by_block_id.setdefault(block_id, []).append(
            {
                "smiles": molecule_smiles.strip(),
                "canonical_smiles": canonical_smiles,
                "inchikey": inchikey,
                "backend": result.get("backend") if isinstance(result.get("backend"), str) else None,
                "confidence": _coerce_float(result.get("confidence")),
                "status": (
                    result.get("status")
                    if isinstance(result.get("status"), str)
                    else record_status
                ),
            }
        )
    return by_block_id


def _build_page_figure_records(
    *,
    settings: ServiceSettings,
    document: LocalDocumentEvidence,
    page_index: int,
    payload: dict[str, Any],
) -> list[dict[str, Any]]:
    normalized = _load_json(str(document.extraction_dir / "normalized.json"))
    block_index_by_id: dict[str, dict[str, Any]] = {}
    for block in normalized.get("block_index", []):
        if not isinstance(block, dict):
            continue
        block_id = block.get("block_id")
        if isinstance(block_id, str):
            block_index_by_id[block_id] = block

    images_by_asset_id: dict[str, dict[str, Any]] = {}
    for image in normalized.get("images", []):
        if not isinstance(image, dict):
            continue
        asset_id = image.get("asset_id")
        if isinstance(asset_id, str):
            images_by_asset_id[asset_id] = image

    smiles_by_block_id = _load_smiles_by_block_id(str(document.extraction_dir))
    title = document.title
    if not isinstance(title, str) or not title.strip():
        metadata = normalized.get("metadata") if isinstance(normalized.get("metadata"), dict) else {}
        title = metadata.get("title") if isinstance(metadata.get("title"), str) else None
    if not isinstance(title, str) or not title.strip():
        title = payload.get("document_title") if isinstance(payload.get("document_title"), str) else payload.get("title")

    records: list[dict[str, Any]] = []
    for figure in normalized.get("figures", []):
        if not isinstance(figure, dict):
            continue
        block_id = figure.get("block_id") if isinstance(figure.get("block_id"), str) else None
        block = block_index_by_id.get(block_id or "", {})
        figure_page_index = _coerce_int(figure.get("page_index"))
        if figure_page_index is None:
            figure_page_index = _coerce_int(block.get("page_index"))
        if figure_page_index != page_index:
            continue

        image_asset_id = figure.get("image_asset_id") if isinstance(figure.get("image_asset_id"), str) else None
        image_record = images_by_asset_id.get(image_asset_id or "", {})
        storage_uri = image_record.get("storage_uri") if isinstance(image_record.get("storage_uri"), str) else None
        bbox_raw = figure.get("bbox") if isinstance(figure.get("bbox"), list) else block.get("bbox")
        bbox = [_coerce_float(value) for value in bbox_raw] if isinstance(bbox_raw, list) else []
        normalized_bbox = [value for value in bbox if value is not None]
        smiles_records = list(smiles_by_block_id.get(block_id or "", []))
        smiles = [record["smiles"] for record in smiles_records if isinstance(record.get("smiles"), str)]
        canonical_smiles = [
            record["canonical_smiles"]
            for record in smiles_records
            if isinstance(record.get("canonical_smiles"), str) and record["canonical_smiles"]
        ]
        smiles_backends = [
            record["backend"]
            for record in smiles_records
            if isinstance(record.get("backend"), str) and record["backend"]
        ]
        smiles_confidences = [
            record["confidence"]
            for record in smiles_records
            if isinstance(record.get("confidence"), float)
        ]
        smiles_review_status = [
            record["status"]
            for record in smiles_records
            if isinstance(record.get("status"), str) and record["status"]
        ]
        block_type = block.get("block_type") if isinstance(block.get("block_type"), str) else None
        page_number = figure_page_index + 1
        records.append(
            {
                "figure_id": figure.get("figure_id") if isinstance(figure.get("figure_id"), str) else None,
                "document_id": document.doc_id,
                "title": title,
                "page_number": page_number,
                "page_index": figure_page_index,
                "block_id": block_id,
                "block_type": block_type,
                "evidence_type": _classify_evidence_type(block_id),
                "figure_kind": block_type,
                "caption_text": figure.get("caption_text") if isinstance(figure.get("caption_text"), str) else None,
                "text": (
                    block.get("text")
                    if isinstance(block.get("text"), str) and block.get("text")
                    else _html_to_text(block.get("html") if isinstance(block.get("html"), str) else None)
                ),
                "bbox": normalized_bbox,
                "image_asset_id": image_asset_id,
                "image_uri": storage_uri,
                "image_path": _resolve_image_path(settings.repo_root, storage_uri),
                "has_smiles": bool(smiles_records),
                "smiles": smiles,
                "canonical_smiles": canonical_smiles,
                "smiles_backends": smiles_backends,
                "smiles_confidences": smiles_confidences,
                "smiles_review_status": smiles_review_status,
                "smiles_records": smiles_records,
            }
        )
    return records


@dataclass
class PointUpdate:
    point_id: Any
    payload: dict[str, Any]
    doc_id: str
    page_number: int | None
    figure_record_count: int
    title: str | None


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Promote extraction-derived figure and SMILES evidence into page-level Qdrant payloads.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write updates to Qdrant. Without this flag the command stays in dry-run mode.",
    )
    parser.add_argument(
        "--scroll-batch-size",
        type=int,
        default=256,
        help="Number of Qdrant points fetched per scroll page.",
    )
    parser.add_argument(
        "--preview-limit",
        type=int,
        default=5,
        help="Maximum number of payload updates to include in the JSON preview.",
    )
    return parser


def main() -> int:
    args = build_argument_parser().parse_args()
    service_settings = ServiceSettings.load()
    client = QdrantClient(
        url=service_settings.qdrant_url,
        api_key=service_settings.qdrant_api_key,
    )
    local_store = load_local_evidence_store()

    updates: list[PointUpdate] = []
    matched_documents = 0
    points_without_document = 0
    points_without_page = 0
    points_with_figures = 0
    scanned_points = 0

    offset: Any = None
    while True:
        points, offset = client.scroll(
            collection_name=service_settings.rag_collection_name,
            scroll_filter=build_page_scroll_filter(models),
            limit=args.scroll_batch_size,
            offset=offset,
            with_payload=True,
            with_vectors=False,
        )
        if not points:
            break

        for point in points:
            scanned_points += 1
            payload = dict(getattr(point, "payload", {}) or {})
            document = _resolve_document(local_store, payload)
            if document is None:
                points_without_document += 1
                continue
            matched_documents += 1

            page_index = _resolve_page_index(payload)
            if page_index is None:
                points_without_page += 1
                continue

            figure_records = _build_page_figure_records(
                settings=service_settings,
                document=document,
                page_index=page_index,
                payload=payload,
            )
            if figure_records:
                points_with_figures += 1

            update_payload = {
                "has_figure_records": bool(figure_records),
                "figure_record_count": len(figure_records),
                "figure_records": figure_records,
            }
            if (
                payload.get("has_figure_records") == update_payload["has_figure_records"]
                and payload.get("figure_record_count") == update_payload["figure_record_count"]
                and payload.get("figure_records") == update_payload["figure_records"]
            ):
                continue

            updates.append(
                PointUpdate(
                    point_id=getattr(point, "id", None),
                    payload=update_payload,
                    doc_id=document.doc_id,
                    page_number=_coerce_int(payload.get("page_number")) or (_coerce_int(payload.get("page"))),
                    figure_record_count=len(figure_records),
                    title=document.title,
                )
            )

        if offset is None:
            break

    summary: dict[str, Any] = {
        "mode": "apply" if args.apply else "dry_run",
        "collection": service_settings.rag_collection_name,
        "scanned_points": scanned_points,
        "matched_documents": matched_documents,
        "points_without_document": points_without_document,
        "points_without_page": points_without_page,
        "points_with_figures": points_with_figures,
        "points_needing_update": len(updates),
        "preview": [
            {
                "point_id": update.point_id,
                "doc_id": update.doc_id,
                "page_number": update.page_number,
                "figure_record_count": update.figure_record_count,
                "title": update.title,
            }
            for update in updates[: args.preview_limit]
        ],
    }

    if not args.apply:
        print(json.dumps(summary, indent=2))
        return 0

    write_calls = 0
    for update in updates:
        client.set_payload(
            collection_name=service_settings.rag_collection_name,
            payload=update.payload,
            points=[update.point_id],
        )
        write_calls += 1
    summary["write_calls"] = write_calls
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
