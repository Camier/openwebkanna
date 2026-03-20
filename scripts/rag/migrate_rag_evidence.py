#!/usr/bin/env python3
"""Run the one-collection `rag_evidence` migration as one coherent operation.

This script intentionally orchestrates the existing canonical migration stages
instead of reimplementing them:

1. copy legacy page dense+sparse vectors into `rag_evidence`
2. index figure multivectors with official ColModernVBERT

It also emits a compact post-migration collection summary so operators can see
whether the target collection is structurally populated without switching
between multiple commands.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from html import unescape
from pathlib import Path
from typing import Any

from qdrant_client import QdrantClient
from qdrant_client.http import models


def make_page_id(document_id: str, page_number: int) -> str:
    return f"{document_id}#page:{page_number}"


def make_point_id(prefix: str, value: str) -> str:
    return f"{prefix}#{value}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Materialization manifest JSON used to resolve figure JSONL.",
    )
    parser.add_argument(
        "--figures-jsonl",
        type=Path,
        help="Materialized figure JSONL. Overrides --manifest for stage 2.",
    )
    parser.add_argument(
        "--extractions-root",
        type=Path,
        default=None,
        help="Extraction root containing */normalized.json. Defaults to datalab_extraction when present.",
    )
    parser.add_argument(
        "--qdrant-url", default=os.environ.get("QDRANT_URL", "http://127.0.0.1:6335")
    )
    parser.add_argument("--qdrant-api-key", default=os.environ.get("QDRANT_API_KEY"))
    parser.add_argument(
        "--source-collection",
        default=os.environ.get("QDRANT_COLLECTION", "pdf_nemotron_hybrid"),
    )
    parser.add_argument("--target-collection", default="rag_evidence")
    parser.add_argument(
        "--source-dense-vector",
        default=os.environ.get("QDRANT_VECTOR_NAME") or "dense_prefetch",
    )
    parser.add_argument(
        "--source-sparse-vector",
        default=os.environ.get("QDRANT_SPARSE_VECTOR_NAME", "sparse_text"),
    )
    parser.add_argument("--target-dense-vector", default="text_dense")
    parser.add_argument("--target-sparse-vector", default="text_sparse")
    parser.add_argument("--vision-vector-name", default="vision_li")
    parser.add_argument("--vision-li-dim", type=int, default=128)
    parser.add_argument("--text-batch-size", type=int, default=128)
    parser.add_argument("--text-limit", type=int, default=0)
    parser.add_argument("--figure-batch-size", type=int, default=4)
    parser.add_argument("--figure-upsert-batch-size", type=int, default=8)
    parser.add_argument("--model-name", default="ModernVBERT/colmodernvbert")
    parser.add_argument("--device", default=None)
    parser.add_argument("--torch-dtype", default="float32")
    parser.add_argument("--attn-implementation", default=None)
    parser.add_argument("--output-report", type=Path, default=None)
    parser.add_argument("--skip-text", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    extractions_root = resolve_extractions_root(repo_root, args.extractions_root)
    client = QdrantClient(url=args.qdrant_url, api_key=args.qdrant_api_key)
    text_dense_dim = resolve_dense_dimension(
        client, args.source_collection, args.source_dense_vector
    )
    generated_figures_jsonl: Path | None = None

    stage_reports: dict[str, Any] = {}
    if not args.skip_text:
        stage_reports["text_migration"] = run_stage(
            repo_root / "scripts/rag/migrate_legacy_text_points.py",
            [
                "--qdrant-url",
                args.qdrant_url,
                "--source-collection",
                args.source_collection,
                "--target-collection",
                args.target_collection,
                "--source-dense-vector",
                args.source_dense_vector,
                "--source-sparse-vector",
                args.source_sparse_vector,
                "--target-dense-vector",
                args.target_dense_vector,
                "--target-sparse-vector",
                args.target_sparse_vector,
                "--vision-li-dim",
                str(args.vision_li_dim),
                "--batch-size",
                str(args.text_batch_size),
                "--limit",
                str(args.text_limit),
            ],
            qdrant_api_key=args.qdrant_api_key,
        )

    if not args.skip_figures:
        figure_stage_input = args.figures_jsonl
        if figure_stage_input is None and args.manifest is None:
            generated_figures_jsonl, materialization_report = (
                materialize_figures_from_extractions(
                    repo_root=repo_root,
                    extractions_root=extractions_root,
                )
            )
            figure_stage_input = generated_figures_jsonl
            stage_reports["figure_materialization"] = materialization_report
        if figure_stage_input is None and args.manifest is None:
            raise SystemExit(
                "Figure indexing requires --manifest, --figures-jsonl, or a usable extraction root."
            )
        figure_args = [
            "--qdrant-url",
            args.qdrant_url,
            "--collection-name",
            args.target_collection,
            "--vision-vector-name",
            args.vision_vector_name,
            "--text-dense-dim",
            str(text_dense_dim),
            "--model-name",
            args.model_name,
            "--batch-size",
            str(args.figure_batch_size),
            "--upsert-batch-size",
            str(args.figure_upsert_batch_size),
            "--torch-dtype",
            args.torch_dtype,
        ]
        if args.device:
            figure_args.extend(["--device", args.device])
        if args.attn_implementation:
            figure_args.extend(["--attn-implementation", args.attn_implementation])
        if figure_stage_input is not None:
            figure_args.extend(["--figures-jsonl", str(figure_stage_input)])
        else:
            figure_args.extend(["--manifest", str(args.manifest)])
        stage_reports["figure_indexing"] = run_stage(
            repo_root / "scripts/rag/index_colmodernvbert_figures.py",
            figure_args,
            qdrant_api_key=args.qdrant_api_key,
        )

    summary = {
        "target_collection": args.target_collection,
        "extractions_root": str(extractions_root),
        "generated_figures_jsonl": str(generated_figures_jsonl)
        if generated_figures_jsonl is not None
        else None,
        "stage_reports": stage_reports,
        "collection_summary": build_collection_summary(client, args.target_collection),
    }
    if args.output_report is not None:
        args.output_report.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


def run_stage(
    script_path: Path,
    script_args: list[str],
    *,
    qdrant_api_key: str | None,
    qdrant_url: str | None = None,
    env_overrides: dict[str, str] | None = None,
) -> Any:
    env = os.environ.copy()
    repo_root = script_path.parents[2]
    existing_pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        f"{repo_root}{os.pathsep}{existing_pythonpath}" if existing_pythonpath else str(repo_root)
    )
    if qdrant_api_key:
        env["QDRANT_API_KEY"] = qdrant_api_key
    if qdrant_url:
        env["QDRANT_URL"] = qdrant_url
    if env_overrides:
        env.update(env_overrides)
    completed = subprocess.run(
        [sys.executable, str(script_path), *script_args],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    stdout = completed.stdout.strip()
    if not stdout:
        return {"ok": True, "stdout": ""}
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        return {"ok": True, "stdout": stdout}


def resolve_dense_dimension(
    client: QdrantClient, collection_name: str, dense_vector_name: str
) -> int:
    collection = client.get_collection(collection_name)
    vectors = getattr(collection.config.params, "vectors", None) or {}
    if isinstance(vectors, dict) and dense_vector_name in vectors:
        return int(vectors[dense_vector_name].size)
    raise SystemExit(
        f"Could not resolve dense vector dimension for {dense_vector_name!r} in {collection_name!r}."
    )


def resolve_extractions_root(repo_root: Path, explicit_root: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit_root is not None:
        candidates.append(explicit_root)
    candidates.extend(
        [
            repo_root / "datalab_extraction",
            repo_root / "data" / "extractions",
        ]
    )
    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved.exists():
            return resolved
    raise SystemExit(
        "Could not resolve extraction root. Pass --extractions-root or create datalab_extraction/."
    )


def materialize_figures_from_extractions(
    *, repo_root: Path, extractions_root: Path
) -> tuple[Path, dict[str, Any]]:
    output_dir = repo_root / "artifacts" / "rag" / "materialized"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "generated.figures.jsonl"

    documents_scanned = 0
    figures_written = 0
    figures_skipped_missing_image = 0
    with output_path.open("w", encoding="utf-8") as handle:
        for normalized_path in sorted(extractions_root.glob("*/normalized.json")):
            documents_scanned += 1
            extraction_dir = normalized_path.parent
            normalized = json.loads(normalized_path.read_text(encoding="utf-8"))
            figures = (
                normalized.get("figures")
                if isinstance(normalized.get("figures"), list)
                else []
            )
            block_index = (
                normalized.get("block_index")
                if isinstance(normalized.get("block_index"), list)
                else []
            )
            blocks_by_id = {
                block.get("block_id"): block
                for block in block_index
                if isinstance(block, dict) and isinstance(block.get("block_id"), str)
            }
            images = (
                normalized.get("images")
                if isinstance(normalized.get("images"), list)
                else []
            )
            images_by_asset_id = {
                image.get("asset_id"): image
                for image in images
                if isinstance(image, dict) and isinstance(image.get("asset_id"), str)
            }
            metadata = (
                normalized.get("metadata")
                if isinstance(normalized.get("metadata"), dict)
                else {}
            )
            document_id = str(normalized.get("doc_id") or extraction_dir.name)
            document_title = str(metadata.get("title") or document_id)
            source_uri = resolve_source_uri(repo_root, extraction_dir)

            for figure_index, figure in enumerate(figures, start=1):
                if not isinstance(figure, dict):
                    continue
                block_id = (
                    figure.get("block_id")
                    if isinstance(figure.get("block_id"), str)
                    else None
                )
                block = blocks_by_id.get(block_id or "", {})
                page_index = _coerce_int(figure.get("page_index"))
                if page_index is None:
                    page_index = _coerce_int(block.get("page_index"))
                if page_index is None:
                    continue
                page_number = page_index + 1
                page_id = make_page_id(document_id, page_number)
                figure_id = (
                    figure.get("figure_id")
                    if isinstance(figure.get("figure_id"), str)
                    and figure.get("figure_id")
                    else f"fig_{document_id}_{page_number:04d}_{figure_index:02d}"
                )
                image_asset_id = (
                    figure.get("image_asset_id")
                    if isinstance(figure.get("image_asset_id"), str)
                    else None
                )
                image_record = images_by_asset_id.get(image_asset_id or "", {})
                storage_uri = (
                    image_record.get("storage_uri")
                    if isinstance(image_record.get("storage_uri"), str)
                    else None
                )
                image_uri = resolve_image_uri(repo_root, extraction_dir, storage_uri)
                if image_uri is None:
                    figures_skipped_missing_image += 1
                    continue
                bbox = normalize_bbox(
                    figure.get("bbox")
                    if isinstance(figure.get("bbox"), list)
                    else block.get("bbox")
                )
                payload = {
                    "point_id": make_point_id("figure", figure_id),
                    "object_type": "figure",
                    "document_id": document_id,
                    "document_title": document_title,
                    "page_id": page_id,
                    "page_number": page_number,
                    "parent_id": page_id,
                    "figure_id": figure_id,
                    "block_id": block_id,
                    "molecule_id": None,
                    "source_uri": source_uri,
                    "pipeline_version": "rag-migration-v1",
                    "created_at": datetime.now(timezone.utc).isoformat(),
                    "caption_text": figure.get("caption_text")
                    if isinstance(figure.get("caption_text"), str)
                    else None,
                    "ocr_text": (
                        block.get("text")
                        if isinstance(block.get("text"), str) and block.get("text")
                        else html_to_text(
                            block.get("html")
                            if isinstance(block.get("html"), str)
                            else None
                        )
                    ),
                    "image_uri": image_uri,
                    "bbox": bbox,
                    "figure_kind": (
                        block.get("block_type")
                        if isinstance(block.get("block_type"), str)
                        and block.get("block_type")
                        else "other"
                    ),
                    "linked_molecule_ids": [],
                }
                record = {
                    "point_id": payload["point_id"],
                    "image_uri": image_uri,
                    "payload": payload,
                }
                handle.write(json.dumps(record, ensure_ascii=True) + "\n")
                figures_written += 1

    return output_path, {
        "extractions_root": str(extractions_root),
        "figures_jsonl": str(output_path),
        "documents_scanned": documents_scanned,
        "figures_written": figures_written,
        "figures_skipped_missing_image": figures_skipped_missing_image,
    }


def resolve_source_uri(repo_root: Path, extraction_dir: Path) -> str:
    pdf_candidate = repo_root / "pdfs" / f"{extraction_dir.name}.pdf"
    if pdf_candidate.exists():
        return str(pdf_candidate.resolve())
    return str((extraction_dir / "normalized.json").resolve())


def resolve_image_uri(
    repo_root: Path, extraction_dir: Path, storage_uri: str | None
) -> str | None:
    if not storage_uri:
        return None
    candidate = Path(storage_uri)
    candidates = [candidate]
    if not candidate.is_absolute():
        candidates.extend([extraction_dir / candidate, repo_root / candidate])
        storage_text = storage_uri.replace("\\", "/")
        asset_marker = "/_assets/"
        if asset_marker in storage_text:
            asset_suffix = storage_text.split(asset_marker, 1)[1]
            candidates.append(extraction_dir / "_assets" / asset_suffix)
    for item in candidates:
        resolved = item.expanduser().resolve()
        if resolved.exists():
            return str(resolved)
    return None


def normalize_bbox(raw_bbox: Any) -> list[float]:
    if not isinstance(raw_bbox, list):
        return []
    values = [_coerce_float(value) for value in raw_bbox]
    return [value for value in values if value is not None]


def html_to_text(value: str | None) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    text = "".join(" " if ch == "\n" else ch for ch in value)
    text = unescape(text)
    text = " ".join(text.replace("<", " <").replace(">", "> ").split())
    text = " ".join(
        part
        for part in text.split()
        if not (part.startswith("<") and part.endswith(">"))
    )
    return text or None


def _coerce_int(value: Any) -> int | None:
    try:
        return int(value) if value is not None and value != "" else None
    except (TypeError, ValueError):
        return None


def _coerce_float(value: Any) -> float | None:
    try:
        return float(value) if value is not None and value != "" else None
    except (TypeError, ValueError):
        return None


def build_collection_summary(
    client: QdrantClient, collection_name: str
) -> dict[str, Any]:
    return {
        "points_total": count_points(client, collection_name),
        "pages": count_points(client, collection_name, object_type="page"),
        "figures": count_points(client, collection_name, object_type="figure"),
        "pages_with_figure_records": count_points(
            client,
            collection_name,
            object_type="page",
            extra_conditions=[
                models.FieldCondition(
                    key="has_figure_records",
                    match=models.MatchValue(value=True),
                )
            ],
        ),
    }


def count_points(
    client: QdrantClient,
    collection_name: str,
    *,
    object_type: str | None = None,
    extra_conditions: list[models.FieldCondition] | None = None,
) -> int:
    conditions: list[models.FieldCondition] = []
    if object_type is not None:
        conditions.append(
            models.FieldCondition(
                key="object_type", match=models.MatchValue(value=object_type)
            )
        )
    if extra_conditions:
        conditions.extend(extra_conditions)
    count_filter = models.Filter(must=conditions) if conditions else None
    response = client.count(
        collection_name=collection_name, count_filter=count_filter, exact=True
    )
    return int(getattr(response, "count", 0))


if __name__ == "__main__":
    raise SystemExit(main())
