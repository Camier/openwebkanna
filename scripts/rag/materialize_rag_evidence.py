#!/usr/bin/env python3
"""Offline materialization entrypoint for the one-collection RAG.

This script reads a small structured JSON spec, materializes atomic
`page` and `figure` evidence records, writes them as JSONL,
and emits an ingestion manifest. It does not embed or upsert anything.

Input spec shape:
{
  "run_id": "2026-03-14T18-00-00Z_rag_v1",
  "pipeline_version": "rag-v1.0.0",
  "collection_name": "rag_evidence",
  "source_corpus_version": "corpus-2026-03-14",
  "embedding_models": {
    "text_dense": {"name": "your-text-model", "dim": 1024, "enabled": true},
    "vision_li": {"name": "vidore/colpali-...", "dim": 128, "enabled": true},
    "text_sparse": {"name": "bm25", "enabled": true}
  },
  "document": {
    "document_id": "doc_example_001",
    "document_title": "Example document",
    "source_uri": "file:///tmp/example.pdf",
    "doi": "10.0000/example",
    "authors": ["Author A"],
    "year": 2026
  },
  "pages": [...],
  "figures": [...],
}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from services.rag import (
    EmbeddingModelSpec,
    FigureInput,
    PageInput,
    SourceDocument,
    build_manifest,
    materialize_figures,
    materialize_pages,
    write_jsonl,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-json",
        required=True,
        help="Path to a structured JSON input spec.",
    )
    parser.add_argument(
        "--output-dir",
        default="artifacts/rag",
        help="Base output directory for manifests and materialized JSONL.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_json)
    output_dir = Path(args.output_dir)

    spec = json.loads(input_path.read_text(encoding="utf-8"))
    run_id = spec["run_id"]
    pipeline_version = spec["pipeline_version"]
    collection_name = spec.get("collection_name", "rag_evidence")
    source_corpus_version = spec.get("source_corpus_version", "unknown-corpus")

    document = SourceDocument(**spec["document"])
    pages = [PageInput(**item) for item in spec.get("pages", [])]
    figures = [FigureInput(**item) for item in spec.get("figures", [])]

    page_records = materialize_pages(
        document,
        pages,
        pipeline_version=pipeline_version,
    )
    figure_records = materialize_figures(
        document,
        figures,
        pipeline_version=pipeline_version,
    )

    materialized_dir = output_dir / "materialized"
    manifests_dir = output_dir / "manifests"

    page_path = write_jsonl(page_records, materialized_dir / f"{run_id}.pages.jsonl")
    figure_path = write_jsonl(figure_records, materialized_dir / f"{run_id}.figures.jsonl")

    embedding_models = {
        key: EmbeddingModelSpec(**value)
        for key, value in spec.get("embedding_models", {}).items()
    }
    manifest = build_manifest(
        run_id=run_id,
        pipeline_version=pipeline_version,
        collection_name=collection_name,
        source_corpus_version=source_corpus_version,
        embedding_models=embedding_models,
        artifacts={
            "source_spec": str(input_path),
            "page_records": str(page_path),
            "figure_records": str(figure_path),
        },
        page_records=page_records,
        figure_records=figure_records,
        points_total=len(page_records) + len(figure_records),
    )
    manifest_path = manifest.write_json(manifests_dir / f"{run_id}.json")

    print(
        json.dumps(
            {
                "run_id": run_id,
                "page_records": len(page_records),
                "figure_records": len(figure_records),
                "manifest_path": str(manifest_path),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
