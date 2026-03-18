"""Manifest and JSONL export helpers for one-collection RAG ingestion."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

from .materialize_evidence import MaterializedRecord


@dataclass(frozen=True)
class EmbeddingModelSpec:
    """Exact model identity for one embedding family."""

    name: str
    dim: int | None = None
    enabled: bool = True


@dataclass
class IngestionManifest:
    """Rebuild ledger for a single offline materialization run."""

    run_id: str
    pipeline_version: str
    collection_name: str
    source_corpus_version: str
    embedding_models: dict[str, EmbeddingModelSpec]
    counts: dict[str, int]
    artifacts: dict[str, str]
    collection_schema: dict[str, Any]
    status: dict[str, bool]
    errors: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["embedding_models"] = {
            key: asdict(value) for key, value in self.embedding_models.items()
        }
        return payload

    def write_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return output_path


def write_jsonl(records: Iterable[MaterializedRecord], path: str | Path) -> Path:
    """Write materialized records as JSONL."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record.as_json(), ensure_ascii=True) + "\n")
    return output_path


def build_manifest(
    *,
    run_id: str,
    pipeline_version: str,
    collection_name: str,
    source_corpus_version: str,
    embedding_models: dict[str, EmbeddingModelSpec],
    artifacts: dict[str, str],
    page_records: list[MaterializedRecord],
    figure_records: list[MaterializedRecord],
    molecule_records: list[MaterializedRecord],
    points_total: int,
    errors: list[str] | None = None,
    notes: list[str] | None = None,
    upsert_completed: bool = False,
    payload_indexes_created: bool = False,
    validation_completed: bool = False,
) -> IngestionManifest:
    """Assemble a manifest from materialized record sets."""

    molecules_total = len(molecule_records)
    molecules_parsed = sum(
        1 for record in molecule_records if record.payload.get("review_status") == "parsed"
    )
    molecules_invalid = molecules_total - molecules_parsed

    return IngestionManifest(
        run_id=run_id,
        pipeline_version=pipeline_version,
        collection_name=collection_name,
        source_corpus_version=source_corpus_version,
        embedding_models=embedding_models,
        counts={
            "documents": len({record.payload["document_id"] for record in page_records + figure_records + molecule_records}),
            "pages": len(page_records),
            "figures": len(figure_records),
            "molecules_total": molecules_total,
            "molecules_parsed": molecules_parsed,
            "molecules_invalid": molecules_invalid,
            "points_total": points_total,
        },
        artifacts=artifacts,
        collection_schema={
            "name": collection_name,
            "point_types": ["page", "figure", "molecule"],
            "vectors": ["text_dense", "text_sparse", "vision_li", "chem_dense"],
        },
        status={
            "upsert_completed": upsert_completed,
            "payload_indexes_created": payload_indexes_created,
            "validation_completed": validation_completed,
        },
        errors=errors or [],
        notes=notes or [],
    )
