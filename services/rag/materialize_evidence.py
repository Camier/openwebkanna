"""Materialize atomic evidence records for the one-collection RAG.

This module converts upstream extraction outputs into deterministic,
typed evidence objects before any embedding or Qdrant upsert step.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Iterable
from uuid import NAMESPACE_URL, uuid5

from .normalize_smiles import NormalizedSmiles, normalize_smiles


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


@dataclass(frozen=True)
class SourceDocument:
    """Stable document identity and metadata."""

    document_id: str
    document_title: str
    source_uri: str
    doi: str | None = None
    authors: list[str] = field(default_factory=list)
    year: int | None = None


@dataclass(frozen=True)
class PageInput:
    """Normalized source input for page materialization."""

    page_number: int
    chunk_text: str
    section: str | None = None
    linked_figure_ids: list[str] = field(default_factory=list)
    linked_molecule_ids: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class FigureInput:
    """Normalized source input for figure materialization."""

    page_number: int
    figure_seq: int
    image_uri: str
    caption_text: str | None = None
    ocr_text: str | None = None
    bbox: list[int] | None = None
    figure_kind: str = "other"
    linked_molecule_ids: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class MoleculeInput:
    """Normalized source input for molecule materialization."""

    page_number: int
    figure_seq: int
    molecule_seq: int
    raw_smiles: str
    source_text: str | None = None
    backend: str | list[str] | None = None
    confidence: float | list[float] | None = None
    image_uri: str | None = None
    bbox: list[int] | None = None
    formula: str | None = None
    figure_kind: str | None = None


@dataclass(frozen=True)
class MaterializedRecord:
    """Generic materialized record ready for embedding."""

    point_id: str
    payload: dict[str, Any]
    text_for_embedding: str | None = None
    image_uri: str | None = None
    smiles_for_embedding: str | None = None

    def as_json(self) -> dict[str, Any]:
        return asdict(self)


def make_page_id(document_id: str, page_number: int) -> str:
    return f"{document_id}_p{page_number:04d}"


def make_figure_id(document_id: str, page_number: int, figure_seq: int) -> str:
    return f"fig_{document_id}_{page_number:04d}_{figure_seq:02d}"


def make_molecule_id(document_id: str, page_number: int, figure_seq: int, molecule_seq: int) -> str:
    return f"mol_{document_id}_{page_number:04d}_{figure_seq:02d}_{molecule_seq:02d}"


def make_point_id(object_type: str, local_id: str) -> str:
    return f"{object_type}::{local_id}"


def make_qdrant_point_id(point_id: str) -> str:
    """Return a deterministic Qdrant-valid UUID for a logical point identifier."""

    return str(uuid5(NAMESPACE_URL, f"rag_evidence:{point_id}"))


def materialize_pages(
    document: SourceDocument,
    pages: Iterable[PageInput],
    *,
    pipeline_version: str,
    created_at: str | None = None,
) -> list[MaterializedRecord]:
    """Materialize page evidence objects."""

    created_at = created_at or utc_now_iso()
    records: list[MaterializedRecord] = []
    for page in pages:
        page_id = make_page_id(document.document_id, page.page_number)
        point_id = make_point_id("page", page_id)
        payload = {
            "point_id": point_id,
            "object_type": "page",
            "document_id": document.document_id,
            "document_title": document.document_title,
            "page_id": page_id,
            "page_number": page.page_number,
            "parent_id": None,
            "figure_id": None,
            "block_id": None,
            "molecule_id": None,
            "source_uri": document.source_uri,
            "pipeline_version": pipeline_version,
            "created_at": created_at,
            "chunk_text": page.chunk_text,
            "section": page.section,
            "doi": document.doi,
            "authors": document.authors,
            "year": document.year,
            "has_figures": bool(page.linked_figure_ids),
            "linked_figure_ids": page.linked_figure_ids,
            "linked_molecule_ids": page.linked_molecule_ids,
        }
        text_for_embedding = _build_page_text(document, page)
        records.append(
            MaterializedRecord(
                point_id=point_id,
                payload=payload,
                text_for_embedding=text_for_embedding,
            )
        )
    return records


def materialize_figures(
    document: SourceDocument,
    figures: Iterable[FigureInput],
    *,
    pipeline_version: str,
    created_at: str | None = None,
) -> list[MaterializedRecord]:
    """Materialize figure evidence objects."""

    created_at = created_at or utc_now_iso()
    records: list[MaterializedRecord] = []
    for figure in figures:
        page_id = make_page_id(document.document_id, figure.page_number)
        figure_id = make_figure_id(document.document_id, figure.page_number, figure.figure_seq)
        point_id = make_point_id("figure", figure_id)
        payload = {
            "point_id": point_id,
            "object_type": "figure",
            "document_id": document.document_id,
            "document_title": document.document_title,
            "page_id": page_id,
            "page_number": figure.page_number,
            "parent_id": page_id,
            "figure_id": figure_id,
            "block_id": None,
            "molecule_id": None,
            "source_uri": document.source_uri,
            "pipeline_version": pipeline_version,
            "created_at": created_at,
            "caption_text": figure.caption_text,
            "ocr_text": figure.ocr_text,
            "image_uri": figure.image_uri,
            "bbox": figure.bbox,
            "figure_kind": figure.figure_kind,
            "has_smiles": bool(figure.linked_molecule_ids),
            "linked_molecule_ids": figure.linked_molecule_ids,
        }
        text_for_embedding = _build_figure_text(figure)
        records.append(
            MaterializedRecord(
                point_id=point_id,
                payload=payload,
                text_for_embedding=text_for_embedding,
                image_uri=figure.image_uri,
            )
        )
    return records


def materialize_molecules(
    document: SourceDocument,
    molecules: Iterable[MoleculeInput],
    *,
    pipeline_version: str,
    created_at: str | None = None,
    dedupe_within_figure: bool = True,
) -> list[MaterializedRecord]:
    """Materialize molecule evidence objects with RDKit normalization."""

    created_at = created_at or utc_now_iso()
    records: list[MaterializedRecord] = []
    seen_keys: set[tuple[str, str]] = set()

    for molecule in molecules:
        normalized = normalize_smiles(molecule.raw_smiles)
        page_id = make_page_id(document.document_id, molecule.page_number)
        figure_id = make_figure_id(document.document_id, molecule.page_number, molecule.figure_seq)
        molecule_id = make_molecule_id(
            document.document_id,
            molecule.page_number,
            molecule.figure_seq,
            molecule.molecule_seq,
        )
        point_id = make_point_id("molecule", molecule_id)

        if dedupe_within_figure and normalized.inchikey:
            dedupe_key = (figure_id, normalized.inchikey)
            if dedupe_key in seen_keys:
                continue
            seen_keys.add(dedupe_key)

        payload = {
            "point_id": point_id,
            "object_type": "molecule",
            "document_id": document.document_id,
            "document_title": document.document_title,
            "page_id": page_id,
            "page_number": molecule.page_number,
            "parent_id": figure_id,
            "figure_id": figure_id,
            "block_id": None,
            "molecule_id": molecule_id,
            "source_uri": document.source_uri,
            "pipeline_version": pipeline_version,
            "created_at": created_at,
            "raw_smiles": molecule.raw_smiles,
            "canonical_smiles": normalized.canonical_smiles,
            "inchikey": normalized.inchikey,
            "formula": molecule.formula,
            "backend": molecule.backend,
            "confidence": molecule.confidence,
            "review_status": normalized.review_status,
            "source_text": molecule.source_text,
            "image_uri": molecule.image_uri,
            "bbox": molecule.bbox,
            "has_smiles": normalized.valid,
            "figure_kind": molecule.figure_kind,
            "normalization_error": normalized.error,
        }
        text_for_embedding = _build_molecule_text(document, molecule, normalized)
        records.append(
            MaterializedRecord(
                point_id=point_id,
                payload=payload,
                text_for_embedding=text_for_embedding,
                image_uri=molecule.image_uri,
                smiles_for_embedding=normalized.canonical_smiles,
            )
        )
    return records


def _build_page_text(document: SourceDocument, page: PageInput) -> str:
    lines = [
        f"Title: {document.document_title}",
        f"Page: {page.page_number}",
    ]
    if page.section:
        lines.append(f"Section: {page.section}")
    lines.append(page.chunk_text.strip())
    return "\n".join(line for line in lines if line)


def _build_figure_text(figure: FigureInput) -> str:
    lines = []
    if figure.caption_text:
        lines.append(f"Caption: {figure.caption_text.strip()}")
    if figure.ocr_text:
        lines.append(f"OCR: {figure.ocr_text.strip()}")
    if figure.figure_kind:
        lines.append(f"Type: {figure.figure_kind}")
    return "\n".join(lines)


def _build_molecule_text(
    document: SourceDocument,
    molecule: MoleculeInput,
    normalized: NormalizedSmiles,
) -> str:
    lines = [f"Document: {document.document_title}"]
    if molecule.source_text:
        lines.append(f"Label: {molecule.source_text.strip()}")
    if normalized.canonical_smiles:
        lines.append(f"Canonical SMILES: {normalized.canonical_smiles}")
    if normalized.inchikey:
        lines.append(f"InChIKey: {normalized.inchikey}")
    if molecule.formula:
        lines.append(f"Formula: {molecule.formula}")
    if molecule.figure_kind:
        lines.append(f"Figure kind: {molecule.figure_kind}")
    return "\n".join(lines)
