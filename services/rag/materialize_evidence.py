"""Materialized evidence dataclasses and helpers for one-collection RAG."""

from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from typing import Any


@dataclass
class MaterializedRecord:
    """Single atomic evidence record for one-collection RAG."""

    point_id: str
    payload: dict[str, Any]
    text_for_embedding: str | None = None
    image_uri: str | None = None

    def as_json(self) -> dict[str, Any]:
        """Serialize to JSON-compatible dict."""
        return {
            "point_id": self.point_id,
            "payload": self.payload,
            "text_for_embedding": self.text_for_embedding,
            "image_uri": self.image_uri,
        }


def make_point_id(object_type: str, local_id: str) -> str:
    """Create a stable point ID from object type and local identifier."""
    return f"{object_type}::{local_id}"


def make_qdrant_point_id(point_id: str) -> str:
    """Convert a point_id to a stable UUID5 for Qdrant."""
    # Use UUID5 with a fixed namespace to generate deterministic UUIDs from point_id
    # This ensures the same point_id always maps to the same Qdrant point ID
    return str(uuid.uuid5(uuid.NAMESPACE_URL, point_id))


# Input dataclasses for script compatibility


@dataclass
class SourceDocument:
    """Source document metadata for materialization."""

    document_id: str
    document_title: str
    source_uri: str
    doi: str | None = None
    authors: list[str] = field(default_factory=list)
    year: int | None = None


@dataclass
class PageInput:
    """Page input for materialization."""

    page_id: str
    page_number: int
    text: str
    document_id: str | None = None


@dataclass
class FigureInput:
    """Figure input for materialization."""

    figure_id: str
    caption: str
    image_uri: str
    page_id: str | None = None
    document_id: str | None = None


# Materialization functions


def materialize_pages(
    document: SourceDocument,
    pages: list[PageInput],
    *,
    pipeline_version: str,
) -> list[MaterializedRecord]:
    """Materialize page records from source document and page inputs."""
    records: list[MaterializedRecord] = []
    for page in pages:
        point_id = make_point_id("page", page.page_id)
        payload = {
            "object_type": "page",
            "document_id": document.document_id,
            "document_title": document.document_title,
            "source_uri": document.source_uri,
            "doi": document.doi,
            "authors": document.authors,
            "year": document.year,
            "page_id": page.page_id,
            "page_number": page.page_number,
            "text": page.text,
            "pipeline_version": pipeline_version,
        }
        records.append(
            MaterializedRecord(
                point_id=point_id,
                payload=payload,
                text_for_embedding=page.text,
            )
        )
    return records


def materialize_figures(
    document: SourceDocument,
    figures: list[FigureInput],
    *,
    pipeline_version: str,
) -> list[MaterializedRecord]:
    """Materialize figure records from source document and figure inputs."""
    records: list[MaterializedRecord] = []
    for figure in figures:
        point_id = make_point_id("figure", figure.figure_id)
        # Combine caption and figure_id for embedding text
        text_for_embedding = f"Figure {figure.figure_id}: {figure.caption}"
        payload = {
            "object_type": "figure",
            "document_id": document.document_id,
            "document_title": document.document_title,
            "source_uri": document.source_uri,
            "doi": document.doi,
            "authors": document.authors,
            "year": document.year,
            "figure_id": figure.figure_id,
            "caption": figure.caption,
            "image_uri": figure.image_uri,
            "page_id": figure.page_id,
            "pipeline_version": pipeline_version,
        }
        records.append(
            MaterializedRecord(
                point_id=point_id,
                payload=payload,
                text_for_embedding=text_for_embedding,
                image_uri=figure.image_uri,
            )
        )
    return records
