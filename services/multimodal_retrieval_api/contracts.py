from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field


RetrievalMode = Literal[
    "hybrid",
    "dense_only",
    "lexical_fallback",
    "hybrid_no_late",
    "hybrid_with_late",
]
EvidenceType = Literal["page", "figure", "chemical_block"]
EvidenceSource = Literal["qdrant", "local_extraction"]


class RetrieveRequest(BaseModel):
    query: str = Field(min_length=1, description="Free-text retrieval query.")
    top_k: int = Field(default=10, ge=1, le=100, description="Maximum number of reranked hits to return.")
    mode: RetrievalMode | None = Field(
        default=None,
        description="Override thesis_graph retrieval mode. Defaults to runtime settings.",
    )


class RetrievalEvidenceHit(BaseModel):
    point_id: str
    score: float
    rank: int | None = None
    stage: str
    doc_id: str | None = None
    page_number: int | None = None
    title: str | None = None
    citation_key: str | None = None
    group_id: str | None = None
    source_pdf_sha256: str | None = None
    retrieval_mode: str | None = None
    vector_name: str | None = None
    base_score: float | None = None
    metadata_bonus: float | None = None
    payload: dict[str, Any] = Field(default_factory=dict)


class RetrievalBackendInfo(BaseModel):
    service: str
    qdrant_url: str
    qdrant_collection: str
    dense_vector_name: str | None = None
    sparse_vector_name: str | None = None
    late_vector_name: str | None = None


class RetrievalEvidenceObject(BaseModel):
    evidence_id: str
    evidence_type: EvidenceType
    source: EvidenceSource
    doc_id: str
    title: str | None = None
    page_number: int | None = None
    page_index: int | None = None
    point_id: str | None = None
    parent_point_id: str | None = None
    rank: int | None = None
    score: float | None = None
    source_pdf_sha256: str | None = None
    citation_key: str | None = None
    figure_id: str | None = None
    block_id: str | None = None
    block_type: str | None = None
    caption_text: str | None = None
    text: str | None = None
    image_asset_id: str | None = None
    image_path: str | None = None
    bbox: list[float] = Field(default_factory=list)
    figure_kind: str | None = None


class RetrieveResponse(BaseModel):
    query: str
    mode: str
    top_k: int
    backend: RetrievalBackendInfo
    candidate_hits: list[RetrievalEvidenceHit] = Field(default_factory=list)
    reranked_hits: list[RetrievalEvidenceHit] = Field(default_factory=list)
    evidence_objects: list[RetrievalEvidenceObject] = Field(default_factory=list)
    diagnostics: dict[str, Any] = Field(default_factory=dict)


class HealthResponse(BaseModel):
    status: Literal["ok"]
    service: str


class ReadyResponse(BaseModel):
    status: Literal["ready", "degraded"]
    service: str
    backend: RetrievalBackendInfo | None = None
    qdrant: dict[str, Any] | None = None
    runtime: dict[str, Any] | None = None
    warnings: list[str] = Field(default_factory=list)
