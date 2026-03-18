from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import json
from pathlib import Path
import re

from services.multimodal_retrieval_api.settings import ServiceSettings


@dataclass(frozen=True)
class LocalSmilesEvidence:
    smiles: str
    canonical_smiles: str | None = None
    backend: str | None = None
    confidence: float | None = None
    status: str | None = None


@dataclass(frozen=True)
class LocalFigureEvidence:
    evidence_id: str
    evidence_type: str
    doc_id: str
    title: str | None
    page_index: int | None
    figure_id: str | None
    block_id: str | None
    block_type: str | None
    caption_text: str | None
    text: str | None
    html: str | None
    image_asset_id: str | None
    image_path: str | None
    bbox: list[float] = field(default_factory=list)
    smiles_records: list[LocalSmilesEvidence] = field(default_factory=list)


@dataclass
class LocalDocumentEvidence:
    doc_id: str
    title: str | None
    normalized_title: str | None
    doi: str | None
    publication_year: int | None
    extraction_dir: Path
    figures_by_page: dict[int, list[LocalFigureEvidence]] = field(default_factory=dict)

    def get_page_figures(self, page_index: int | None) -> list[LocalFigureEvidence]:
        if page_index is None:
            return []
        return list(self.figures_by_page.get(page_index, []))

    def all_figures(self) -> list[LocalFigureEvidence]:
        figures: list[LocalFigureEvidence] = []
        for page_index in sorted(self.figures_by_page):
            figures.extend(self.figures_by_page[page_index])
        return figures


@dataclass
class LocalEvidenceStore:
    extractions_root: Path
    documents_by_doc_id: dict[str, LocalDocumentEvidence]
    documents_by_normalized_title: dict[str, list[LocalDocumentEvidence]]
    documents_by_doi: dict[str, list[LocalDocumentEvidence]]

    def get_document(self, doc_id: str | None) -> LocalDocumentEvidence | None:
        if not doc_id:
            return None
        return self.documents_by_doc_id.get(doc_id)

    def get_documents_by_normalized_title(self, title: str | None) -> list[LocalDocumentEvidence]:
        normalized = _normalize_title(title)
        if normalized is None:
            return []
        return list(self.documents_by_normalized_title.get(normalized, []))

    def get_documents_by_doi(self, doi: str | None) -> list[LocalDocumentEvidence]:
        normalized = _normalize_doi(doi)
        if normalized is None:
            return []
        return list(self.documents_by_doi.get(normalized, []))


def _normalize_title(value: str | None) -> str | None:
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


def _coerce_int(value: object) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _coerce_float(value: object) -> float | None:
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


def _load_smiles_by_block_id(smiles_path: Path) -> dict[str, list[LocalSmilesEvidence]]:
    if not smiles_path.exists():
        return {}
    payload = json.loads(smiles_path.read_text(encoding="utf-8"))
    by_block_id: dict[str, list[LocalSmilesEvidence]] = {}
    seen: set[tuple[str, str]] = set()
    for result in payload.get("results", []):
        if not isinstance(result, dict):
            continue
        block_id = result.get("block_id")
        molecule_smiles = result.get("molecule_smiles")
        if not isinstance(block_id, str) or not isinstance(molecule_smiles, str) or not molecule_smiles.strip():
            continue
        raw_candidates = result.get("raw_candidates") if isinstance(result.get("raw_candidates"), list) else []
        canonical_smiles = (
            result.get("canonical_smiles")
            if isinstance(result.get("canonical_smiles"), str) and result.get("canonical_smiles").strip()
            else None
        )
        if canonical_smiles is None:
            for candidate in raw_candidates:
                if not isinstance(candidate, dict):
                    continue
                candidate_smiles = candidate.get("canonical_smiles")
                if isinstance(candidate_smiles, str) and candidate_smiles.strip():
                    canonical_smiles = candidate_smiles.strip()
                    break
        dedupe_key = (block_id, molecule_smiles)
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        by_block_id.setdefault(block_id, []).append(
            LocalSmilesEvidence(
                smiles=molecule_smiles,
                backend=result.get("backend") if isinstance(result.get("backend"), str) else None,
                confidence=_coerce_float(result.get("confidence")),
                canonical_smiles=canonical_smiles,
                status=result.get("status") if isinstance(result.get("status"), str) else None,
            )
        )
    return by_block_id


def _classify_evidence_type(block_id: str | None) -> str:
    if isinstance(block_id, str) and "/ChemicalBlock/" in block_id:
        return "chemical_block"
    return "figure"


def _load_document_evidence(
    *,
    normalized_path: Path,
    repo_root: Path,
) -> LocalDocumentEvidence | None:
    payload = json.loads(normalized_path.read_text(encoding="utf-8"))
    doc_id = payload.get("doc_id")
    if not isinstance(doc_id, str) or not doc_id:
        return None

    metadata = payload.get("metadata") if isinstance(payload.get("metadata"), dict) else {}
    title = metadata.get("title") if isinstance(metadata.get("title"), str) else None
    doi = _normalize_doi(metadata.get("doi")) if isinstance(metadata.get("doi"), str) else None
    publication_year = _coerce_int(metadata.get("publication_year"))
    extraction_dir = normalized_path.parent

    block_index_by_id: dict[str, dict] = {}
    for block in payload.get("block_index", []):
        if not isinstance(block, dict):
            continue
        block_id = block.get("block_id")
        if isinstance(block_id, str):
            block_index_by_id[block_id] = block

    images_by_asset_id: dict[str, dict] = {}
    for image in payload.get("images", []):
        if not isinstance(image, dict):
            continue
        asset_id = image.get("asset_id")
        if isinstance(asset_id, str):
            images_by_asset_id[asset_id] = image

    smiles_by_block_id = _load_smiles_by_block_id(extraction_dir / "smiles_extracted.json")
    figures_by_page: dict[int, list[LocalFigureEvidence]] = {}

    for figure in payload.get("figures", []):
        if not isinstance(figure, dict):
            continue
        block_id = figure.get("block_id") if isinstance(figure.get("block_id"), str) else None
        block = block_index_by_id.get(block_id or "", {})
        page_index = _coerce_int(figure.get("page_index"))
        if page_index is None:
            page_index = _coerce_int(block.get("page_index"))
        if page_index is None:
            continue
        image_asset_id = figure.get("image_asset_id") if isinstance(figure.get("image_asset_id"), str) else None
        image_record = images_by_asset_id.get(image_asset_id or "", {})
        bbox = figure.get("bbox") if isinstance(figure.get("bbox"), list) else block.get("bbox")
        normalized_bbox = [float(value) for value in bbox] if isinstance(bbox, list) else []
        smiles_records = list(smiles_by_block_id.get(block_id or "", []))
        entry = LocalFigureEvidence(
            evidence_id=str(figure.get("figure_id") or block_id or f"{doc_id}:page:{page_index}:figure"),
            evidence_type=_classify_evidence_type(block_id),
            doc_id=doc_id,
            title=title,
            page_index=page_index,
            figure_id=figure.get("figure_id") if isinstance(figure.get("figure_id"), str) else None,
            block_id=block_id,
            block_type=block.get("block_type") if isinstance(block.get("block_type"), str) else None,
            caption_text=figure.get("caption_text") if isinstance(figure.get("caption_text"), str) else None,
            text=block.get("text") if isinstance(block.get("text"), str) else None,
            html=block.get("html") if isinstance(block.get("html"), str) else None,
            image_asset_id=image_asset_id,
            image_path=_resolve_image_path(
                repo_root,
                image_record.get("storage_uri") if isinstance(image_record.get("storage_uri"), str) else None,
            ),
            bbox=normalized_bbox,
            smiles_records=smiles_records,
        )
        figures_by_page.setdefault(page_index, []).append(entry)

    return LocalDocumentEvidence(
        doc_id=doc_id,
        title=title,
        normalized_title=_normalize_title(title),
        doi=doi,
        publication_year=publication_year,
        extraction_dir=extraction_dir,
        figures_by_page=figures_by_page,
    )


@lru_cache(maxsize=1)
def load_local_evidence_store() -> LocalEvidenceStore:
    settings = ServiceSettings.load()
    documents: dict[str, LocalDocumentEvidence] = {}
    documents_by_normalized_title: dict[str, list[LocalDocumentEvidence]] = {}
    documents_by_doi: dict[str, list[LocalDocumentEvidence]] = {}
    if settings.extractions_root.exists():
        for normalized_path in sorted(settings.extractions_root.glob("*/normalized.json")):
            try:
                document = _load_document_evidence(
                    normalized_path=normalized_path,
                    repo_root=settings.repo_root,
                )
            except Exception:
                continue
            if document is None:
                continue
            documents[document.doc_id] = document
            if document.normalized_title:
                documents_by_normalized_title.setdefault(document.normalized_title, []).append(document)
            if document.doi:
                documents_by_doi.setdefault(document.doi, []).append(document)
    return LocalEvidenceStore(
        extractions_root=settings.extractions_root,
        documents_by_doc_id=documents,
        documents_by_normalized_title=documents_by_normalized_title,
        documents_by_doi=documents_by_doi,
    )
