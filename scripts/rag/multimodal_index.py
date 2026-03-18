#!/usr/bin/env python3

from __future__ import annotations

import base64
import binascii
import html
import io
import json
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image


FIGURE_TYPES = {"Figure", "Picture", "ChemicalBlock"}
CAPTION_LOOKAROUND = (-1, 1, 2, -2)


@dataclass
class CorpusDoc:
    paper_id: str
    title: str
    doi: str
    citekey: str
    doc_id: str
    extraction_path: Path
    aliases: set[str]


@dataclass
class FigureRecord:
    figure_id: str
    paper_id: str
    doc_id: str
    title: str
    doi: str
    citekey: str
    extraction_path: str
    chunk_path: str
    block_id: str
    block_type: str
    page: int | None
    image_name: str
    description: str
    caption: str
    retrieval_text: str
    molecule_smiles: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def slugify(value: str) -> str:
    slug = re.sub(r"[^0-9A-Za-z._-]+", "-", value).strip("-")
    return slug or "item"


def model_slug(model_id: str) -> str:
    return slugify(model_id.replace("/", "-"))


def normalize_text(value: str) -> str:
    lowered = value.casefold()
    lowered = re.sub(r"[^0-9a-z]+", " ", lowered)
    return " ".join(lowered.split())


def tokenize(value: str) -> set[str]:
    return {token for token in normalize_text(value).split() if len(token) >= 3}


def strip_tags(value: str) -> str:
    text = re.sub(r"<[^>]+>", " ", value)
    return " ".join(html.unescape(text).split())


def extract_img_alt(value: str) -> str:
    match = re.search(r'<img[^>]*\salt="([^"]+)"', value, flags=re.IGNORECASE)
    return html.unescape(match.group(1)).strip() if match else ""


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


@lru_cache(maxsize=None)
def load_smiles_records(smiles_path_str: str) -> tuple[dict[str, Any], ...]:
    smiles_path = Path(smiles_path_str)
    if not smiles_path.exists():
        return ()
    payload = load_json(smiles_path)
    picture_data = payload.get("picture_data") or []
    if not isinstance(picture_data, list):
        return ()
    return tuple(item for item in picture_data if isinstance(item, dict))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=True) + "\n")


def l2_normalize(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms = np.where(norms == 0.0, 1.0, norms)
    return (matrix / norms).astype(np.float32, copy=False)


def decode_image_bytes(value: str) -> bytes:
    cleaned = value.strip()
    if "," in cleaned and cleaned.lower().startswith("data:image/"):
        cleaned = cleaned.split(",", 1)[1]
    return base64.b64decode(cleaned)


def build_corpus_docs(repo_root: Path) -> list[CorpusDoc]:
    docs: list[CorpusDoc] = []
    for row in load_jsonl(repo_root / "data/corpus/biblio_corpus.curated.jsonl"):
        extraction_raw = str(row.get("extraction_path") or "").strip()
        paper_id = str(row.get("paper_id") or "").strip()
        if not extraction_raw or not paper_id:
            continue
        extraction_path = repo_root / extraction_raw
        aliases: set[str] = set()
        for key in ("paper_id", "title", "doi", "citekey", "doc_id"):
            value = str(row.get(key) or "").strip()
            normalized = normalize_text(value)
            if normalized:
                aliases.add(normalized)
        aliases.add(normalize_text(extraction_path.name))
        aliases.add(normalize_text(extraction_path.name.rsplit("_", 1)[0]))
        docs.append(
            CorpusDoc(
                paper_id=paper_id,
                title=str(row.get("title") or "").strip(),
                doi=str(row.get("doi") or "").strip(),
                citekey=str(row.get("citekey") or "").strip(),
                doc_id=str(row.get("doc_id") or "").strip(),
                extraction_path=extraction_path,
                aliases=aliases,
            )
        )
    return docs


def candidate_caption(blocks: list[dict[str, Any]], index: int) -> str:
    figure_block = blocks[index]
    figure_page = figure_block.get("page")
    captions: list[str] = []
    seen: set[str] = set()
    for offset in CAPTION_LOOKAROUND:
        probe = index + offset
        if probe < 0 or probe >= len(blocks):
            continue
        block = blocks[probe]
        if block.get("block_type") != "Caption":
            continue
        if figure_page is not None and block.get("page") != figure_page:
            continue
        caption = strip_tags(str(block.get("html") or ""))
        normalized = normalize_text(caption)
        if caption and normalized and normalized not in seen:
            captions.append(caption)
            seen.add(normalized)
    return " ".join(captions).strip()


def figure_description(block: dict[str, Any], caption: str = "") -> str:
    html_value = str(block.get("html") or "")
    parts = [extract_img_alt(html_value), strip_tags(html_value), caption]
    deduped: list[str] = []
    seen: set[str] = set()
    for part in parts:
        normalized = normalize_text(part)
        if part and normalized and normalized not in seen:
            deduped.append(part)
            seen.add(normalized)
    return " ".join(deduped).strip()


def build_retrieval_text(*parts: str) -> str:
    deduped: list[str] = []
    seen: set[str] = set()
    for part in parts:
        part = part.strip()
        normalized = normalize_text(part)
        if part and normalized and normalized not in seen:
            deduped.append(part)
            seen.add(normalized)
    return " ".join(deduped).strip()


def resolve_molecule_smiles(doc: CorpusDoc, *, block_id: str, image_name: str) -> list[str]:
    smiles: list[str] = []
    seen: set[str] = set()
    smiles_path = doc.extraction_path / "smiles_extracted.json"
    for record in load_smiles_records(str(smiles_path)):
        record_block_id = str(record.get("block_id") or "").strip()
        record_image = Path(str(record.get("image") or "")).name
        record_source_image = Path(str(record.get("source_image") or "")).name
        block_match = bool(record_block_id and record_block_id == block_id)
        image_match = bool(image_name and image_name in {record_image, record_source_image})
        if not (block_match or image_match):
            continue
        value = str(record.get("molecule_smiles") or "").strip()
        if value and value not in seen:
            smiles.append(value)
            seen.add(value)
    return smiles


def iter_figure_records(repo_root: Path) -> list[FigureRecord]:
    records: list[FigureRecord] = []
    for doc in build_corpus_docs(repo_root):
        chunk_path = doc.extraction_path / "rag/chunks.json"
        if not chunk_path.exists():
            continue
        payload = load_json(chunk_path)
        blocks = payload.get("blocks") or []
        for index, block in enumerate(blocks):
            if block.get("block_type") not in FIGURE_TYPES:
                continue
            images = block.get("images") or {}
            if not isinstance(images, dict) or not images:
                continue
            caption = candidate_caption(blocks, index)
            description = figure_description(block, caption)
            if not description:
                continue
            image_name = str(next(iter(images.keys())))
            block_id = str(block.get("id") or image_name)
            molecule_smiles = resolve_molecule_smiles(doc, block_id=block_id, image_name=image_name)
            retrieval_text = build_retrieval_text(
                doc.title,
                description,
                caption,
                " ".join(molecule_smiles),
                doc.paper_id,
                doc.doi,
                doc.citekey,
            )
            records.append(
                FigureRecord(
                    figure_id=f"{doc.doc_id}:{block_id}",
                    paper_id=doc.paper_id,
                    doc_id=doc.doc_id,
                    title=doc.title,
                    doi=doc.doi,
                    citekey=doc.citekey,
                    extraction_path=str(doc.extraction_path.relative_to(repo_root)),
                    chunk_path=str(chunk_path.relative_to(repo_root)),
                    block_id=block_id,
                    block_type=str(block.get("block_type") or ""),
                    page=block.get("page") if isinstance(block.get("page"), int) else None,
                    image_name=image_name,
                    description=description,
                    caption=caption,
                    molecule_smiles=molecule_smiles,
                    retrieval_text=retrieval_text,
                )
            )
    return records


def load_figure_image_bytes(repo_root: Path, record: FigureRecord) -> bytes:
    payload = load_json(repo_root / record.chunk_path)
    blocks = payload.get("blocks") or []
    for block in blocks:
        block_id = str(block.get("id") or "")
        if block_id != record.block_id:
            continue
        images = block.get("images") or {}
        if not isinstance(images, dict) or not images:
            break
        if record.image_name in images:
            return decode_image_bytes(str(images[record.image_name]))
        image_value = next(iter(images.values()))
        return decode_image_bytes(str(image_value))
    raise FileNotFoundError(f"Unable to resolve image payload for {record.figure_id}")


def load_figure_image(repo_root: Path, record: FigureRecord) -> Image.Image:
    return Image.open(io.BytesIO(load_figure_image_bytes(repo_root, record))).convert("RGB")


def load_clip_bundle(model_id: str, device: str = "auto") -> tuple[Any, Any, str]:
    import torch
    from transformers import AutoProcessor, CLIPModel

    resolved_device = device
    if resolved_device == "auto":
        resolved_device = "cuda" if torch.cuda.is_available() else "cpu"
    processor = AutoProcessor.from_pretrained(model_id)
    model = CLIPModel.from_pretrained(model_id)
    model.eval()
    model.to(resolved_device)
    return processor, model, resolved_device


def embed_texts(texts: list[str], *, processor: Any, model: Any, device: str, batch_size: int = 32) -> np.ndarray:
    import torch

    vectors: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(texts), batch_size):
            batch = texts[start : start + batch_size]
            inputs = processor(text=batch, return_tensors="pt", padding=True, truncation=True)
            inputs = {key: value.to(device) for key, value in inputs.items()}
            features = model.get_text_features(**inputs)
            vectors.append(features.cpu().numpy().astype(np.float32))
    return l2_normalize(np.vstack(vectors)) if vectors else np.zeros((0, 512), dtype=np.float32)


def embed_images(images: list[Image.Image], *, processor: Any, model: Any, device: str, batch_size: int = 8) -> np.ndarray:
    import torch

    vectors: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(images), batch_size):
            batch = images[start : start + batch_size]
            inputs = processor(images=batch, return_tensors="pt")
            inputs = {key: value.to(device) for key, value in inputs.items()}
            features = model.get_image_features(**inputs)
            vectors.append(features.cpu().numpy().astype(np.float32))
    return l2_normalize(np.vstack(vectors)) if vectors else np.zeros((0, 512), dtype=np.float32)


def save_index(
    index_dir: Path,
    *,
    model_id: str,
    records: list[FigureRecord],
    text_embeddings: np.ndarray,
    image_embeddings: np.ndarray,
) -> None:
    index_dir.mkdir(parents=True, exist_ok=True)
    write_json(
        index_dir / "manifest.json",
        {
            "generated_at": utc_now_iso(),
            "model_id": model_id,
            "record_count": len(records),
            "embedding_dim": int(text_embeddings.shape[1]) if text_embeddings.size else 0,
            "artifacts": {
                "records": "records.jsonl",
                "embeddings": "embeddings.npz",
            },
        },
    )
    write_jsonl(index_dir / "records.jsonl", [record.to_dict() for record in records])
    np.savez_compressed(index_dir / "embeddings.npz", text=text_embeddings, image=image_embeddings)


def load_index(index_dir: Path) -> tuple[dict[str, Any], list[FigureRecord], np.ndarray, np.ndarray]:
    manifest = load_json(index_dir / "manifest.json")
    records = [FigureRecord(**row) for row in load_jsonl(index_dir / "records.jsonl")]
    embeddings = np.load(index_dir / "embeddings.npz")
    text_embeddings = embeddings["text"].astype(np.float32)
    image_embeddings = embeddings["image"].astype(np.float32)
    return manifest, records, text_embeddings, image_embeddings


def cosine_scores(query_embedding: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    return matrix @ query_embedding.astype(np.float32)


def safe_load_figure_image(repo_root: Path, record: FigureRecord) -> Image.Image | None:
    try:
        return load_figure_image(repo_root, record)
    except (FileNotFoundError, OSError, ValueError, binascii.Error):
        return None
