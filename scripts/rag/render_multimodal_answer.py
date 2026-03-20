#!/usr/bin/env python3

from __future__ import annotations

import argparse
import base64
import binascii
import html
import json
import os
import re
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any


FIGURE_TYPES = {"Figure", "Picture", "ChemicalBlock"}
BOILERPLATE_PATTERNS = (
    "elsevier logo",
    "crossmark",
    "cover image",
    "journal cover",
    "publisher logo",
)
CHEMISTRY_TERMS = {
    "alkaloid",
    "alkaloids",
    "chemical",
    "chemistry",
    "compound",
    "compounds",
    "mesembrine",
    "mesembrane",
    "mesembrenone",
    "joubertiamine",
    "structure",
    "stereochemistry",
}


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def normalize_text(value: str) -> str:
    lowered = value.casefold()
    lowered = re.sub(r"[^0-9a-z]+", " ", lowered)
    return " ".join(lowered.split())


def slugify(value: str) -> str:
    slug = re.sub(r"[^0-9A-Za-z._-]+", "-", value).strip("-")
    return slug or "item"


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


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def strip_tags(value: str) -> str:
    text = re.sub(r"<[^>]+>", " ", value)
    return " ".join(html.unescape(text).split())


def extract_img_alt(value: str) -> str:
    match = re.search(r'<img[^>]*\salt="([^"]+)"', value, flags=re.IGNORECASE)
    return html.unescape(match.group(1)).strip() if match else ""


def tokenize(value: str) -> set[str]:
    return {token for token in normalize_text(value).split() if len(token) >= 3}


def request_json(url: str, *, token: str, method: str = "GET", payload: Any | None = None) -> Any:
    body = None
    headers: dict[str, str] = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    if payload is not None:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"HTTP {exc.code} from {url}: {detail[:400]}") from exc


def normalize_knowledge_items(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict) and isinstance(payload.get("items"), list):
        return payload["items"]
    return []


def extract_hits(payload: dict[str, Any]) -> list[dict[str, Any]]:
    ids = payload.get("ids") or []
    docs = payload.get("documents") or []
    metas = payload.get("metadatas") or []
    dists = payload.get("distances") or []

    id_row = ids[0] if ids and isinstance(ids[0], list) else []
    doc_row = docs[0] if docs and isinstance(docs[0], list) else []
    meta_row = metas[0] if metas and isinstance(metas[0], list) else []
    dist_row = dists[0] if dists and isinstance(dists[0], list) else []

    count = max(len(id_row), len(doc_row), len(meta_row), len(dist_row))
    hits: list[dict[str, Any]] = []
    for index in range(count):
        metadata = meta_row[index] if index < len(meta_row) and isinstance(meta_row[index], dict) else {}
        hits.append(
            {
                "id": id_row[index] if index < len(id_row) else None,
                "document": doc_row[index] if index < len(doc_row) else "",
                "metadata": metadata,
                "distance": dist_row[index] if index < len(dist_row) else None,
                "rank": index + 1,
            }
        )
    return hits


def first_model_id(models_payload: Any) -> str | None:
    data = []
    if isinstance(models_payload, dict):
        data = models_payload.get("data") or []
    if not isinstance(data, list):
        return None
    for item in data:
        if not isinstance(item, dict):
            continue
        model_id = str(item.get("id") or "").strip()
        if model_id and model_id != "arena-model" and "*" not in model_id and not model_id.startswith("ollama-cloud/"):
            return model_id
    return None


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
class FigureMatch:
    doc: CorpusDoc
    block_id: str
    block_type: str
    page: int | None
    description: str
    image_name: str
    image_bytes: bytes
    score: float
    source_hit_rank: int
    source_hit_page: int | None
    source_hit_title: str
    molecule_smiles: list[str] = field(default_factory=list)
    molecule_records: list[dict[str, Any]] = field(default_factory=list)


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


def _slim_molecule_record(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "block_id": record.get("block_id"),
        "image": record.get("image"),
        "source_image": record.get("source_image"),
        "crop_index": record.get("crop_index"),
        "backend": record.get("backend"),
        "detector": record.get("detector"),
        "status": record.get("status"),
        "confidence": record.get("confidence"),
        "molecule_smiles": record.get("molecule_smiles"),
        "raw_smiles": record.get("raw_smiles"),
    }


@lru_cache(maxsize=None)
def load_molecule_records_for_path(smiles_path_str: str) -> tuple[dict[str, Any], ...]:
    smiles_path = Path(smiles_path_str)
    if not smiles_path.exists():
        return ()
    payload = load_json(smiles_path)
    picture_data = payload.get("picture_data") or []
    if not isinstance(picture_data, list):
        return ()
    return tuple(item for item in picture_data if isinstance(item, dict))


def load_molecule_records(doc: CorpusDoc) -> tuple[dict[str, Any], ...]:
    return load_molecule_records_for_path(str(doc.extraction_path / "smiles_extracted.json"))


def resolve_figure_molecules(doc: CorpusDoc, *, block_id: str, image_name: str) -> tuple[list[str], list[dict[str, Any]]]:
    valid_records: list[dict[str, Any]] = []
    fallback_records: list[dict[str, Any]] = []

    for record in load_molecule_records(doc):
        record_block_id = str(record.get("block_id") or "").strip()
        record_image = Path(str(record.get("image") or "")).name
        record_source_image = Path(str(record.get("source_image") or "")).name

        block_match = bool(record_block_id and record_block_id == block_id)
        image_match = bool(image_name and image_name in {record_image, record_source_image})
        if not (block_match or image_match):
            continue

        slim = _slim_molecule_record(record)
        if slim.get("molecule_smiles"):
            valid_records.append(slim)
        else:
            fallback_records.append(slim)

    chosen_records = valid_records or fallback_records
    smiles: list[str] = []
    seen_smiles: set[str] = set()
    for record in valid_records:
        value = str(record.get("molecule_smiles") or "").strip()
        if value and value not in seen_smiles:
            smiles.append(value)
            seen_smiles.add(value)
    return smiles, chosen_records


def match_hit_to_doc(hit: dict[str, Any], corpus_docs: list[CorpusDoc]) -> tuple[CorpusDoc | None, float]:
    metadata = hit.get("metadata") or {}
    fields = [
        str(metadata.get("title") or ""),
        str(metadata.get("doi") or ""),
        str(metadata.get("name") or ""),
        str(metadata.get("source") or ""),
        str(hit.get("document") or "")[:300],
    ]
    hit_aliases = {normalize_text(field) for field in fields if normalize_text(field)}
    if not hit_aliases:
        return None, 0.0

    best_doc: CorpusDoc | None = None
    best_score = 0.0
    title_tokens = tokenize(str(metadata.get("title") or ""))
    file_tokens = tokenize(str(metadata.get("name") or "") + " " + str(metadata.get("source") or ""))

    for doc in corpus_docs:
        score = 0.0
        if hit_aliases & doc.aliases:
            score += 100.0
        if doc.doi and normalize_text(doc.doi) in hit_aliases:
            score += 60.0
        doc_title_tokens = tokenize(doc.title)
        score += len(title_tokens & doc_title_tokens) * 3.0
        score += len(file_tokens & doc.aliases) * 1.0
        if score > best_score:
            best_doc = doc
            best_score = score

    if best_score < 8.0:
        return None, best_score
    return best_doc, best_score


def decode_image_bytes(value: str) -> bytes:
    cleaned = value.strip()
    if "," in cleaned and cleaned.lower().startswith("data:image/"):
        cleaned = cleaned.split(",", 1)[1]
    return base64.b64decode(cleaned)


def figure_description(block: dict[str, Any]) -> str:
    html_value = str(block.get("html") or "")
    parts = [extract_img_alt(html_value), strip_tags(html_value)]
    deduped: list[str] = []
    for part in parts:
        normalized = normalize_text(part)
        if part and normalized and normalized not in {normalize_text(item) for item in deduped}:
            deduped.append(part)
    return " ".join(deduped).strip()


def is_boilerplate_figure(description: str) -> bool:
    normalized = normalize_text(description)
    return any(pattern in normalized for pattern in BOILERPLATE_PATTERNS)


def figure_score(
    *,
    query_tokens: set[str],
    block: dict[str, Any],
    description: str,
    hit_rank: int,
    hit_page: int | None,
) -> float:
    score = max(0.0, 30.0 - (hit_rank - 1) * 3.0)
    block_page = block.get("page")
    if isinstance(block_page, int) and hit_page is not None:
        page_delta = abs(block_page - hit_page)
        if page_delta == 0:
            score += 30.0
        elif page_delta == 1:
            score += 18.0
        elif page_delta == 2:
            score += 8.0
    score += 6.0 if block.get("block_type") == "ChemicalBlock" else 0.0
    figure_tokens = tokenize(description)
    overlap = query_tokens & figure_tokens
    score += float(len(overlap) * 6)
    if query_tokens & CHEMISTRY_TERMS and figure_tokens & CHEMISTRY_TERMS:
        score += 10.0
    return score


def select_figures(
    *,
    hits: list[dict[str, Any]],
    corpus_docs: list[CorpusDoc],
    query: str,
    max_figures: int,
    max_docs: int,
) -> list[FigureMatch]:
    selected: list[FigureMatch] = []
    seen_blocks: set[tuple[str, str]] = set()
    seen_docs: set[str] = set()
    query_tokens = tokenize(query)

    for hit in hits:
        matched_doc, match_score = match_hit_to_doc(hit, corpus_docs)
        if matched_doc is None or matched_doc.paper_id in seen_docs:
            continue

        chunks_path = matched_doc.extraction_path / "rag/chunks.json"
        if not chunks_path.exists():
            continue

        payload = load_json(chunks_path)
        blocks = payload.get("blocks") or []
        hit_page = hit.get("metadata", {}).get("page")
        if not isinstance(hit_page, int):
            try:
                hit_page = int(hit_page)
            except (TypeError, ValueError):
                hit_page = None

        for block in blocks:
            if block.get("block_type") not in FIGURE_TYPES:
                continue
            images = block.get("images") or {}
            if not isinstance(images, dict) or not images:
                continue
            description = figure_description(block)
            if not description or is_boilerplate_figure(description):
                continue
            score = match_score + figure_score(
                query_tokens=query_tokens,
                block=block,
                description=description,
                hit_rank=int(hit.get("rank") or 999),
                hit_page=hit_page,
            )
            if score <= 0:
                continue
            image_name, image_value = next(iter(images.items()))
            block_key = (matched_doc.paper_id, str(block.get("id") or image_name))
            if block_key in seen_blocks:
                continue
            try:
                image_bytes = decode_image_bytes(str(image_value))
            except (ValueError, binascii.Error):
                continue
            molecule_smiles, molecule_records = resolve_figure_molecules(
                matched_doc,
                block_id=str(block.get("id") or ""),
                image_name=str(image_name),
            )
            selected.append(
                FigureMatch(
                    doc=matched_doc,
                    block_id=str(block.get("id") or ""),
                    block_type=str(block.get("block_type") or ""),
                    page=block.get("page") if isinstance(block.get("page"), int) else None,
                    description=description,
                    image_name=str(image_name),
                    image_bytes=image_bytes,
                    score=score,
                    source_hit_rank=int(hit.get("rank") or 999),
                    source_hit_page=hit_page,
                    source_hit_title=str((hit.get("metadata") or {}).get("title") or ""),
                    molecule_smiles=molecule_smiles,
                    molecule_records=molecule_records,
                )
            )
            seen_blocks.add(block_key)

        seen_docs.add(matched_doc.paper_id)
        if len(seen_docs) >= max_docs:
            break

    selected.sort(key=lambda item: item.score, reverse=True)
    return selected[:max_figures]


def resolve_collection_id(url: str, token: str, *, kb_name: str | None, collection_id: str | None) -> tuple[str, str | None]:
    if collection_id:
        return collection_id, kb_name
    if not kb_name:
        raise SystemExit("Pass either --kb-name or --collection-id")
    payload = request_json(f"{url.rstrip('/')}/api/v1/knowledge/", token=token)
    for item in normalize_knowledge_items(payload):
        if str(item.get("name") or "") == kb_name:
            return str(item.get("id")), kb_name
    raise SystemExit(f"Knowledge base not found: {kb_name}")


def resolve_retrieval_target(
    *,
    retrieval_backend: str,
    openwebui_url: str,
    token: str,
    kb_name: str | None,
    collection_id: str | None,
) -> tuple[str | None, str | None]:
    if retrieval_backend == "canonical":
        return None, kb_name
    return resolve_collection_id(
        openwebui_url,
        token,
        kb_name=kb_name,
        collection_id=collection_id,
    )


def canonical_hit_to_retrieval_hit(hit: dict[str, Any]) -> dict[str, Any]:
    payload = hit.get("payload") or {}
    title = hit.get("title") or payload.get("document_title") or payload.get("title")
    return {
        "id": hit.get("point_id"),
        "document": payload.get("chunk_text")
        or payload.get("text")
        or payload.get("content")
        or payload.get("caption_text")
        or "",
        "metadata": {
            "title": title,
            "name": title,
            "source": payload.get("source") or payload.get("document_title") or title,
            "doc_id": hit.get("doc_id") or payload.get("document_id") or payload.get("doc_id"),
            "paper_id": hit.get("doc_id") or payload.get("document_id") or payload.get("doc_id"),
            "doi": payload.get("doi"),
            "citekey": hit.get("citation_key") or payload.get("citation_key"),
            "page": hit.get("page_number") or payload.get("page_number"),
        },
        "distance": None,
        "rank": int(hit.get("rank") or 0) or None,
    }


def extract_canonical_hits(payload: dict[str, Any]) -> list[dict[str, Any]]:
    hits = payload.get("reranked_hits")
    if not isinstance(hits, list):
        return []
    return [canonical_hit_to_retrieval_hit(hit) for hit in hits if isinstance(hit, dict)]


def retrieve_hits(
    openwebui_url: str,
    token: str,
    *,
    collection_id: str | None,
    query: str,
    k: int,
    retrieval_backend: str,
    multimodal_retrieval_api_url: str,
) -> list[dict[str, Any]]:
    if retrieval_backend == "canonical":
        payload = request_json(
            f"{multimodal_retrieval_api_url.rstrip('/')}/api/v1/retrieve",
            token="",
            method="POST",
            payload={"query": query, "top_k": k},
        )
        return extract_canonical_hits(payload)

    if not collection_id:
        raise SystemExit("collection_id is required for legacy-openwebui retrieval")

    payload = request_json(
        f"{openwebui_url.rstrip('/')}/api/v1/retrieval/query/doc",
        token=token,
        method="POST",
        payload={"collection_name": collection_id, "query": query, "k": k},
    )
    return extract_hits(payload)


def build_answer_context(query: str, hits: list[dict[str, Any]], figures: list[FigureMatch]) -> str:
    hit_lines = []
    for hit in hits[:5]:
        metadata = hit.get("metadata") or {}
        snippet = " ".join(str(hit.get("document") or "").split())
        hit_lines.append(
            f"- rank={hit.get('rank')} title={metadata.get('title') or metadata.get('name') or 'untitled'} "
            f"page={metadata.get('page')} snippet={snippet[:900]}"
        )

    figure_lines = []
    for figure in figures[:6]:
        figure_lines.append(
            f"- doc={figure.doc.title or figure.doc.paper_id} page={figure.page} "
            f"block_type={figure.block_type} smiles={', '.join(figure.molecule_smiles[:3]) or 'none'} "
            f"description={figure.description[:900]}"
        )

    return (
        "User question:\n"
        f"{query}\n\n"
        "Retrieved text evidence:\n"
        + ("\n".join(hit_lines) if hit_lines else "- none")
        + "\n\nRetrieved figure evidence:\n"
        + ("\n".join(figure_lines) if figure_lines else "- none")
    )


def generate_answer(url: str, token: str, *, model: str, query: str, hits: list[dict[str, Any]], figures: list[FigureMatch]) -> str:
    context = build_answer_context(query, hits, figures)
    payload = {
        "model": model,
        "messages": [
            {
                "role": "system",
                "content": (
                    "Answer using only the provided retrieval evidence. "
                    "If the evidence is insufficient, say so. "
                    "When figures are relevant, mention them explicitly as figure evidence."
                ),
            },
            {"role": "user", "content": context},
        ],
        "temperature": 0.2,
        "max_tokens": 500,
    }
    response = request_json(f"{url.rstrip('/')}/api/chat/completions", token=token, method="POST", payload=payload)
    choices = response.get("choices") or []
    if not choices:
        return ""
    message = choices[0].get("message") or {}
    return str(message.get("content") or "").strip()


def write_markdown(
    path: Path,
    *,
    query: str,
    kb_name: str | None,
    collection_id: str | None,
    answer: str,
    hits: list[dict[str, Any]],
    figures: list[tuple[FigureMatch, Path]],
) -> None:
    lines = [
        "# Multimodal Retrieval Answer",
        "",
        f"- Generated: {utc_now_iso()}",
        f"- Query: {query}",
        f"- Knowledge base: {kb_name or '(collection id only)'}",
        f"- Collection id: {collection_id or 'n/a'}",
        "",
        "## Answer",
        "",
        answer or "_No answer was generated._",
        "",
        "## Retrieval Hits",
        "",
    ]

    if hits:
        for hit in hits:
            metadata = hit.get("metadata") or {}
            title = metadata.get("title") or metadata.get("name") or "untitled"
            snippet = " ".join(str(hit.get("document") or "").split())
            lines.append(
                f"- Rank {hit.get('rank')}: {title} (page {metadata.get('page')}, distance {hit.get('distance')})"
            )
            lines.append(f"  Snippet: {snippet[:320]}")
    else:
        lines.append("- No hits returned.")

    lines.extend(["", "## Figures", ""])
    if figures:
        for figure, image_path in figures:
            lines.append(
                f"### {figure.doc.title or figure.doc.paper_id} | page {figure.page} | {figure.block_type}"
            )
            lines.append("")
            lines.append(f"Score: {figure.score:.1f} | Source hit rank: {figure.source_hit_rank}")
            lines.append("")
            if figure.molecule_smiles:
                lines.append(f"SMILES: `{', '.join(figure.molecule_smiles)}`")
                lines.append("")
            lines.append(figure.description)
            lines.append("")
            lines.append(f"![{figure.description[:120]}]({image_path.as_posix()})")
            lines.append("")
    else:
        lines.append("- No relevant figures were recovered from local migration artifacts.")
        lines.append("")

    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render canonical multimodal retrieval hits with migration-artifact figure attachments."
    )
    parser.add_argument("--query", required=True, help="User query to send to retrieval.")
    parser.add_argument("--kb-name", help="Legacy OpenWebUI knowledge base name.")
    parser.add_argument("--collection-id", help="Legacy OpenWebUI collection UUID.")
    parser.add_argument("--k", type=int, default=8, help="Top-k retrieval hits to request (default: 8).")
    parser.add_argument("--max-docs", type=int, default=3, help="Maximum matched corpus docs to mine for figures.")
    parser.add_argument("--max-figures", type=int, default=4, help="Maximum figures to export.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Default: artifacts/multimodal-answer/<timestamp>-<slug>",
    )
    parser.add_argument("--openwebui-url", default=os.environ.get("OPENWEBUI_URL", "http://localhost:3000"))
    parser.add_argument(
        "--multimodal-retrieval-api-url",
        default=os.environ.get("MULTIMODAL_RETRIEVAL_API_URL", "http://127.0.0.1:8510"),
    )
    parser.add_argument(
        "--retrieval-backend",
        choices=("canonical", "legacy-openwebui"),
        default=os.environ.get("RETRIEVAL_BACKEND", "canonical"),
        help="Retrieval backend to query.",
    )
    parser.add_argument("--chat-model", default=None, help="Model for answer synthesis. Default: first available model.")
    parser.add_argument("--no-answer", action="store_true", help="Skip answer synthesis and export retrieval evidence only.")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    token = os.environ.get("OPENWEBUI_MULTIMODAL_TOKEN") or os.environ.get("OPENWEBUI_API_KEY") or os.environ.get("API_KEY")
    if (args.retrieval_backend == "legacy-openwebui" or not args.no_answer) and not token:
        raise SystemExit("OPENWEBUI_MULTIMODAL_TOKEN or OPENWEBUI_API_KEY is required")

    repo_root = Path(__file__).resolve().parents[2]
    output_dir = args.out_dir or (
        repo_root
        / "artifacts"
        / "multimodal-answer"
        / f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{slugify(args.query)[:48]}"
    )
    images_dir = output_dir / "images"
    images_dir.mkdir(parents=True, exist_ok=True)

    collection_id, resolved_kb_name = resolve_retrieval_target(
        retrieval_backend=args.retrieval_backend,
        openwebui_url=args.openwebui_url,
        token=token or "",
        kb_name=args.kb_name,
        collection_id=args.collection_id,
    )
    hits = retrieve_hits(
        args.openwebui_url,
        token or "",
        collection_id=collection_id,
        query=args.query,
        k=args.k,
        retrieval_backend=args.retrieval_backend,
        multimodal_retrieval_api_url=args.multimodal_retrieval_api_url,
    )
    corpus_docs = build_corpus_docs(repo_root)
    figures = select_figures(
        hits=hits,
        corpus_docs=corpus_docs,
        query=args.query,
        max_figures=args.max_figures,
        max_docs=args.max_docs,
    )

    exported_figures: list[tuple[FigureMatch, Path]] = []
    json_figures: list[dict[str, Any]] = []
    for index, figure in enumerate(figures, start=1):
        extension = Path(figure.image_name).suffix or ".jpg"
        image_filename = f"{index:02d}-{slugify(figure.doc.paper_id)[:80]}{extension}"
        image_path = images_dir / image_filename
        image_path.write_bytes(figure.image_bytes)
        relative_path = image_path.relative_to(output_dir)
        exported_figures.append((figure, relative_path))
        json_figures.append(
            {
                "paper_id": figure.doc.paper_id,
                "title": figure.doc.title,
                "page": figure.page,
                "block_id": figure.block_id,
                "block_type": figure.block_type,
                "description": figure.description,
                "molecule_smiles": figure.molecule_smiles,
                "molecule_records": figure.molecule_records,
                "score": figure.score,
                "source_hit_rank": figure.source_hit_rank,
                "source_hit_page": figure.source_hit_page,
                "image_path": relative_path.as_posix(),
            }
        )

    answer = ""
    answer_model = ""
    if not args.no_answer:
        answer_model = args.chat_model or first_model_id(
            request_json(f"{args.openwebui_url.rstrip('/')}/api/models", token=token or "")
        )
        if answer_model:
            answer = generate_answer(
                args.openwebui_url,
                token or "",
                model=answer_model,
                query=args.query,
                hits=hits,
                figures=figures,
            )

    json_payload = {
        "generated_at": utc_now_iso(),
        "retrieval_backend": args.retrieval_backend,
        "multimodal_retrieval_api_url": args.multimodal_retrieval_api_url,
        "query": args.query,
        "knowledge_base": resolved_kb_name,
        "collection_id": collection_id,
        "chat_model": answer_model or None,
        "answer": answer,
        "hits": hits,
        "figures": json_figures,
    }

    write_json(output_dir / "answer.json", json_payload)
    write_markdown(
        output_dir / "answer.md",
        query=args.query,
        kb_name=resolved_kb_name,
        collection_id=collection_id,
        answer=answer,
        hits=hits,
        figures=exported_figures,
    )

    print(json.dumps({"out_dir": str(output_dir), "hits": len(hits), "figures": len(exported_figures)}, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
