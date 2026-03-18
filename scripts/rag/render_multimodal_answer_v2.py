#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from multimodal_index import (
    cosine_scores,
    embed_texts,
    load_clip_bundle,
    load_figure_image_bytes,
    load_index,
    model_slug,
    tokenize,
)
from render_multimodal_answer import (
    FigureMatch,
    build_corpus_docs,
    first_model_id,
    generate_answer,
    is_boilerplate_figure,
    match_hit_to_doc,
    request_json,
    resolve_figure_molecules,
    resolve_collection_id,
    retrieve_hits,
    slugify,
    utc_now_iso,
    write_json,
    write_markdown,
)


@dataclass
class HitMatch:
    match_score: float
    hit_rank: int
    hit_page: int | None
    hit_title: str


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render answers with a CLIP-backed figure index fused with OpenWebUI text retrieval."
    )
    parser.add_argument("--query", required=True, help="User query to send to multimodal retrieval.")
    parser.add_argument("--kb-name", help="OpenWebUI knowledge base name.")
    parser.add_argument("--collection-id", help="OpenWebUI collection UUID.")
    parser.add_argument("--k", type=int, default=8, help="Top-k text retrieval hits to request (default: 8).")
    parser.add_argument("--max-figures", type=int, default=4, help="Maximum figures to export.")
    parser.add_argument(
        "--figure-index-dir",
        type=Path,
        default=None,
        help="Figure index directory. Default: artifacts/multimodal-index/<clip-model-slug>",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Default: artifacts/multimodal-answer-v2/<timestamp>-<slug>",
    )
    parser.add_argument("--openwebui-url", default=os.environ.get("OPENWEBUI_URL", "http://localhost:3000"))
    parser.add_argument("--chat-model", default=None, help="Model for answer synthesis. Default: first available model.")
    parser.add_argument("--no-answer", action="store_true", help="Skip answer synthesis and export retrieval evidence only.")
    return parser.parse_args(argv)


def normalize_hit_page(value: Any) -> int | None:
    if isinstance(value, int):
        return value
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def build_hit_map(hits: list[dict[str, Any]], corpus_docs: list[Any]) -> dict[str, HitMatch]:
    matched: dict[str, HitMatch] = {}
    for hit in hits:
        doc, match_score = match_hit_to_doc(hit, corpus_docs)
        if doc is None:
            continue
        metadata = hit.get("metadata") or {}
        candidate = HitMatch(
            match_score=match_score,
            hit_rank=int(hit.get("rank") or 999),
            hit_page=normalize_hit_page(metadata.get("page")),
            hit_title=str(metadata.get("title") or metadata.get("name") or ""),
        )
        previous = matched.get(doc.paper_id)
        if previous is None or candidate.hit_rank < previous.hit_rank or candidate.match_score > previous.match_score:
            matched[doc.paper_id] = candidate
    return matched


def overlap_score(query: str, text: str) -> float:
    return float(len(tokenize(query) & tokenize(text)) * 2)


def hit_rank_boost(hit_rank: int | None) -> float:
    if hit_rank is None:
        return 0.0
    return max(0.0, 18.0 - (hit_rank - 1) * 3.0)


def page_boost(figure_page: int | None, hit_page: int | None) -> float:
    if figure_page is None or hit_page is None:
        return 0.0
    delta = abs(figure_page - hit_page)
    if delta == 0:
        return 12.0
    if delta == 1:
        return 6.0
    if delta == 2:
        return 2.5
    return 0.0

def main(argv: list[str]) -> int:
    args = parse_args(argv)
    token = (
        os.environ.get("OPENWEBUI_MULTIMODAL_TOKEN")
        or os.environ.get("OPENWEBUI_API_KEY")
        or os.environ.get("API_KEY")
    )
    if not token:
        raise SystemExit("OPENWEBUI_MULTIMODAL_TOKEN or OPENWEBUI_API_KEY is required")

    repo_root = Path(__file__).resolve().parents[2]
    default_index_dir = repo_root / "artifacts" / "multimodal-index" / model_slug("openai/clip-vit-base-patch32")
    figure_index_dir = args.figure_index_dir or default_index_dir
    if not figure_index_dir.exists():
        raise SystemExit(
            f"Figure index not found: {figure_index_dir}. Run ./scripts/rag/build-multimodal-index.sh first."
        )

    output_dir = args.out_dir or (
        repo_root
        / "artifacts"
        / "multimodal-answer-v2"
        / f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{slugify(args.query)[:48]}"
    )
    images_dir = output_dir / "images"
    images_dir.mkdir(parents=True, exist_ok=True)

    collection_id, resolved_kb_name = resolve_collection_id(
        args.openwebui_url,
        token,
        kb_name=args.kb_name,
        collection_id=args.collection_id,
    )
    hits = retrieve_hits(args.openwebui_url, token, collection_id=collection_id, query=args.query, k=args.k)
    corpus_docs = build_corpus_docs(repo_root)
    corpus_doc_map = {doc.paper_id: doc for doc in corpus_docs}
    hit_map = build_hit_map(hits, corpus_docs)

    manifest, records, text_embeddings, image_embeddings = load_index(figure_index_dir)
    processor, model, device = load_clip_bundle(str(manifest.get("model_id") or "openai/clip-vit-base-patch32"))
    query_embedding = embed_texts([args.query], processor=processor, model=model, device=device)[0]

    selected = []
    text_scores = cosine_scores(query_embedding, text_embeddings)
    image_scores = cosine_scores(query_embedding, image_embeddings)
    ranked: list[tuple[float, int, dict[str, Any]]] = []
    for index, record in enumerate(records):
        if is_boilerplate_figure(record.description):
            continue
        multimodal_score = float(image_scores[index] * 60.0 + text_scores[index] * 40.0)
        hit_match = hit_map.get(record.paper_id)
        fused_score = multimodal_score + overlap_score(args.query, record.retrieval_text)
        if hit_match is not None:
            fused_score += 8.0
            fused_score += hit_rank_boost(hit_match.hit_rank)
            fused_score += page_boost(record.page, hit_match.hit_page)
        ranked.append(
            (
                fused_score,
                index,
                {
                    "multimodal_score": multimodal_score,
                    "image_score": float(image_scores[index]),
                    "text_score": float(text_scores[index]),
                },
            )
        )
    ranked.sort(key=lambda item: item[0], reverse=True)

    seen_blocks: set[str] = set()
    per_doc: dict[str, int] = {}
    exported_figures: list[tuple[FigureMatch, Path]] = []
    json_figures: list[dict[str, Any]] = []
    selected_matches: list[FigureMatch] = []

    for fused_score, index, score_breakdown in ranked:
        record = records[index]
        if record.figure_id in seen_blocks:
            continue
        if per_doc.get(record.paper_id, 0) >= 2:
            continue
        try:
            image_bytes = load_figure_image_bytes(repo_root, record)
        except Exception:
            continue
        matched_doc = corpus_doc_map.get(record.paper_id)
        if matched_doc is None:
            continue
        hit_match = hit_map.get(record.paper_id)
        molecule_smiles, molecule_records = resolve_figure_molecules(
            matched_doc,
            block_id=record.block_id,
            image_name=record.image_name,
        )
        figure = FigureMatch(
            doc=matched_doc,
            block_id=record.block_id,
            block_type=record.block_type,
            page=record.page,
            description=record.description,
            image_name=record.image_name,
            image_bytes=image_bytes,
            score=fused_score,
            source_hit_rank=hit_match.hit_rank if hit_match else 999,
            source_hit_page=hit_match.hit_page if hit_match else None,
            source_hit_title=hit_match.hit_title if hit_match else record.title,
            molecule_smiles=molecule_smiles,
            molecule_records=molecule_records,
        )
        image_filename = f"{len(exported_figures)+1:02d}-{slugify(record.paper_id)[:80]}{Path(record.image_name).suffix or '.jpg'}"
        image_path = images_dir / image_filename
        image_path.write_bytes(image_bytes)
        relative_path = image_path.relative_to(output_dir)

        exported_figures.append((figure, relative_path))
        selected_matches.append(figure)
        json_figures.append(
            {
                "figure_id": record.figure_id,
                "paper_id": record.paper_id,
                "title": record.title,
                "page": record.page,
                "block_id": record.block_id,
                "block_type": record.block_type,
                "description": record.description,
                "caption": record.caption,
                "molecule_smiles": figure.molecule_smiles,
                "molecule_records": figure.molecule_records,
                "score": fused_score,
                "source_hit_rank": figure.source_hit_rank,
                "source_hit_page": figure.source_hit_page,
                "image_path": relative_path.as_posix(),
                "scores": score_breakdown,
            }
        )
        seen_blocks.add(record.figure_id)
        per_doc[record.paper_id] = per_doc.get(record.paper_id, 0) + 1
        if len(exported_figures) >= args.max_figures:
            break

    answer = ""
    answer_model = ""
    if not args.no_answer:
        answer_model = args.chat_model or first_model_id(request_json(f"{args.openwebui_url.rstrip('/')}/api/models", token=token))
        if answer_model:
            answer = generate_answer(
                args.openwebui_url,
                token,
                model=answer_model,
                query=args.query,
                hits=hits,
                figures=selected_matches,
            )

    json_payload = {
        "generated_at": utc_now_iso(),
        "query": args.query,
        "knowledge_base": resolved_kb_name,
        "collection_id": collection_id,
        "figure_index_dir": str(figure_index_dir.relative_to(repo_root) if figure_index_dir.is_relative_to(repo_root) else figure_index_dir),
        "figure_index_model": manifest.get("model_id"),
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
    print(f"[v2] output_dir={output_dir}")
    print(f"[v2] hits={len(hits)}")
    print(f"[v2] figures={len(exported_figures)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
