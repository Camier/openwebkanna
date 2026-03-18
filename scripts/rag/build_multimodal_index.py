#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from multimodal_index import (
    embed_images,
    embed_texts,
    iter_figure_records,
    load_clip_bundle,
    model_slug,
    safe_load_figure_image,
    save_index,
)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a local CLIP-backed figure index from extraction artifacts.")
    parser.add_argument(
        "--model-id",
        default="openai/clip-vit-base-patch32",
        help="Hugging Face model id for CLIP embeddings.",
    )
    parser.add_argument("--batch-size-text", type=int, default=32, help="Batch size for text embedding.")
    parser.add_argument("--batch-size-image", type=int, default=8, help="Batch size for image embedding.")
    parser.add_argument("--device", default="auto", help="Embedding device: auto, cpu, or cuda.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Default: artifacts/multimodal-index/<model-slug>",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    repo_root = Path(__file__).resolve().parents[2]
    out_dir = args.out_dir or repo_root / "artifacts" / "multimodal-index" / model_slug(args.model_id)

    print(f"[build] repo_root={repo_root}")
    print(f"[build] model_id={args.model_id}")
    print(f"[build] out_dir={out_dir}")

    records = iter_figure_records(repo_root)
    if not records:
        raise SystemExit("No figure records found under data/extractions/*/rag/chunks.json")
    print(f"[build] discovered_records={len(records)}")

    processor, model, device = load_clip_bundle(args.model_id, device=args.device)
    print(f"[build] device={device}")

    valid_records = []
    images = []
    for record in records:
        image = safe_load_figure_image(repo_root, record)
        if image is None:
            continue
        valid_records.append(record)
        images.append(image)
    print(f"[build] embeddable_records={len(valid_records)}")

    retrieval_texts = [record.retrieval_text for record in valid_records]
    text_embeddings = embed_texts(
        retrieval_texts,
        processor=processor,
        model=model,
        device=device,
        batch_size=args.batch_size_text,
    )
    image_embeddings = embed_images(
        images,
        processor=processor,
        model=model,
        device=device,
        batch_size=args.batch_size_image,
    )

    save_index(
        out_dir,
        model_id=args.model_id,
        records=valid_records,
        text_embeddings=text_embeddings,
        image_embeddings=image_embeddings,
    )
    print(f"[build] wrote manifest={out_dir / 'manifest.json'}")
    print(f"[build] wrote records={out_dir / 'records.jsonl'}")
    print(f"[build] wrote embeddings={out_dir / 'embeddings.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
