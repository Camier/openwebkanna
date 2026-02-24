#!/usr/bin/env python3
"""
Transform prod_max to embed base64 images for multimodal RAG.
Uses normalized.json which already contains base64 images.
"""

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

PROD_MAX_DIR = Path("/LAB/@thesis/openwebui/prod_max")
OUTPUT_DIR = Path("/LAB/@thesis/openwebui/prod_max_multimodal")


def _looks_like_base64_image(value: object) -> bool:
    if not isinstance(value, str):
        return False
    return value.startswith("/9j/") or value.startswith("iVBOR")


def _as_data_uri(value: str) -> str:
    if value.startswith("data:image/"):
        return value
    if value.startswith("/9j/"):
        return f"data:image/jpeg;base64,{value}"
    if value.startswith("iVBOR"):
        return f"data:image/png;base64,{value}"
    return value


def _load_molecule_map(paper_dir: Path) -> Tuple[Dict[str, Dict], List[Dict]]:
    molecules_file = paper_dir / "smiles_extracted.json"
    if not molecules_file.exists():
        return {}, []

    try:
        with open(molecules_file) as f:
            data = json.load(f)
    except Exception:
        return {}, []

    mapping: Dict[str, Dict] = {}
    ordered: List[Dict] = []
    for item in data.get("picture_data", []):
        hint = {
            "molecule_smiles": item.get("molecule_smiles"),
            "status": item.get("status"),
            "description": item.get("description"),
            "raw_candidates": item.get("raw_candidates", []),
        }
        ordered.append(hint)
        image_path = item.get("image", "")
        if not isinstance(image_path, str) or not image_path:
            continue
        image_name = Path(image_path).name
        mapping[image_name] = hint

    return mapping, ordered


def process_paper(
    paper_dir: Path,
    dry_run: bool = False,
    include_molecule_metadata: bool = False,
) -> Dict:
    """Process a single paper using normalized.json."""
    paper_name = paper_dir.name

    norm_file = paper_dir / "normalized.json"
    if not norm_file.exists():
        return {"status": "skipped", "reason": "no normalized.json"}

    with open(norm_file) as f:
        data = json.load(f)

    blocks = data.get("raw", {}).get("chunks", {}).get("blocks", [])
    molecule_map, ordered_hints = (
        _load_molecule_map(paper_dir) if include_molecule_metadata else ({}, [])
    )
    hint_index = 0

    chunks_with_images = 0
    transformed_chunks = []

    for block in blocks:
        if block.get("images") and len(block.get("images", {})) > 0:
            images = block.get("images", {})
            has_base64 = any(_looks_like_base64_image(v) for v in images.values())
            molecule_hints = {}

            if has_base64:
                html = block.get("html", "")
                for filename, b64_data in images.items():
                    b64_data = _as_data_uri(b64_data)
                    html = re.sub(
                        rf'src=["\']?{re.escape(filename)}["\']?',
                        f'src="{b64_data}"',
                        html,
                    )
                    if include_molecule_metadata:
                        if filename in molecule_map:
                            molecule_hints[filename] = molecule_map[filename]
                        elif hint_index < len(ordered_hints):
                            molecule_hints[filename] = ordered_hints[hint_index]
                            hint_index += 1
                block["html"] = html
                chunks_with_images += 1
                if molecule_hints:
                    block["molecule_hints"] = molecule_hints

        if block.get("block_type") in [
            "Text",
            "SectionHeader",
            "Caption",
            "Figure",
            "Picture",
            "ChemicalBlock",
        ]:
            chunk = {
                "id": block.get("id", ""),
                "block_type": block.get("block_type", "Text"),
                "html": block.get("html", ""),
                "text": block.get("text", ""),
                "page": block.get("page", 0),
            }
            if include_molecule_metadata and block.get("molecule_hints"):
                chunk["molecule_hints"] = block.get("molecule_hints")
            transformed_chunks.append(chunk)

    if not dry_run:
        output_dir = OUTPUT_DIR / paper_name
        output_dir.mkdir(parents=True, exist_ok=True)

        with open(output_dir / "chunks.jsonl", "w") as f:
            for chunk in transformed_chunks:
                f.write(json.dumps(chunk, ensure_ascii=False) + "\n")

    return {
        "status": "success",
        "paper": paper_name,
        "total_chunks": len(transformed_chunks),
        "images_embedded": chunks_with_images,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--paper", type=str)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--include-molecule-metadata", action="store_true")
    args = parser.parse_args()

    paper_dirs = sorted(PROD_MAX_DIR.iterdir())

    if args.paper:
        paper_dirs = [d for d in paper_dirs if args.paper.lower() in d.name.lower()]
    if args.limit:
        paper_dirs = paper_dirs[: args.limit]

    print(f"Processing {len(paper_dirs)} papers...")

    results = []
    for i, paper_dir in enumerate(paper_dirs):
        if not paper_dir.is_dir():
            continue

        print(f"[{i + 1}/{len(paper_dirs)}] {paper_dir.name[:50]}...", end=" ")

        result = process_paper(
            paper_dir,
            dry_run=args.dry_run,
            include_molecule_metadata=args.include_molecule_metadata,
        )

        if result["status"] == "success":
            print(
                f"✓ {result['total_chunks']} chunks, {result['images_embedded']} images"
            )
        else:
            print(f"⊘ {result.get('reason')}")

        results.append(result)

    total_images = sum(r.get("images_embedded", 0) for r in results)
    print(f"\nTotal images embedded: {total_images}")

    if not args.dry_run:
        print(f"Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
