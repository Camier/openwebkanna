#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert SMILES JSONL rows into retrieval-friendly markdown"
    )
    parser.add_argument("--input", required=True, help="Path to input JSONL file")
    parser.add_argument("--output", required=True, help="Path to output markdown file")
    parser.add_argument(
        "--title",
        default="SMILES High Confidence Molecules",
        help="Top-level markdown heading",
    )
    parser.add_argument(
        "--paper-id",
        default="",
        help="Optional paper identifier to include in the header",
    )
    return parser.parse_args()


def iter_rows(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(payload, dict):
                yield payload


def row_md(index: int, row: dict) -> str:
    smiles = str(row.get("smiles", "")).strip() or "N/A"
    backend = str(row.get("backend", "")).strip() or "unknown"
    confidence = row.get("confidence")
    if confidence in (None, ""):
        confidence_text = "n/a"
    else:
        confidence_text = str(confidence)
    image = str(row.get("image", "")).strip() or "N/A"
    quality = str(row.get("quality", "")).strip() or "unspecified"
    source_file = str(row.get("source_file", "")).strip() or "N/A"
    notes = str(row.get("notes", "")).strip()

    lines = [
        f"### Molecule {index}",
        f"- SMILES: `{smiles}`",
        f"- Backend: `{backend}`",
        f"- Confidence: `{confidence_text}`",
        f"- Quality: `{quality}`",
        f"- Image: `{image}`",
        f"- Source file: `{source_file}`",
    ]

    if notes:
        lines.append(f"- Notes: {notes}")

    for key in sorted(row.keys()):
        if key in {
            "smiles",
            "backend",
            "confidence",
            "quality",
            "image",
            "source_file",
            "notes",
        }:
            continue
        value = row.get(key)
        if value in (None, ""):
            continue
        if isinstance(value, (dict, list)):
            value_text = json.dumps(value, ensure_ascii=False)
        else:
            value_text = str(value)
        lines.append(f"- {key}: `{value_text}`")

    return "\n".join(lines)


def build_markdown(title: str, paper_id: str, rows: list[dict]) -> str:
    parts = [f"# {title}", ""]

    if paper_id:
        parts.extend([f"- Paper: `{paper_id}`", ""])

    parts.extend([f"- Molecule count: `{len(rows)}`", ""])

    if not rows:
        parts.append("No valid molecule rows were found in the source JSONL.")
        return "\n".join(parts) + "\n"

    parts.extend(["## Molecules", ""])
    for idx, row in enumerate(rows, start=1):
        parts.append(row_md(idx, row))
        parts.append("")

    return "\n".join(parts).rstrip() + "\n"


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists() or not input_path.is_file():
        raise SystemExit(f"Input file not found: {input_path}")

    rows = list(iter_rows(input_path))
    markdown = build_markdown(args.title, args.paper_id, rows)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(markdown, encoding="utf-8")


if __name__ == "__main__":
    main()
