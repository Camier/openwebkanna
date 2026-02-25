#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


PROD_MAX_DIR = Path("/LAB/@thesis/openwebui/prod_max")
MULTIMODAL_DIR = Path("/LAB/@thesis/openwebui/prod_max_multimodal")


def _iter_papers(root: Path, paper_filter: str | None, limit: int | None):
    papers = [item for item in sorted(root.iterdir()) if item.is_dir()]
    if paper_filter:
        token = paper_filter.lower()
        papers = [item for item in papers if token in item.name.lower()]
    if limit:
        papers = papers[:limit]
    return papers


def _count_jsonl_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def _count_hint_rows(path: Path) -> int:
    if not path.exists():
        return 0
    count = 0
    with path.open() as handle:
        for line in handle:
            if '"molecule_hints"' in line:
                count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate SMILES QA summary")
    parser.add_argument("--root", type=str, default=str(PROD_MAX_DIR))
    parser.add_argument("--multimodal-root", type=str, default=str(MULTIMODAL_DIR))
    parser.add_argument("--paper", type=str)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=str, default="smiles_qa_summary.json")
    args = parser.parse_args()

    root = Path(args.root)
    mm_root = Path(args.multimodal_root)
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"Invalid root: {root}")

    papers = _iter_papers(root, args.paper, args.limit)
    if not papers:
        raise SystemExit("No papers matched")

    report = {
        "root": str(root),
        "multimodal_root": str(mm_root),
        "paper_count": len(papers),
        "totals": {
            "images_processed": 0,
            "pictures_with_smiles": 0,
            "smiles_candidates": 0,
            "valid_smiles": 0,
            "invalid_smiles": 0,
            "molecules_rows": 0,
            "high_conf_rows": 0,
            "chunks_with_hints": 0,
        },
        "papers": [],
    }

    for paper_dir in papers:
        smiles_file = paper_dir / "smiles_extracted.json"
        molecules_file = paper_dir / "molecules.jsonl"
        high_conf_file = paper_dir / "molecules_high_confidence.jsonl"
        chunks_file = mm_root / paper_dir.name / "chunks.jsonl"

        if smiles_file.exists():
            with smiles_file.open() as handle:
                smiles_payload = json.load(handle)
        else:
            smiles_payload = {}

        paper_metrics = {
            "paper": paper_dir.name,
            "images_processed": int(smiles_payload.get("pictures_processed", 0) or 0),
            "pictures_with_smiles": int(
                smiles_payload.get("pictures_with_smiles", 0) or 0
            ),
            "smiles_candidates": len(smiles_payload.get("smiles_extracted", []) or []),
            "valid_smiles": len(smiles_payload.get("valid_smiles", []) or []),
            "invalid_smiles": len(smiles_payload.get("invalid_smiles", []) or []),
            "backend": smiles_payload.get("backend"),
            "molecules_rows": _count_jsonl_rows(molecules_file),
            "high_conf_rows": _count_jsonl_rows(high_conf_file),
            "chunks_with_hints": _count_hint_rows(chunks_file),
        }
        report["papers"].append(paper_metrics)

        for key in report["totals"]:
            report["totals"][key] += int(paper_metrics[key])

    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = root / output_path
    with output_path.open("w") as handle:
        json.dump(report, handle, indent=2)

    print(f"Wrote QA summary: {output_path}")
    print(json.dumps(report["totals"], indent=2))


if __name__ == "__main__":
    main()
