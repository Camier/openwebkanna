#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


DEFAULT_ROOT = Path("/LAB/@thesis/openwebui/prod_max")


def _load_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def _iter_papers(root: Path, paper_filter: str | None, limit: int | None):
    papers = [item for item in sorted(root.iterdir()) if item.is_dir()]
    if paper_filter:
        token = paper_filter.lower()
        papers = [item for item in papers if token in item.name.lower()]
    if limit:
        papers = papers[:limit]
    return papers


def _validate_smiles_json(path: Path):
    required = {
        "paper",
        "paper_dir",
        "pictures_found",
        "pictures_processed",
        "pictures_with_smiles",
        "rdkit_available",
        "results",
        "picture_data",
        "smiles_extracted",
        "valid_smiles",
        "invalid_smiles",
    }
    payload = _load_json(path)
    errors: list[str] = []

    missing = [key for key in sorted(required) if key not in payload]
    if missing:
        errors.append(f"missing_keys={','.join(missing)}")

    for key in (
        "smiles_extracted",
        "valid_smiles",
        "invalid_smiles",
        "results",
        "picture_data",
    ):
        if key in payload and not isinstance(payload[key], list):
            errors.append(f"{key}_not_list")

    if "rdkit_available" in payload and not isinstance(
        payload["rdkit_available"], bool
    ):
        errors.append("rdkit_available_not_bool")

    if "pictures_processed" in payload and "results" in payload:
        if isinstance(payload["pictures_processed"], int) and isinstance(
            payload["results"], list
        ):
            if payload["pictures_processed"] != len(payload["results"]):
                errors.append("pictures_processed_mismatch_results")

    return payload, errors


def _validate_molecules_jsonl(path: Path):
    if not path.exists():
        return 0, ["missing_molecules_jsonl"]

    required = {
        "paper",
        "paper_dir",
        "image",
        "page",
        "bbox",
        "status",
        "description",
        "smiles",
    }
    count = 0
    errors: list[str] = []

    with path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            count += 1
            try:
                payload = json.loads(text)
            except Exception:
                errors.append(f"invalid_jsonl_line_{line_no}")
                continue

            missing = [key for key in sorted(required) if key not in payload]
            if missing:
                errors.append(f"line_{line_no}_missing_keys={','.join(missing)}")
            if "smiles" in payload and not isinstance(payload["smiles"], str):
                errors.append(f"line_{line_no}_smiles_not_str")

    return count, errors


def main():
    parser = argparse.ArgumentParser(description="Validate smiles extraction outputs")
    parser.add_argument("--root", type=str, default=str(DEFAULT_ROOT))
    parser.add_argument("--paper", type=str)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    root = Path(args.root)
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"Invalid root: {root}")

    papers = _iter_papers(root, args.paper, args.limit)
    if not papers:
        raise SystemExit("No papers matched for validation")

    total = len(papers)
    ok = 0
    failures = 0

    for paper_dir in papers:
        smiles_path = paper_dir / "smiles_extracted.json"
        molecules_path = paper_dir / "molecules.jsonl"

        issues: list[str] = []
        smiles_payload = None

        if not smiles_path.exists():
            issues.append("missing_smiles_extracted_json")
        else:
            try:
                smiles_payload, errs = _validate_smiles_json(smiles_path)
                issues.extend(errs)
            except Exception as error:
                issues.append(f"smiles_json_read_error={error}")

        mol_count = 0
        if smiles_payload is not None:
            try:
                mol_count, errs = _validate_molecules_jsonl(molecules_path)
                issues.extend(errs)
                valid_smiles_count = len(smiles_payload.get("valid_smiles", []))
                if mol_count < valid_smiles_count:
                    issues.append("molecules_rows_less_than_valid_smiles")
            except Exception as error:
                issues.append(f"molecules_jsonl_read_error={error}")

        if issues:
            failures += 1
            print(f"FAIL {paper_dir.name}: {'; '.join(issues)}")
        else:
            ok += 1
            print(f"OK   {paper_dir.name}: validated")

    print(f"\nValidation summary: total={total}, ok={ok}, failed={failures}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
