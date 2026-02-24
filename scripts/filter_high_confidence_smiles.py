#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


DEFAULT_ROOT = Path("/LAB/@thesis/openwebui/prod_max")


def _iter_papers(root: Path, paper_filter: str | None, limit: int | None):
    papers = [item for item in sorted(root.iterdir()) if item.is_dir()]
    if paper_filter:
        token = paper_filter.lower()
        papers = [item for item in papers if token in item.name.lower()]
    if limit:
        papers = papers[:limit]
    return papers


def _safe_int(value, default: int) -> int:
    try:
        return int(value)
    except Exception:
        return default


def _load_rdkit():
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    return Chem, Descriptors


def _classify_smiles(smiles: str, chem, descriptors, min_heavy_atoms: int):
    text = (smiles or "").strip()
    if not text:
        return False, "empty"
    if "*" in text:
        return False, "wildcard"

    try:
        mol = chem.MolFromSmiles(text)
    except Exception:
        return False, "rdkit_parse_error"

    if mol is None:
        return False, "invalid"

    heavy_atoms = _safe_int(mol.GetNumHeavyAtoms(), 0)
    if heavy_atoms < min_heavy_atoms:
        return False, f"too_small_{heavy_atoms}"

    mw = float(descriptors.MolWt(mol))
    if mw < 50.0:
        return False, "low_mw"

    canonical = chem.MolToSmiles(mol)
    return True, canonical


def main() -> None:
    parser = argparse.ArgumentParser(description="Create high-confidence SMILES JSONL")
    parser.add_argument("--root", type=str, default=str(DEFAULT_ROOT))
    parser.add_argument("--paper", type=str)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--min-heavy-atoms", type=int, default=6)
    args = parser.parse_args()

    root = Path(args.root)
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"Invalid root: {root}")

    chem, descriptors = _load_rdkit()
    papers = _iter_papers(root, args.paper, args.limit)
    if not papers:
        raise SystemExit("No matching papers")

    total_rows = 0
    total_kept = 0

    for paper_dir in papers:
        source = paper_dir / "molecules.jsonl"
        target = paper_dir / "molecules_high_confidence.jsonl"
        if not source.exists():
            target.write_text("")
            print(f"{paper_dir.name}: source missing, wrote empty high-confidence file")
            continue

        kept_rows = []
        with source.open() as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                total_rows += 1
                try:
                    payload = json.loads(line)
                except Exception:
                    continue

                smiles = str(payload.get("smiles", ""))
                keep, detail = _classify_smiles(
                    smiles,
                    chem,
                    descriptors,
                    args.min_heavy_atoms,
                )
                if not keep:
                    continue

                payload["smiles"] = detail
                payload["quality"] = "high_confidence"
                kept_rows.append(payload)
                total_kept += 1

        with target.open("w") as handle:
            for row in kept_rows:
                handle.write(json.dumps(row, ensure_ascii=False) + "\n")

        print(f"{paper_dir.name}: kept {len(kept_rows)} high-confidence rows")

    print(f"\nHigh-confidence summary: rows_in={total_rows}, rows_kept={total_kept}")


if __name__ == "__main__":
    main()
