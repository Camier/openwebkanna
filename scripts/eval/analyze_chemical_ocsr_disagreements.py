#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT_DIR = REPO_ROOT / "artifacts" / "chemical-ocsr-eval" / "latest"
PAIR_SIMILARITY_THRESHOLD = 0.35


@dataclass(frozen=True)
class MoleculeRecord:
    backend: str
    smiles: str
    nonisomeric_smiles: str
    scaffold_smiles: str
    formula: str
    exact_mw: float
    heavy_atoms: int
    ring_count: int
    aromatic_ring_count: int
    confidence: float | None
    detector_confidence: float | None
    block_id: str
    crop_index: int
    description: str
    flags: tuple[str, ...]
    fingerprint: Any


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def canonical_nonisomeric_smiles(mol: Chem.Mol) -> str:
    clone = Chem.Mol(mol)
    Chem.RemoveStereochemistry(clone)
    return Chem.MolToSmiles(clone, canonical=True)


def scaffold_smiles(mol: Chem.Mol) -> str:
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        return ""
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return ""
    return Chem.MolToSmiles(scaffold, canonical=True)


def molecule_flags(smiles: str, mol: Chem.Mol) -> tuple[str, ...]:
    flags: set[str] = set()
    if "*" in smiles:
        flags.add("dummy_atom")
    if any(token in smiles for token in ("[Rf]", "[R", "[1*]", "[2*]", "[3*]")):
        flags.add("placeholder_group")
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        flags.add("dummy_atom")
    if any(atom.GetNumRadicalElectrons() > 0 for atom in mol.GetAtoms()):
        flags.add("radical")
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        flags.add("charged")
    if any(atom.IsInRing() for atom in mol.GetAtoms()):
        flags.add("ring_system")
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        flags.add("aromatic")
    return tuple(sorted(flags))


def build_molecule_record(backend: str, row: dict[str, Any]) -> MoleculeRecord:
    smiles = str(row.get("molecule_smiles") or "")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse validated SMILES: {smiles}")
    return MoleculeRecord(
        backend=backend,
        smiles=smiles,
        nonisomeric_smiles=canonical_nonisomeric_smiles(mol),
        scaffold_smiles=scaffold_smiles(mol),
        formula=rdMolDescriptors.CalcMolFormula(mol),
        exact_mw=float(Descriptors.ExactMolWt(mol)),
        heavy_atoms=int(mol.GetNumHeavyAtoms()),
        ring_count=int(rdMolDescriptors.CalcNumRings(mol)),
        aromatic_ring_count=int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        confidence=float(row["confidence"]) if isinstance(row.get("confidence"), (int, float)) else None,
        detector_confidence=float(row["detector_confidence"])
        if isinstance(row.get("detector_confidence"), (int, float))
        else None,
        block_id=str(row.get("block_id") or ""),
        crop_index=int(row.get("crop_index") or 0),
        description=str(row.get("description") or ""),
        flags=molecule_flags(smiles, mol),
        fingerprint=AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048),
    )


def index_validated_records(payload: dict[str, Any], backend: str) -> dict[str, MoleculeRecord]:
    indexed: dict[str, MoleculeRecord] = {}
    for row in payload.get("results") or payload.get("picture_data") or []:
        if not isinstance(row, dict):
            continue
        if str(row.get("status") or "") != "validated":
            continue
        smiles = str(row.get("molecule_smiles") or "")
        if not smiles or smiles in indexed:
            continue
        indexed[smiles] = build_molecule_record(backend, row)
    return indexed


def tanimoto(left: MoleculeRecord, right: MoleculeRecord) -> float:
    return float(DataStructs.TanimotoSimilarity(left.fingerprint, right.fingerprint))


def classify_pair(left: MoleculeRecord, right: MoleculeRecord, similarity: float) -> str:
    if left.nonisomeric_smiles == right.nonisomeric_smiles and left.smiles != right.smiles:
        return "stereo_only"
    if "dummy_atom" in left.flags or "dummy_atom" in right.flags or "placeholder_group" in left.flags or "placeholder_group" in right.flags:
        return "placeholder_or_attachment"
    if left.scaffold_smiles and left.scaffold_smiles == right.scaffold_smiles:
        if similarity >= 0.8:
            return "substituent_drift"
        return "same_scaffold_variant"
    if similarity >= 0.55:
        return "close_analog"
    return "scaffold_shift"


def classify_single(record: MoleculeRecord) -> str:
    if "dummy_atom" in record.flags or "placeholder_group" in record.flags:
        return "placeholder_or_attachment"
    if "radical" in record.flags:
        return "radical_or_query_atom"
    return "coverage_gap"


def greedy_pairs(
    molgrapher_only: list[MoleculeRecord],
    molscribe_only: list[MoleculeRecord],
) -> tuple[list[dict[str, Any]], list[MoleculeRecord], list[MoleculeRecord]]:
    scored_pairs: list[tuple[float, int, int]] = []
    for left_idx, left in enumerate(molgrapher_only):
        for right_idx, right in enumerate(molscribe_only):
            scored_pairs.append((tanimoto(left, right), left_idx, right_idx))
    scored_pairs.sort(reverse=True)

    used_left: set[int] = set()
    used_right: set[int] = set()
    matched: list[dict[str, Any]] = []
    for similarity, left_idx, right_idx in scored_pairs:
        if similarity < PAIR_SIMILARITY_THRESHOLD:
            break
        if left_idx in used_left or right_idx in used_right:
            continue
        left = molgrapher_only[left_idx]
        right = molscribe_only[right_idx]
        category = classify_pair(left, right, similarity)
        matched.append(
            {
                "category": category,
                "similarity": similarity,
                "molgrapher": serialize_record(left),
                "molscribe": serialize_record(right),
            }
        )
        used_left.add(left_idx)
        used_right.add(right_idx)

    unmatched_left = [record for idx, record in enumerate(molgrapher_only) if idx not in used_left]
    unmatched_right = [record for idx, record in enumerate(molscribe_only) if idx not in used_right]
    return matched, unmatched_left, unmatched_right


def serialize_record(record: MoleculeRecord) -> dict[str, Any]:
    return {
        "backend": record.backend,
        "smiles": record.smiles,
        "nonisomeric_smiles": record.nonisomeric_smiles,
        "scaffold_smiles": record.scaffold_smiles,
        "formula": record.formula,
        "exact_mw": round(record.exact_mw, 4),
        "heavy_atoms": record.heavy_atoms,
        "ring_count": record.ring_count,
        "aromatic_ring_count": record.aromatic_ring_count,
        "confidence": record.confidence,
        "detector_confidence": record.detector_confidence,
        "block_id": record.block_id,
        "crop_index": record.crop_index,
        "flags": list(record.flags),
        "description": record.description,
    }


def build_paper_analysis(paper_run: dict[str, Any]) -> dict[str, Any]:
    mg_payload = load_json(REPO_ROOT / paper_run["backend_runs"]["molgrapher"]["output_path"])
    ms_payload = load_json(REPO_ROOT / paper_run["backend_runs"]["molscribe"]["output_path"])
    mg_records = index_validated_records(mg_payload, "molgrapher")
    ms_records = index_validated_records(ms_payload, "molscribe")
    overlap = paper_run["comparison"]["valid_smiles_overlap"]

    mg_only = [mg_records[smiles] for smiles in overlap["left_only"] if smiles in mg_records]
    ms_only = [ms_records[smiles] for smiles in overlap["right_only"] if smiles in ms_records]
    matched_pairs, unmatched_mg, unmatched_ms = greedy_pairs(mg_only, ms_only)

    unmatched_rows: list[dict[str, Any]] = []
    for record in unmatched_mg:
        unmatched_rows.append(
            {
                "category": classify_single(record),
                "backend": "molgrapher",
                "record": serialize_record(record),
            }
        )
    for record in unmatched_ms:
        unmatched_rows.append(
            {
                "category": classify_single(record),
                "backend": "molscribe",
                "record": serialize_record(record),
            }
        )

    category_counts = Counter(row["category"] for row in matched_pairs)
    category_counts.update(row["category"] for row in unmatched_rows)
    return {
        "paper_name": paper_run["paper_name"],
        "paper_slug": paper_run["paper_slug"],
        "jaccard": overlap["jaccard"],
        "counts": {
            "shared_valid_smiles": len(overlap["intersection"]),
            "molgrapher_only_valid_smiles": len(overlap["left_only"]),
            "molscribe_only_valid_smiles": len(overlap["right_only"]),
            "matched_disagreement_pairs": len(matched_pairs),
            "unmatched_disagreements": len(unmatched_rows),
        },
        "categories": dict(sorted(category_counts.items())),
        "matched_pairs": matched_pairs,
        "unmatched": unmatched_rows,
    }


def render_report(run_payload: dict[str, Any], analysis_payload: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# Chemical OCSR Disagreement Analysis",
        "",
        f"- Source run: `{analysis_payload['source_run_json']}`",
        f"- Detector: `{run_payload['config']['detector']}`",
        f"- Papers: `{len(analysis_payload['per_paper'])}`",
        "",
        "## Aggregate Categories",
        "",
        "| Category | Count |",
        "| --- | ---: |",
    ]
    for category, count in analysis_payload["aggregate"]["categories"].items():
        lines.append(f"| {category} | {count} |")

    lines.extend(["", "## Per Paper", ""])
    for paper in analysis_payload["per_paper"]:
        lines.append(f"### {paper['paper_slug']}")
        lines.append("")
        lines.append(
            "- Shared `{shared}` | MolGrapher-only `{mg}` | MolScribe-only `{ms}` | Matched pairs `{pairs}` | Unmatched `{unmatched}` | Jaccard `{jaccard:.3f}`".format(
                shared=paper["counts"]["shared_valid_smiles"],
                mg=paper["counts"]["molgrapher_only_valid_smiles"],
                ms=paper["counts"]["molscribe_only_valid_smiles"],
                pairs=paper["counts"]["matched_disagreement_pairs"],
                unmatched=paper["counts"]["unmatched_disagreements"],
                jaccard=paper["jaccard"],
            )
        )
        if paper["categories"]:
            lines.append(
                "- Categories: "
                + ", ".join(f"`{name}={count}`" for name, count in paper["categories"].items())
            )
        else:
            lines.append("- Categories: `none`")
        lines.append("")

        for pair in paper["matched_pairs"][:5]:
            mg = pair["molgrapher"]
            ms = pair["molscribe"]
            lines.append(
                "- Pair `{category}` sim=`{similarity:.3f}`: MolGrapher `{mg_smiles}` vs MolScribe `{ms_smiles}`".format(
                    category=pair["category"],
                    similarity=pair["similarity"],
                    mg_smiles=mg["smiles"],
                    ms_smiles=ms["smiles"],
                )
            )
        for row in paper["unmatched"][:5]:
            record = row["record"]
            lines.append(
                "- Unmatched `{category}` from `{backend}`: `{smiles}`".format(
                    category=row["category"],
                    backend=row["backend"],
                    smiles=record["smiles"],
                )
            )
        lines.append("")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify MolGrapher vs MolScribe disagreement sets with RDKit-backed structure diagnostics.")
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Artifact directory containing run.json from run-chemical-ocsr-eval.sh. Default: artifacts/chemical-ocsr-eval/latest",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO_ROOT / out_dir
    run_json_path = out_dir / "run.json"
    run_payload = load_json(run_json_path)

    per_paper = [build_paper_analysis(paper_run) for paper_run in run_payload["per_paper"]]
    aggregate_categories = Counter()
    for paper in per_paper:
        aggregate_categories.update(paper["categories"])
    analysis_payload = {
        "source_run_json": str(run_json_path.relative_to(REPO_ROOT)),
        "generated_from_run_at_utc": run_payload["run_at_utc"],
        "per_paper": per_paper,
        "aggregate": {
            "categories": dict(sorted(aggregate_categories.items())),
        },
    }

    write_json(out_dir / "disagreement_analysis.json", analysis_payload)
    render_report(run_payload, analysis_payload, out_dir / "disagreement_analysis.md")
    print(f"[analysis] json={out_dir / 'disagreement_analysis.json'}")
    print(f"[analysis] report={out_dir / 'disagreement_analysis.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
