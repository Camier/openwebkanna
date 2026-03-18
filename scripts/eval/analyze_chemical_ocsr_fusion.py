#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT_DIR = REPO_ROOT / "artifacts" / "chemical-ocsr-eval" / "latest"
RISKY_CATEGORIES = {"placeholder_or_attachment", "radical_or_query_atom"}
REVIEW_CATEGORIES = {"stereo_only", "same_scaffold_variant", "substituent_drift", "close_analog", "scaffold_shift"}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def safe_slug(value: str) -> str:
    lowered = value.casefold()
    chars = [char if char.isalnum() else "-" for char in lowered]
    slug = "".join(chars)
    while "--" in slug:
        slug = slug.replace("--", "-")
    return slug.strip("-") or "paper"


def index_validated_rows(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for row in payload.get("results") or payload.get("picture_data") or []:
        if not isinstance(row, dict):
            continue
        if str(row.get("status") or "") != "validated":
            continue
        smiles = str(row.get("molecule_smiles") or "")
        if smiles and smiles not in indexed:
            indexed[smiles] = row
    return indexed


def tier_for_disagreement_category(category: str, matched_pair: bool) -> str:
    if category in RISKY_CATEGORIES:
        return "risky"
    if matched_pair or category in REVIEW_CATEGORIES:
        return "review_required"
    if category == "coverage_gap":
        return "complementary"
    return "review_required"


def normalize_record(record: dict[str, Any], backend: str) -> dict[str, Any]:
    return {
        "backend": backend,
        "smiles": str(record.get("molecule_smiles") or record.get("smiles") or ""),
        "confidence": record.get("confidence"),
        "detector_confidence": record.get("detector_confidence"),
        "block_id": str(record.get("block_id") or ""),
        "crop_index": int(record.get("crop_index") or 0),
        "description": str(record.get("description") or ""),
        "page": record.get("page"),
        "source_image": record.get("source_image"),
    }


def shared_candidates(
    shared_smiles: list[str],
    mg_rows: dict[str, dict[str, Any]],
    ms_rows: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for smiles in shared_smiles:
        mg_row = mg_rows.get(smiles)
        ms_row = ms_rows.get(smiles)
        if not mg_row or not ms_row:
            continue
        rows.append(
            {
                "tier": "consensus",
                "reason": "shared_exact_smiles",
                "smiles": smiles,
                "backends": ["molgrapher", "molscribe"],
                "molgrapher": normalize_record(mg_row, "molgrapher"),
                "molscribe": normalize_record(ms_row, "molscribe"),
            }
        )
    return rows


def build_fusion_payload(run_payload: dict[str, Any], disagreement_payload: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    disagreement_by_slug = {row["paper_slug"]: row for row in disagreement_payload["per_paper"]}
    per_paper_payload: list[dict[str, Any]] = []
    aggregate_tiers: Counter[str] = Counter()
    aggregate_reasons: Counter[str] = Counter()
    aggregate_backends: Counter[str] = Counter()
    upper_bound_union = 0

    for paper_run in run_payload["per_paper"]:
        paper_slug = paper_run["paper_slug"]
        disagreement_row = disagreement_by_slug[paper_slug]
        mg_payload = load_json(REPO_ROOT / paper_run["backend_runs"]["molgrapher"]["output_path"])
        ms_payload = load_json(REPO_ROOT / paper_run["backend_runs"]["molscribe"]["output_path"])
        mg_rows = index_validated_rows(mg_payload)
        ms_rows = index_validated_rows(ms_payload)
        overlap = paper_run["comparison"]["valid_smiles_overlap"]

        candidates: list[dict[str, Any]] = []
        candidates.extend(shared_candidates(overlap["intersection"], mg_rows, ms_rows))

        for pair in disagreement_row["matched_pairs"]:
            category = str(pair["category"])
            tier = tier_for_disagreement_category(category, matched_pair=True)
            for backend_name in ("molgrapher", "molscribe"):
                candidate_record = normalize_record(pair[backend_name], backend_name)
                candidates.append(
                    {
                        "tier": tier,
                        "reason": category,
                        "smiles": candidate_record["smiles"],
                        "backends": [backend_name],
                        "paired_with_backend": "molscribe" if backend_name == "molgrapher" else "molgrapher",
                        "paired_with_smiles": pair["molscribe"]["smiles"] if backend_name == "molgrapher" else pair["molgrapher"]["smiles"],
                        "similarity": pair["similarity"],
                        backend_name: candidate_record,
                    }
                )

        for row in disagreement_row["unmatched"]:
            category = str(row["category"])
            backend_name = str(row["backend"])
            tier = tier_for_disagreement_category(category, matched_pair=False)
            candidate_record = normalize_record(row["record"], backend_name)
            candidates.append(
                {
                    "tier": tier,
                    "reason": category,
                    "smiles": candidate_record["smiles"],
                    "backends": [backend_name],
                    backend_name: candidate_record,
                }
            )

        by_tier = Counter(candidate["tier"] for candidate in candidates)
        by_reason = Counter(candidate["reason"] for candidate in candidates)
        by_backend = Counter()
        for candidate in candidates:
            for backend_name in candidate["backends"]:
                by_backend[backend_name] += 1

        upper_bound = len(overlap["intersection"]) + len(overlap["left_only"]) + len(overlap["right_only"])
        upper_bound_union += upper_bound
        aggregate_tiers.update(by_tier)
        aggregate_reasons.update(by_reason)
        aggregate_backends.update(by_backend)
        per_paper_payload.append(
            {
                "paper_name": paper_run["paper_name"],
                "paper_slug": paper_slug,
                "counts": {
                    "upper_bound_union": upper_bound,
                    "consensus": by_tier["consensus"],
                    "complementary": by_tier["complementary"],
                    "review_required": by_tier["review_required"],
                    "risky": by_tier["risky"],
                    "auto_keep": by_tier["consensus"] + by_tier["complementary"],
                },
                "reasons": dict(sorted(by_reason.items())),
                "backend_contributions": dict(sorted(by_backend.items())),
                "candidates": candidates,
            }
        )

    return {
        "source_run_json": str((out_dir / "run.json").relative_to(REPO_ROOT)),
        "generated_from_run_at_utc": run_payload["run_at_utc"],
        "aggregate": {
            "upper_bound_union": upper_bound_union,
            "consensus": aggregate_tiers["consensus"],
            "complementary": aggregate_tiers["complementary"],
            "review_required": aggregate_tiers["review_required"],
            "risky": aggregate_tiers["risky"],
            "auto_keep": aggregate_tiers["consensus"] + aggregate_tiers["complementary"],
            "reasons": dict(sorted(aggregate_reasons.items())),
            "backend_contributions": dict(sorted(aggregate_backends.items())),
        },
        "per_paper": per_paper_payload,
    }


def render_report(run_payload: dict[str, Any], fusion_payload: dict[str, Any], out_path: Path) -> None:
    agg = fusion_payload["aggregate"]
    lines = [
        "# Chemical OCSR Fusion Analysis",
        "",
        f"- Source run: `{run_payload['run_at_utc']}`",
        f"- Papers: `{len(run_payload['per_paper'])}`",
        f"- Union upper bound: `{agg['upper_bound_union']}` candidate SMILES",
        f"- Auto-keep (`consensus + complementary`): `{agg['auto_keep']}`",
        f"- Review required: `{agg['review_required']}`",
        f"- Risky: `{agg['risky']}`",
        "",
        "## Aggregate Tier Summary",
        "",
        "| Tier | Candidate Count | Meaning |",
        "| --- | ---: | --- |",
        f"| consensus | {agg['consensus']} | exact agreement across both backends |",
        f"| complementary | {agg['complementary']} | backend-unique validated candidate without a competing analog |",
        f"| review_required | {agg['review_required']} | competing analog/variant disagreement between backends |",
        f"| risky | {agg['risky']} | placeholder, attachment, or radical/query-style artifact |",
        "",
        "## Aggregate Reasons",
        "",
        "| Reason | Candidate Count |",
        "| --- | ---: |",
    ]
    for reason, count in agg["reasons"].items():
        lines.append(f"| {reason} | {count} |")

    lines.extend(
        [
            "",
            "## Per Paper",
            "",
            "| Paper | Union Upper Bound | Auto-Keep | Consensus | Complementary | Review Required | Risky |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for paper in fusion_payload["per_paper"]:
        counts = paper["counts"]
        lines.append(
            "| {paper} | {upper_bound_union} | {auto_keep} | {consensus} | {complementary} | {review_required} | {risky} |".format(
                paper=paper["paper_slug"],
                **counts,
            )
        )

    for paper in fusion_payload["per_paper"]:
        counts = paper["counts"]
        lines.extend(
            [
                "",
                f"### {paper['paper_slug']}",
                "",
                "- Union upper bound `{upper_bound_union}` | Auto-keep `{auto_keep}` | Review required `{review_required}` | Risky `{risky}`".format(
                    **counts
                ),
                f"- Backend contributions: `{paper['backend_contributions']}`",
                f"- Reason mix: `{paper['reasons']}`",
                "",
            ]
        )
        by_tier: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for candidate in paper["candidates"]:
            by_tier[candidate["tier"]].append(candidate)
        for tier_name in ("consensus", "complementary", "review_required", "risky"):
            examples = by_tier[tier_name][:3]
            if not examples:
                continue
            lines.append(f"- {tier_name} examples:")
            for candidate in examples:
                notes = [candidate["smiles"], f"reason={candidate['reason']}"]
                if "similarity" in candidate:
                    notes.append(f"sim={candidate['similarity']:.3f}")
                lines.append(f"  - `{', '.join(notes)}`")
        lines.append("")

    lines.extend(
        [
            "## Notes",
            "",
            "- `Union upper bound` is the deduplicated RDKit-valid union across both backends. It is a coverage ceiling, not a correctness guarantee.",
            "- `Auto-keep` keeps exact backend agreement plus backend-unique candidates without a competing analog in the other backend.",
            "- `Review required` means both backends found structurally related but non-identical molecules on the same benchmark slice.",
            "- `Risky` flags placeholder, attachment-point, or radical/query-style artifacts that should stay out of automated downstream ingestion.",
        ]
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_review_queue_payload(fusion_payload: dict[str, Any]) -> dict[str, Any]:
    index_rows: list[dict[str, Any]] = []
    per_paper_rows: list[dict[str, Any]] = []
    aggregate = Counter()

    for paper in fusion_payload["per_paper"]:
        grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for candidate in paper["candidates"]:
            grouped[candidate["tier"]].append(candidate)
        row = {
            "paper_name": paper["paper_name"],
            "paper_slug": paper["paper_slug"],
            "paper_file_slug": safe_slug(paper["paper_slug"]),
            "counts": paper["counts"],
            "reasons": paper["reasons"],
            "backend_contributions": paper["backend_contributions"],
            "auto_keep": grouped["consensus"] + grouped["complementary"],
            "review_required": grouped["review_required"],
            "risky": grouped["risky"],
        }
        per_paper_rows.append(row)
        index_rows.append(
            {
                "paper_name": paper["paper_name"],
                "paper_slug": paper["paper_slug"],
                "paper_file_slug": row["paper_file_slug"],
                "counts": paper["counts"],
            }
        )
        aggregate.update(paper["counts"])

    return {
        "generated_from_run_at_utc": fusion_payload["generated_from_run_at_utc"],
        "source_run_json": fusion_payload["source_run_json"],
        "aggregate": {
            "paper_count": len(per_paper_rows),
            "upper_bound_union": aggregate["upper_bound_union"],
            "auto_keep": aggregate["auto_keep"],
            "review_required": aggregate["review_required"],
            "risky": aggregate["risky"],
        },
        "index": index_rows,
        "per_paper": per_paper_rows,
    }


def render_review_queue_markdown(review_queue_payload: dict[str, Any], out_path: Path) -> None:
    aggregate = review_queue_payload["aggregate"]
    lines = [
        "# Chemical OCSR Review Queue",
        "",
        f"- Source run: `{review_queue_payload['generated_from_run_at_utc']}`",
        f"- Papers: `{aggregate['paper_count']}`",
        f"- Auto-keep: `{aggregate['auto_keep']}`",
        f"- Review required: `{aggregate['review_required']}`",
        f"- Risky: `{aggregate['risky']}`",
        "",
        "## Papers",
        "",
        "| Paper | Auto-Keep | Review Required | Risky | Queue File |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for row in review_queue_payload["index"]:
        counts = row["counts"]
        lines.append(
            "| {paper} | {auto_keep} | {review_required} | {risky} | `review_queue/{file_slug}.md` |".format(
                paper=row["paper_slug"],
                file_slug=row["paper_file_slug"],
                **counts,
            )
        )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_paper_queue_markdown(paper_payload: dict[str, Any], out_path: Path) -> None:
    counts = paper_payload["counts"]
    lines = [
        f"# Review Queue: {paper_payload['paper_slug']}",
        "",
        f"- Auto-keep: `{counts['auto_keep']}`",
        f"- Review required: `{counts['review_required']}`",
        f"- Risky: `{counts['risky']}`",
        f"- Backend contributions: `{paper_payload['backend_contributions']}`",
        f"- Reason mix: `{paper_payload['reasons']}`",
        "",
    ]

    sections = (
        ("auto_keep", "Auto-Keep"),
        ("review_required", "Review Required"),
        ("risky", "Risky"),
    )
    for key, title in sections:
        rows = paper_payload[key]
        lines.extend([f"## {title}", ""])
        if not rows:
            lines.append("- None")
            lines.append("")
            continue
        for row in rows:
            detail = [row["smiles"], f"reason={row['reason']}"]
            if "similarity" in row:
                detail.append(f"sim={row['similarity']:.3f}")
            detail.append(f"backend={','.join(row['backends'])}")
            lines.append(f"- `{', '.join(detail)}`")
        lines.append("")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def default_adjudication_action(tier: str) -> str:
    if tier == "consensus":
        return "accept"
    if tier == "complementary":
        return "accept_review_optional"
    if tier == "review_required":
        return "needs_review"
    return "reject_pending_review"


def flatten_candidate_row(paper_payload: dict[str, Any], candidate: dict[str, Any], candidate_index: int) -> dict[str, Any]:
    molgrapher = candidate.get("molgrapher") or {}
    molscribe = candidate.get("molscribe") or {}
    return {
        "queue_id": f"{paper_payload['paper_file_slug']}:{candidate['tier']}:{candidate_index:03d}",
        "paper_name": paper_payload["paper_name"],
        "paper_slug": paper_payload["paper_slug"],
        "paper_file_slug": paper_payload["paper_file_slug"],
        "tier": candidate["tier"],
        "reason": candidate["reason"],
        "default_action": default_adjudication_action(candidate["tier"]),
        "smiles": candidate["smiles"],
        "backends": ",".join(candidate["backends"]),
        "paired_with_backend": candidate.get("paired_with_backend", ""),
        "paired_with_smiles": candidate.get("paired_with_smiles", ""),
        "similarity": candidate.get("similarity", ""),
        "molgrapher_smiles": molgrapher.get("smiles", ""),
        "molgrapher_confidence": molgrapher.get("confidence", ""),
        "molgrapher_detector_confidence": molgrapher.get("detector_confidence", ""),
        "molgrapher_page": molgrapher.get("page", ""),
        "molgrapher_block_id": molgrapher.get("block_id", ""),
        "molgrapher_crop_index": molgrapher.get("crop_index", ""),
        "molgrapher_description": molgrapher.get("description", ""),
        "molgrapher_source_image": molgrapher.get("source_image", ""),
        "molscribe_smiles": molscribe.get("smiles", ""),
        "molscribe_confidence": molscribe.get("confidence", ""),
        "molscribe_detector_confidence": molscribe.get("detector_confidence", ""),
        "molscribe_page": molscribe.get("page", ""),
        "molscribe_block_id": molscribe.get("block_id", ""),
        "molscribe_crop_index": molscribe.get("crop_index", ""),
        "molscribe_description": molscribe.get("description", ""),
        "molscribe_source_image": molscribe.get("source_image", ""),
        "review_status": "pending",
        "review_decision": "",
        "reviewed_smiles": "",
        "reviewer": "",
        "review_notes": "",
    }


def build_adjudication_exports(review_queue_payload: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    flat_rows: list[dict[str, Any]] = []
    per_paper_rows: dict[str, list[dict[str, Any]]] = {}
    for paper_payload in review_queue_payload["per_paper"]:
        paper_rows: list[dict[str, Any]] = []
        candidate_index = 1
        for key in ("auto_keep", "review_required", "risky"):
            for candidate in paper_payload[key]:
                row = flatten_candidate_row(paper_payload, candidate, candidate_index)
                flat_rows.append(row)
                paper_rows.append(row)
                candidate_index += 1
        per_paper_rows[paper_payload["paper_file_slug"]] = paper_rows
    return flat_rows, per_paper_rows


def write_adjudication_exports(review_queue_payload: dict[str, Any], out_dir: Path, queue_dir: Path) -> None:
    fieldnames = [
        "queue_id",
        "paper_name",
        "paper_slug",
        "paper_file_slug",
        "tier",
        "reason",
        "default_action",
        "smiles",
        "backends",
        "paired_with_backend",
        "paired_with_smiles",
        "similarity",
        "molgrapher_smiles",
        "molgrapher_confidence",
        "molgrapher_detector_confidence",
        "molgrapher_page",
        "molgrapher_block_id",
        "molgrapher_crop_index",
        "molgrapher_description",
        "molgrapher_source_image",
        "molscribe_smiles",
        "molscribe_confidence",
        "molscribe_detector_confidence",
        "molscribe_page",
        "molscribe_block_id",
        "molscribe_crop_index",
        "molscribe_description",
        "molscribe_source_image",
        "review_status",
        "review_decision",
        "reviewed_smiles",
        "reviewer",
        "review_notes",
    ]
    flat_rows, per_paper_rows = build_adjudication_exports(review_queue_payload)
    write_jsonl(out_dir / "review_queue.jsonl", flat_rows)
    write_csv(out_dir / "review_queue.csv", flat_rows, fieldnames)
    for paper_payload in review_queue_payload["per_paper"]:
        file_slug = paper_payload["paper_file_slug"]
        paper_rows = per_paper_rows[file_slug]
        write_jsonl(queue_dir / f"{file_slug}.jsonl", paper_rows)
        write_csv(queue_dir / f"{file_slug}.csv", paper_rows, fieldnames)


def write_review_queue(fusion_payload: dict[str, Any], out_dir: Path) -> None:
    queue_dir = out_dir / "review_queue"
    queue_dir.mkdir(parents=True, exist_ok=True)
    review_queue_payload = build_review_queue_payload(fusion_payload)
    write_json(out_dir / "review_queue.json", review_queue_payload)
    render_review_queue_markdown(review_queue_payload, out_dir / "review_queue.md")
    write_adjudication_exports(review_queue_payload, out_dir, queue_dir)
    for paper_payload in review_queue_payload["per_paper"]:
        file_slug = paper_payload["paper_file_slug"]
        write_json(queue_dir / f"{file_slug}.json", paper_payload)
        render_paper_queue_markdown(paper_payload, queue_dir / f"{file_slug}.md")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize a conservative fusion candidate set from MolGrapher vs MolScribe OCSR benchmark artifacts.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help=f"Artifact directory containing run.json and disagreement_analysis.json. Default: {DEFAULT_OUT_DIR}")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO_ROOT / out_dir

    run_payload = load_json(out_dir / "run.json")
    disagreement_payload = load_json(out_dir / "disagreement_analysis.json")
    fusion_payload = build_fusion_payload(run_payload, disagreement_payload, out_dir)
    write_json(out_dir / "fusion_candidates.json", fusion_payload)
    render_report(run_payload, fusion_payload, out_dir / "fusion_report.md")
    write_review_queue(fusion_payload, out_dir)
    print(f"[fusion] json={out_dir / 'fusion_candidates.json'}")
    print(f"[fusion] report={out_dir / 'fusion_report.md'}")
    print(f"[fusion] review_queue={out_dir / 'review_queue.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
