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
ACCEPT_DECISIONS = {"accept", "accept_as_is", "accept_review_optional"}
REJECT_DECISIONS = {"reject", "reject_pending_review"}


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            rows.append(json.loads(stripped))
    return rows


def load_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


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


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Materialize accepted, rejected, and pending OCSR candidates from the review queue export.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help=f"Artifact directory containing review_queue.csv or review_queue.jsonl. Default: {DEFAULT_OUT_DIR}")
    parser.add_argument("--input", help="Optional explicit queue file. Supports .csv and .jsonl. Defaults to review_queue.csv, then review_queue.jsonl.")
    parser.add_argument(
        "--apply-default-actions",
        action="store_true",
        help="Apply default_action to rows that have not been explicitly reviewed. Without this flag, pending rows stay pending until review_decision is filled.",
    )
    return parser.parse_args(argv)


def resolve_input_path(out_dir: Path, explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit)
        return path if path.is_absolute() else REPO_ROOT / path
    for candidate in (out_dir / "review_queue.csv", out_dir / "review_queue.jsonl"):
        if candidate.exists():
            return candidate
    raise SystemExit(f"No review queue input found in {out_dir}; expected review_queue.csv or review_queue.jsonl")


def load_rows(path: Path) -> list[dict[str, Any]]:
    if path.suffix == ".csv":
        return load_csv(path)
    if path.suffix == ".jsonl":
        return load_jsonl(path)
    raise SystemExit(f"Unsupported input format: {path.name}")


def normalize_value(value: Any) -> str:
    return str(value or "").strip()


def parse_backend_list(value: Any) -> list[str]:
    raw = normalize_value(value)
    if not raw:
        return []
    return [part.strip() for part in raw.split(",") if part.strip()]


def adjudication_state(row: dict[str, Any], apply_default_actions: bool) -> tuple[str, str]:
    review_status = normalize_value(row.get("review_status")).casefold()
    review_decision = normalize_value(row.get("review_decision")).casefold()
    default_action = normalize_value(row.get("default_action")).casefold()

    decision = review_decision
    if not decision and apply_default_actions:
        decision = default_action

    if review_status == "reviewed" and decision in ACCEPT_DECISIONS:
        return "accepted", decision
    if review_status == "reviewed" and decision in REJECT_DECISIONS:
        return "rejected", decision
    if review_status == "reviewed":
        return "pending", decision or "reviewed_unknown"

    if not apply_default_actions:
        return "pending", "pending_review"

    if decision in ACCEPT_DECISIONS:
        return "accepted", decision
    if decision in REJECT_DECISIONS:
        return "rejected", decision
    return "pending", decision or "pending_review"


def materialize_smiles(row: dict[str, Any], status: str) -> str:
    if status != "accepted":
        return ""
    reviewed_smiles = normalize_value(row.get("reviewed_smiles"))
    if reviewed_smiles:
        return reviewed_smiles
    return normalize_value(row.get("smiles"))


def annotate_rows(rows: list[dict[str, Any]], apply_default_actions: bool) -> list[dict[str, Any]]:
    annotated: list[dict[str, Any]] = []
    for row in rows:
        status, decision = adjudication_state(row, apply_default_actions)
        enriched = dict(row)
        enriched["adjudication_status"] = status
        enriched["effective_decision"] = decision
        enriched["effective_smiles"] = materialize_smiles(enriched, status)
        annotated.append(enriched)
    return annotated


def backend_snapshot(row: dict[str, Any], backend: str) -> dict[str, str]:
    prefix = f"{backend}_"
    return {
        "backend": backend,
        "smiles": normalize_value(row.get(f"{prefix}smiles")),
        "confidence": normalize_value(row.get(f"{prefix}confidence")),
        "detector_confidence": normalize_value(row.get(f"{prefix}detector_confidence")),
        "page": normalize_value(row.get(f"{prefix}page")),
        "block_id": normalize_value(row.get(f"{prefix}block_id")),
        "crop_index": normalize_value(row.get(f"{prefix}crop_index")),
        "description": normalize_value(row.get(f"{prefix}description")),
        "source_image": normalize_value(row.get(f"{prefix}source_image")),
    }


def build_accepted_catalog_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    catalog_rows: list[dict[str, Any]] = []
    for row in rows:
        effective_smiles = normalize_value(row.get("effective_smiles"))
        if not effective_smiles:
            continue

        source_backends = parse_backend_list(row.get("backends"))
        provenance = [backend_snapshot(row, backend) for backend in source_backends]
        provenance = [entry for entry in provenance if any(entry.values())]
        matched_backends = [
            entry["backend"] for entry in provenance if entry.get("smiles") == effective_smiles
        ]
        primary_backend = matched_backends[0] if matched_backends else (source_backends[0] if source_backends else "")
        primary = next((entry for entry in provenance if entry["backend"] == primary_backend), {})
        secondary = next((entry for entry in provenance if entry["backend"] != primary_backend), {})

        catalog_rows.append(
            {
                "queue_id": normalize_value(row.get("queue_id")),
                "paper_name": normalize_value(row.get("paper_name")),
                "paper_slug": normalize_value(row.get("paper_slug")),
                "paper_file_slug": normalize_value(row.get("paper_file_slug")),
                "tier": normalize_value(row.get("tier")),
                "reason": normalize_value(row.get("reason")),
                "effective_decision": normalize_value(row.get("effective_decision")),
                "effective_smiles": effective_smiles,
                "source_backends": ",".join(source_backends),
                "matched_backends": ",".join(matched_backends),
                "primary_backend": normalize_value(primary.get("backend")),
                "primary_smiles": normalize_value(primary.get("smiles")),
                "primary_confidence": normalize_value(primary.get("confidence")),
                "primary_detector_confidence": normalize_value(primary.get("detector_confidence")),
                "primary_page": normalize_value(primary.get("page")),
                "primary_block_id": normalize_value(primary.get("block_id")),
                "primary_crop_index": normalize_value(primary.get("crop_index")),
                "primary_source_image": normalize_value(primary.get("source_image")),
                "secondary_backend": normalize_value(secondary.get("backend")),
                "secondary_smiles": normalize_value(secondary.get("smiles")),
                "secondary_confidence": normalize_value(secondary.get("confidence")),
                "secondary_detector_confidence": normalize_value(secondary.get("detector_confidence")),
                "secondary_page": normalize_value(secondary.get("page")),
                "secondary_block_id": normalize_value(secondary.get("block_id")),
                "secondary_crop_index": normalize_value(secondary.get("crop_index")),
                "secondary_source_image": normalize_value(secondary.get("source_image")),
                "provenance": provenance,
            }
        )
    return catalog_rows


def render_report(summary: dict[str, Any], out_path: Path) -> None:
    aggregate = summary["aggregate"]
    lines = [
        "# Chemical OCSR Adjudication Export",
        "",
        f"- Input: `{summary['input_path']}`",
        f"- Default actions applied: `{summary['apply_default_actions']}`",
        f"- Rows: `{aggregate['rows_total']}`",
        f"- Accepted: `{aggregate['accepted']}`",
        f"- Rejected: `{aggregate['rejected']}`",
        f"- Pending: `{aggregate['pending']}`",
        "",
        "## Per Paper",
        "",
        "| Paper | Accepted | Rejected | Pending |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in summary["per_paper"]:
        lines.append(
            "| {paper_slug} | {accepted} | {rejected} | {pending} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `accepted` rows populate `effective_smiles`, using `reviewed_smiles` when present and falling back to the candidate `smiles` otherwise.",
            "- `rejected` rows stay out of downstream ingestion exports.",
            "- `pending` rows remain unresolved; if `--apply-default-actions` is enabled, only rows whose default action is non-final remain pending.",
        ]
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def mode_prefix(apply_default_actions: bool) -> str:
    return "adjudication_defaulted" if apply_default_actions else "adjudication_manual"


def build_summary(rows: list[dict[str, Any]], input_path: Path, apply_default_actions: bool) -> dict[str, Any]:
    aggregate = Counter()
    by_paper: dict[str, Counter[str]] = defaultdict(Counter)

    for row in rows:
        status = row["adjudication_status"]
        paper_slug = normalize_value(row.get("paper_slug")) or normalize_value(row.get("paper_file_slug")) or "paper"
        aggregate[status] += 1
        aggregate["rows_total"] += 1
        by_paper[paper_slug][status] += 1

    per_paper = []
    for paper_slug in sorted(by_paper):
        counts = by_paper[paper_slug]
        per_paper.append(
            {
                "paper_slug": paper_slug,
                "paper_file_slug": safe_slug(paper_slug),
                "accepted": counts["accepted"],
                "rejected": counts["rejected"],
                "pending": counts["pending"],
            }
        )

    return {
        "input_path": str(input_path.relative_to(REPO_ROOT) if input_path.is_absolute() else input_path),
        "apply_default_actions": apply_default_actions,
        "aggregate": {
            "rows_total": aggregate["rows_total"],
            "accepted": aggregate["accepted"],
            "rejected": aggregate["rejected"],
            "pending": aggregate["pending"],
        },
        "per_paper": per_paper,
    }


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO_ROOT / out_dir
    input_path = resolve_input_path(out_dir, args.input)
    rows = load_rows(input_path)
    annotated_rows = annotate_rows(rows, args.apply_default_actions)
    summary = build_summary(annotated_rows, input_path, args.apply_default_actions)

    fieldnames = list(annotated_rows[0].keys()) if annotated_rows else [
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
        "adjudication_status",
        "effective_decision",
        "effective_smiles",
    ]

    by_status: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in annotated_rows:
        by_status[row["adjudication_status"]].append(row)

    accepted_rows = by_status["accepted"]
    rejected_rows = by_status["rejected"]
    pending_rows = by_status["pending"]

    prefix = mode_prefix(args.apply_default_actions)
    write_json(out_dir / f"{prefix}_summary.json", summary)
    render_report(summary, out_dir / f"{prefix}_report.md")
    write_jsonl(out_dir / f"{prefix}_accepted_candidates.jsonl", accepted_rows)
    write_jsonl(out_dir / f"{prefix}_rejected_candidates.jsonl", rejected_rows)
    write_jsonl(out_dir / f"{prefix}_pending_candidates.jsonl", pending_rows)
    write_csv(out_dir / f"{prefix}_accepted_candidates.csv", accepted_rows, fieldnames)
    write_csv(out_dir / f"{prefix}_rejected_candidates.csv", rejected_rows, fieldnames)
    write_csv(out_dir / f"{prefix}_pending_candidates.csv", pending_rows, fieldnames)

    accepted_catalog_rows = build_accepted_catalog_rows(accepted_rows)
    accepted_catalog_fieldnames = [
        "queue_id",
        "paper_name",
        "paper_slug",
        "paper_file_slug",
        "tier",
        "reason",
        "effective_decision",
        "effective_smiles",
        "source_backends",
        "matched_backends",
        "primary_backend",
        "primary_smiles",
        "primary_confidence",
        "primary_detector_confidence",
        "primary_page",
        "primary_block_id",
        "primary_crop_index",
        "primary_source_image",
        "secondary_backend",
        "secondary_smiles",
        "secondary_confidence",
        "secondary_detector_confidence",
        "secondary_page",
        "secondary_block_id",
        "secondary_crop_index",
        "secondary_source_image",
    ]
    write_jsonl(out_dir / f"{prefix}_accepted_catalog.jsonl", accepted_catalog_rows)
    write_csv(out_dir / f"{prefix}_accepted_catalog.csv", accepted_catalog_rows, accepted_catalog_fieldnames)

    print(f"[adjudication] summary={out_dir / f'{prefix}_summary.json'}")
    print(f"[adjudication] report={out_dir / f'{prefix}_report.md'}")
    print(f"[adjudication] accepted={out_dir / f'{prefix}_accepted_candidates.jsonl'}")
    print(f"[adjudication] rejected={out_dir / f'{prefix}_rejected_candidates.jsonl'}")
    print(f"[adjudication] pending={out_dir / f'{prefix}_pending_candidates.jsonl'}")
    print(f"[adjudication] accepted_catalog={out_dir / f'{prefix}_accepted_catalog.jsonl'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
