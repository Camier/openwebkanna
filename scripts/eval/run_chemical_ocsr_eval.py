#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
EXTRACTIONS_DIR = REPO_ROOT / "data" / "extractions"
DEFAULT_PAPER_FILTERS = (
    "20260131_170728_9279ee63_1977_-_Capps_-_Sceletium_alkaloids_Part_7_Structure_and_abso",
    "20260131_171507_9336cc60_2007_-_Stafford_-_Monoamine_oxidase_inhibition_by_southern_A",
    "20260131_171629_97b0222d_2008_-_Gericke_-_Sceletium_a_review_update",
)
DEFAULT_OUT_ROOT = REPO_ROOT / "artifacts" / "chemical-ocsr-eval"
DEFAULT_DETECTOR = "moldetv2-general"
DEFAULT_OUTPUT_TEMPLATE = "smiles_extracted_{backend}_bench.json"


@dataclass(frozen=True)
class BackendSpec:
    name: str
    python_bin: str


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def safe_slug(value: str) -> str:
    lowered = value.casefold()
    chars = [char if char.isalnum() else "-" for char in lowered]
    slug = "".join(chars)
    while "--" in slug:
        slug = slug.replace("--", "-")
    return slug.strip("-") or "run"


def percentile(values: list[float], pct: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    index = math.ceil((pct / 100.0) * len(ordered)) - 1
    index = max(0, min(index, len(ordered) - 1))
    return float(ordered[index])


def average(values: list[float]) -> float:
    if not values:
        return 0.0
    return float(sum(values) / len(values))


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def iter_extraction_dirs() -> list[Path]:
    return sorted(path for path in EXTRACTIONS_DIR.iterdir() if path.is_dir() and (path / "normalized.json").exists())


def resolve_paper_dirs(filters: list[str]) -> list[Path]:
    candidates = iter_extraction_dirs()
    resolved: list[Path] = []
    for paper_filter in filters:
        matches = [path for path in candidates if paper_filter.casefold() in path.name.casefold()]
        if not matches:
            raise SystemExit(f"No extraction directory matched filter: {paper_filter}")
        if len(matches) > 1:
            joined = ", ".join(path.name for path in matches)
            raise SystemExit(f"Filter is ambiguous: {paper_filter} -> {joined}")
        resolved.append(matches[0])
    return resolved


def resolve_backend_specs(args: argparse.Namespace) -> list[BackendSpec]:
    specs = [
        BackendSpec(
            name="molgrapher",
            python_bin=str(args.molgrapher_python_bin or os.environ.get("MOLGRAPHER_PYTHON_BIN", "")).strip(),
        ),
        BackendSpec(
            name="molscribe",
            python_bin=str(args.molscribe_python_bin or os.environ.get("MOLSCRIBE_PYTHON_BIN", "")).strip(),
        ),
    ]
    missing = [spec.name for spec in specs if not spec.python_bin]
    if missing:
        raise SystemExit(
            "Missing backend python interpreter(s): "
            + ", ".join(missing)
            + ". Pass --molgrapher-python-bin/--molscribe-python-bin or export the matching env vars."
        )
    return specs


def build_out_dir(raw_out_dir: str | None, filters: list[str]) -> Path:
    if raw_out_dir:
        out_dir = Path(raw_out_dir)
        return out_dir if out_dir.is_absolute() else REPO_ROOT / out_dir
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = safe_slug("-".join(filters))
    return DEFAULT_OUT_ROOT / f"{timestamp}-{suffix}"


def run_extract(
    *,
    paper_dir: Path,
    backend: BackendSpec,
    detector: str,
    detector_model_path: str | None,
    detector_model_repo: str,
    detector_confidence: float,
    detector_iou: float,
    detector_padding: float,
    min_confidence: float | None,
    timeout_seconds: int,
    output_name: str,
) -> dict[str, Any]:
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "rag" / "extract_chemical_smiles.py"),
        "--paper",
        paper_dir.name,
        "--limit",
        "1",
        "--backend",
        backend.name,
        "--python-bin",
        backend.python_bin,
        "--detector",
        detector,
        "--detector-model-repo",
        detector_model_repo,
        "--detector-confidence",
        str(detector_confidence),
        "--detector-iou",
        str(detector_iou),
        "--detector-padding",
        str(detector_padding),
        "--timeout-seconds",
        str(timeout_seconds),
        "--output-name",
        output_name,
        "--overwrite",
    ]
    if detector_model_path:
        cmd.extend(["--detector-model-path", detector_model_path])
    if min_confidence is not None:
        cmd.extend(["--min-confidence", str(min_confidence)])

    started = datetime.now(timezone.utc)
    try:
        result = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=max(timeout_seconds + 120, timeout_seconds),
        )
        returncode = result.returncode
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.TimeoutExpired as exc:
        returncode = 124
        stdout = exc.stdout or ""
        stderr = exc.stderr or f"Timed out after {timeout_seconds} seconds"
    completed = datetime.now(timezone.utc)
    output_path = paper_dir / output_name
    payload = load_json(output_path) if output_path.exists() else None
    return {
        "command": cmd,
        "returncode": returncode,
        "stdout": stdout,
        "stderr": stderr,
        "started_at_utc": started.replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "completed_at_utc": completed.replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "duration_seconds": (completed - started).total_seconds(),
        "output_path": str(output_path.relative_to(REPO_ROOT)),
        "payload": payload,
    }


def summarize_payload(payload: dict[str, Any] | None) -> dict[str, Any]:
    if not payload:
        return {
            "pictures_detected_total": 0,
            "pictures_found": 0,
            "pictures_processed": 0,
            "pictures_with_smiles": 0,
            "pictures_skipped_no_detections": 0,
            "valid_smiles_count": 0,
            "invalid_smiles_count": 0,
            "unique_valid_smiles": [],
            "unique_invalid_smiles": [],
            "status_counts": {},
            "confidence_avg": 0.0,
            "confidence_p50": 0.0,
            "confidence_p95": 0.0,
        }

    results = payload.get("results") or payload.get("picture_data") or []
    status_counts = Counter(str(item.get("status") or "unknown") for item in results if isinstance(item, dict))
    confidences = [
        float(item["confidence"])
        for item in results
        if isinstance(item, dict) and isinstance(item.get("confidence"), (int, float))
    ]
    unique_valid = [str(item) for item in payload.get("valid_smiles") or payload.get("smiles_extracted") or [] if item]
    unique_invalid = [str(item) for item in payload.get("invalid_smiles") or [] if item]

    return {
        "pictures_detected_total": int(payload.get("pictures_detected_total") or 0),
        "pictures_found": int(payload.get("pictures_found") or 0),
        "pictures_processed": int(payload.get("pictures_processed") or 0),
        "pictures_with_smiles": int(payload.get("pictures_with_smiles") or 0),
        "pictures_skipped_no_detections": int(payload.get("pictures_skipped_no_detections") or 0),
        "valid_smiles_count": len(unique_valid),
        "invalid_smiles_count": len(unique_invalid),
        "unique_valid_smiles": unique_valid,
        "unique_invalid_smiles": unique_invalid,
        "status_counts": dict(status_counts),
        "confidence_avg": average(confidences),
        "confidence_p50": percentile(confidences, 50),
        "confidence_p95": percentile(confidences, 95),
    }


def set_stats(left: list[str], right: list[str]) -> dict[str, Any]:
    left_set = set(left)
    right_set = set(right)
    union = left_set | right_set
    intersection = left_set & right_set
    return {
        "left_only": sorted(left_set - right_set),
        "right_only": sorted(right_set - left_set),
        "intersection": sorted(intersection),
        "jaccard": (len(intersection) / len(union)) if union else 1.0,
    }


def shorten_list(values: list[str], limit: int = 5) -> list[str]:
    return values[:limit]


def build_disagreements_payload(per_paper: list[dict[str, Any]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for paper_run in per_paper:
        overlap = paper_run["comparison"]["valid_smiles_overlap"]
        rows.append(
            {
                "paper_name": paper_run["paper_name"],
                "paper_dir": paper_run["paper_dir"],
                "molgrapher_output_path": paper_run["backend_runs"]["molgrapher"]["output_path"],
                "molscribe_output_path": paper_run["backend_runs"]["molscribe"]["output_path"],
                "shared_valid_smiles": overlap["intersection"],
                "molgrapher_only_valid_smiles": overlap["left_only"],
                "molscribe_only_valid_smiles": overlap["right_only"],
                "jaccard": overlap["jaccard"],
            }
        )
    return {"per_paper": rows}


def aggregate_backend_summaries(per_paper: list[dict[str, Any]], backend_names: list[str]) -> dict[str, dict[str, Any]]:
    aggregate: dict[str, dict[str, Any]] = {}
    for backend_name in backend_names:
        rows = [item["backend_runs"][backend_name]["summary"] for item in per_paper]
        aggregate[backend_name] = {
            "paper_count": len(rows),
            "pictures_detected_total": sum(int(row["pictures_detected_total"]) for row in rows),
            "pictures_found": sum(int(row["pictures_found"]) for row in rows),
            "pictures_processed": sum(int(row["pictures_processed"]) for row in rows),
            "pictures_with_smiles": sum(int(row["pictures_with_smiles"]) for row in rows),
            "pictures_skipped_no_detections": sum(int(row["pictures_skipped_no_detections"]) for row in rows),
            "valid_smiles_count": sum(int(row["valid_smiles_count"]) for row in rows),
            "invalid_smiles_count": sum(int(row["invalid_smiles_count"]) for row in rows),
            "avg_valid_smiles_per_paper": average([float(row["valid_smiles_count"]) for row in rows]),
            "avg_detector_crops_per_paper": average([float(row["pictures_found"]) for row in rows]),
            "confidence_avg": average([float(row["confidence_avg"]) for row in rows]),
        }
    return aggregate


def render_report(run_payload: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# Chemical OCSR Evaluation Report",
        "",
        f"- Run at: `{run_payload['run_at_utc']}`",
        f"- Detector: `{run_payload['config']['detector']}`",
        f"- Paper filters: `{', '.join(run_payload['paper_filters'])}`",
        f"- Paper count: `{len(run_payload['per_paper'])}`",
        "",
        "## Aggregate Backend Summary",
        "",
        "| Backend | Papers | Source Images | Detector Crops | With SMILES | Unique Valid SMILES | Unique Invalid SMILES | Avg Valid / Paper | Avg Crops / Paper | Avg Confidence |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for backend_name, summary in run_payload["aggregate"]["by_backend"].items():
        lines.append(
            "| {backend} | {paper_count} | {pictures_detected_total} | {pictures_found} | {pictures_with_smiles} | {valid_smiles_count} | {invalid_smiles_count} | {avg_valid_smiles_per_paper:.2f} | {avg_detector_crops_per_paper:.2f} | {confidence_avg:.3f} |".format(
                backend=backend_name,
                **summary,
            )
        )

    lines.extend(
        [
            "",
            "## Per Paper",
            "",
            "| Paper | Backend | Source Images | Detector Crops | With SMILES | Unique Valid | Unique Invalid | P50 Confidence | P95 Confidence |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for paper_run in run_payload["per_paper"]:
        for backend_name, backend_run in paper_run["backend_runs"].items():
            summary = backend_run["summary"]
            lines.append(
                "| {paper} | {backend} | {pictures_detected_total} | {pictures_found} | {pictures_with_smiles} | {valid_smiles_count} | {invalid_smiles_count} | {confidence_p50:.3f} | {confidence_p95:.3f} |".format(
                    paper=paper_run["paper_slug"],
                    backend=backend_name,
                    **summary,
                )
            )

    lines.extend(
        [
            "",
            "## Backend Overlap",
            "",
            "| Paper | Shared Valid SMILES | MolGrapher Only | MolScribe Only | Jaccard |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for paper_run in run_payload["per_paper"]:
        overlap = paper_run["comparison"]["valid_smiles_overlap"]
        lines.append(
            "| {paper} | {shared} | {left_only} | {right_only} | {jaccard:.3f} |".format(
                paper=paper_run["paper_slug"],
                shared=len(overlap["intersection"]),
                left_only=len(overlap["left_only"]),
                right_only=len(overlap["right_only"]),
                jaccard=float(overlap["jaccard"]),
            )
        )

    lines.extend(["", "## Example Disagreements", ""])
    for paper_run in run_payload["per_paper"]:
        overlap = paper_run["comparison"]["valid_smiles_overlap"]
        lines.append(f"### {paper_run['paper_slug']}")
        lines.append("")
        lines.append(
            f"- Shared examples: `{', '.join(shorten_list(overlap['intersection'])) or '-'}`"
        )
        lines.append(
            f"- MolGrapher only: `{', '.join(shorten_list(overlap['left_only'])) or '-'}`"
        )
        lines.append(
            f"- MolScribe only: `{', '.join(shorten_list(overlap['right_only'])) or '-'}`"
        )
        lines.append("")

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- This benchmark compares backend yield on the same MolDetv2 crops. It is not a ground-truth chemical accuracy score.",
            "- `Unique Valid SMILES` counts canonicalized RDKit-valid strings written by the extraction pipeline.",
            "- Disagreements between backends are useful for routing and review, not as an automatic correctness label.",
        ]
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare MolGrapher and MolScribe on the same MolDetv2 detector crops.")
    parser.add_argument(
        "--paper",
        action="append",
        dest="papers",
        help="Substring filter for one extraction directory. Repeat to compare multiple papers.",
    )
    parser.add_argument("--out-dir", help="Artifact output directory. Defaults to artifacts/chemical-ocsr-eval/<timestamp>-<papers>.")
    parser.add_argument(
        "--detector",
        choices=("none", "moldetv2-general", "moldetv2-doc"),
        default=DEFAULT_DETECTOR,
        help=f"Detection stage applied before both backends. Default: {DEFAULT_DETECTOR}",
    )
    parser.add_argument(
        "--detector-model-path",
        default=os.environ.get("MOLDETV2_MODEL_PATH", "").strip() or None,
        help="Optional local MolDetv2 ONNX path. If omitted, extract_chemical_smiles.py resolves from Hugging Face.",
    )
    parser.add_argument(
        "--detector-model-repo",
        default=os.environ.get("MOLDETV2_MODEL_REPO", "UniParser/MolDetv2"),
        help="Hugging Face repo for MolDetv2 weights.",
    )
    parser.add_argument("--detector-confidence", type=float, default=0.5, help="MolDetv2 confidence threshold.")
    parser.add_argument("--detector-iou", type=float, default=0.5, help="MolDetv2 IoU threshold.")
    parser.add_argument("--detector-padding", type=float, default=0.03, help="Relative padding around detector boxes.")
    parser.add_argument("--min-confidence", type=float, default=0.5, help="Backend prediction confidence threshold.")
    parser.add_argument("--timeout-seconds", type=int, default=1800, help="Per-paper backend timeout.")
    parser.add_argument(
        "--output-template",
        default=DEFAULT_OUTPUT_TEMPLATE,
        help=f"Per-backend output filename template. Must contain {{backend}}. Default: {DEFAULT_OUTPUT_TEMPLATE}",
    )
    parser.add_argument("--molgrapher-python-bin", default=None, help="Python interpreter for the MolGrapher backend.")
    parser.add_argument("--molscribe-python-bin", default=None, help="Python interpreter for the MolScribe backend.")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if "{backend}" not in args.output_template:
        raise SystemExit("--output-template must contain {backend}")

    paper_filters = list(args.papers or DEFAULT_PAPER_FILTERS)
    paper_dirs = resolve_paper_dirs(paper_filters)
    backend_specs = resolve_backend_specs(args)
    out_dir = build_out_dir(args.out_dir, paper_filters)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[eval] out_dir={out_dir}")
    print(f"[eval] detector={args.detector}")
    print(f"[eval] papers={len(paper_dirs)}")
    print(f"[eval] backends={','.join(spec.name for spec in backend_specs)}")

    per_paper: list[dict[str, Any]] = []
    failures = 0
    for paper_dir in paper_dirs:
        print(f"[eval] paper={paper_dir.name}")
        backend_runs: dict[str, Any] = {}
        for backend_spec in backend_specs:
            output_name = args.output_template.format(backend=backend_spec.name)
            print(f"  - backend={backend_spec.name}", flush=True)
            run_result = run_extract(
                paper_dir=paper_dir,
                backend=backend_spec,
                detector=args.detector,
                detector_model_path=args.detector_model_path,
                detector_model_repo=args.detector_model_repo,
                detector_confidence=args.detector_confidence,
                detector_iou=args.detector_iou,
                detector_padding=args.detector_padding,
                min_confidence=args.min_confidence,
                timeout_seconds=args.timeout_seconds,
                output_name=output_name,
            )
            run_result["summary"] = summarize_payload(run_result["payload"])
            if run_result["returncode"] != 0:
                failures += 1
                print(f"    returncode={run_result['returncode']}")
            else:
                print(
                    "    crops={pictures_found} valid={valid_smiles_count} invalid={invalid_smiles_count}".format(
                        **run_result["summary"]
                    )
                )
            backend_runs[backend_spec.name] = run_result

        molgrapher_valid = backend_runs["molgrapher"]["summary"]["unique_valid_smiles"]
        molscribe_valid = backend_runs["molscribe"]["summary"]["unique_valid_smiles"]
        comparison = {"valid_smiles_overlap": set_stats(molgrapher_valid, molscribe_valid)}
        per_paper.append(
            {
                "paper_dir": str(paper_dir.relative_to(REPO_ROOT)),
                "paper_name": paper_dir.name,
                "paper_slug": paper_dir.name.split("_-_", 1)[1] if "_-_" in paper_dir.name else paper_dir.name,
                "backend_runs": backend_runs,
                "comparison": comparison,
            }
        )

    aggregate = {
        "by_backend": aggregate_backend_summaries(per_paper, [spec.name for spec in backend_specs]),
        "overlap": {
            "jaccard_avg": average(
                [float(item["comparison"]["valid_smiles_overlap"]["jaccard"]) for item in per_paper]
            ),
            "jaccard_p50": percentile(
                [float(item["comparison"]["valid_smiles_overlap"]["jaccard"]) for item in per_paper],
                50,
            ),
        },
    }

    run_payload = {
        "run_at_utc": utc_now_iso(),
        "paper_filters": paper_filters,
        "config": {
            "detector": args.detector,
            "detector_model_path": args.detector_model_path,
            "detector_model_repo": args.detector_model_repo,
            "detector_confidence": args.detector_confidence,
            "detector_iou": args.detector_iou,
            "detector_padding": args.detector_padding,
            "min_confidence": args.min_confidence,
            "timeout_seconds": args.timeout_seconds,
            "output_template": args.output_template,
            "molgrapher_python_bin": backend_specs[0].python_bin,
            "molscribe_python_bin": backend_specs[1].python_bin,
        },
        "per_paper": per_paper,
        "aggregate": aggregate,
        "failures": failures,
    }
    write_json(out_dir / "run.json", run_payload)
    write_json(out_dir / "disagreements.json", build_disagreements_payload(per_paper))
    render_report(run_payload, out_dir / "report.md")
    print(f"[eval] report={out_dir / 'report.md'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
