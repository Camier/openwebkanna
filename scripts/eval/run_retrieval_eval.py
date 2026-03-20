#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


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


def parse_k_values(raw: str) -> list[int]:
    values: list[int] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        values.append(int(item))
    if not values:
        raise ValueError("At least one k value is required")
    return sorted(set(values))


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def normalize_text(value: str) -> str:
    lowered = value.casefold()
    lowered = re.sub(r"[^0-9a-z]+", " ", lowered)
    return " ".join(lowered.split())


def register_alias(alias_map: dict[str, set[str]], alias: str, canonical_id: str) -> None:
    normalized = normalize_text(alias)
    if not normalized:
        return
    alias_map.setdefault(normalized, set()).add(canonical_id)


def build_paper_id_aliases(repo_root: Path) -> dict[str, set[str]]:
    alias_map: dict[str, set[str]] = {}
    sources = [
        repo_root / "data/corpus/biblio_corpus.curated.jsonl",
        repo_root / "data/metadata/paper_catalog.jsonl",
    ]
    for path in sources:
        for row in load_jsonl(path):
            canonical_id = str(row.get("paper_id") or "").strip()
            if not canonical_id:
                continue
            register_alias(alias_map, canonical_id, canonical_id)
            for key in ("title", "citekey", "doi", "doc_id"):
                value = str(row.get(key) or "").strip()
                if value:
                    register_alias(alias_map, value, canonical_id)
    return alias_map


def resolve_alias_ids(
    values: list[str], alias_map: dict[str, set[str]], *, include_unresolved: bool = False
) -> set[str]:
    resolved: set[str] = set()
    for value in values:
        normalized = normalize_text(value)
        if not normalized:
            continue
        alias_matches = alias_map.get(normalized, set())
        if alias_matches:
            resolved.update(alias_matches)
        elif include_unresolved:
            resolved.add(value)
    return resolved


def request_json(url: str, *, token: str, method: str = "GET", payload: Any | None = None) -> Any:
    body = None
    headers: dict[str, str] = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    if payload is not None:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=body, headers=headers, method=method)
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def normalize_knowledge_items(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        items = payload.get("items")
        if isinstance(items, list):
            return items
    return []


def extract_hits(payload: dict[str, Any]) -> list[dict[str, Any]]:
    ids = payload.get("ids") or []
    docs = payload.get("documents") or []
    metas = payload.get("metadatas") or []
    dists = payload.get("distances") or []

    id_row = ids[0] if ids and isinstance(ids[0], list) else []
    doc_row = docs[0] if docs and isinstance(docs[0], list) else []
    meta_row = metas[0] if metas and isinstance(metas[0], list) else []
    dist_row = dists[0] if dists and isinstance(dists[0], list) else []

    count = max(len(id_row), len(doc_row), len(meta_row), len(dist_row))
    hits: list[dict[str, Any]] = []
    for idx in range(count):
        hits.append(
            {
                "id": id_row[idx] if idx < len(id_row) else None,
                "document": doc_row[idx] if idx < len(doc_row) else "",
                "metadata": meta_row[idx] if idx < len(meta_row) else {},
                "distance": dist_row[idx] if idx < len(dist_row) else None,
                "rank": idx + 1,
            }
        )
    return hits


def lower_list(values: list[str]) -> list[str]:
    return [value.lower() for value in values if value]


def hit_text_fields(hit: dict[str, Any]) -> list[str]:
    metadata = hit.get("metadata") or {}
    return [
        str(metadata.get("title") or ""),
        str(metadata.get("name") or ""),
        str(metadata.get("source") or ""),
        str(metadata.get("paper_id") or ""),
        str(metadata.get("doi") or ""),
        str(metadata.get("citekey") or ""),
        str(metadata.get("doc_id") or ""),
        str(hit.get("document") or ""),
    ]


def hit_alias_ids(hit: dict[str, Any], alias_map: dict[str, set[str]]) -> set[str]:
    return resolve_alias_ids(hit_text_fields(hit)[:7], alias_map)


def match_hit(hit: dict[str, Any], query_spec: dict[str, Any], ctx: "EvalContext") -> tuple[bool, list[str]]:
    fields = hit_text_fields(hit)
    haystack_titles = " ".join(field for field in fields[:7] if field).lower()
    haystack_document = fields[-1].lower()
    haystack_full = " ".join(field for field in fields if field).lower()

    title_needles = lower_list(query_spec.get("expected_title_substrings") or [])
    source_needles = lower_list(query_spec.get("expected_source_substrings") or [])
    text_needles = lower_list(query_spec.get("expected_text_substrings") or [])
    required_needles = lower_list(query_spec.get("must_include") or [])
    expected_ids = [str(item) for item in (query_spec.get("expected_paper_ids") or []) if item]
    resolved_expected_ids = resolve_alias_ids(expected_ids, ctx.paper_id_aliases, include_unresolved=True)
    matched_ids: set[str] = set()
    is_match = False

    if title_needles and any(needle in haystack_titles for needle in title_needles):
        is_match = True
    if not is_match and source_needles and any(needle in haystack_titles for needle in source_needles):
        is_match = True
    if not is_match and text_needles and any(needle in haystack_document for needle in text_needles):
        is_match = True
    if resolved_expected_ids:
        matched_ids = hit_alias_ids(hit, ctx.paper_id_aliases)
        if matched_ids.intersection(resolved_expected_ids):
            is_match = True
    if is_match and required_needles and not all(needle in haystack_full for needle in required_needles):
        return False, sorted(matched_ids)
    return is_match, sorted(matched_ids)


def reciprocal_rank(rank: int | None, k: int) -> float:
    if rank is None or rank > k or rank <= 0:
        return 0.0
    return 1.0 / float(rank)


@dataclass
class EvalContext:
    dataset_path: Path
    out_dir: Path
    url: str
    backend: str
    token: str | None
    k_values: list[int]
    baseline_path: Path | None
    max_recall_drop: float
    max_mrr_drop: float
    enforce_regression_gate: bool
    collection_name_override: str | None
    knowledge_name_override: str | None
    paper_id_aliases: dict[str, set[str]]


def resolve_collection_name(ctx: EvalContext, dataset: dict[str, Any]) -> tuple[str | None, str | None]:
    if ctx.backend == "canonical":
        return None, dataset.get("knowledge_name")

    if ctx.collection_name_override:
        return ctx.collection_name_override, dataset.get("knowledge_name")

    if dataset.get("collection_name"):
        return str(dataset["collection_name"]), dataset.get("knowledge_name")

    knowledge_name = ctx.knowledge_name_override or dataset.get("knowledge_name")
    if not knowledge_name:
        raise SystemExit("Dataset must define knowledge_name or collection_name, or pass an override")

    payload = request_json(f"{ctx.url.rstrip('/')}/api/v1/knowledge/", token=ctx.token)
    for item in normalize_knowledge_items(payload):
        if str(item.get("name")) == str(knowledge_name):
            return str(item["id"]), str(knowledge_name)

    raise SystemExit(f"Knowledge base not found: {knowledge_name}")


def canonical_hit_to_legacy_shape(hit: dict[str, Any]) -> dict[str, Any]:
    payload = hit.get("payload") or {}
    title = hit.get("title") or payload.get("document_title") or payload.get("title")
    document = (
        payload.get("chunk_text")
        or payload.get("text")
        or payload.get("content")
        or payload.get("caption_text")
        or ""
    )
    metadata = {
        "title": title,
        "name": title,
        "source": payload.get("source") or payload.get("document_title") or title,
        "paper_id": hit.get("doc_id") or payload.get("document_id") or payload.get("doc_id"),
        "doi": payload.get("doi"),
        "citekey": hit.get("citation_key") or payload.get("citation_key"),
        "doc_id": hit.get("doc_id") or payload.get("document_id") or payload.get("doc_id"),
        "page": hit.get("page_number") or payload.get("page_number"),
    }
    return {
        "id": hit.get("point_id"),
        "document": str(document or ""),
        "metadata": metadata,
        "distance": None,
        "rank": int(hit.get("rank") or 0) or None,
    }


def canonical_extract_hits(payload: dict[str, Any]) -> list[dict[str, Any]]:
    hits = payload.get("reranked_hits")
    if not isinstance(hits, list):
        return []
    return [canonical_hit_to_legacy_shape(hit) for hit in hits if isinstance(hit, dict)]


def evaluate_query(
    ctx: EvalContext,
    collection_name: str | None,
    query_spec: dict[str, Any],
) -> dict[str, Any]:
    query = str(query_spec["query"])
    modality = str(query_spec.get("modality") or "text")
    max_k = max(ctx.k_values)
    if ctx.backend == "canonical":
        payload = {"query": query, "top_k": max_k}
    else:
        payload = {"collection_name": collection_name, "query": query, "k": max_k}

    started = time.perf_counter()
    try:
        if ctx.backend == "canonical":
            response = request_json(
                f"{ctx.url.rstrip('/')}/api/v1/retrieve",
                token="",
                method="POST",
                payload=payload,
            )
        else:
            response = request_json(
                f"{ctx.url.rstrip('/')}/api/v1/retrieval/query/doc",
                token=ctx.token or "",
                method="POST",
                payload=payload,
            )
        error = None
    except urllib.error.HTTPError as exc:
        response = {"error": f"HTTP {exc.code}"}
        error = f"HTTP {exc.code}"
    except Exception as exc:  # noqa: BLE001
        response = {"error": str(exc)}
        error = str(exc)
    latency_ms = (time.perf_counter() - started) * 1000.0

    if error is None:
        hits = canonical_extract_hits(response) if ctx.backend == "canonical" else extract_hits(response)
    else:
        hits = []
    relevant_rank = None
    matched_title = None
    matched_paper_ids: list[str] = []
    for hit in hits:
        is_match, candidate_paper_ids = match_hit(hit, query_spec, ctx)
        if is_match:
            relevant_rank = hit["rank"]
            metadata = hit.get("metadata") or {}
            matched_title = metadata.get("title") or metadata.get("name") or metadata.get("source")
            matched_paper_ids = candidate_paper_ids
            break

    metrics: dict[str, float] = {}
    metrics_by_k: list[dict[str, Any]] = []
    for k in ctx.k_values:
        recall = 1.0 if relevant_rank is not None and relevant_rank <= k else 0.0
        mrr = reciprocal_rank(relevant_rank, k)
        metrics[f"pre_recall@{k}"] = recall
        metrics[f"post_recall@{k}"] = recall
        metrics[f"pre_mrr@{k}"] = mrr
        metrics[f"post_mrr@{k}"] = mrr
        metrics_by_k.append(
            {
                "k": k,
                "pre_recall": recall,
                "post_recall": recall,
                "pre_mrr": mrr,
                "post_mrr": mrr,
                "relevant_found_pre": recall > 0.0,
                "relevant_found_post": recall > 0.0,
                "relevant_rank_pre": relevant_rank if recall > 0.0 else None,
                "relevant_rank_post": relevant_rank if recall > 0.0 else None,
            }
        )

    return {
        "query_id": query_spec["query_id"],
        "query_text": query,
        "modality": modality,
        "latency_ms": latency_ms,
        "retrieved_count": len(hits),
        "relevant_rank": relevant_rank,
        "matched_title": matched_title,
        "matched_paper_ids": matched_paper_ids,
        "error": error,
        "metrics": metrics,
        "metrics_by_k": metrics_by_k,
        "top_hits": [
            {
                "rank": hit["rank"],
                "title": (hit.get("metadata") or {}).get("title")
                or (hit.get("metadata") or {}).get("name")
                or (hit.get("metadata") or {}).get("source"),
                "page": (hit.get("metadata") or {}).get("page"),
                "distance": hit.get("distance"),
            }
            for hit in hits[:5]
        ],
    }


def build_summary_rows(per_query: list[dict[str, Any]], k_values: list[int]) -> list[dict[str, Any]]:
    modalities = sorted({str(item.get("modality") or "text") for item in per_query})
    rows: list[dict[str, Any]] = []
    for modality in modalities:
        modality_queries = [item for item in per_query if str(item.get("modality") or "text") == modality]
        latencies = [float(item["latency_ms"]) for item in modality_queries]
        for k in k_values:
            pre_recalls = [float(item["metrics"][f"pre_recall@{k}"]) for item in modality_queries]
            post_recalls = [float(item["metrics"][f"post_recall@{k}"]) for item in modality_queries]
            pre_mrrs = [float(item["metrics"][f"pre_mrr@{k}"]) for item in modality_queries]
            post_mrrs = [float(item["metrics"][f"post_mrr@{k}"]) for item in modality_queries]
            rows.append(
                {
                    "modality": modality,
                    "k": k,
                    "query_count": len(modality_queries),
                    "pre_recall": average(pre_recalls),
                    "post_recall": average(post_recalls),
                    "pre_mrr": average(pre_mrrs),
                    "post_mrr": average(post_mrrs),
                    "p95_latency_ms": percentile(latencies, 95),
                    "p50_latency_ms": percentile(latencies, 50),
                }
            )
    return rows


def compare_to_baseline(
    summary_rows: list[dict[str, Any]],
    baseline_payload: dict[str, Any] | None,
    max_recall_drop: float,
    max_mrr_drop: float,
) -> tuple[str | None, dict[str, Any] | None]:
    if baseline_payload is None:
        return None, None

    baseline_rows = {
        (str(row["modality"]), int(row["k"])): row for row in baseline_payload.get("summary_rows", [])
    }
    comparisons: list[dict[str, Any]] = []
    regressions: list[dict[str, Any]] = []
    recall_drops: list[float] = []
    mrr_drops: list[float] = []

    for row in summary_rows:
        key = (str(row["modality"]), int(row["k"]))
        baseline = baseline_rows.get(key)
        current_recall = float(row["post_recall"])
        current_mrr = float(row["post_mrr"])
        if baseline is None:
            comparison = {
                "modality": row["modality"],
                "k": row["k"],
                "baseline_found": False,
                "baseline_post_recall": None,
                "current_post_recall": current_recall,
                "delta_post_recall": None,
                "baseline_post_mrr": None,
                "current_post_mrr": current_mrr,
                "delta_post_mrr": None,
                "recall_regression": False,
                "mrr_regression": False,
                "is_regression": False,
            }
        else:
            baseline_recall = float(baseline["post_recall"])
            baseline_mrr = float(baseline["post_mrr"])
            delta_recall = current_recall - baseline_recall
            delta_mrr = current_mrr - baseline_mrr
            recall_regression = delta_recall < -max_recall_drop
            mrr_regression = delta_mrr < -max_mrr_drop
            comparison = {
                "modality": row["modality"],
                "k": row["k"],
                "baseline_found": True,
                "baseline_post_recall": baseline_recall,
                "current_post_recall": current_recall,
                "delta_post_recall": delta_recall,
                "baseline_post_mrr": baseline_mrr,
                "current_post_mrr": current_mrr,
                "delta_post_mrr": delta_mrr,
                "recall_regression": recall_regression,
                "mrr_regression": mrr_regression,
                "is_regression": recall_regression or mrr_regression,
            }
            recall_drops.append(max(0.0, baseline_recall - current_recall))
            mrr_drops.append(max(0.0, baseline_mrr - current_mrr))
            if comparison["is_regression"]:
                regressions.append(comparison)
        comparisons.append(comparison)

    return str(baseline_payload.get("_path")) if baseline_payload.get("_path") else None, {
        "max_recall_drop": max(recall_drops) if recall_drops else 0.0,
        "max_mrr_drop": max(mrr_drops) if mrr_drops else 0.0,
        "comparisons": comparisons,
        "regressions": regressions,
    }


def render_report(
    dataset: dict[str, Any],
    dataset_path: Path,
    run_payload: dict[str, Any],
    out_path: Path,
    knowledge_name: str | None,
    collection_name: str | None,
) -> None:
    query_count = len(run_payload["per_query"])
    k_values = ", ".join(str(k) for k in run_payload["k_values"])
    lines = [
        "# Retrieval Evaluation Report",
        "",
        f"- Dataset version: `{run_payload['dataset_version']}`",
        f"- Source file: `{dataset_path}`",
        f"- Run at: `{run_payload['run_at_utc']}`",
        f"- K values: `{k_values}`",
        f"- Query count: `{query_count}`",
        f"- Backend: `{run_payload['backend']}`",
        f"- Knowledge base: `{knowledge_name or 'n/a'}`",
        f"- Collection name: `{collection_name or 'n/a'}`",
        "",
        "## Summary",
        "",
        "| Modality | k | Query Count | Pre Recall | Post Recall | Pre MRR | Post MRR | P50 Latency (ms) | P95 Latency (ms) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in run_payload["summary_rows"]:
        lines.append(
            "| {modality} | {k} | {query_count} | {pre_recall:.3f} | {post_recall:.3f} | {pre_mrr:.3f} | {post_mrr:.3f} | {p50_latency_ms:.2f} | {p95_latency_ms:.2f} |".format(
                **row
            )
        )

    lines.extend(["", "## Per Query", "", "| Query ID | Modality | Relevant Rank | Latency (ms) | Matched Title | Error |", "| --- | --- | ---: | ---: | --- | --- |"])
    for item in run_payload["per_query"]:
        rank = item["relevant_rank"] if item["relevant_rank"] is not None else "-"
        matched = item.get("matched_title") or "-"
        error = item.get("error") or "-"
        lines.append(
            f"| {item['query_id']} | {item['modality']} | {rank} | {item['latency_ms']:.2f} | {matched} | {error} |"
        )

    baseline_comparison = run_payload.get("baseline_comparison")
    if baseline_comparison:
        lines.extend(
            [
                "",
                "## Baseline Comparison",
                "",
                f"- Baseline path: `{run_payload.get('baseline_path')}`",
                f"- Max recall drop: `{baseline_comparison['max_recall_drop']:.3f}`",
                f"- Max MRR drop: `{baseline_comparison['max_mrr_drop']:.3f}`",
                "",
                "| Modality | k | Baseline Found | Baseline Post Recall | Current Post Recall | Delta Recall | Baseline Post MRR | Current Post MRR | Delta MRR | Regression |",
                "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
            ]
        )
        for item in baseline_comparison["comparisons"]:
            lines.append(
                "| {modality} | {k} | {baseline_found} | {baseline_post_recall} | {current_post_recall} | {delta_post_recall} | {baseline_post_mrr} | {current_post_mrr} | {delta_post_mrr} | {is_regression} |".format(
                    **item
                )
            )

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_results_csv(run_payload: dict[str, Any], out_path: Path) -> None:
    fieldnames = [
        "dataset_version",
        "run_at_utc",
        "modality",
        "k",
        "query_count",
        "pre_recall",
        "post_recall",
        "delta_recall",
        "pre_mrr",
        "post_mrr",
        "delta_mrr",
        "p50_latency_ms",
        "p95_latency_ms",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in run_payload["summary_rows"]:
            writer.writerow(
                {
                    "dataset_version": run_payload["dataset_version"],
                    "run_at_utc": run_payload["run_at_utc"],
                    "modality": row["modality"],
                    "k": row["k"],
                    "query_count": row["query_count"],
                    "pre_recall": row["pre_recall"],
                    "post_recall": row["post_recall"],
                    "delta_recall": row["post_recall"] - row["pre_recall"],
                    "pre_mrr": row["pre_mrr"],
                    "post_mrr": row["post_mrr"],
                    "delta_mrr": row["post_mrr"] - row["pre_mrr"],
                    "p50_latency_ms": row["p50_latency_ms"],
                    "p95_latency_ms": row["p95_latency_ms"],
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a retrieval evaluation against the canonical multimodal service or legacy OpenWebUI retrieval."
    )
    parser.add_argument(
        "--dataset",
        default="data/eval/retrieval/sceletium-retrieval-sample-v1.json",
        help="Path to the retrieval evaluation dataset JSON.",
    )
    parser.add_argument(
        "--out-dir",
        default="artifacts/retrieval-evaluation-harness",
        help="Directory where run.json, report.md, and results.csv will be written.",
    )
    parser.add_argument("--k-values", default="1,5,10", help="Comma-separated k values.")
    parser.add_argument("--baseline", default=None, help="Optional baseline run.json to compare against.")
    parser.add_argument("--max-recall-drop", type=float, default=0.0, help="Allowed recall drop before gating fails.")
    parser.add_argument("--max-mrr-drop", type=float, default=0.0, help="Allowed MRR drop before gating fails.")
    parser.add_argument(
        "--enforce-regression-gate",
        action="store_true",
        help="Exit non-zero when the baseline comparison exceeds allowed drops.",
    )
    parser.add_argument("--collection-name", default=None, help="Optional OpenWebUI collection/knowledge UUID override.")
    parser.add_argument("--knowledge-name", default=None, help="Optional knowledge base name override.")
    parser.add_argument("--openwebui-url", default=os.environ.get("OPENWEBUI_URL", "http://localhost:3000"))
    parser.add_argument(
        "--multimodal-retrieval-api-url",
        default=os.environ.get("MULTIMODAL_RETRIEVAL_API_URL", "http://127.0.0.1:8510"),
    )
    parser.add_argument(
        "--backend",
        choices=("canonical", "legacy-openwebui"),
        default=os.environ.get("RETRIEVAL_EVAL_BACKEND", "canonical"),
        help="Retrieval backend to evaluate.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    token = os.environ.get("OPENWEBUI_EVAL_TOKEN") or os.environ.get("OPENWEBUI_API_KEY")
    if args.backend == "legacy-openwebui" and not token:
        print("OPENWEBUI_EVAL_TOKEN or OPENWEBUI_API_KEY is required for legacy-openwebui backend", file=sys.stderr)
        return 1

    dataset_path = Path(args.dataset)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ctx = EvalContext(
        dataset_path=dataset_path,
        out_dir=out_dir,
        url=args.multimodal_retrieval_api_url if args.backend == "canonical" else args.openwebui_url,
        backend=args.backend,
        token=token,
        k_values=parse_k_values(args.k_values),
        baseline_path=Path(args.baseline) if args.baseline else None,
        max_recall_drop=args.max_recall_drop,
        max_mrr_drop=args.max_mrr_drop,
        enforce_regression_gate=args.enforce_regression_gate,
        collection_name_override=args.collection_name,
        knowledge_name_override=args.knowledge_name,
        paper_id_aliases=build_paper_id_aliases(Path(__file__).resolve().parents[2]),
    )

    dataset = load_json(dataset_path)
    collection_name, knowledge_name = resolve_collection_name(ctx, dataset)

    per_query = [evaluate_query(ctx, collection_name, query_spec) for query_spec in dataset["queries"]]
    summary_rows = build_summary_rows(per_query, ctx.k_values)

    baseline_payload = None
    if ctx.baseline_path:
        baseline_payload = load_json(ctx.baseline_path)
        if isinstance(baseline_payload, dict):
            baseline_payload["_path"] = str(ctx.baseline_path)

    baseline_path, baseline_comparison = compare_to_baseline(
        summary_rows, baseline_payload, ctx.max_recall_drop, ctx.max_mrr_drop
    )
    regression_failed = bool(baseline_comparison and baseline_comparison.get("regressions"))
    run_payload: dict[str, Any] = {
        "schema_version": "retrieval-eval-run/v1",
        "backend": args.backend,
        "service_url": ctx.url,
        "dataset_version": dataset["dataset_version"],
        "dataset_path": str(dataset_path),
        "run_at_utc": utc_now_iso(),
        "k_values": ctx.k_values,
        "query_count": len(per_query),
        "rrf_k": int(dataset.get("rrf_k", 60)),
        "weights": dataset.get("weights", {}),
        "enable_smiles_embedding": bool(dataset.get("enable_smiles_embedding", False)),
        "knowledge_name": knowledge_name,
        "collection_name": collection_name,
        "summary_rows": summary_rows,
        "per_query": per_query,
        "baseline_path": baseline_path,
        "baseline_comparison": baseline_comparison,
        "regression_gate": {
            "enabled": bool(ctx.enforce_regression_gate and baseline_payload is not None),
            "failed": bool(ctx.enforce_regression_gate and regression_failed),
            "thresholds": {
                "max_recall_drop": ctx.max_recall_drop,
                "max_mrr_drop": ctx.max_mrr_drop,
            },
        },
    }

    write_json(out_dir / "run.json", run_payload)
    write_results_csv(run_payload, out_dir / "results.csv")
    render_report(dataset, dataset_path, run_payload, out_dir / "report.md", knowledge_name, collection_name)

    if ctx.enforce_regression_gate and regression_failed:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
