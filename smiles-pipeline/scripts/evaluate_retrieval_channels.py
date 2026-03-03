#!/usr/bin/env python3
"""
Evaluate retrieval channels separately for text and SMILES queries.

Input JSON schema:
{
  "queries": [
    {
      "query_id": "q1",
      "modality": "text" | "smiles",
      "relevant_doc_ids": ["doc-1", "doc-9"],
      "channels": {
        "text_dense": [{"doc_id": "doc-1", "score": 0.95}],
        "bm25": [{"doc_id": "doc-2", "score": 12.1}],
        "smiles_fingerprint": [{"doc_id": "doc-9", "score": 0.88}],
        "smiles_embedding": [{"doc_id": "doc-7", "score": 0.66}]
      }
    }
  ]
}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from retrieval.rrf_fusion import compute_recall_at_k, fuse_rankings_rrf, normalize_ranked_docs


def parse_k_values(raw: str) -> List[int]:
    values = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        values.append(int(part))
    if not values:
        raise ValueError("At least one k value is required")
    return sorted(set(values))


def parse_weights(raw: str) -> Dict[str, float]:
    if not raw.strip():
        return {}
    result: Dict[str, float] = {}
    for part in raw.split(","):
        channel, _, value = part.partition("=")
        if not channel or not value:
            raise ValueError(f"Invalid weight fragment: {part!r}")
        result[channel.strip()] = float(value.strip())
    return result


def ranked_ids(channel_hits: Iterable[Dict[str, Any]]) -> List[str]:
    return [row["doc_id"] for row in normalize_ranked_docs(list(channel_hits))]


def compute_mean(values: List[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def evaluate_queries(
    queries: List[Dict[str, Any]],
    k_values: List[int],
    rrf_k: int,
    weights: Dict[str, float],
    enable_smiles_embedding: bool,
) -> Dict[str, Any]:
    by_modality = {
        "text": {"pre": {k: [] for k in k_values}, "post": {k: [] for k in k_values}},
        "smiles": {"pre": {k: [] for k in k_values}, "post": {k: [] for k in k_values}},
    }
    per_query: List[Dict[str, Any]] = []

    for query in queries:
        modality = str(query.get("modality", "text")).strip().lower()
        if modality not in by_modality:
            continue
        channels = query.get("channels") or {}
        relevant = [str(doc_id) for doc_id in (query.get("relevant_doc_ids") or [])]

        text_dense_ids = ranked_ids(channels.get("text_dense", []))
        bm25_ids = ranked_ids(channels.get("bm25", []))
        smiles_fp_ids = ranked_ids(channels.get("smiles_fingerprint", []))
        smiles_emb_ids = ranked_ids(channels.get("smiles_embedding", []))

        if modality == "text":
            pre_ids = text_dense_ids
        else:
            pre_ids = smiles_fp_ids

        channel_rankings = {
            "text_dense": channels.get("text_dense", []),
            "bm25": channels.get("bm25", []),
            "smiles_fingerprint": channels.get("smiles_fingerprint", []),
        }
        if enable_smiles_embedding:
            channel_rankings["smiles_embedding"] = channels.get("smiles_embedding", [])

        fused = fuse_rankings_rrf(
            channel_rankings,
            top_k=max(k_values),
            rrf_k=rrf_k,
            channel_weights=weights,
        )
        post_ids = [row["doc_id"] for row in fused]

        query_metrics: Dict[str, Any] = {
            "query_id": query.get("query_id"),
            "modality": modality,
            "recall": {},
        }

        for k in k_values:
            pre_recall = compute_recall_at_k(pre_ids, relevant, k)
            post_recall = compute_recall_at_k(post_ids, relevant, k)
            by_modality[modality]["pre"][k].append(pre_recall)
            by_modality[modality]["post"][k].append(post_recall)
            query_metrics["recall"][f"pre@{k}"] = pre_recall
            query_metrics["recall"][f"post@{k}"] = post_recall

        query_metrics["channels"] = {
            "text_dense_top": text_dense_ids[:5],
            "bm25_top": bm25_ids[:5],
            "smiles_fingerprint_top": smiles_fp_ids[:5],
            "smiles_embedding_top": smiles_emb_ids[:5],
            "fused_top": post_ids[:5],
        }
        per_query.append(query_metrics)

    summary: Dict[str, Any] = {"text": {}, "smiles": {}}
    for modality in ("text", "smiles"):
        for k in k_values:
            summary[modality][f"Recall@{k}"] = {
                "pre_fusion": compute_mean(by_modality[modality]["pre"][k]),
                "post_fusion": compute_mean(by_modality[modality]["post"][k]),
            }

    return {
        "config": {
            "k_values": k_values,
            "rrf_k": rrf_k,
            "weights": weights,
            "enable_smiles_embedding": enable_smiles_embedding,
        },
        "summary": summary,
        "per_query": per_query,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate Recall@k for text and SMILES channels before/after RRF fusion."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to retrieval evaluation JSON input",
    )
    parser.add_argument(
        "--k",
        default="5,10",
        help="Comma-separated Recall@k values (default: 5,10)",
    )
    parser.add_argument(
        "--rrf-k",
        type=int,
        default=60,
        help="RRF rank constant (default: 60)",
    )
    parser.add_argument(
        "--weights",
        default="",
        help="Optional channel weights, e.g. text_dense=1.0,bm25=0.8,smiles_fingerprint=1.2",
    )
    parser.add_argument(
        "--enable-smiles-embedding",
        action="store_true",
        help="Include smiles_embedding channel in fused ranking",
    )
    parser.add_argument(
        "--output",
        help="Optional output path for JSON report",
    )
    args = parser.parse_args()

    k_values = parse_k_values(args.k)
    weights = parse_weights(args.weights)

    input_path = Path(args.input)
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    queries = payload.get("queries", [])
    if not isinstance(queries, list):
        raise ValueError("Input file must contain a list at key 'queries'")

    report = evaluate_queries(
        queries=queries,
        k_values=k_values,
        rrf_k=args.rrf_k,
        weights=weights,
        enable_smiles_embedding=args.enable_smiles_embedding,
    )

    rendered = json.dumps(report, indent=2, ensure_ascii=True)
    if args.output:
        Path(args.output).write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
