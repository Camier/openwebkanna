#!/bin/bash

###############################################################################
# audit-data-quality.sh - Deterministic corpus data quality audit
#
# Usage:
#   ./scripts/audit-data-quality.sh [--output-dir <path>] [--chunks-file <path>] [--biblio-file <path>]
#
# Outputs:
#   <output-dir>/data_quality_audit.json
#   <output-dir>/data_quality_audit.md
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

OUTPUT_DIR="${PROJECT_ROOT}/artifacts/data-quality"
CHUNKS_FILE="${PROJECT_ROOT}/data/corpus/chunks_corpus.jsonl"
BIBLIO_FILE="${PROJECT_ROOT}/data/corpus/biblio_corpus.jsonl"
CHUNKS_FILE_SET=0
BIBLIO_FILE_SET=0

if [[ -n ${CORPUS_CHUNKS_FILE:-} ]]; then
    CHUNKS_FILE="${CORPUS_CHUNKS_FILE}"
    CHUNKS_FILE_SET=1
fi
if [[ -n ${CORPUS_BIBLIO_FILE:-} ]]; then
    BIBLIO_FILE="${CORPUS_BIBLIO_FILE}"
    BIBLIO_FILE_SET=1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --chunks-file)
            CHUNKS_FILE="$2"
            CHUNKS_FILE_SET=1
            shift 2
            ;;
        --biblio-file)
            BIBLIO_FILE="$2"
            BIBLIO_FILE_SET=1
            shift 2
            ;;
        --help | -h)
            echo "Usage: $0 [--output-dir <path>] [--chunks-file <path>] [--biblio-file <path>]"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

if [[ ${CHUNKS_FILE_SET} -eq 0 && -f "${PROJECT_ROOT}/data/corpus/chunks_corpus.curated.jsonl" ]]; then
    CHUNKS_FILE="${PROJECT_ROOT}/data/corpus/chunks_corpus.curated.jsonl"
fi
if [[ ${BIBLIO_FILE_SET} -eq 0 && -f "${PROJECT_ROOT}/data/corpus/biblio_corpus.curated.jsonl" ]]; then
    BIBLIO_FILE="${PROJECT_ROOT}/data/corpus/biblio_corpus.curated.jsonl"
fi

mkdir -p "${OUTPUT_DIR}"

python - "${PROJECT_ROOT}" "${OUTPUT_DIR}" "${CHUNKS_FILE}" "${BIBLIO_FILE}" <<'PY'
import json
import re
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import sys

project_root = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
chunks_file_arg = Path(sys.argv[3])
biblio_file_arg = Path(sys.argv[4])

biblio_path = biblio_file_arg if biblio_file_arg.is_absolute() else (project_root / biblio_file_arg)
chunks_path = chunks_file_arg if chunks_file_arg.is_absolute() else (project_root / chunks_file_arg)
catalog_path = project_root / "data/metadata/paper_catalog.jsonl"
pdfs_dir = project_root / "data/pdfs"


def load_jsonl(path: Path):
    rows = []
    invalid = 0
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                invalid += 1
    return rows, invalid


biblio, biblio_invalid = load_jsonl(biblio_path)
chunks, chunks_invalid = load_jsonl(chunks_path)
catalog, catalog_invalid = load_jsonl(catalog_path)

biblio_doc_ids = {r.get("doc_id") for r in biblio if r.get("doc_id")}
biblio_doc_to_paper_id = {
    r.get("doc_id"): r.get("paper_id")
    for r in biblio
    if r.get("doc_id")
}
biblio_paper_ids = {r.get("paper_id") for r in biblio if r.get("paper_id")}
catalog_paper_ids = {r.get("paper_id") for r in catalog if r.get("paper_id")}

chunk_doc_ids = Counter()
chunk_key_counter = Counter()
missing_paper_id = 0
missing_source_page = 0
empty_text = 0
text_lengths = []
normalized_texts = Counter()
doc_total_rows = Counter()
doc_short_rows = Counter()
doc_long_rows = Counter()
doc_normalized_texts = Counter()

for r in chunks:
    doc_id = r.get("doc_id")
    chunk_id = r.get("chunk_id")
    if doc_id:
        chunk_doc_ids[doc_id] += 1
        doc_total_rows[doc_id] += 1
    if doc_id and chunk_id:
        chunk_key_counter[(doc_id, chunk_id)] += 1

    if not r.get("paper_id"):
        missing_paper_id += 1
    if r.get("source_page") in (None, ""):
        missing_source_page += 1

    text = r.get("text", "")
    if not isinstance(text, str):
        text = str(text)
    if not text.strip():
        empty_text += 1
    l = len(text)
    text_lengths.append(l)
    if doc_id and l < 200:
        doc_short_rows[doc_id] += 1
    if doc_id and l > 2000:
        doc_long_rows[doc_id] += 1

    norm = re.sub(r"\s+", " ", text).strip().lower()
    if norm:
        normalized_texts[norm] += 1
        if doc_id:
            doc_normalized_texts[(doc_id, norm)] += 1

unknown_doc_rows = sum(1 for r in chunks if r.get("doc_id") not in biblio_doc_ids)
duplicate_chunk_key_rows = sum(c - 1 for c in chunk_key_counter.values() if c > 1)
paper_ids_in_biblio_not_catalog = sorted(list(biblio_paper_ids - catalog_paper_ids))
paper_ids_in_catalog_not_biblio = sorted(list(catalog_paper_ids - biblio_paper_ids))

chunk_count_mismatches = 0
for r in biblio:
    doc_id = r.get("doc_id")
    declared = r.get("chunk_count")
    observed = chunk_doc_ids.get(doc_id, 0)
    if isinstance(declared, int) and declared != observed:
        chunk_count_mismatches += 1

missing_extraction_path = 0
for r in biblio:
    p = r.get("extraction_path")
    if p and not (project_root / p).exists():
        missing_extraction_path += 1

pdf_non_pdf_entries = []
if pdfs_dir.exists():
    for p in sorted(pdfs_dir.iterdir()):
        if p.is_file() and p.suffix.lower() != ".pdf":
            pdf_non_pdf_entries.append(p.name)

total_chunks = len(chunks)
total_norm_texts = sum(normalized_texts.values())
duplicate_text_rows = sum(c - 1 for c in normalized_texts.values() if c > 1)
duplicate_text_pct = round((duplicate_text_rows / total_norm_texts) * 100, 2) if total_norm_texts else 0.0

doc_duplicate_rows = Counter()
for (doc_id, _txt), count in doc_normalized_texts.items():
    if count > 1:
        doc_duplicate_rows[doc_id] += count - 1

if text_lengths:
    sorted_lengths = sorted(text_lengths)
    def pct(p):
        idx = min(len(sorted_lengths) - 1, int(p * len(sorted_lengths)))
        return sorted_lengths[idx]

    short_chunks = sum(1 for l in text_lengths if l < 200)
    long_chunks = sum(1 for l in text_lengths if l > 2000)
    summary_lengths = {
        "min": min(text_lengths),
        "p50": pct(0.50),
        "p90": pct(0.90),
        "p95": pct(0.95),
        "p99": pct(0.99),
        "max": max(text_lengths),
        "mean": round(statistics.mean(text_lengths), 2),
        "chunks_lt_200_pct": round((short_chunks / len(text_lengths)) * 100, 2),
        "chunks_gt_2000_pct": round((long_chunks / len(text_lengths)) * 100, 2),
    }
else:
    summary_lengths = {
        "min": 0,
        "p50": 0,
        "p90": 0,
        "p95": 0,
        "p99": 0,
        "max": 0,
        "mean": 0,
        "chunks_lt_200_pct": 0,
        "chunks_gt_2000_pct": 0,
    }

top_duplicate_texts = []
for txt, count in normalized_texts.most_common(10):
    if count < 2:
        break
    top_duplicate_texts.append({"count": count, "sample": txt[:160]})


def pct_or_zero(num, den):
    if den <= 0:
        return 0.0
    return round((num / den) * 100, 2)


doc_hotspots = []
for doc_id, total in doc_total_rows.items():
    if total <= 0:
        continue
    doc_hotspots.append(
        {
            "doc_id": doc_id,
            "paper_id": biblio_doc_to_paper_id.get(doc_id),
            "chunk_rows": total,
            "short_lt_200_rows": int(doc_short_rows.get(doc_id, 0)),
            "short_lt_200_pct": pct_or_zero(doc_short_rows.get(doc_id, 0), total),
            "long_gt_2000_rows": int(doc_long_rows.get(doc_id, 0)),
            "long_gt_2000_pct": pct_or_zero(doc_long_rows.get(doc_id, 0), total),
            "duplicate_text_rows": int(doc_duplicate_rows.get(doc_id, 0)),
            "duplicate_text_pct": pct_or_zero(doc_duplicate_rows.get(doc_id, 0), total),
        }
    )

top_short_hotspots = sorted(
    doc_hotspots,
    key=lambda x: (x["short_lt_200_pct"], x["chunk_rows"]),
    reverse=True,
)[:15]
top_duplicate_hotspots = sorted(
    [x for x in doc_hotspots if x["duplicate_text_rows"] > 0],
    key=lambda x: (x["duplicate_text_pct"], x["duplicate_text_rows"]),
    reverse=True,
)[:15]
top_long_hotspots = sorted(
    [x for x in doc_hotspots if x["long_gt_2000_rows"] > 0],
    key=lambda x: (x["long_gt_2000_pct"], x["long_gt_2000_rows"]),
    reverse=True,
)[:15]

report = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "paths": {
        "biblio": str(biblio_path.relative_to(project_root)),
        "chunks": str(chunks_path.relative_to(project_root)),
        "paper_catalog": str(catalog_path.relative_to(project_root)),
        "pdfs_dir": str(pdfs_dir.relative_to(project_root)),
    },
    "totals": {
        "biblio_rows": len(biblio),
        "chunks_rows": len(chunks),
        "paper_catalog_rows": len(catalog),
    },
    "integrity": {
        "biblio_invalid_json": biblio_invalid,
        "chunks_invalid_json": chunks_invalid,
        "paper_catalog_invalid_json": catalog_invalid,
        "unknown_chunk_doc_rows": unknown_doc_rows,
        "duplicate_chunk_key_rows": duplicate_chunk_key_rows,
        "empty_chunk_text_rows": empty_text,
        "missing_chunk_paper_id_rows": missing_paper_id,
        "missing_chunk_source_page_rows": missing_source_page,
        "biblio_extraction_path_missing_rows": missing_extraction_path,
        "biblio_chunk_count_mismatch_docs": chunk_count_mismatches,
        "paper_ids_in_biblio_not_catalog_count": len(paper_ids_in_biblio_not_catalog),
        "paper_ids_in_catalog_not_biblio_count": len(paper_ids_in_catalog_not_biblio),
        "pdf_non_pdf_entry_count": len(pdf_non_pdf_entries),
    },
    "distribution": summary_lengths,
    "duplicates": {
        "duplicate_text_rows": duplicate_text_rows,
        "duplicate_text_pct": duplicate_text_pct,
        "top_duplicate_texts": top_duplicate_texts,
    },
    "hotspots": {
        "top_short_chunks_pct": top_short_hotspots,
        "top_duplicate_text_pct": top_duplicate_hotspots,
        "top_long_chunks_pct": top_long_hotspots,
    },
    "samples": {
        "paper_ids_in_biblio_not_catalog_head": paper_ids_in_biblio_not_catalog[:10],
        "paper_ids_in_catalog_not_biblio_head": paper_ids_in_catalog_not_biblio[:10],
        "pdf_non_pdf_entries": pdf_non_pdf_entries,
    },
}

json_path = output_dir / "data_quality_audit.json"
md_path = output_dir / "data_quality_audit.md"
hotspots_json_path = output_dir / "data_quality_hotspots.json"
hotspots_md_path = output_dir / "data_quality_hotspots.md"
json_path.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
hotspots_json_path.write_text(
    json.dumps(report["hotspots"], indent=2, ensure_ascii=True) + "\n",
    encoding="utf-8",
)

lines = []
lines.append("# Data Quality Audit")
lines.append("")
lines.append(f"- Generated (UTC): `{report['generated_at_utc']}`")
lines.append(f"- Biblio rows: `{report['totals']['biblio_rows']}`")
lines.append(f"- Chunk rows: `{report['totals']['chunks_rows']}`")
lines.append(f"- Catalog rows: `{report['totals']['paper_catalog_rows']}`")
lines.append("")
lines.append("## Integrity")
for k, v in report["integrity"].items():
    lines.append(f"- `{k}`: `{v}`")
lines.append("")
lines.append("## Length Distribution")
for k, v in report["distribution"].items():
    lines.append(f"- `{k}`: `{v}`")
lines.append("")
lines.append("## Duplicate Text (Top)")
for item in report["duplicates"]["top_duplicate_texts"]:
    lines.append(f"- `{item['count']}`: `{item['sample']}`")
lines.append("")

lines.append("## Doc Hotspots (Short Chunks)")
for item in report["hotspots"]["top_short_chunks_pct"][:10]:
    lines.append(
        "- "
        f"`{item.get('paper_id') or item.get('doc_id')}`: "
        f"short<200={item['short_lt_200_rows']}/{item['chunk_rows']} ({item['short_lt_200_pct']}%)"
    )
lines.append("")

lines.append("## Doc Hotspots (Duplicate Text)")
for item in report["hotspots"]["top_duplicate_text_pct"][:10]:
    lines.append(
        "- "
        f"`{item.get('paper_id') or item.get('doc_id')}`: "
        f"duplicate={item['duplicate_text_rows']}/{item['chunk_rows']} ({item['duplicate_text_pct']}%)"
    )
lines.append("")

md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

hotspot_lines = []
hotspot_lines.append("# Data Quality Hotspots")
hotspot_lines.append("")
hotspot_lines.append("## Top Short Chunk Documents")
for item in report["hotspots"]["top_short_chunks_pct"]:
    hotspot_lines.append(
        "- "
        f"`{item.get('paper_id') or item.get('doc_id')}` "
        f"(doc_id={item['doc_id']}): short<200={item['short_lt_200_rows']}/{item['chunk_rows']} "
        f"({item['short_lt_200_pct']}%), long>2000={item['long_gt_2000_rows']}/{item['chunk_rows']} "
        f"({item['long_gt_2000_pct']}%), duplicate={item['duplicate_text_rows']}/{item['chunk_rows']} "
        f"({item['duplicate_text_pct']}%)"
    )
hotspot_lines.append("")
hotspot_lines.append("## Top Duplicate Text Documents")
for item in report["hotspots"]["top_duplicate_text_pct"]:
    hotspot_lines.append(
        "- "
        f"`{item.get('paper_id') or item.get('doc_id')}` "
        f"(doc_id={item['doc_id']}): duplicate={item['duplicate_text_rows']}/{item['chunk_rows']} "
        f"({item['duplicate_text_pct']}%), short<200={item['short_lt_200_rows']}/{item['chunk_rows']} "
        f"({item['short_lt_200_pct']}%)"
    )
hotspot_lines.append("")
hotspot_lines.append("## Top Long Chunk Documents")
for item in report["hotspots"]["top_long_chunks_pct"]:
    hotspot_lines.append(
        "- "
        f"`{item.get('paper_id') or item.get('doc_id')}` "
        f"(doc_id={item['doc_id']}): long>2000={item['long_gt_2000_rows']}/{item['chunk_rows']} "
        f"({item['long_gt_2000_pct']}%), max contribution to retrieval latency risk."
    )
hotspot_lines.append("")
hotspots_md_path.write_text("\n".join(hotspot_lines) + "\n", encoding="utf-8")

print(json_path)
print(md_path)
print(hotspots_json_path)
print(hotspots_md_path)
PY

echo "Audit artifacts written to ${OUTPUT_DIR}"
