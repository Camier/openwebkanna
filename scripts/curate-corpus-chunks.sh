#!/bin/bash

###############################################################################
# curate-corpus-chunks.sh - Enrich and curate chunk corpus for retrieval usage
#
# Actions:
#  1) Enrich data/corpus/chunks_corpus.jsonl in place with:
#     - paper_id (mapped from biblio by doc_id)
#     - source_page (derived from page, fallback chunk_id /page/<n>/...)
#  2) Build curated retrieval corpus at:
#     data/corpus/chunks_corpus.curated.jsonl
#  3) Quarantine removed rows with reasons and emit stats artifact
#
# Usage:
#   ./scripts/curate-corpus-chunks.sh
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BIBLIO_PATH="${PROJECT_ROOT}/data/corpus/biblio_corpus.jsonl"
RAW_CHUNKS_PATH="${PROJECT_ROOT}/data/corpus/chunks_corpus.jsonl"
CURATED_CHUNKS_PATH="${PROJECT_ROOT}/data/corpus/chunks_corpus.curated.jsonl"
BIBLIO_CURATED_PATH="${PROJECT_ROOT}/data/corpus/biblio_corpus.curated.jsonl"

QUARANTINE_DIR="${PROJECT_ROOT}/data/metadata/quarantine/chunks"
BACKUP_DIR="${PROJECT_ROOT}/data/metadata/quarantine/backups"
ARTIFACT_DIR="${PROJECT_ROOT}/artifacts/data-quality"

mkdir -p "${QUARANTINE_DIR}" "${BACKUP_DIR}" "${ARTIFACT_DIR}"

TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RAW_BACKUP="${BACKUP_DIR}/chunks_corpus.jsonl.pre_curate.${TIMESTAMP}.bak"
REMOVED_PATH="${QUARANTINE_DIR}/chunks_removed.${TIMESTAMP}.jsonl"
STATS_PATH="${ARTIFACT_DIR}/curation_stats.${TIMESTAMP}.json"

cp "${RAW_CHUNKS_PATH}" "${RAW_BACKUP}"

python - "${BIBLIO_PATH}" "${RAW_CHUNKS_PATH}" "${CURATED_CHUNKS_PATH}" "${BIBLIO_CURATED_PATH}" "${REMOVED_PATH}" "${STATS_PATH}" <<'PY'
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
import sys

biblio_path = Path(sys.argv[1])
raw_chunks_path = Path(sys.argv[2])
curated_chunks_path = Path(sys.argv[3])
biblio_curated_path = Path(sys.argv[4])
removed_path = Path(sys.argv[5])
stats_path = Path(sys.argv[6])

tmp_enriched = raw_chunks_path.with_suffix(raw_chunks_path.suffix + ".tmp.enriched")

doc_meta = {}
with biblio_path.open("r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        obj = json.loads(line)
        doc_id = obj.get("doc_id")
        if doc_id:
            doc_meta[doc_id] = {
                "paper_id": obj.get("paper_id"),
                "title": obj.get("title"),
            }

boilerplate_exact = {
    "author manuscript",
    "references",
    "notes and references",
    "acknowledgements",
    "introduction",
    "1. introduction",
    "article info",
    "elsevier logo featuring a tree and the word elsevier",
}

boilerplate_regexes = [
    ("page_furniture_blank", re.compile(r"^\s*this page intentionally left blank\s*$", re.IGNORECASE)),
    ("page_furniture_updates", re.compile(r"^\s*check for updates (icon|button)\s*$", re.IGNORECASE)),
    ("page_furniture_continued", re.compile(r"^\s*\(continued on following page\)\s*$", re.IGNORECASE)),
    ("page_furniture_cover", re.compile(r"^\s*cover image of the journal of ethnopharmacology\s*$", re.IGNORECASE)),
]

drop_block_types = {
    "ComplexRegion",
    "Equation",
    "TableOfContents",
    "Footnote",
}

max_chunk_len = 4000
target_len = 2500
overlap = 250
min_semantic_chunk_len = 180
max_merged_chunk_len = 900

source_page_re = re.compile(r"^/page/(\d+)/")


def normalize_text(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip()).lower()


def split_long_text(text: str):
    text = text.strip()
    if len(text) <= max_chunk_len:
        return [text]
    out = []
    i = 0
    while i < len(text):
        j = min(i + target_len, len(text))
        out.append(text[i:j])
        if j == len(text):
            break
        i = max(i + target_len - overlap, i + 1)
    return out


def match_boilerplate_regex(text: str):
    for name, pattern in boilerplate_regexes:
        if pattern.match(text):
            return name
    return None


enrich_stats = Counter()
curation_stats = Counter()
curated_counts_by_doc = Counter()
seen_short_text_by_doc = defaultdict(set)

# Pass 1: Enrich raw chunks in place (via temp file)
with raw_chunks_path.open("r", encoding="utf-8") as in_f, tmp_enriched.open("w", encoding="utf-8") as out_f:
    for line in in_f:
        line = line.strip()
        if not line:
            continue
        obj = json.loads(line)
        enrich_stats["rows_in"] += 1
        doc_id = obj.get("doc_id")
        meta = doc_meta.get(doc_id, {})

        if not obj.get("paper_id") and meta.get("paper_id"):
            obj["paper_id"] = meta["paper_id"]
            enrich_stats["paper_id_backfilled"] += 1

        if not obj.get("paper_title") and meta.get("title"):
            obj["paper_title"] = meta["title"]
            enrich_stats["paper_title_backfilled"] += 1

        if obj.get("source_page") in (None, ""):
            page = obj.get("page")
            if isinstance(page, int):
                obj["source_page"] = page
                enrich_stats["source_page_from_page"] += 1
            else:
                chunk_id = obj.get("chunk_id", "")
                m = source_page_re.match(chunk_id)
                if m:
                    obj["source_page"] = int(m.group(1))
                    enrich_stats["source_page_from_chunk_id"] += 1

        out_f.write(json.dumps(obj, ensure_ascii=True) + "\n")
        enrich_stats["rows_out"] += 1

tmp_enriched.replace(raw_chunks_path)

# Pass 2: Build curated corpus and quarantine removed chunks
with raw_chunks_path.open("r", encoding="utf-8") as in_f, \
        curated_chunks_path.open("w", encoding="utf-8") as curated_f, \
        removed_path.open("w", encoding="utf-8") as removed_f:
    state = {"pending_chunk": None}

    def write_curated(out_obj):
        merge_count = int(out_obj.pop("__merge_count", 1))
        if merge_count > 1:
            out_obj["merged_chunks"] = merge_count
        curated_f.write(json.dumps(out_obj, ensure_ascii=True) + "\n")
        curation_stats["rows_out"] += 1
        if out_obj.get("doc_id"):
            curated_counts_by_doc[out_obj["doc_id"]] += 1

    def flush_pending():
        if state["pending_chunk"] is None:
            return
        write_curated(state["pending_chunk"])
        state["pending_chunk"] = None

    def queue_or_merge(next_chunk):
        if state["pending_chunk"] is None:
            state["pending_chunk"] = next_chunk
            return

        same_doc = state["pending_chunk"].get("doc_id") == next_chunk.get("doc_id")
        same_page = state["pending_chunk"].get("source_page") == next_chunk.get("source_page")
        pending_text = state["pending_chunk"].get("text", "")
        next_text = next_chunk.get("text", "")
        can_merge = (
            same_doc
            and same_page
            and len(pending_text) < min_semantic_chunk_len
            and (len(pending_text) + 2 + len(next_text)) <= max_merged_chunk_len
        )

        if can_merge:
            state["pending_chunk"]["text"] = f"{pending_text}\n\n{next_text}"
            state["pending_chunk"]["chunk_id"] = f"{state['pending_chunk'].get('chunk_id')}+{next_chunk.get('chunk_id')}"
            state["pending_chunk"]["__merge_count"] = int(state["pending_chunk"].get("__merge_count", 1)) + int(next_chunk.get("__merge_count", 1))
            curation_stats["rows_merged_pairs"] += 1
            curation_stats["rows_merged_source_rows"] += int(next_chunk.get("__merge_count", 1))
            return

        flush_pending()
        state["pending_chunk"] = next_chunk

    for line in in_f:
        line = line.strip()
        if not line:
            continue
        obj = json.loads(line)
        curation_stats["rows_in"] += 1

        text = (obj.get("text") or "").strip()
        norm = normalize_text(text) if text else ""
        block_type = obj.get("block_type", "")

        drop_reason = None
        if not text:
            drop_reason = "empty_text"
        elif norm in boilerplate_exact:
            drop_reason = "boilerplate_exact"
        else:
            regex_name = match_boilerplate_regex(text)
            if regex_name:
                drop_reason = f"boilerplate_regex:{regex_name}"

        if not drop_reason:
            doc_id = obj.get("doc_id") or "_unknown"
            if len(text) <= 250 and norm in seen_short_text_by_doc[doc_id]:
                drop_reason = "duplicate_short_text_in_doc"
            elif len(text) <= 250 and norm:
                seen_short_text_by_doc[doc_id].add(norm)

        if not drop_reason:
            if block_type in drop_block_types:
                drop_reason = f"block_type:{block_type}"
            elif block_type == "SectionHeader" and len(text) < 40:
                drop_reason = "short_section_header_lt40"
            elif block_type == "Text" and len(text) < 20:
                drop_reason = "short_text_lt20"
            elif block_type == "Caption" and len(text) < 20:
                drop_reason = "short_caption_lt20"

        if drop_reason:
            removed = {
                "doc_id": obj.get("doc_id"),
                "paper_id": obj.get("paper_id"),
                "chunk_id": obj.get("chunk_id"),
                "block_type": block_type,
                "text_len": len(text),
                "reason": drop_reason,
                "text_preview": text[:200],
            }
            removed_f.write(json.dumps(removed, ensure_ascii=True) + "\n")
            curation_stats["rows_removed"] += 1
            curation_stats[f"removed:{drop_reason}"] += 1
            continue

        split_parts = split_long_text(text)
        total_parts = len(split_parts)

        for idx, part in enumerate(split_parts, start=1):
            out = {
                "doc_id": obj.get("doc_id"),
                "paper_id": obj.get("paper_id"),
                "paper_title": obj.get("paper_title"),
                "chunk_id": obj.get("chunk_id") if total_parts == 1 else f"{obj.get('chunk_id')}#seg-{idx}",
                "source_page": obj.get("source_page"),
                "page": obj.get("page"),
                "block_type": obj.get("block_type"),
                "text": part,
                "__merge_count": 1,
            }
            if total_parts > 1:
                out["split_part"] = idx
                out["split_parts_total"] = total_parts
                curation_stats["rows_split_generated"] += 1

            queue_or_merge(out)

    flush_pending()

# Build curated biblio with chunk_count aligned to curated chunks
with biblio_path.open("r", encoding="utf-8") as in_f, biblio_curated_path.open("w", encoding="utf-8") as out_f:
    for line in in_f:
        line = line.strip()
        if not line:
            continue
        obj = json.loads(line)
        doc_id = obj.get("doc_id")
        obj["chunk_count_raw"] = obj.get("chunk_count", 0)
        obj["chunk_count"] = int(curated_counts_by_doc.get(doc_id, 0))
        out_f.write(json.dumps(obj, ensure_ascii=True) + "\n")

# Summary quality on curated output
dup_counter = Counter()
lengths = []
with curated_chunks_path.open("r", encoding="utf-8") as f:
    for line in f:
        obj = json.loads(line)
        text = (obj.get("text") or "").strip()
        lengths.append(len(text))
        n = normalize_text(text)
        if n:
            dup_counter[n] += 1

duplicate_rows = sum(c - 1 for c in dup_counter.values() if c > 1)
short_rows = sum(1 for l in lengths if l < 200)

stats = {
    "enrichment": dict(enrich_stats),
    "curation": dict(curation_stats),
    "curated_metrics": {
        "rows": len(lengths),
        "duplicate_text_rows": duplicate_rows,
        "duplicate_text_pct": round((duplicate_rows / len(lengths)) * 100, 2) if lengths else 0.0,
        "short_lt_200_rows": short_rows,
        "short_lt_200_pct": round((short_rows / len(lengths)) * 100, 2) if lengths else 0.0,
        "max_text_len": max(lengths) if lengths else 0,
    },
    "paths": {
        "raw_chunks": str(raw_chunks_path),
        "curated_chunks": str(curated_chunks_path),
        "curated_biblio": str(biblio_curated_path),
        "removed_chunks": str(removed_path),
    },
}

stats_path.write_text(json.dumps(stats, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
print(json.dumps(stats, ensure_ascii=True))
PY

echo "Raw chunks backup: ${RAW_BACKUP}"
echo "Curated corpus: ${CURATED_CHUNKS_PATH}"
echo "Removed quarantine: ${REMOVED_PATH}"
echo "Stats: ${STATS_PATH}"
