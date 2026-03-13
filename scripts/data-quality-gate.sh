#!/bin/bash

###############################################################################
# data-quality-gate.sh - Enforce data quality thresholds on corpus artifacts
#
# Defaults:
#   - Uses curated corpus files when present
#   - Fails if curated files are missing (unless --allow-raw-fallback or
#     --skip-if-missing is provided)
#
# Usage:
#   ./scripts/data-quality-gate.sh
#   ./scripts/data-quality-gate.sh --max-duplicate-text-pct 3.5
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

OUTPUT_DIR="${PROJECT_ROOT}/artifacts/data-quality/gate"

CHUNKS_FILE="${PROJECT_ROOT}/data/corpus/chunks_corpus.curated.jsonl"
BIBLIO_FILE="${PROJECT_ROOT}/data/corpus/biblio_corpus.curated.jsonl"

REQUIRE_CURATED=1
SKIP_IF_MISSING=0
PROFILE="default"

MAX_DUPLICATE_TEXT_PCT="3.0"
MAX_SHORT_CHUNKS_PCT="35.0"
MAX_CHUNK_LENGTH=4000
MAX_MISSING_PAPER_ID=0
MAX_MISSING_SOURCE_PAGE=0
MAX_BIBLIO_MISSING_PATH=0
WARN_DUPLICATE_TEXT_PCT=""
WARN_SHORT_CHUNKS_PCT=""
WARN_LONG_CHUNKS_PCT=""

if [[ -n ${CORPUS_CHUNKS_FILE:-} ]]; then
    CHUNKS_FILE="${CORPUS_CHUNKS_FILE}"
fi
if [[ -n ${CORPUS_BIBLIO_FILE:-} ]]; then
    BIBLIO_FILE="${CORPUS_BIBLIO_FILE}"
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --chunks-file)
            CHUNKS_FILE="$2"
            shift 2
            ;;
        --biblio-file)
            BIBLIO_FILE="$2"
            shift 2
            ;;
        --max-duplicate-text-pct)
            MAX_DUPLICATE_TEXT_PCT="$2"
            shift 2
            ;;
        --max-short-chunks-pct)
            MAX_SHORT_CHUNKS_PCT="$2"
            shift 2
            ;;
        --max-chunk-length)
            MAX_CHUNK_LENGTH="$2"
            shift 2
            ;;
        --max-missing-paper-id)
            MAX_MISSING_PAPER_ID="$2"
            shift 2
            ;;
        --max-missing-source-page)
            MAX_MISSING_SOURCE_PAGE="$2"
            shift 2
            ;;
        --max-biblio-missing-paths)
            MAX_BIBLIO_MISSING_PATH="$2"
            shift 2
            ;;
        --allow-raw-fallback)
            REQUIRE_CURATED=0
            shift
            ;;
        --skip-if-missing)
            SKIP_IF_MISSING=1
            shift
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --help | -h)
            cat <<'EOF'
Usage: ./scripts/data-quality-gate.sh [options]

Options:
  --output-dir <path>               Audit artifact output directory
  --chunks-file <path>              Override chunks file
  --biblio-file <path>              Override biblio file
  --max-duplicate-text-pct <float>  Threshold for duplicate text percentage
  --max-short-chunks-pct <float>    Threshold for chunks shorter than 200 chars
  --max-chunk-length <int>          Threshold for max chunk length
  --max-missing-paper-id <int>      Threshold for missing paper_id rows
  --max-missing-source-page <int>   Threshold for missing source_page rows
  --max-biblio-missing-paths <int>  Threshold for missing extraction_path rows
  --allow-raw-fallback              Allow fallback to non-curated corpus files
  --skip-if-missing                 Exit 0 when required corpus files are missing
  --profile <default|strict>        strict adds warning thresholds (non-fatal)
EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

if [[ ! -f ${CHUNKS_FILE} || ! -f ${BIBLIO_FILE} ]]; then
    if [[ ${REQUIRE_CURATED} -eq 0 ]]; then
        [[ -f ${CHUNKS_FILE} ]] || CHUNKS_FILE="${PROJECT_ROOT}/data/corpus/chunks_corpus.jsonl"
        [[ -f ${BIBLIO_FILE} ]] || BIBLIO_FILE="${PROJECT_ROOT}/data/corpus/biblio_corpus.jsonl"
    fi
fi

if [[ ! -f ${CHUNKS_FILE} || ! -f ${BIBLIO_FILE} ]]; then
    if [[ ${SKIP_IF_MISSING} -eq 1 ]]; then
        echo "Data quality gate skipped: corpus files missing."
        echo "  chunks: ${CHUNKS_FILE}"
        echo "  biblio: ${BIBLIO_FILE}"
        exit 0
    fi
    echo "Data quality gate failed: required corpus files missing." >&2
    echo "  chunks: ${CHUNKS_FILE}" >&2
    echo "  biblio: ${BIBLIO_FILE}" >&2
    exit 2
fi

case "${PROFILE}" in
    default) ;;
    strict)
        WARN_DUPLICATE_TEXT_PCT="1.0"
        WARN_SHORT_CHUNKS_PCT="5.0"
        WARN_LONG_CHUNKS_PCT="5.0"
        ;;
    *)
        echo "Unknown profile: ${PROFILE}. Expected: default|strict" >&2
        exit 2
        ;;
esac

mkdir -p "${OUTPUT_DIR}"

"${PROJECT_ROOT}/scripts/audit-data-quality.sh" \
    --output-dir "${OUTPUT_DIR}" \
    --chunks-file "${CHUNKS_FILE}" \
    --biblio-file "${BIBLIO_FILE}" >/dev/null

AUDIT_JSON="${OUTPUT_DIR}/data_quality_audit.json"

python - \
    "${AUDIT_JSON}" \
    "${MAX_DUPLICATE_TEXT_PCT}" \
    "${MAX_SHORT_CHUNKS_PCT}" \
    "${MAX_CHUNK_LENGTH}" \
    "${MAX_MISSING_PAPER_ID}" \
    "${MAX_MISSING_SOURCE_PAGE}" \
    "${MAX_BIBLIO_MISSING_PATH}" \
    "${PROFILE}" \
    "${WARN_DUPLICATE_TEXT_PCT}" \
    "${WARN_SHORT_CHUNKS_PCT}" \
    "${WARN_LONG_CHUNKS_PCT}" <<'PY'
import json
import sys
from pathlib import Path

audit_path = Path(sys.argv[1])
max_duplicate_text_pct = float(sys.argv[2])
max_short_chunks_pct = float(sys.argv[3])
max_chunk_length = int(sys.argv[4])
max_missing_paper_id = int(sys.argv[5])
max_missing_source_page = int(sys.argv[6])
max_biblio_missing_paths = int(sys.argv[7])
profile = sys.argv[8]
warn_duplicate_text_pct = float(sys.argv[9]) if sys.argv[9] else None
warn_short_chunks_pct = float(sys.argv[10]) if sys.argv[10] else None
warn_long_chunks_pct = float(sys.argv[11]) if sys.argv[11] else None

audit = json.loads(audit_path.read_text(encoding="utf-8"))
integrity = audit["integrity"]
distribution = audit["distribution"]
duplicate_pct = float(audit["duplicates"]["duplicate_text_pct"])
short_pct = float(distribution["chunks_lt_200_pct"])
long_pct = float(distribution.get("chunks_gt_2000_pct", 0))
max_len = int(distribution["max"])

checks = [
    (
        "missing_chunk_paper_id_rows",
        int(integrity["missing_chunk_paper_id_rows"]),
        max_missing_paper_id,
        "le",
    ),
    (
        "missing_chunk_source_page_rows",
        int(integrity["missing_chunk_source_page_rows"]),
        max_missing_source_page,
        "le",
    ),
    (
        "biblio_extraction_path_missing_rows",
        int(integrity["biblio_extraction_path_missing_rows"]),
        max_biblio_missing_paths,
        "le",
    ),
    (
        "duplicate_text_pct",
        duplicate_pct,
        max_duplicate_text_pct,
        "le",
    ),
    (
        "chunks_lt_200_pct",
        short_pct,
        max_short_chunks_pct,
        "le",
    ),
    (
        "max_chunk_length",
        max_len,
        max_chunk_length,
        "le",
    ),
]

failures = []
for name, actual, threshold, op in checks:
    if op == "le" and actual > threshold:
        failures.append((name, actual, threshold))

print(f"Data quality gate summary (profile={profile}):")
for name, actual, threshold, _ in checks:
    print(f"  {name}: actual={actual} threshold<={threshold}")

if failures:
    print("Data quality gate failed:")
    for name, actual, threshold in failures:
        print(f"  - {name}: actual={actual} exceeds threshold={threshold}")
    sys.exit(1)

warnings = []
if warn_duplicate_text_pct is not None and duplicate_pct > warn_duplicate_text_pct:
    warnings.append(("duplicate_text_pct", duplicate_pct, warn_duplicate_text_pct))
if warn_short_chunks_pct is not None and short_pct > warn_short_chunks_pct:
    warnings.append(("chunks_lt_200_pct", short_pct, warn_short_chunks_pct))
if warn_long_chunks_pct is not None and long_pct > warn_long_chunks_pct:
    warnings.append(("chunks_gt_2000_pct", long_pct, warn_long_chunks_pct))

if warnings:
    print("Data quality strict warnings:")
    for name, actual, threshold in warnings:
        print(f"  - {name}: actual={actual} exceeds warn-threshold={threshold}")

print("Data quality gate passed.")
PY
