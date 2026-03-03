#!/bin/bash

set -e
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$REPO_ROOT"

INPUT_DIR="data/extractions"
OUTPUT_DIR="/tmp/smiles-smoke/one-command"
BACKEND_ORDER="molscribe,decimer"
MIN_VALIDATION_SCORE=""
OPENWEBUI_URL="http://localhost:${WEBUI_PORT:-3000}"
LIMIT_ARGS=()
GPU_ARGS=(--no-gpu)
KB_NAME="SMILES One-Command $(date +%Y%m%d-%H%M%S)"
KB_DESCRIPTION="One-command extraction/import smoke run"

show_help() {
    cat <<'EOF'
Usage: ./smiles-pipeline/run_extract_and_import.sh [OPTIONS]

Runs extraction and KB import in one command.

Options:
  --input-dir DIR            Input directory (default: data/extractions)
  --output-dir DIR           Output directory (default: /tmp/smiles-smoke/one-command)
  --backend-order LIST       Backend order (default: molscribe,decimer)
  --min-validation-score N   Validation threshold override (default: config file)
  --limit N                  Process first N images
  --gpu                      Use GPU
  --no-gpu                   Force CPU (default)
  --url URL                  OpenWebUI URL (default: http://localhost:${WEBUI_PORT:-3000})
  --kb-name NAME             Knowledge base name
  --kb-description TEXT      Knowledge base description
  -h, --help                 Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --backend-order)
            BACKEND_ORDER="$2"
            shift 2
            ;;
        --min-validation-score)
            MIN_VALIDATION_SCORE="$2"
            shift 2
            ;;
        --limit)
            LIMIT_ARGS=(--limit "$2")
            shift 2
            ;;
        --gpu)
            GPU_ARGS=(--gpu)
            shift
            ;;
        --no-gpu)
            GPU_ARGS=(--no-gpu)
            shift
            ;;
        --url)
            OPENWEBUI_URL="$2"
            shift 2
            ;;
        --kb-name)
            KB_NAME="$2"
            shift 2
            ;;
        --kb-description)
            KB_DESCRIPTION="$2"
            shift 2
            ;;
        -h | --help)
            show_help
            exit 0
            ;;
        *)
            printf 'Unknown option: %s\n' "$1" >&2
            exit 1
            ;;
    esac
done

mkdir -p "$OUTPUT_DIR"

PY_CMD=(
    python
    "smiles-pipeline/scripts/extract_smiles_pipeline.py"
    --input-dir "$INPUT_DIR"
    --output-dir "$OUTPUT_DIR"
    --backend-order "$BACKEND_ORDER"
)
if [ -n "$MIN_VALIDATION_SCORE" ]; then
    PY_CMD+=(--min-validation-score "$MIN_VALIDATION_SCORE")
fi
PY_CMD+=("${GPU_ARGS[@]}")
if [ ${#LIMIT_ARGS[@]} -gt 0 ]; then
    PY_CMD+=("${LIMIT_ARGS[@]}")
fi

if command -v micromamba >/dev/null 2>&1; then
    if ! micromamba env list | awk '{print $1}' | grep -qx "smiles-extraction"; then
        printf 'ERROR: micromamba env smiles-extraction not found\n' >&2
        exit 1
    fi
    micromamba run -n smiles-extraction "${PY_CMD[@]}"
elif command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)"/etc/profile.d/conda.sh
    if ! conda env list | awk '{print $1}' | grep -qx "smiles-extraction"; then
        printf 'ERROR: conda env smiles-extraction not found\n' >&2
        exit 1
    fi
    conda activate smiles-extraction
    "${PY_CMD[@]}"
else
    printf 'ERROR: micromamba/conda not available\n' >&2
    exit 1
fi

"${SCRIPT_DIR}/import-smiles-to-kb.sh" \
    --source-root "$OUTPUT_DIR" \
    --find-min-depth 1 \
    --find-max-depth 1 \
    --limit 1 \
    --url "$OPENWEBUI_URL" \
    --kb-name "$KB_NAME" \
    --kb-description "$KB_DESCRIPTION"
