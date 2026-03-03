#!/bin/bash
# SMILES Pipeline - Quick Start Script
#
# Usage:
#   ./run_extraction.sh --limit 10    # Test on 10 images
#   ./run_extraction.sh --gpu          # Full run with GPU
#
# Defaults:
#   - Backend order: molscribe,decimer
#   - Min validation score: 70
#   - Output: smiles-pipeline/data/validated/

set -e
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

# Default values
INPUT_DIR="data/extractions"
OUTPUT_DIR="smiles-pipeline/data/validated"
LIMIT_ARGS=()
GPU_ARGS=()
BACKEND_ORDER="molscribe,decimer"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
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
        -h | --help)
            echo "SMILES Pipeline - Extract molecules from chemical structure images"
            echo ""
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --limit N              Process only first N images (for testing)"
            echo "  --gpu                  Enable GPU acceleration (default: auto)"
            echo "  --no-gpu              Force CPU mode"
            echo "  --input-dir DIR       Input directory with images (default: data/extractions)"
            echo "  --output-dir DIR      Output directory (default: smiles-pipeline/data/validated)"
            echo "  --backend-order LIST  Comma-separated backend order (default: molscribe,decimer)"
            echo "  -h, --help            Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Run pipeline
echo "Starting SMILES extraction pipeline..."
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Backends: $BACKEND_ORDER"
[ ${#GPU_ARGS[@]} -gt 0 ] && echo "GPU flag: ${GPU_ARGS[*]}"
[ ${#LIMIT_ARGS[@]} -gt 0 ] && echo "Limit: ${LIMIT_ARGS[1]}"
echo ""

PY_CMD=(
    python
    smiles-pipeline/scripts/extract_smiles_pipeline.py
    --input-dir "$INPUT_DIR"
    --output-dir "$OUTPUT_DIR"
    --backend-order "$BACKEND_ORDER"
)
if [ ${#GPU_ARGS[@]} -gt 0 ]; then
    PY_CMD+=("${GPU_ARGS[@]}")
fi
if [ ${#LIMIT_ARGS[@]} -gt 0 ]; then
    PY_CMD+=("${LIMIT_ARGS[@]}")
fi

# Guard runtime so extraction always runs in the dedicated smiles-extraction environment.
if command -v micromamba >/dev/null 2>&1; then
    if micromamba env list | awk '{print $1}' | grep -qx "smiles-extraction"; then
        micromamba run -n smiles-extraction "${PY_CMD[@]}"
    else
        echo "ERROR: micromamba env 'smiles-extraction' not found."
        echo "Create it with:"
        echo "  micromamba env create -f smiles-pipeline/envs/environment-smiles-extraction.yml --yes"
        exit 1
    fi
elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)"/etc/profile.d/conda.sh
    if ! conda env list | awk '{print $1}' | grep -qx "smiles-extraction"; then
        echo "ERROR: conda env 'smiles-extraction' not found."
        echo "Create it with:"
        echo "  conda env create -f smiles-pipeline/envs/environment-smiles-extraction.yml"
        exit 1
    fi
    conda activate smiles-extraction
    "${PY_CMD[@]}"
else
    echo "ERROR: Neither micromamba nor conda is available."
    echo "Install micromamba (recommended) or conda, then create smiles-extraction env."
    exit 1
fi

# Show results
echo ""
echo "Pipeline complete!"
echo ""
if [ -f "$OUTPUT_DIR/molecules.jsonl" ]; then
    COUNT=$(wc -l <"$OUTPUT_DIR/molecules.jsonl")
    echo "Results: $OUTPUT_DIR/molecules.jsonl ($COUNT molecules)"
    echo ""
    echo "Next steps:"
    echo "  1. Review: head -n 5 $OUTPUT_DIR/molecules.jsonl | jq"
    echo "  2. Review manual queue: head -n 5 $OUTPUT_DIR/molecules_manual_review.jsonl | jq"
    echo "  3. Import READY set: ./smiles-pipeline/import-smiles-to-kb.sh --source-root $OUTPUT_DIR --find-min-depth 1 --find-max-depth 1 --file-name molecules_high_confidence.jsonl"
else
    echo "No molecules extracted. Check logs for errors."
fi
