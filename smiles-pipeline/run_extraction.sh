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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

# Default values
INPUT_DIR="data/extractions"
OUTPUT_DIR="smiles-pipeline/data/validated"
LIMIT=""
GPU_FLAG=""
BACKEND_ORDER="molscribe,decimer"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --limit)
            LIMIT="--limit $2"
            shift 2
            ;;
        --gpu)
            GPU_FLAG="--gpu"
            shift
            ;;
        --no-gpu)
            GPU_FLAG="--no-gpu"
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

# Check if conda environment exists
if ! conda env list | grep -q "smiles-extraction"; then
    echo "Creating conda environment..."
    conda env create -f smiles-pipeline/envs/environment-smiles-extraction.yml
fi

# Activate environment
source "$(conda info --base)"/etc/profile.d/conda.sh
conda activate smiles-extraction

# Run pipeline
echo "Starting SMILES extraction pipeline..."
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Backends: $BACKEND_ORDER"
[ -n "$GPU_FLAG" ] && echo "GPU: enabled"
[ -n "$LIMIT" ] && echo "Limit: $LIMIT"
echo ""

python smiles-pipeline/scripts/extract_smiles_pipeline.py \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --backend-order "$BACKEND_ORDER" \
    $GPU_FLAG \
    $LIMIT

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
    echo "  2. Convert to KB: python smiles-pipeline/scripts/jsonl_to_markdown.py --input $OUTPUT_DIR/molecules.jsonl"
    echo "  3. Import: ./smiles-pipeline/import-smiles-to-kb.sh --input-dir $OUTPUT_DIR"
else
    echo "No molecules extracted. Check logs for errors."
fi
