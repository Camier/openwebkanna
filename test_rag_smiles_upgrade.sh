#!/usr/bin/env bash
# RAG & SMILES Architecture Upgrade Tests
# Run inside micromamba: micromamba run -n smiles-extraction python test_rag_smiles_upgrade.py

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "${SCRIPT_DIR}"

# Check if micromamba env exists
if command -v micromamba >/dev/null 2>&1 && micromamba env list | grep -q "smiles-extraction"; then
    echo "Running inside micromamba smiles-extraction environment..."
    micromamba run -n smiles-extraction python test_rag_smiles_upgrade.py
else
    echo "ERROR: micromamba smiles-extraction environment not found."
    echo "Please create it first:"
    echo "  cd smiles-pipeline/envs"
    echo "  micromamba env create -f environment-smiles-extraction.yml"
    exit 1
fi
