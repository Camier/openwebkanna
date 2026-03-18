#!/bin/bash

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
cd "$SCRIPT_DIR"

exec python3 "${SELF_DIR}/extract_chemical_smiles.py" "$@"
