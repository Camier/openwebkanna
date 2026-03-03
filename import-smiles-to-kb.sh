#!/bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_SCRIPT="${SCRIPT_DIR}/smiles-pipeline/import-smiles-to-kb.sh"

if [ ! -f "$TARGET_SCRIPT" ]; then
    printf '[ERROR] Missing target script: %s\n' "$TARGET_SCRIPT" >&2
    exit 1
fi

exec bash "$TARGET_SCRIPT" "$@"
