#!/usr/bin/env bash
set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SELF_DIR}/../.." && pwd)"

cd "${REPO_ROOT}"

exec python3 "${SELF_DIR}/analyze_chemical_ocsr_fusion.py" "$@"
