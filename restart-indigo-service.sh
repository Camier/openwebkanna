#!/bin/bash

###############################################################################
# Restart Indigo Service sidecar.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

"${SCRIPT_DIR}/stop-indigo-service.sh"
"${SCRIPT_DIR}/start-indigo-service.sh"
