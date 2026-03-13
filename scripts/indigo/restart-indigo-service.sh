#!/bin/bash

###############################################################################
# Restart Indigo Service sidecar.
###############################################################################

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"${SELF_DIR}/stop-indigo-service.sh"
"${SELF_DIR}/start-indigo-service.sh"
