#!/bin/bash

###############################################################################
# Check Indigo Service health and basic API readiness
###############################################################################

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

INDIGO_SERVICE_HOST="${INDIGO_SERVICE_HOST:-${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}}"
INDIGO_SERVICE_PORT="${INDIGO_SERVICE_PORT:-8012}"
INDIGO_SERVICE_BASE_URL="${INDIGO_SERVICE_BASE_URL:-http://${INDIGO_SERVICE_HOST}:${INDIGO_SERVICE_PORT}}"
INDIGO_SERVICE_TIMEOUT="${INDIGO_SERVICE_TIMEOUT:-8}"
QUIET=false

while [ $# -gt 0 ]; do
    case "$1" in
        -q | --quiet)
            QUIET=true
            shift
            ;;
        -h | --help)
            echo "Usage: $0 [--quiet]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [ "$QUIET" = false ]; then
    print_header "Indigo Service Check"
fi

tmp_file="$(mktemp /tmp/indigo_service_info.XXXXXX.json)"
endpoint_hit=""
trap 'rm -f "$tmp_file"' EXIT

probe_endpoint() {
    local url="$1"
    local code
    code="$(curl -sS -m "$INDIGO_SERVICE_TIMEOUT" -o "$tmp_file" -w "%{http_code}" "$url" 2>/dev/null || true)"
    if [ "$code" = "200" ]; then
        endpoint_hit="$url"
        return 0
    fi
    return 1
}

if ! probe_endpoint "${INDIGO_SERVICE_BASE_URL}/v2/indigo/info" &&
    ! probe_endpoint "${INDIGO_SERVICE_BASE_URL}/indigo/info" &&
    ! probe_endpoint "${INDIGO_SERVICE_BASE_URL}/v2/info"; then
    if [ "$QUIET" = false ]; then
        print_error "Indigo Service endpoint probe failed at ${INDIGO_SERVICE_BASE_URL}"
    fi
    exit 1
fi

if [ "$QUIET" = false ]; then
    print_success "Indigo Service healthy"
    print_info "Endpoint: ${endpoint_hit}"
    if command_exists jq; then
        print_info "Version: $(jq -r '.Indigo.version // .indigo_version // .version // .indigoVersion // .server_version // "unknown"' "$tmp_file" 2>/dev/null || echo "unknown")"
    fi
fi

exit 0
