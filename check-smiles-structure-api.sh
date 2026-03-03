#!/bin/bash

###############################################################################
# Check SMILES Structure Search API health
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

SMILES_STRUCTURE_API_DOCKER_CONTAINER="${SMILES_STRUCTURE_API_DOCKER_CONTAINER:-smiles_structure_api}"
SMILES_STRUCTURE_API_HOST="${SMILES_STRUCTURE_API_HOST:-127.0.0.1}"
SMILES_STRUCTURE_API_PORT="${SMILES_STRUCTURE_API_PORT:-8011}"
SMILES_STRUCTURE_API_BASE_URL="${SMILES_STRUCTURE_API_BASE_URL:-http://${SMILES_STRUCTURE_API_HOST}:${SMILES_STRUCTURE_API_PORT}}"
SMILES_STRUCTURE_API_HEALTH_URL="${SMILES_STRUCTURE_API_HEALTH_URL:-${SMILES_STRUCTURE_API_BASE_URL}/health}"
SMILES_STRUCTURE_API_STATS_URL="${SMILES_STRUCTURE_API_STATS_URL:-${SMILES_STRUCTURE_API_BASE_URL}/v1/structure/stats}"
SMILES_STRUCTURE_API_TIMEOUT="${SMILES_STRUCTURE_API_TIMEOUT:-8}"
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
    print_header "SMILES Structure API Check"
    print_step "Checking ${SMILES_STRUCTURE_API_HEALTH_URL}"
fi

health_code="$(curl -sS -m "$SMILES_STRUCTURE_API_TIMEOUT" -o /tmp/smiles_api_health.json -w "%{http_code}" "$SMILES_STRUCTURE_API_HEALTH_URL" 2>/dev/null || true)"

if [ "$health_code" != "200" ] && command_exists docker && docker ps --format '{{.Names}}' | grep -qx "$SMILES_STRUCTURE_API_DOCKER_CONTAINER"; then
    # Fallback: probe from inside container if host port mapping is unavailable.
    health_code="$(docker exec "$SMILES_STRUCTURE_API_DOCKER_CONTAINER" sh -lc "curl -sS -m ${SMILES_STRUCTURE_API_TIMEOUT} -o /tmp/smiles_api_health.json -w '%{http_code}' http://127.0.0.1:8011/health" || true)"
fi

if [ "$health_code" != "200" ]; then
    if [ "$QUIET" = false ]; then
        print_error "Health endpoint failed (HTTP ${health_code})"
    fi
    exit 1
fi

stats_code="$(curl -sS -m "$SMILES_STRUCTURE_API_TIMEOUT" -o /tmp/smiles_api_stats.json -w "%{http_code}" "$SMILES_STRUCTURE_API_STATS_URL" 2>/dev/null || true)"

if [ "$stats_code" != "200" ] && command_exists docker && docker ps --format '{{.Names}}' | grep -qx "$SMILES_STRUCTURE_API_DOCKER_CONTAINER"; then
    stats_code="$(docker exec "$SMILES_STRUCTURE_API_DOCKER_CONTAINER" sh -lc "curl -sS -m ${SMILES_STRUCTURE_API_TIMEOUT} -o /tmp/smiles_api_stats.json -w '%{http_code}' http://127.0.0.1:8011/v1/structure/stats" || true)"
    if [ "$stats_code" = "200" ]; then
        docker cp "${SMILES_STRUCTURE_API_DOCKER_CONTAINER}:/tmp/smiles_api_stats.json" /tmp/smiles_api_stats.json >/dev/null 2>&1 || true
    fi
fi

if [ "$stats_code" != "200" ]; then
    if [ "$QUIET" = false ]; then
        print_warning "Stats endpoint failed (HTTP ${stats_code})"
        print_success "Health endpoint OK"
    fi
    exit 0
fi

if [ "$QUIET" = false ]; then
    print_success "SMILES API healthy"
    if command_exists jq; then
        print_info "Coverage: $(jq -r '.coverage // 0' /tmp/smiles_api_stats.json 2>/dev/null || echo "n/a")"
        print_info "Documents with fingerprints: $(jq -r '.documents_with_fingerprints // 0' /tmp/smiles_api_stats.json 2>/dev/null || echo "0")"
    fi
fi

exit 0
