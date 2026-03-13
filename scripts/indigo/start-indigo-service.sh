#!/bin/bash

###############################################################################
# Start Indigo Service sidecar (Docker Compose profile).
###############################################################################

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

INDIGO_SERVICE_DOCKER_SERVICE="${INDIGO_SERVICE_DOCKER_SERVICE:-indigo-service}"
INDIGO_SERVICE_PORT="${INDIGO_SERVICE_PORT:-8012}"
INDIGO_SERVICE_BIND_ADDRESS="${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}"
INDIGO_SERVICE_BASE_URL="${INDIGO_SERVICE_BASE_URL:-http://${INDIGO_SERVICE_BIND_ADDRESS}:${INDIGO_SERVICE_PORT}}"
INDIGO_SERVICE_STARTUP_RETRIES="${INDIGO_SERVICE_STARTUP_RETRIES:-60}"

main() {
    print_header "Start Indigo Service"

    if ! command_exists docker; then
        print_error "Docker is required"
        exit 1
    fi

    if ! docker info >/dev/null 2>&1; then
        print_error "Docker daemon is not running"
        exit 1
    fi

    if ! init_compose_cmd; then
        exit 1
    fi

    if ! is_true "${INDIGO_SERVICE_ENABLED:-false}"; then
        print_warning "INDIGO_SERVICE_ENABLED=false in .env (service will still be started for this session)"
    fi

    print_step "Starting ${INDIGO_SERVICE_DOCKER_SERVICE} via Compose profile"
    docker_compose --profile indigo-service up -d "$INDIGO_SERVICE_DOCKER_SERVICE"

    print_step "Waiting for Indigo Service health"
    local attempt=1
    while [ "$attempt" -le "$INDIGO_SERVICE_STARTUP_RETRIES" ]; do
        if "${SELF_DIR}/check-indigo-service.sh" --quiet >/dev/null 2>&1; then
            print_success "Indigo Service is healthy"
            print_info "URL: ${INDIGO_SERVICE_BASE_URL%/}"
            return 0
        fi
        sleep 1
        attempt=$((attempt + 1))
    done

    print_error "Indigo Service did not become healthy"
    docker_compose --profile indigo-service logs --tail 120 "$INDIGO_SERVICE_DOCKER_SERVICE" || true
    exit 1
}

main "$@"
