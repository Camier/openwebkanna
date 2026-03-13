#!/bin/bash

###############################################################################
# Stop Indigo Service sidecar (Docker Compose profile).
###############################################################################

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

INDIGO_SERVICE_DOCKER_SERVICE="${INDIGO_SERVICE_DOCKER_SERVICE:-indigo-service}"

main() {
    print_header "Stop Indigo Service"

    if ! command_exists docker; then
        print_info "Docker not installed; nothing to stop"
        return 0
    fi

    if ! docker info >/dev/null 2>&1; then
        print_info "Docker daemon not running; nothing to stop"
        return 0
    fi

    if ! init_compose_cmd; then
        exit 1
    fi

    print_step "Stopping ${INDIGO_SERVICE_DOCKER_SERVICE}"
    docker_compose --profile indigo-service stop "$INDIGO_SERVICE_DOCKER_SERVICE" >/dev/null 2>&1 || true
    docker_compose --profile indigo-service rm -f -s "$INDIGO_SERVICE_DOCKER_SERVICE" >/dev/null 2>&1 || true

    print_success "Indigo Service stopped"
}

main "$@"
