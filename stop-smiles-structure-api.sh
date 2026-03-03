#!/bin/bash

###############################################################################
# Stop SMILES Structure Search API
# Stops Docker Compose service when available, then host fallback process.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

SMILES_STRUCTURE_API_DOCKER_SERVICE="${SMILES_STRUCTURE_API_DOCKER_SERVICE:-smiles-structure-api}"
SMILES_STRUCTURE_API_DOCKER_CONTAINER="${SMILES_STRUCTURE_API_DOCKER_CONTAINER:-smiles_structure_api}"
SMILES_STRUCTURE_API_PID_FILE="${SMILES_STRUCTURE_API_PID_FILE:-${SCRIPT_DIR}/.smiles-structure-api.pid}"
SMILES_STRUCTURE_API_PORT="${SMILES_STRUCTURE_API_PORT:-8011}"
SMILES_STRUCTURE_API_STOP_TIMEOUT="${SMILES_STRUCTURE_API_STOP_TIMEOUT:-20}"

pid_is_running() {
    local pid="$1"
    kill -0 "$pid" 2>/dev/null
}

stop_pid() {
    local pid="$1"
    local waited=0

    kill "$pid" 2>/dev/null || true

    while pid_is_running "$pid" && [ "$waited" -lt "$SMILES_STRUCTURE_API_STOP_TIMEOUT" ]; do
        sleep 1
        waited=$((waited + 1))
    done

    if pid_is_running "$pid"; then
        print_warning "Process still running after ${SMILES_STRUCTURE_API_STOP_TIMEOUT}s; sending SIGKILL"
        kill -9 "$pid" 2>/dev/null || true
    fi
}

stop_docker_service_if_present() {
    if ! command_exists docker; then
        return 0
    fi
    if ! init_compose_cmd; then
        return 0
    fi

    if ! docker_compose config --services 2>/dev/null | grep -qx "$SMILES_STRUCTURE_API_DOCKER_SERVICE"; then
        return 0
    fi

    if docker ps --format '{{.Names}}' | grep -qx "$SMILES_STRUCTURE_API_DOCKER_CONTAINER"; then
        print_step "Stopping Docker service ${SMILES_STRUCTURE_API_DOCKER_SERVICE}"
        docker_compose stop "$SMILES_STRUCTURE_API_DOCKER_SERVICE" >/dev/null || true
        print_success "Docker SMILES API stopped"
    fi
}

main() {
    print_header "Stop SMILES Structure API"

    local stopped=false

    stop_docker_service_if_present || true

    if [ -f "$SMILES_STRUCTURE_API_PID_FILE" ]; then
        local pid
        pid="$(cat "$SMILES_STRUCTURE_API_PID_FILE" 2>/dev/null || true)"
        if [[ $pid =~ ^[0-9]+$ ]] && pid_is_running "$pid"; then
            print_step "Stopping PID ${pid}"
            stop_pid "$pid"
            stopped=true
        fi
        rm -f "$SMILES_STRUCTURE_API_PID_FILE"
    fi

    if [ "$stopped" = false ] && command_exists lsof; then
        local port_pid
        port_pid="$(lsof -nP -iTCP:"$SMILES_STRUCTURE_API_PORT" -sTCP:LISTEN -t 2>/dev/null | head -n 1 || true)"
        if [[ $port_pid =~ ^[0-9]+$ ]]; then
            print_step "Stopping listener on port ${SMILES_STRUCTURE_API_PORT} (PID ${port_pid})"
            stop_pid "$port_pid"
            stopped=true
        fi
    fi

    if [ "$stopped" = true ]; then
        print_success "Host SMILES API stopped"
    else
        print_info "No host SMILES API process running"
    fi
}

main "$@"
