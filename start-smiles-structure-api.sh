#!/bin/bash

###############################################################################
# Start SMILES Structure Search API
# Mode auto:
# - Prefer Docker Compose service `smiles-structure-api` when available
# - Fallback to host process (micromamba/conda/system python)
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

SMILES_STRUCTURE_API_MODE="${SMILES_STRUCTURE_API_MODE:-auto}"
SMILES_STRUCTURE_API_DOCKER_SERVICE="${SMILES_STRUCTURE_API_DOCKER_SERVICE:-smiles-structure-api}"
SMILES_STRUCTURE_API_DOCKER_CONTAINER="${SMILES_STRUCTURE_API_DOCKER_CONTAINER:-smiles_structure_api}"

SMILES_STRUCTURE_API_HOST="${SMILES_STRUCTURE_API_HOST:-0.0.0.0}"
SMILES_STRUCTURE_API_PORT="${SMILES_STRUCTURE_API_PORT:-8011}"
SMILES_STRUCTURE_API_WORKDIR="${SMILES_STRUCTURE_API_WORKDIR:-${SCRIPT_DIR}/smiles-pipeline/src}"
SMILES_STRUCTURE_API_LOG_FILE="${SMILES_STRUCTURE_API_LOG_FILE:-${SCRIPT_DIR}/logs/smiles-structure-api.log}"
SMILES_STRUCTURE_API_PID_FILE="${SMILES_STRUCTURE_API_PID_FILE:-${SCRIPT_DIR}/.smiles-structure-api.pid}"
SMILES_STRUCTURE_API_DB_URL="${SMILES_STRUCTURE_API_DB_URL:-}"
SMILES_STRUCTURE_API_DB_HOST="${SMILES_STRUCTURE_API_DB_HOST:-}"
SMILES_STRUCTURE_API_DB_PORT="${SMILES_STRUCTURE_API_DB_PORT:-}"
SMILES_STRUCTURE_API_POSTGRES_CONTAINER="${SMILES_STRUCTURE_API_POSTGRES_CONTAINER:-openwebui_postgres}"
SMILES_STRUCTURE_API_HEALTH_URL="${SMILES_STRUCTURE_API_HEALTH_URL:-http://127.0.0.1:${SMILES_STRUCTURE_API_PORT}/health}"
SMILES_STRUCTURE_API_STARTUP_RETRIES="${SMILES_STRUCTURE_API_STARTUP_RETRIES:-90}"

require_dep() {
    local cmd="$1"
    if ! command_exists "$cmd"; then
        print_error "Missing dependency: $cmd"
        exit 1
    fi
}

pid_is_running() {
    local pid="$1"
    kill -0 "$pid" 2>/dev/null
}

is_port_in_use() {
    local port="$1"
    if command_exists lsof; then
        lsof -nP -iTCP:"$port" -sTCP:LISTEN >/dev/null 2>&1
        return $?
    fi
    return 1
}

resolve_postgres_endpoint() {
    if [ -n "$SMILES_STRUCTURE_API_DB_HOST" ]; then
        local host="$SMILES_STRUCTURE_API_DB_HOST"
        local port="${SMILES_STRUCTURE_API_DB_PORT:-5432}"
        echo "${host}:${port}"
        return
    fi

    if command_exists docker && docker ps --format '{{.Names}}' | grep -qx "$SMILES_STRUCTURE_API_POSTGRES_CONTAINER"; then
        local pg_ip
        pg_ip="$(docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}} {{end}}' "$SMILES_STRUCTURE_API_POSTGRES_CONTAINER" 2>/dev/null | awk '{print $1}')"
        if [ -n "$pg_ip" ]; then
            echo "${pg_ip}:5432"
            return
        fi
    fi

    echo "127.0.0.1:${POSTGRES_PORT:-5432}"
}

build_db_url() {
    if [ -n "$SMILES_STRUCTURE_API_DB_URL" ]; then
        echo "$SMILES_STRUCTURE_API_DB_URL"
        return
    fi

    local endpoint db_host db_port
    endpoint="$(resolve_postgres_endpoint)"
    db_host="${endpoint%:*}"
    db_port="${endpoint##*:}"
    echo "postgresql://${POSTGRES_USER:-openwebui}:${POSTGRES_PASSWORD:-openwebui_pgvector_pass_2026}@${db_host}:${db_port}/${POSTGRES_DB:-openwebui}"
}

has_docker_service() {
    if ! command_exists docker; then
        return 1
    fi
    if ! init_compose_cmd; then
        return 1
    fi
    docker_compose config --services 2>/dev/null | grep -qx "$SMILES_STRUCTURE_API_DOCKER_SERVICE"
}

start_docker_mode() {
    require_dep docker

    print_step "Starting Docker service: ${SMILES_STRUCTURE_API_DOCKER_SERVICE}"
    docker_compose up -d "$SMILES_STRUCTURE_API_DOCKER_SERVICE" >/dev/null

    local attempt=1
    while [ "$attempt" -le "$SMILES_STRUCTURE_API_STARTUP_RETRIES" ]; do
        local cid
        cid="$(docker ps --filter "name=^/${SMILES_STRUCTURE_API_DOCKER_CONTAINER}$" --format '{{.ID}}' | head -n 1)"
        if [ -n "$cid" ] && docker exec "$SMILES_STRUCTURE_API_DOCKER_CONTAINER" sh -lc "curl -fsS -m 3 http://127.0.0.1:8011/health >/dev/null" >/dev/null 2>&1; then
            print_success "SMILES API running via Docker service (${SMILES_STRUCTURE_API_DOCKER_CONTAINER})"
            print_info "Health endpoint: ${SMILES_STRUCTURE_API_HEALTH_URL}"
            return
        fi
        sleep 1
        attempt=$((attempt + 1))
    done

    print_error "SMILES Docker service failed to become healthy"
    docker_compose logs --tail 80 "$SMILES_STRUCTURE_API_DOCKER_SERVICE" || true
    exit 1
}

start_host_mode() {
    require_dep curl
    require_dep python

    if [ -f "$SMILES_STRUCTURE_API_PID_FILE" ]; then
        local existing_pid
        existing_pid="$(cat "$SMILES_STRUCTURE_API_PID_FILE" 2>/dev/null || true)"
        if [[ $existing_pid =~ ^[0-9]+$ ]] && pid_is_running "$existing_pid"; then
            print_info "SMILES API already running (PID: $existing_pid)"
            print_info "Health: $SMILES_STRUCTURE_API_HEALTH_URL"
            return
        fi
        rm -f "$SMILES_STRUCTURE_API_PID_FILE"
    fi

    if is_port_in_use "$SMILES_STRUCTURE_API_PORT"; then
        print_error "Port ${SMILES_STRUCTURE_API_PORT} already in use"
        exit 1
    fi

    mkdir -p "$(dirname "$SMILES_STRUCTURE_API_LOG_FILE")"
    print_step "Starting SMILES API on ${SMILES_STRUCTURE_API_HOST}:${SMILES_STRUCTURE_API_PORT} (host mode)"

    local runner=""
    if command_exists micromamba; then
        runner="micromamba run -n smiles-extraction python -m uvicorn api.structure_search:app --host ${SMILES_STRUCTURE_API_HOST} --port ${SMILES_STRUCTURE_API_PORT}"
    elif command_exists conda; then
        runner="bash -lc 'source \"$(conda info --base)/etc/profile.d/conda.sh\" && conda activate smiles-extraction && python -m uvicorn api.structure_search:app --host ${SMILES_STRUCTURE_API_HOST} --port ${SMILES_STRUCTURE_API_PORT}'"
    else
        print_warning "micromamba/conda not found; using system python"
        runner="python -m uvicorn api.structure_search:app --host ${SMILES_STRUCTURE_API_HOST} --port ${SMILES_STRUCTURE_API_PORT}"
    fi

    local db_url db_url_masked
    db_url="$(build_db_url)"
    db_url_masked="$(printf '%s' "$db_url" | sed -E 's#(://[^:]+:)[^@]+@#\1***@#')"
    print_info "Using DATABASE_URL=${db_url_masked}"

    (
        cd "$SMILES_STRUCTURE_API_WORKDIR"
        export DATABASE_URL="$db_url"
        if [[ $runner == bash\ -lc* ]]; then
            nohup bash -lc "${runner#bash -lc }" >>"$SMILES_STRUCTURE_API_LOG_FILE" 2>&1 &
        else
            nohup bash -lc "$runner" >>"$SMILES_STRUCTURE_API_LOG_FILE" 2>&1 &
        fi
        echo $! >"$SMILES_STRUCTURE_API_PID_FILE"
    )

    local attempt=1
    while [ "$attempt" -le "$SMILES_STRUCTURE_API_STARTUP_RETRIES" ]; do
        if curl -fsS -m 3 "$SMILES_STRUCTURE_API_HEALTH_URL" >/dev/null 2>&1; then
            local pid
            pid="$(cat "$SMILES_STRUCTURE_API_PID_FILE" 2>/dev/null || true)"
            print_success "SMILES API running (PID: ${pid:-unknown})"
            print_info "Health endpoint: $SMILES_STRUCTURE_API_HEALTH_URL"
            return
        fi
        sleep 1
        attempt=$((attempt + 1))
    done

    print_error "SMILES API failed to become healthy (host mode)"
    print_info "Last log lines:"
    tail -n 40 "$SMILES_STRUCTURE_API_LOG_FILE" || true
    exit 1
}

main() {
    print_header "Start SMILES Structure API"

    local mode="$SMILES_STRUCTURE_API_MODE"
    if [ "$mode" = "auto" ]; then
        if has_docker_service; then
            mode="docker"
        else
            mode="host"
        fi
    fi

    if [ "$mode" = "docker" ]; then
        start_docker_mode
        return
    fi

    if [ "$mode" = "host" ]; then
        start_host_mode
        return
    fi

    print_error "Invalid SMILES_STRUCTURE_API_MODE='${SMILES_STRUCTURE_API_MODE}' (expected auto|docker|host)"
    exit 2
}

main "$@"
