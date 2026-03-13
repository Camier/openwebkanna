#!/usr/bin/env bash

###############################################################################
# Docker Compose Helper Library
# Provides reusable functions for detecting and invoking docker compose.
#
# Usage:
#   source /path/to/lib/docker-helpers.sh
#   init_compose_cmd  # Optional; docker_compose calls this automatically
#   docker_compose up -d
###############################################################################

# Global array to hold the compose command (e.g., "docker compose" or "docker-compose")
# Scripts should declare this as global if they want to use it after sourcing.
if [ -z "${COMPOSE_CMD+x}" ]; then
    COMPOSE_CMD=()
fi

# Default compose file (can be overridden by sourcing scripts)
if [ -z "${COMPOSE_FILE+x}" ]; then
    COMPOSE_FILE="${COMPOSE_FILE:-docker-compose.yml}"
fi

###############################################################################
# init_compose_cmd
# Detects and initializes the COMPOSE_CMD array.
# Prefers "docker compose" (v2 plugin) over "docker-compose" (v1 standalone).
#
# Globals:
#   COMPOSE_CMD - Array set to the detected compose command
#
# Returns:
#   0 on success, 1 if docker compose is not available
###############################################################################
init_compose_cmd() {
    # Already initialized
    if [ ${#COMPOSE_CMD[@]} -gt 0 ]; then
        return 0
    fi

    # Prefer Docker Compose v2 (plugin)
    if docker compose version >/dev/null 2>&1; then
        COMPOSE_CMD=(docker compose)
        return 0
    fi

    # Fallback to Docker Compose v1 (standalone)
    if command -v docker-compose >/dev/null 2>&1; then
        COMPOSE_CMD=(docker-compose)
        return 0
    fi

    echo "Error: Docker Compose is not available (plugin or docker-compose binary)" >&2
    return 1
}

###############################################################################
# docker_compose
# Wrapper function that invokes the detected compose command with the
# configured compose file.
#
# Globals:
#   COMPOSE_CMD  - Array containing the compose command
#   COMPOSE_FILE - Path to the compose file (default: docker-compose.yml)
#
# Arguments:
#   $@ - All arguments passed through to docker compose
#
# Returns:
#   Exit code from docker compose, or 1 if init fails
#
# Example:
#   docker_compose up -d
#   docker_compose logs -f openwebui
#   docker_compose exec -T openwebui python --version
###############################################################################
docker_compose() {
    init_compose_cmd || return 1
    "${COMPOSE_CMD[@]}" -f "$COMPOSE_FILE" "$@"
}

###############################################################################
# compose_service_running
# Check whether a compose service currently has a running container.
#
# Arguments:
#   $1 - Compose service name
#
# Returns:
#   0 if the service is running, 1 otherwise
###############################################################################
compose_service_running() {
    local service="$1"
    local services=""

    services="$(docker_compose ps --status running --services 2>/dev/null || true)"
    printf "%s\n" "$services" | grep -qx "$service"
}

###############################################################################
# resolve_container_id
# Resolve a running container ID for a compose service, with an optional
# fallback to a literal Docker container name for compatibility copies.
#
# Arguments:
#   $1 - Compose service name
#   $2 - Optional fallback container name (defaults to the service name)
#
# Output:
#   Container ID to stdout
#
# Returns:
#   0 if a running container ID was found, 1 otherwise
###############################################################################
resolve_container_id() {
    local service="$1"
    local fallback_name="${2:-$service}"
    local container_id=""

    container_id="$(docker_compose ps -q "$service" 2>/dev/null | tr -d '\r' | head -n 1 || true)"
    if [ -n "$container_id" ]; then
        printf "%s" "$container_id"
        return 0
    fi

    container_id="$(docker ps -q --filter "name=^/${fallback_name}$" 2>/dev/null | tr -d '\r' | head -n 1 || true)"
    if [ -n "$container_id" ]; then
        printf "%s" "$container_id"
        return 0
    fi

    return 1
}

###############################################################################
# resolve_host_python
# Resolve a usable Python interpreter on the host, preferring python3.
#
# Output:
#   Interpreter name to stdout
#
# Returns:
#   0 if a supported interpreter was found, 1 otherwise
###############################################################################
resolve_host_python() {
    if command -v python3 >/dev/null 2>&1; then
        printf "python3"
        return 0
    fi

    if command -v python >/dev/null 2>&1; then
        printf "python"
        return 0
    fi

    return 1
}

###############################################################################
# resolve_container_python
# Resolve a usable Python interpreter inside a running container, preferring
# python3 over python.
#
# Arguments:
#   $1 - Running container ID or name
#
# Output:
#   Interpreter name to stdout
#
# Returns:
#   0 if a supported interpreter was found, 1 otherwise
###############################################################################
resolve_container_python() {
    local container_id="$1"
    local python_bin=""

    python_bin="$(docker exec "$container_id" sh -lc 'if command -v python3 >/dev/null 2>&1; then printf python3; elif command -v python >/dev/null 2>&1; then printf python; fi' 2>/dev/null || true)"
    if [ -z "$python_bin" ]; then
        return 1
    fi

    printf "%s" "$python_bin"
}
