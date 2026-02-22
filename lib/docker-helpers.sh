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
