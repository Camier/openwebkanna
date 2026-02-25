#!/bin/bash

###############################################################################
# refresh-models.sh - Refresh OpenWebUI model cache
#
# OpenWebUI caches models in memory (app.state.MODELS) on startup.
# This script gracefully restarts the container to pick up new models
# from connected providers (LiteLLM, Ollama, etc.)
#
# Usage:
#   ./refresh-models.sh [--wait]
#
# Options:
#   --wait    Wait for all services to be healthy before exiting
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/print-utils.sh"
source "${SCRIPT_DIR}/lib/docker-helpers.sh"
source "${SCRIPT_DIR}/lib/env-loader.sh"

load_env_defaults
cd "$SCRIPT_DIR"

# Configuration
WEBUI_PORT="${WEBUI_PORT:-3000}"
MAX_WAIT_SECONDS="${MAX_WAIT_SECONDS:-120}"
WAIT_FOR_HEALTH="${WAIT_FOR_HEALTH:-false}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-test@example.com}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --wait)
            WAIT_FOR_HEALTH=true
            shift
            ;;
        -h | --help)
            echo "Usage: $0 [--wait]"
            echo ""
            echo "Refreshes the OpenWebUI model cache by restarting the container."
            echo ""
            echo "Options:"
            echo "  --wait    Wait for all services to be healthy before exiting"
            exit 0
            ;;
        *)
            print_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

print_header "Refresh OpenWebUI Models"

###############################################################################
# get_model_count
# Returns the current number of models available via API.
###############################################################################
get_model_count() {
    local token count
    token=$(curl -sf "http://127.0.0.1:${WEBUI_PORT}/api/v1/auths/signin" \
        -H "Content-Type: application/json" \
        -d "{\"email\":\"${OPENWEBUI_SIGNIN_EMAIL}\",\"password\":\"${OPENWEBUI_SIGNIN_PASSWORD}\"}" 2>/dev/null | jq -r '.token' 2>/dev/null)

    if [[ -z $token || $token == "null" ]]; then
        echo "0"
        return
    fi

    count=$(curl -sf "http://127.0.0.1:${WEBUI_PORT}/api/v1/models" \
        -H "Authorization: Bearer $token" 2>/dev/null | jq '.data | length' 2>/dev/null)
    echo "${count:-0}"
}

###############################################################################
# wait_for_healthy
# Waits for OpenWebUI to become healthy.
###############################################################################
wait_for_healthy() {
    local elapsed=0
    local interval=3

    print_step "Waiting for OpenWebUI to be healthy..."

    while [[ $elapsed -lt $MAX_WAIT_SECONDS ]]; do
        if curl -sf "http://127.0.0.1:${WEBUI_PORT}/health" >/dev/null 2>&1; then
            print_success "OpenWebUI is healthy"
            return 0
        fi
        sleep $interval
        elapsed=$((elapsed + interval))
        echo -n "."
    done

    echo ""
    print_error "Timeout waiting for OpenWebUI to become healthy"
    return 1
}

###############################################################################
# main
###############################################################################
main() {
    local before_count after_count

    # Get current model count
    print_step "Checking current model count..."
    before_count=$(get_model_count)
    print_info "Current models: $before_count"

    # Restart OpenWebUI container
    print_step "Restarting OpenWebUI container..."
    if ! docker_compose restart openwebui; then
        print_error "Failed to restart OpenWebUI"
        exit 1
    fi
    print_success "Container restarted"

    # Wait for health if requested
    if [[ $WAIT_FOR_HEALTH == "true" ]]; then
        wait_for_healthy
    else
        # Brief wait for container to start
        sleep 5
    fi

    # Get new model count
    print_step "Checking new model count..."
    sleep 5 # Additional wait for models to load
    after_count=$(get_model_count)
    print_info "New models: $after_count"

    # Summary
    echo ""
    print_section "Summary"
    echo "  Before: $before_count models"
    echo "  After:  $after_count models"

    if [[ $after_count -gt $before_count ]]; then
        print_success "Gained $((after_count - before_count)) models"
    elif [[ $after_count -lt $before_count ]]; then
        print_warning "Lost $((before_count - after_count)) models"
    else
        print_info "Model count unchanged"
    fi

    print_success "Model cache refresh complete"
}

main "$@"
