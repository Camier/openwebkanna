#!/bin/bash

###############################################################################
# OpenWebUI RAG Deployment Script
# Deploys the complete RAG stack: OpenWebUI + Postgres + Jupyter + MCPO.
# LiteLLM proxy is the upstream OpenAI-compatible endpoint (external container).
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_PORT="${WEBUI_PORT:-3000}"
COMPOSE_FILE="docker-compose.yml"
MCPO_BASE_URL="${MCPO_BASE_URL:-http://127.0.0.1:${MCPO_PORT:-8000}}"
PLUGIN_AUDIT_ENABLED="${PLUGIN_AUDIT_ENABLED:-true}"
PLUGIN_AUDIT_SCRIPT="${PLUGIN_AUDIT_SCRIPT:-./audit-openwebui-plugins.sh}"
PLUGIN_AUDIT_FOCUS="${PLUGIN_AUDIT_FOCUS:-all}"
# shellcheck disable=SC2034
COMPOSE_CMD=()

###############################################################################
# Helpers
###############################################################################

check_command() {
    if ! command -v "$1" &>/dev/null; then
        print_error "$1 is not installed."
        exit 1
    fi
}

start_docker_compose() {
    print_step "Starting Docker Compose services"

    if [ ! -f "$COMPOSE_FILE" ]; then
        print_error "Docker Compose file not found: $COMPOSE_FILE"
        exit 1
    fi

    if ! docker info &>/dev/null; then
        print_error "Docker is not running."
        exit 1
    fi

    docker_compose up -d
    print_success "Docker Compose services started"
}

wait_for_openwebui() {
    local url="http://localhost:${OPENWEBUI_PORT}"
    print_step "Waiting for OpenWebUI to be ready"

    local attempts=60 i=1
    while [ $i -le $attempts ]; do
        if curl -sf "$url/health" &>/dev/null || curl -sf "$url" &>/dev/null; then
            print_success "OpenWebUI is ready"
            return 0
        fi
        echo -n "."
        sleep 2
        i=$((i + 1))
    done

    echo
    print_warning "OpenWebUI may still be starting. Check: ./logs.sh"
    return 1
}

run_plugin_audit_gate() {
    if [ "$PLUGIN_AUDIT_ENABLED" != "true" ]; then
        print_info "Plugin audit skipped (PLUGIN_AUDIT_ENABLED=${PLUGIN_AUDIT_ENABLED})"
        return 0
    fi

    if [ ! -x "$PLUGIN_AUDIT_SCRIPT" ]; then
        print_warning "Plugin audit script not found: $PLUGIN_AUDIT_SCRIPT"
        return 0
    fi

    print_step "Running OpenWebUI plugin syntax audit"
    if "$PLUGIN_AUDIT_SCRIPT" --focus "$PLUGIN_AUDIT_FOCUS"; then
        print_success "Plugin syntax audit passed"
        return 0
    fi

    print_error "Plugin syntax audit failed — fix tool/function syntax in webui.db before proceeding"
    return 1
}

show_access_info() {
    echo
    print_success "Deployment complete"
    echo
    echo -e "  ${CYAN}OpenWebUI:${NC}  http://localhost:${OPENWEBUI_PORT}"
    echo -e "  ${CYAN}MCPO:${NC}       ${MCPO_BASE_URL}"
    echo
    echo -e "  ${YELLOW}Logs:${NC}    ./logs.sh"
    echo -e "  ${YELLOW}Status:${NC}  ./status.sh"
    echo -e "  ${YELLOW}Stop:${NC}    ./cleanup.sh"
    echo
}

###############################################################################
# Main
###############################################################################

main() {
    print_header "OpenWebUI RAG Deployment"

    SKIP_DOCKER=false
    SHOW_LOGS=true

    while [[ $# -gt 0 ]]; do
        case $1 in
            --skip-docker)
                SKIP_DOCKER=true
                shift
                ;;
            --no-logs)
                SHOW_LOGS=false
                shift
                ;;
            -h | --help)
                echo "Usage: $0 [--skip-docker] [--no-logs]"
                echo
                echo "Options:"
                echo "  --skip-docker   Skip docker compose up (use if already running)"
                echo "  --no-logs       Don't show logs after deployment"
                echo
                echo "Env vars:"
                echo "  PLUGIN_AUDIT_ENABLED   Gate deploy on plugin syntax audit (default: true)"
                echo "  PLUGIN_AUDIT_FOCUS     Audit scope: all|tool|function (default: all)"
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                echo "Use -h for help"
                exit 1
                ;;
        esac
    done

    print_step "Checking prerequisites"
    check_command "docker"
    init_compose_cmd || exit 1
    check_command "curl"
    print_success "Prerequisites met"

    if [ "$SKIP_DOCKER" = false ]; then
        start_docker_compose
    else
        print_info "Skipping docker compose (--skip-docker)"
    fi

    wait_for_openwebui || true
    run_plugin_audit_gate || exit 1

    if [ "$SHOW_LOGS" = true ]; then
        print_step "Recent logs"
        docker_compose logs --tail=20
    fi

    docker_compose ps
    show_access_info
}

main "$@"
