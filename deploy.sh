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
PLUGIN_AUDIT_SCRIPT="${PLUGIN_AUDIT_SCRIPT:-./scripts/testing/audit-openwebui-plugins.sh}"
PLUGIN_AUDIT_FOCUS="${PLUGIN_AUDIT_FOCUS:-all}"
INDIGO_TOOL_AUTOCONFIGURE="${INDIGO_TOOL_AUTOCONFIGURE:-true}"
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

require_nonempty_env() {
    local var_name="$1"
    local help_text="$2"
    local value="${!var_name:-}"
    if [ -n "$value" ]; then
        return 0
    fi

    print_error "${var_name} is required. ${help_text}"
    exit 1
}

validate_install_env() {
    print_step "Validating required environment"

    require_nonempty_env "WEBUI_SECRET_KEY" "Set a stable secret in .env before deploy."
    if [ "${#WEBUI_SECRET_KEY}" -lt 32 ]; then
        print_error "WEBUI_SECRET_KEY must be at least 32 characters for stable sessions and secure cookie signing."
        exit 1
    fi

    require_nonempty_env "JUPYTER_TOKEN" "Set a non-empty token in .env before deploy."
    require_nonempty_env "POSTGRES_PASSWORD" "Replace the placeholder value in .env before deploy."
    require_nonempty_env "OPENAI_API_KEY" "Set the LiteLLM master key in .env before deploy."

    if [ "${POSTGRES_PASSWORD}" = "change-me-strong-password" ]; then
        print_error "POSTGRES_PASSWORD still uses the example placeholder. Replace it in .env before deploy."
        exit 1
    fi

    if [ "${OPENAI_API_KEY}" = "sk-litellm-replace-with-your-master-key" ]; then
        print_error "OPENAI_API_KEY still uses the example placeholder. Replace it with your LiteLLM master key in .env before deploy."
        exit 1
    fi

    if is_true "${ENABLE_CODE_EXECUTION:-true}" &&
        [ "${CODE_EXECUTION_ENGINE:-jupyter}" = "jupyter" ] &&
        [ "${CODE_EXECUTION_JUPYTER_AUTH:-token}" = "token" ] &&
        [ "${CODE_EXECUTION_JUPYTER_AUTH_TOKEN:-}" != "${JUPYTER_TOKEN}" ]; then
        print_error "CODE_EXECUTION_JUPYTER_AUTH_TOKEN must exactly match JUPYTER_TOKEN for Jupyter-backed code execution."
        exit 1
    fi

    if is_true "${ENABLE_CODE_INTERPRETER:-true}" &&
        [ "${CODE_INTERPRETER_ENGINE:-jupyter}" = "jupyter" ] &&
        [ "${CODE_INTERPRETER_JUPYTER_AUTH:-token}" = "token" ] &&
        [ "${CODE_INTERPRETER_JUPYTER_AUTH_TOKEN:-}" != "${JUPYTER_TOKEN}" ]; then
        print_error "CODE_INTERPRETER_JUPYTER_AUTH_TOKEN must exactly match JUPYTER_TOKEN for Jupyter-backed interpreter mode."
        exit 1
    fi

    if ! docker_compose config >/dev/null; then
        print_error "docker compose config failed. Fix the .env values above before deploy."
        exit 1
    fi

    print_success "Environment validation passed"
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

    # CLIProxyAPI is a legacy sidecar behind an explicit flag/profile.
    if is_true "${CLIPROXYAPI_ENABLED:-false}"; then
        print_step "CLIPROXYAPI enabled; starting legacy sidecar profile"
        docker_compose --profile legacy-cliproxy up -d cliproxyapi
    else
        print_info "CLIPROXYAPI disabled; ensuring legacy sidecar is not running"
        docker_compose --profile legacy-cliproxy stop cliproxyapi >/dev/null 2>&1 || true
        docker_compose --profile legacy-cliproxy rm -f -s cliproxyapi >/dev/null 2>&1 || true
    fi

    # SearXNG is optional and only needed when OpenWebUI web search is enabled.
    if is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; then
        print_step "Web search enabled; starting web-search profile"
        docker_compose --profile web-search up -d searxng
    else
        print_info "Web search disabled; ensuring searxng sidecar is not running"
        docker_compose --profile web-search stop searxng >/dev/null 2>&1 || true
        docker_compose --profile web-search rm -f -s searxng >/dev/null 2>&1 || true
    fi

    # Open Terminal sidecar is optional and explicitly opt-in.
    if is_true "${OPEN_TERMINAL_ENABLED:-false}"; then
        if [ -z "${OPEN_TERMINAL_API_KEY:-}" ]; then
            print_error "OPEN_TERMINAL_ENABLED=true but OPEN_TERMINAL_API_KEY is empty. Set it in .env before deploy."
            exit 1
        fi
        print_step "OPEN_TERMINAL enabled; starting open-terminal profile"
        docker_compose --profile open-terminal up -d open-terminal
    else
        print_info "OPEN_TERMINAL disabled; ensuring open-terminal sidecar is not running"
        docker_compose --profile open-terminal stop open-terminal >/dev/null 2>&1 || true
        docker_compose --profile open-terminal rm -f -s open-terminal >/dev/null 2>&1 || true
    fi

    # Indigo Service sidecar is optional and explicitly opt-in.
    if is_true "${INDIGO_SERVICE_ENABLED:-false}"; then
        print_step "INDIGO_SERVICE enabled; starting indigo-service profile"
        docker_compose --profile indigo-service up -d indigo-service
    else
        print_info "INDIGO_SERVICE disabled; ensuring indigo-service sidecar is not running"
        docker_compose --profile indigo-service stop indigo-service >/dev/null 2>&1 || true
        docker_compose --profile indigo-service rm -f -s indigo-service >/dev/null 2>&1 || true
    fi

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

configure_indigo_tool() {
    if ! is_true "${INDIGO_SERVICE_ENABLED:-false}"; then
        print_info "Indigo tool setup skipped (INDIGO_SERVICE_ENABLED=false)"
        return 0
    fi

    if ! is_true "${INDIGO_TOOL_AUTOCONFIGURE:-true}"; then
        print_info "Indigo tool setup skipped (INDIGO_TOOL_AUTOCONFIGURE=false)"
        return 0
    fi

    if [ ! -x "./scripts/indigo/enable-indigo-live.sh" ]; then
        print_error "Indigo tool setup script is missing or not executable: ./scripts/indigo/enable-indigo-live.sh"
        return 1
    fi

    print_step "Provisioning Indigo tool in OpenWebUI"
    ./scripts/indigo/enable-indigo-live.sh
    print_success "Indigo tool provisioned"
}

sync_openwebui_web_search_config() {
    if [ ! -x "./scripts/admin/sync-openwebui-web-search-config.sh" ]; then
        print_warning "Web-search sync script missing or not executable: ./scripts/admin/sync-openwebui-web-search-config.sh"
        return 0
    fi

    print_step "Syncing OpenWebUI retrieval web-search config"
    ./scripts/admin/sync-openwebui-web-search-config.sh
    print_success "OpenWebUI retrieval web-search config synced"
}

show_access_info() {
    local indigo_url="${INDIGO_SERVICE_BASE_URL:-http://${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}:${INDIGO_SERVICE_PORT:-8012}}"

    echo
    print_success "Deployment complete"
    echo
    echo -e "  ${CYAN}OpenWebUI:${NC}  http://localhost:${OPENWEBUI_PORT}"
    echo -e "  ${CYAN}MCPO:${NC}       ${MCPO_BASE_URL}"
    if is_true "${INDIGO_SERVICE_ENABLED:-false}"; then
        echo -e "  ${CYAN}Indigo:${NC}     ${indigo_url}"
        if is_true "${INDIGO_TOOL_AUTOCONFIGURE:-true}"; then
            echo -e "  ${CYAN}Tool:${NC}       indigo_chemistry (auto-provisioned)"
        fi
    fi
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
                echo "  INDIGO_TOOL_AUTOCONFIGURE  Auto-register Indigo tool when INDIGO_SERVICE_ENABLED=true (default: true)"
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
    validate_install_env
    print_success "Prerequisites met"

    if [ "$SKIP_DOCKER" = false ]; then
        start_docker_compose
    else
        print_info "Skipping docker compose (--skip-docker)"
    fi

    local openwebui_ready=false
    if wait_for_openwebui; then
        openwebui_ready=true
        sync_openwebui_web_search_config || exit 1
    fi

    if [ "$openwebui_ready" = false ]; then
        print_warning "Skipping retrieval web-search sync because OpenWebUI is not ready yet"
    fi

    run_plugin_audit_gate || exit 1
    configure_indigo_tool || exit 1

    if [ "$SHOW_LOGS" = true ]; then
        print_step "Recent logs"
        docker_compose logs --tail=20
    fi

    docker_compose ps
    show_access_info
}

main "$@"
