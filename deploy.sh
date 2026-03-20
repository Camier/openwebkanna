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
COMPOSE_FILE="config/compose/docker-compose.yml"
MCPO_BASE_URL="${MCPO_BASE_URL:-http://127.0.0.1:${MCPO_PORT:-8000}}"
MULTIMODAL_RETRIEVAL_API_URL="${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}"
MULTIMODAL_RETRIEVAL_API_HOST="${MULTIMODAL_RETRIEVAL_API_HOST:-127.0.0.1}"
MULTIMODAL_RETRIEVAL_API_PORT="${MULTIMODAL_RETRIEVAL_API_PORT:-8510}"
MULTIMODAL_RETRIEVAL_API_STARTUP_TIMEOUT_SECONDS="${MULTIMODAL_RETRIEVAL_API_STARTUP_TIMEOUT_SECONDS:-240}"
MULTIMODAL_RETRIEVAL_API_STATE_DIR="${MULTIMODAL_RETRIEVAL_API_STATE_DIR:-${SCRIPT_DIR}/.codex-state/multimodal_retrieval_api}"
MULTIMODAL_RETRIEVAL_API_PID_FILE="${MULTIMODAL_RETRIEVAL_API_PID_FILE:-${MULTIMODAL_RETRIEVAL_API_STATE_DIR}/multimodal_retrieval_api.pid}"
MULTIMODAL_RETRIEVAL_API_LOG_FILE="${MULTIMODAL_RETRIEVAL_API_LOG_FILE:-${SCRIPT_DIR}/logs/multimodal_retrieval_api.log}"
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

read_env_value_from_file() {
    local path="$1"
    local key="$2"
    if [ -z "$path" ] || [ ! -f "$path" ]; then
        return 1
    fi

    python3 - "$path" "$key" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
key = sys.argv[2]
for raw_line in path.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip()
    if not line or line.startswith("#"):
        continue
    if line.startswith("export "):
        line = line[len("export ") :].strip()
    current_key, separator, value = line.partition("=")
    if not separator or current_key.strip() != key:
        continue
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1]
    print(value)
    break
PY
}

resolve_multimodal_retrieval_text_query_model_path() {
    local text_query_model_path="${MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH:-}"
    if [ -n "$text_query_model_path" ]; then
        printf '%s\n' "$text_query_model_path"
        return 0
    fi

    local compat_env_file="${MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE:-}"
    if [ -n "$compat_env_file" ]; then
        read_env_value_from_file "$compat_env_file" "NEMOTRON_MODEL_PATH"
        return 0
    fi

    return 0
}

ensure_multimodal_retrieval_runtime_config() {
    local text_query_model_path
    text_query_model_path="$(resolve_multimodal_retrieval_text_query_model_path)"

    if [ -z "$text_query_model_path" ]; then
        print_error "Canonical retrieval startup requires MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH or MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE with NEMOTRON_MODEL_PATH."
        exit 1
    fi

    if [ ! -e "$text_query_model_path" ]; then
        print_error "Canonical retrieval text query model path does not exist: $text_query_model_path"
        exit 1
    fi

    if [ -n "${MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE:-}" ] && [ ! -f "${MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE}" ]; then
        print_error "MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE does not exist: ${MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE}"
        exit 1
    fi

    mkdir -p "${MULTIMODAL_RETRIEVAL_API_STATE_DIR}" "$(dirname "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}")"
}

multimodal_retrieval_health_ok() {
    local base_url="${MULTIMODAL_RETRIEVAL_API_URL%/}"
    local health_response
    health_response="$(curl -sS -m 5 "${base_url}/health" 2>/dev/null || true)"
    echo "$health_response" | jq -e '.status == "ok"' >/dev/null 2>&1
}

multimodal_retrieval_ready_ok() {
    local base_url="${MULTIMODAL_RETRIEVAL_API_URL%/}"
    local ready_payload
    ready_payload="$(mktemp)"
    local ready_http
    ready_http="$(curl -sS -m 8 -o "$ready_payload" -w "%{http_code}" "${base_url}/ready" 2>/dev/null || true)"
    if [ "$ready_http" = "200" ] && jq -e '.status == "ready"' "$ready_payload" >/dev/null 2>&1; then
        rm -f "$ready_payload"
        return 0
    fi
    rm -f "$ready_payload"
    return 1
}

stop_managed_multimodal_retrieval_process() {
    if [ ! -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}" ]; then
        return 0
    fi

    local pid
    pid="$(cat "${MULTIMODAL_RETRIEVAL_API_PID_FILE}" 2>/dev/null || true)"
    if [ -z "$pid" ]; then
        rm -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}"
        return 0
    fi

    if ! kill -0 "$pid" 2>/dev/null; then
        rm -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}"
        return 0
    fi

    print_info "Stopping managed canonical retrieval process (pid=${pid})"
    kill "$pid" 2>/dev/null || true
    for _ in $(seq 1 20); do
        if ! kill -0 "$pid" 2>/dev/null; then
            rm -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}"
            return 0
        fi
        sleep 1
    done

    print_warning "Canonical retrieval process did not exit cleanly; sending SIGKILL"
    kill -9 "$pid" 2>/dev/null || true
    rm -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}"
}

fail_if_unmanaged_multimodal_retrieval_listener() {
    if multimodal_retrieval_ready_ok; then
        print_error "Canonical retrieval endpoint ${MULTIMODAL_RETRIEVAL_API_URL} is already serving /ready from an unmanaged process. Stop that process or set MULTIMODAL_RETRIEVAL_API_URL to the intended endpoint before deploy."
        exit 1
    fi

    if multimodal_retrieval_health_ok; then
        print_error "Canonical retrieval endpoint ${MULTIMODAL_RETRIEVAL_API_URL} is already serving /health from an unmanaged process. Stop that process or set MULTIMODAL_RETRIEVAL_API_URL to the intended endpoint before deploy."
        exit 1
    fi
}

launch_multimodal_retrieval_api() {
    ensure_multimodal_retrieval_runtime_config

    stop_managed_multimodal_retrieval_process
    fail_if_unmanaged_multimodal_retrieval_listener

    print_step "Starting canonical retrieval service"
    setsid env \
        PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}" \
        python3 -m uvicorn \
        --app-dir "$SCRIPT_DIR" \
        services.multimodal_retrieval_api.app:app \
        --host "${MULTIMODAL_RETRIEVAL_API_HOST}" \
        --port "${MULTIMODAL_RETRIEVAL_API_PORT}" \
        </dev/null >>"${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" 2>&1 &
    local pid=$!
    printf '%s\n' "$pid" >"${MULTIMODAL_RETRIEVAL_API_PID_FILE}"
    disown "$pid" 2>/dev/null || true
    print_info "Canonical retrieval log: ${MULTIMODAL_RETRIEVAL_API_LOG_FILE}"
}

wait_for_multimodal_retrieval_api() {
    local base_url="${MULTIMODAL_RETRIEVAL_API_URL%/}"
    local attempts
    attempts=$((MULTIMODAL_RETRIEVAL_API_STARTUP_TIMEOUT_SECONDS / 2))
    if [ "$attempts" -lt 1 ]; then
        attempts=1
    fi

    print_step "Waiting for canonical retrieval service to be ready"
    for _ in $(seq 1 "$attempts"); do
        if multimodal_retrieval_ready_ok; then
            print_success "Canonical retrieval service is ready"
            return 0
        fi

        if [ -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}" ]; then
            local pid
            pid="$(cat "${MULTIMODAL_RETRIEVAL_API_PID_FILE}" 2>/dev/null || true)"
            if [ -n "$pid" ] && ! kill -0 "$pid" 2>/dev/null; then
                print_error "Canonical retrieval service exited during startup. Recent log output:"
                tail -n 40 "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" 2>/dev/null || true
                rm -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}"
                return 1
            fi
        fi

        echo -n "."
        sleep 2
    done

    echo
    print_error "Canonical retrieval service did not become ready at ${base_url}/ready within ${MULTIMODAL_RETRIEVAL_API_STARTUP_TIMEOUT_SECONDS}s"
    if [ -f "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" ]; then
        print_info "Recent canonical retrieval log output:"
        tail -n 40 "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" 2>/dev/null || true
    fi
    return 1
}

check_docker_runtime() {
    print_step "Running Docker runtime preflight"

    if ! docker info &>/dev/null; then
        print_error "Docker daemon is not reachable."
        exit 1
    fi

    local cgroup_driver
    local cgroup_version
    cgroup_driver="$(docker info --format '{{.CgroupDriver}}' 2>/dev/null || true)"
    cgroup_version="$(docker info --format '{{.CgroupVersion}}' 2>/dev/null || true)"

    if [ -z "$cgroup_driver" ] || [ -z "$cgroup_version" ]; then
        print_error "Unable to read Docker cgroup metadata from 'docker info'."
        exit 1
    fi

    print_info "Docker CgroupDriver: ${cgroup_driver}"
    print_info "Docker CgroupVersion: ${cgroup_version}"

    local probe_output
    local probe_name="openwebui-pref-${RANDOM:-0}"
    if ! probe_output="$(docker run --rm --name "$probe_name" --network none busybox:latest true 2>&1)"; then
        if printf '%s\n' "$probe_output" | grep -Eq "unable to apply cgroup configuration|failed to create task|Launch helper"; then
            print_error "Docker runtime preflight failed due to container cgroup setup (hard blocker for OpenWebUI/Jupyter startup)."
            echo
            echo "Recommended host fix:"
            echo "  - set Docker cgroup driver to cgroupfs in /etc/docker/daemon.json:"
            echo '    "exec-opts": ["native.cgroupdriver=cgroupfs"]'
            echo "  - restart docker daemon: sudo systemctl restart docker"
            echo
            print_error "Keep the default runtime value (nvidia/others) unless your CUDA toolchain requires a different value."
            echo
            print_error "$probe_output"
            exit 1
        fi

        print_error "Docker runtime preflight container test failed."
        print_error "$probe_output"
        echo
        print_error "Common causes: corrupted Docker image cache, daemon not healthy, or container runtime pull/start issues."
        exit 1
    fi

    print_success "Docker runtime preflight passed"
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

sync_openwebui_retrieval_config() {
    if [ ! -x "./scripts/admin/sync-openwebui-retrieval-config.sh" ]; then
        print_warning "Retrieval sync script missing or not executable: ./scripts/admin/sync-openwebui-retrieval-config.sh"
        return 0
    fi

    print_step "Syncing OpenWebUI retrieval config"
    ./scripts/admin/sync-openwebui-retrieval-config.sh
    print_success "OpenWebUI retrieval config synced"
}

show_access_info() {
    local indigo_url="${INDIGO_SERVICE_BASE_URL:-http://${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}:${INDIGO_SERVICE_PORT:-8012}}"

    echo
    print_success "Deployment complete"
    echo
    echo -e "  ${CYAN}OpenWebUI:${NC}  http://localhost:${OPENWEBUI_PORT}"
    echo -e "  ${CYAN}MCPO:${NC}       ${MCPO_BASE_URL}"
    echo -e "  ${CYAN}Canonical RAG:${NC} ${MULTIMODAL_RETRIEVAL_API_URL}"
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
    echo -e "  ${YELLOW}Canonical log:${NC} ${MULTIMODAL_RETRIEVAL_API_LOG_FILE}"
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
    check_docker_runtime
    check_command "curl"
    check_command "jq"
    check_command "python3"
    validate_install_env
    print_success "Prerequisites met"

    if [ "$SKIP_DOCKER" = false ]; then
        start_docker_compose
    else
        print_info "Skipping docker compose (--skip-docker)"
    fi

    launch_multimodal_retrieval_api

    local openwebui_ready=false
    if wait_for_openwebui; then
        openwebui_ready=true
        sync_openwebui_retrieval_config || exit 1
        sync_openwebui_web_search_config || exit 1
    fi

    if [ "$openwebui_ready" = false ]; then
        print_warning "Skipping retrieval config sync steps because OpenWebUI is not ready yet"
    fi

    run_plugin_audit_gate || exit 1
    configure_indigo_tool || exit 1
    wait_for_multimodal_retrieval_api || exit 1

    if [ "$SHOW_LOGS" = true ]; then
        print_step "Recent logs"
        docker_compose logs --tail=20
    fi

    docker_compose ps
    show_access_info
}

main "$@"
