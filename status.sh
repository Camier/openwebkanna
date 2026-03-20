#!/bin/bash

###############################################################################
# OpenWebUI RAG System Status
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_PORT="${WEBUI_PORT:-3000}"
# shellcheck disable=SC2034
COMPOSE_FILE="config/compose/docker-compose.yml"
MULTIMODAL_RETRIEVAL_API_URL="${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}"
MCPO_BASE_URL="${MCPO_BASE_URL:-http://127.0.0.1:${MCPO_PORT:-8000}}"
LITELLM_URL="${LITELLM_URL:-http://localhost:4000}"
OPEN_TERMINAL_ENABLED="${OPEN_TERMINAL_ENABLED:-false}"
OPEN_TERMINAL_BIND_ADDRESS="${OPEN_TERMINAL_BIND_ADDRESS:-127.0.0.1}"
OPEN_TERMINAL_PORT="${OPEN_TERMINAL_PORT:-8320}"
OPEN_TERMINAL_BASE_URL="${OPEN_TERMINAL_BASE_URL:-http://${OPEN_TERMINAL_BIND_ADDRESS}:${OPEN_TERMINAL_PORT}}"
INDIGO_SERVICE_ENABLED="${INDIGO_SERVICE_ENABLED:-false}"
INDIGO_SERVICE_BIND_ADDRESS="${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}"
INDIGO_SERVICE_PORT="${INDIGO_SERVICE_PORT:-8012}"
INDIGO_SERVICE_BASE_URL="${INDIGO_SERVICE_BASE_URL:-http://${INDIGO_SERVICE_BIND_ADDRESS}:${INDIGO_SERVICE_PORT}}"
INDIGO_TOOL_ENABLED="${INDIGO_TOOL_ENABLED:-true}"
INDIGO_TOOL_ID="${INDIGO_TOOL_ID:-indigo_chemistry}"
# shellcheck disable=SC2034
COMPOSE_CMD=()
LITELLM_PROBE_SOURCE=""
LITELLM_PROBE_URL=""

###############################################################################
# Section helpers (use lib print functions, no local redefinitions)
###############################################################################

print_svc() {
    local name="$1" state="$2" detail="${3:-}"
    case "$state" in
        running) echo -e "  ${GREEN}●${NC} ${BOLD}${name}${NC}: ${GREEN}Running${NC} ${detail}" ;;
        stopped) echo -e "  ${RED}○${NC} ${BOLD}${name}${NC}: ${RED}Stopped${NC} ${detail}" ;;
        warning) echo -e "  ${YELLOW}◐${NC} ${BOLD}${name}${NC}: ${YELLOW}Warning${NC} ${detail}" ;;
        *) echo -e "  ${YELLOW}?${NC} ${BOLD}${name}${NC}: ${YELLOW}Unknown${NC} ${detail}" ;;
    esac
}

###############################################################################
# Checks
###############################################################################

check_docker() {
    print_section "Docker"
    if ! docker info &>/dev/null 2>&1; then
        print_svc "Docker Daemon" "stopped"
        return 1
    fi
    print_svc "Docker Daemon" "running" "($(docker --version | sed 's/Docker version //'))"
}

check_docker_compose() {
    print_section "Docker Compose Services"

    if ! init_compose_cmd 2>/dev/null; then
        print_svc "Docker Compose" "stopped" "(compose command unavailable)"
        return 1
    fi

    local -a profile_args=()
    if is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; then
        profile_args+=(--profile web-search)
    fi
    if is_true "${OPEN_TERMINAL_ENABLED:-false}"; then
        profile_args+=(--profile open-terminal)
    fi
    if is_true "${INDIGO_SERVICE_ENABLED:-false}"; then
        profile_args+=(--profile indigo-service)
    fi

    local services ps_json
    services=$("${COMPOSE_CMD[@]}" -f "$COMPOSE_FILE" "${profile_args[@]}" config --services 2>/dev/null || echo "")
    ps_json=$(docker_compose ps --format json 2>/dev/null | jq -sr '.' 2>/dev/null || echo "[]")

    if [ -z "$services" ] && [ "$ps_json" = "[]" ]; then
        print_svc "Docker Compose" "stopped" "(no services)"
        return 1
    fi

    echo -e "\n  ${BOLD}Configured Services:${NC}"
    for svc in $services; do
        local cid state health ports started
        cid=$(docker_compose ps -q "$svc" 2>/dev/null || true)
        if [ -z "$cid" ]; then
            echo -e "    ${YELLOW}?${NC} ${BOLD}${svc}${NC} (not created)"
            continue
        fi
        state=$(docker inspect -f '{{.State.Status}}' "$cid" 2>/dev/null || echo "unknown")
        health=$(docker inspect -f '{{.State.Health.Status}}' "$cid" 2>/dev/null || true)
        ports=$(docker port "$cid" 2>/dev/null | head -1 || echo "")
        started=$(docker inspect -f '{{.State.StartedAt}}' "$cid" 2>/dev/null |
            xargs -I{} date -d {} "+%Y-%m-%d %H:%M" 2>/dev/null || echo "")

        if [ "$state" = "running" ]; then
            local health_str=""
            [ -n "$health" ] && health_str=" | health: ${health}"
            echo -e "    ${GREEN}●${NC} ${BOLD}${svc}${NC} (running${health_str})"
            [ -n "$ports" ] && echo -e "       port: ${ports}"
            [ -n "$started" ] && echo -e "       started: ${started}"
        else
            echo -e "    ${RED}○${NC} ${BOLD}${svc}${NC} (${state})"
        fi
    done

    local actual_services legacy_found=false
    actual_services=$(printf '%s' "$ps_json" | jq -r '.[].Service' 2>/dev/null | sort -u)
    for svc in $actual_services; do
        if [ -z "$svc" ]; then
            continue
        fi
        if printf '%s\n' "$services" | grep -qx "$svc"; then
            continue
        fi

        if [ "$legacy_found" = false ]; then
            echo -e "\n  ${BOLD}Legacy / Orphan Project Containers:${NC}"
            legacy_found=true
        fi

        local cid state health ports started
        cid=$(docker_compose ps -q "$svc" 2>/dev/null || true)
        if [ -z "$cid" ]; then
            cid=$(docker ps -aq -f "label=com.docker.compose.project=openwebui" -f "label=com.docker.compose.service=${svc}" | head -1)
        fi

        state=$(docker inspect -f '{{.State.Status}}' "$cid" 2>/dev/null || echo "unknown")
        health=$(docker inspect -f '{{.State.Health.Status}}' "$cid" 2>/dev/null || true)
        ports=$(docker port "$cid" 2>/dev/null | head -1 || echo "")
        started=$(docker inspect -f '{{.State.StartedAt}}' "$cid" 2>/dev/null |
            xargs -I{} date -d {} "+%Y-%m-%d %H:%M" 2>/dev/null || echo "")

        local health_str=""
        [ -n "$health" ] && health_str=" | health: ${health}"
        echo -e "    ${YELLOW}◐${NC} ${BOLD}${svc}${NC} (${state}${health_str} | not in current compose config)"
        [ -n "$ports" ] && echo -e "       port: ${ports}"
        [ -n "$started" ] && echo -e "       started: ${started}"
    done
    echo
}

check_openwebui() {
    print_section "OpenWebUI"
    local url="http://localhost:${OPENWEBUI_PORT}"
    if curl -sf -m 5 "$url/health" &>/dev/null || curl -sf -m 5 "$url" &>/dev/null; then
        local image
        image=$(docker inspect -f '{{.Config.Image}}' openwebui 2>/dev/null || echo "N/A")
        print_svc "OpenWebUI" "running" "(${url})"
        print_info "Image: ${image}"
        return 0
    fi
    print_svc "OpenWebUI" "stopped" "(not accessible at ${url})"
    return 1
}

http_code_is_reachable() {
    local code="$1"
    [ "$code" = "200" ] || [ "$code" = "401" ]
}

probe_litellm_from_host() {
    local models_url="${LITELLM_URL%/}/v1/models"
    local code
    code=$(curl -so /dev/null -w "%{http_code}" -m 5 "$models_url" 2>/dev/null || echo "000")
    if http_code_is_reachable "$code"; then
        LITELLM_PROBE_SOURCE="host"
        LITELLM_PROBE_URL="$models_url"
        return 0
    fi
    return 1
}

probe_litellm_from_openwebui() {
    local is_running
    is_running=$(docker inspect -f '{{.State.Running}}' openwebui 2>/dev/null || echo "false")
    [ "$is_running" = "true" ] || return 1

    local upstream_url
    upstream_url=$(docker inspect -f '{{range .Config.Env}}{{println .}}{{end}}' openwebui 2>/dev/null |
        awk -F= '/^OPENAI_API_BASE_URL=/{print substr($0, index($0, "=") + 1); exit}')
    [ -n "$upstream_url" ] || upstream_url="${OPENAI_API_BASE_URL:-}"
    [ -n "$upstream_url" ] || return 1

    local models_url="${upstream_url%/}/models"
    local code
    code=$(docker exec openwebui sh -lc "curl -so /dev/null -w '%{http_code}' -m 5 \"$models_url\" 2>/dev/null || true" 2>/dev/null | tail -n1)
    if http_code_is_reachable "$code"; then
        LITELLM_PROBE_SOURCE="openwebui"
        LITELLM_PROBE_URL="$models_url"
        return 0
    fi
    return 1
}

resolve_litellm_status() {
    LITELLM_PROBE_SOURCE=""
    LITELLM_PROBE_URL=""

    local openwebui_running
    openwebui_running=$(docker inspect -f '{{.State.Running}}' openwebui 2>/dev/null || echo "false")

    # When OpenWebUI is running, validate the effective live route from the container first.
    if [ "$openwebui_running" = "true" ]; then
        probe_litellm_from_openwebui && return 0

        # Diagnostic fallback: host may still see LiteLLM even if OpenWebUI route is broken.
        if probe_litellm_from_host; then
            LITELLM_PROBE_SOURCE="host_only"
        fi
        return 1
    fi

    # If OpenWebUI is not running, host reachability is the best available signal.
    probe_litellm_from_host && return 0
    return 1
}

check_litellm() {
    print_section "LiteLLM Proxy (upstream)"
    if resolve_litellm_status; then
        if [ "$LITELLM_PROBE_SOURCE" = "openwebui" ]; then
            print_svc "LiteLLM" "running" "(reachable from OpenWebUI via ${LITELLM_PROBE_URL})"
            if probe_litellm_from_host; then
                print_info "Host probe (${LITELLM_URL}) is also reachable"
            else
                print_info "Host probe (${LITELLM_URL}) unavailable, but live OpenWebUI upstream route is healthy"
            fi
        else
            print_svc "LiteLLM" "running" "(${LITELLM_URL})"
        fi
        return 0
    fi
    if [ "$LITELLM_PROBE_SOURCE" = "host_only" ]; then
        print_svc "LiteLLM" "warning" "(host reachable, but OpenWebUI cannot reach its configured upstream)"
        print_info "Check OPENAI_API_BASE_URL from inside container and host-gateway routing"
        return 1
    fi
    print_svc "LiteLLM" "stopped" "(not reachable from host or OpenWebUI container)"
    return 1
}

check_mcpo() {
    print_section "MCPO (MCP Proxy)"
    local url="${MCPO_BASE_URL%/}"
    if curl -sf -m 5 "${url}/docs" &>/dev/null 2>&1; then
        print_svc "MCPO" "running" "(${url})"
        return 0
    fi
    print_svc "MCPO" "stopped" "(not reachable at ${url})"
    return 1
}

check_multimodal_retrieval_api() {
    print_section "Canonical Retrieval Service"

    local base_url="${MULTIMODAL_RETRIEVAL_API_URL%/}"
    local health_response ready_response ready_http

    health_response="$(curl -sS -m 5 "${base_url}/health" 2>/dev/null || true)"
    if ! echo "$health_response" | jq -e '.status == "ok"' >/dev/null 2>&1; then
        print_svc "multimodal_retrieval_api" "stopped" "(health probe failed at ${base_url}/health)"
        return 1
    fi

    ready_response="$(mktemp)"
    ready_http="$(curl -sS -m 8 -o "$ready_response" -w "%{http_code}" "${base_url}/ready" 2>/dev/null || true)"
    if [ "$ready_http" = "200" ] && jq -e '.status == "ready"' "$ready_response" >/dev/null 2>&1; then
        print_svc "multimodal_retrieval_api" "running" "(${base_url})"
        rm -f "$ready_response"
        return 0
    fi

    print_svc "multimodal_retrieval_api" "warning" "(${base_url} | /ready HTTP ${ready_http:-000})"
    if [ -s "$ready_response" ]; then
        print_info "Ready payload: $(tr '\n' ' ' <"$ready_response" | head -c 220)"
    fi
    rm -f "$ready_response"
    return 1
}

check_open_terminal() {
    print_section "Open Terminal (optional sidecar)"

    if ! is_true "$OPEN_TERMINAL_ENABLED"; then
        print_info "Open Terminal disabled (OPEN_TERMINAL_ENABLED=false)"
        return 0
    fi

    local url="${OPEN_TERMINAL_BASE_URL%/}"
    if curl -sf -m 5 "${url}/health" &>/dev/null 2>&1; then
        print_svc "Open Terminal" "running" "(${url})"
        return 0
    fi

    print_svc "Open Terminal" "warning" "(not reachable at ${url})"
    return 1
}

check_indigo_service() {
    print_section "Indigo Service (optional sidecar)"

    if ! is_true "$INDIGO_SERVICE_ENABLED"; then
        print_info "Indigo Service disabled (INDIGO_SERVICE_ENABLED=false)"
        return 0
    fi

    local url="${INDIGO_SERVICE_BASE_URL%/}"
    local check_script="${SCRIPT_DIR}/scripts/indigo/check-indigo-service.sh"
    if [ -x "$check_script" ] && "$check_script" --quiet &>/dev/null 2>&1; then
        print_svc "Indigo Service" "running" "(${url})"
        return 0
    fi

    print_svc "Indigo Service" "warning" "(not reachable at ${url})"
    return 1
}

check_indigo_tool_registration() {
    print_section "Indigo Tool (OpenWebUI)"

    if ! is_true "$INDIGO_SERVICE_ENABLED"; then
        print_info "Indigo tool check skipped (INDIGO_SERVICE_ENABLED=false)"
        return 0
    fi

    if ! is_true "$INDIGO_TOOL_ENABLED"; then
        print_info "Indigo tool check skipped (INDIGO_TOOL_ENABLED=false)"
        return 0
    fi

    local running
    running="$(docker inspect -f '{{.State.Running}}' openwebui 2>/dev/null || echo "false")"
    if [ "$running" != "true" ]; then
        print_svc "Indigo Tool" "warning" "(openwebui container is not running)"
        return 1
    fi

    local output
    if ! output="$(
        docker exec -i openwebui python - "$INDIGO_TOOL_ID" <<'PY'
import sqlite3
import sys

tool_id = sys.argv[1]
required = {
    "info",
    "check",
    "convert",
    "calculate",
    "render",
    "aromatize",
    "dearomatize",
    "layout",
    "clean",
}

conn = sqlite3.connect("/app/backend/data/webui.db")
cur = conn.cursor()
cur.execute("SELECT id, name, is_active, specs FROM tool WHERE id = ?", (tool_id,))
row = cur.fetchone()
conn.close()
if row is None:
    raise SystemExit(f"missing:{tool_id}")

specs_raw = row[3] or "[]"
try:
    import json
    specs = json.loads(specs_raw)
except Exception:
    specs = []

actual = {item.get("name") for item in specs if isinstance(item, dict) and isinstance(item.get("name"), str)}
missing = sorted(required - actual)
if missing:
    raise SystemExit(f"missing_specs:{','.join(missing)}")

print(f"id:{row[0]} name:{row[1]} active:{row[2]} specs:{','.join(sorted(actual))}")
PY
    )"; then
        print_svc "Indigo Tool" "warning" "(not registered or incomplete: ${INDIGO_TOOL_ID})"
        return 1
    fi

    print_svc "Indigo Tool" "running" "(${INDIGO_TOOL_ID})"
    print_info "$output"
    return 0
}

check_rag_profile() {
    print_section "RAG Profile"
    print_info "MULTIMODAL_RETRIEVAL_API_URL=${MULTIMODAL_RETRIEVAL_API_URL}"
    print_info "RAG_EMBEDDING_MODEL=${RAG_EMBEDDING_MODEL:-<unset>}"
    print_info "RAG_TOP_K=${RAG_TOP_K:-<unset>} | CHUNK_SIZE=${CHUNK_SIZE:-<unset>}"
    print_info "VECTOR_DB=${VECTOR_DB:-<unset>}"

    local web_search=false
    { is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; } && web_search=true

    if [ "$web_search" = true ]; then
        print_info "Web search: enabled (${SEARXNG_QUERY_URL:-<unset>})"
        local probe_url="${SEARXNG_QUERY_URL//\{query\}/status+probe}"
        local host_probe="$probe_url"
        if echo "$host_probe" | grep -q "host.docker.internal"; then
            host_probe="${host_probe//host.docker.internal/127.0.0.1}"
        fi

        if curl -sf -m 8 "$host_probe" 2>/dev/null | grep -q "results"; then
            print_svc "SearXNG (host)" "running"
        elif docker ps -q -f "name=^openwebui$" | grep -q . &&
            docker exec openwebui sh -lc "curl -sf -m 8 \"$probe_url\" 2>/dev/null | grep -q results"; then
            print_svc "SearXNG (container route)" "running"
        else
            print_svc "SearXNG" "warning" "(probe failed from host and OpenWebUI container route)"
        fi
    else
        print_info "Web search: disabled"
    fi
}

check_resources() {
    print_section "System Resources"
    print_info "Disk: $(df -h . | tail -1 | awk '{print $5 " used, " $4 " free"}')"
    command -v free &>/dev/null && print_info "Memory: $(free -h | awk '/^Mem:/{print $3"/"$2" used"}')"
}

print_summary() {
    print_section "Summary"

    local all_good=true

    # OpenWebUI
    if curl -sf -m 5 "http://localhost:${OPENWEBUI_PORT}/health" &>/dev/null ||
        curl -sf -m 5 "http://localhost:${OPENWEBUI_PORT}" &>/dev/null; then
        echo -e "  ${GREEN}●${NC} OpenWebUI:  ${GREEN}OK${NC}"
    else
        echo -e "  ${RED}○${NC} OpenWebUI:  ${RED}Not accessible${NC}"
        all_good=false
    fi

    # LiteLLM (host or effective OpenWebUI upstream route)
    if resolve_litellm_status; then
        if [ "$LITELLM_PROBE_SOURCE" = "openwebui" ]; then
            echo -e "  ${GREEN}●${NC} LiteLLM:    ${GREEN}OK${NC} (via OpenWebUI upstream route)"
        else
            echo -e "  ${GREEN}●${NC} LiteLLM:    ${GREEN}OK${NC}"
        fi
    elif [ "$LITELLM_PROBE_SOURCE" = "host_only" ]; then
        echo -e "  ${YELLOW}◐${NC} LiteLLM:    ${YELLOW}Host OK / OpenWebUI route broken${NC}"
        all_good=false
    else
        echo -e "  ${RED}○${NC} LiteLLM:    ${RED}Not reachable${NC}"
        all_good=false
    fi

    # Docker
    if docker info &>/dev/null 2>&1; then
        echo -e "  ${GREEN}●${NC} Docker:     ${GREEN}OK${NC}"
    else
        echo -e "  ${RED}○${NC} Docker:     ${RED}Not running${NC}"
        all_good=false
    fi

    # MCPO
    if curl -sf -m 5 "${MCPO_BASE_URL%/}/docs" &>/dev/null 2>&1; then
        echo -e "  ${GREEN}●${NC} MCPO:       ${GREEN}OK${NC}"
    else
        echo -e "  ${YELLOW}◐${NC} MCPO:       ${YELLOW}Not reachable${NC}"
    fi

    # Canonical retrieval service
    local canonical_url="${MULTIMODAL_RETRIEVAL_API_URL%/}"
    local canonical_health canonical_ready_http canonical_ready_file
    canonical_health="$(curl -sS -m 5 "${canonical_url}/health" 2>/dev/null || true)"
    canonical_ready_file="$(mktemp)"
    canonical_ready_http="$(curl -sS -m 8 -o "$canonical_ready_file" -w "%{http_code}" "${canonical_url}/ready" 2>/dev/null || true)"
    if echo "$canonical_health" | jq -e '.status == "ok"' >/dev/null 2>&1 &&
        [ "$canonical_ready_http" = "200" ] &&
        jq -e '.status == "ready"' "$canonical_ready_file" >/dev/null 2>&1; then
        echo -e "  ${GREEN}●${NC} CanonicalRAG:${GREEN} OK${NC}"
    elif echo "$canonical_health" | jq -e '.status == "ok"' >/dev/null 2>&1; then
        echo -e "  ${YELLOW}◐${NC} CanonicalRAG:${YELLOW} Health OK / Ready degraded${NC}"
        all_good=false
    else
        echo -e "  ${RED}○${NC} CanonicalRAG:${RED} Not reachable${NC}"
        all_good=false
    fi
    rm -f "$canonical_ready_file"

    # Open Terminal (optional)
    if is_true "$OPEN_TERMINAL_ENABLED"; then
        if curl -sf -m 5 "${OPEN_TERMINAL_BASE_URL%/}/health" &>/dev/null 2>&1; then
            echo -e "  ${GREEN}●${NC} OpenTerm:   ${GREEN}OK${NC}"
        else
            echo -e "  ${YELLOW}◐${NC} OpenTerm:   ${YELLOW}Not reachable${NC}"
        fi
    fi

    # Indigo Service (optional)
    if is_true "$INDIGO_SERVICE_ENABLED"; then
        if [ -x "${SCRIPT_DIR}/scripts/indigo/check-indigo-service.sh" ] && "${SCRIPT_DIR}/scripts/indigo/check-indigo-service.sh" --quiet &>/dev/null 2>&1; then
            echo -e "  ${GREEN}●${NC} Indigo:     ${GREEN}OK${NC}"
        else
            echo -e "  ${YELLOW}◐${NC} Indigo:     ${YELLOW}Not reachable${NC}"
        fi

        if is_true "$INDIGO_TOOL_ENABLED"; then
            if docker exec -i openwebui python - "$INDIGO_TOOL_ID" <<'PY' >/dev/null 2>&1; then
import sqlite3
import sys
tool_id = sys.argv[1]
conn = sqlite3.connect("/app/backend/data/webui.db")
cur = conn.cursor()
cur.execute("SELECT 1 FROM tool WHERE id = ? AND is_active = 1", (tool_id,))
row = cur.fetchone()
conn.close()
raise SystemExit(0 if row else 1)
PY
                echo -e "  ${GREEN}●${NC} IndigoTool: ${GREEN}OK${NC}"
            else
                echo -e "  ${YELLOW}◐${NC} IndigoTool: ${YELLOW}Missing in OpenWebUI${NC}"
            fi
        fi
    fi

    echo
    if [ "$all_good" = true ]; then
        print_success "All core services operational — http://localhost:${OPENWEBUI_PORT}"
    else
        print_warning "Some services are down — run ./deploy.sh to start"
    fi
    echo
}

###############################################################################
# Main
###############################################################################

main() {
    print_header "OpenWebUI RAG System Status"

    WATCH_MODE=false
    QUIET_MODE=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            -w | --watch)
                WATCH_MODE=true
                shift
                ;;
            -q | --quiet)
                QUIET_MODE=true
                shift
                ;;
            -h | --help)
                echo "Usage: $0 [-w|--watch] [-q|--quiet]"
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                exit 1
                ;;
        esac
    done

    run_checks() {
        if [ "$QUIET_MODE" = false ]; then
            check_docker || true
            check_docker_compose || true
            check_openwebui || true
            check_litellm || true
            check_mcpo || true
            check_multimodal_retrieval_api || true
            check_open_terminal || true
            check_indigo_service || true
            check_indigo_tool_registration || true
            check_rag_profile || true
            check_resources || true
        fi
        print_summary
    }

    if [ "$WATCH_MODE" = true ]; then
        while true; do
            clear
            print_header "OpenWebUI RAG System Status"
            run_checks
            echo -e "${CYAN}Refreshing every 5s — Ctrl+C to exit${NC}"
            sleep 5
        done
    else
        run_checks
    fi
}

main "$@"
