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
COMPOSE_FILE="docker-compose.yml"
MCPO_BASE_URL="${MCPO_BASE_URL:-http://127.0.0.1:${MCPO_PORT:-8000}}"
LITELLM_URL="${LITELLM_URL:-http://localhost:4000}"
SMILES_STRUCTURE_API_HOST="${SMILES_STRUCTURE_API_HOST:-127.0.0.1}"
SMILES_STRUCTURE_API_PORT="${SMILES_STRUCTURE_API_PORT:-8011}"
SMILES_STRUCTURE_API_BASE_URL="${SMILES_STRUCTURE_API_BASE_URL:-http://${SMILES_STRUCTURE_API_HOST}:${SMILES_STRUCTURE_API_PORT}}"
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

    local services
    services=$(docker_compose ps --services 2>/dev/null || echo "")
    if [ -z "$services" ]; then
        print_svc "Docker Compose" "stopped" "(no services)"
        return 1
    fi

    echo -e "\n  ${BOLD}Services:${NC}"
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

check_smiles_structure_api() {
    print_section "SMILES Structure API"
    local url="${SMILES_STRUCTURE_API_BASE_URL%/}"
    local check_script="${SCRIPT_DIR}/check-smiles-structure-api.sh"
    if curl -sf -m 5 "${url}/health" &>/dev/null 2>&1; then
        print_svc "SMILES API" "running" "(${url})"
        return 0
    fi
    if [ -x "$check_script" ] && "$check_script" --quiet &>/dev/null 2>&1; then
        print_svc "SMILES API" "running" "(reachable via Docker network)"
        return 0
    fi
    print_svc "SMILES API" "warning" "(not reachable at ${url})"
    return 1
}

check_cliproxyapi() {
    print_section "CLIProxyAPI (optional sidecar)"
    # CLIProxyAPI is no longer the primary upstream; report informational status only.
    local check_script="${SCRIPT_DIR}/check-cliproxyapi.sh"
    if [ "${CLIPROXYAPI_ENABLED:-false}" = "false" ]; then
        print_info "CLIProxyAPI disabled (CLIPROXYAPI_ENABLED=false)"
        return 0
    fi
    if [ -x "$check_script" ] && CLIPROXYAPI_ENABLED=true "$check_script" --quiet &>/dev/null 2>&1; then
        print_svc "CLIProxyAPI" "running" "(optional sidecar — not the primary upstream)"
        return 0
    fi
    print_svc "CLIProxyAPI" "warning" "(not responding — this is OK; LiteLLM is the primary upstream)"
    return 0
}

check_rag_profile() {
    print_section "RAG Profile"
    print_info "RAG_EMBEDDING_MODEL=${RAG_EMBEDDING_MODEL:-<unset>}"
    print_info "RAG_TOP_K=${RAG_TOP_K:-<unset>} | CHUNK_SIZE=${CHUNK_SIZE:-<unset>}"
    print_info "VECTOR_DB=${VECTOR_DB:-<unset>}"

    local web_search=false
    { is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; } && web_search=true

    if [ "$web_search" = true ]; then
        print_info "Web search: enabled (${SEARXNG_QUERY_URL:-<unset>})"
        local probe_url="${SEARXNG_QUERY_URL//\{query\}/status+probe}"
        local host_probe="${probe_url//host.docker.internal/127.0.0.1}"
        if curl -sf -m 8 "$host_probe" 2>/dev/null | grep -q "results"; then
            print_svc "SearXNG (host)" "running"
        else
            print_svc "SearXNG (host)" "warning" "(probe failed — check port 8888)"
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

    # SMILES API
    if curl -sf -m 5 "${SMILES_STRUCTURE_API_BASE_URL%/}/health" &>/dev/null 2>&1 ||
        [ -x "${SCRIPT_DIR}/check-smiles-structure-api.sh" ] && "${SCRIPT_DIR}/check-smiles-structure-api.sh" --quiet &>/dev/null 2>&1; then
        echo -e "  ${GREEN}●${NC} SMILES API: ${GREEN}OK${NC}"
    else
        echo -e "  ${YELLOW}◐${NC} SMILES API: ${YELLOW}Not reachable${NC}"
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
            check_smiles_structure_api || true
            check_cliproxyapi || true
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
