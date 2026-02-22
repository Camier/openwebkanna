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
# shellcheck disable=SC2034
COMPOSE_CMD=()

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

check_litellm() {
    print_section "LiteLLM Proxy (upstream)"
    local health_url="${LITELLM_URL}/health"
    if curl -sf -m 5 "$health_url" &>/dev/null 2>&1 ||
        curl -sf -m 5 "${LITELLM_URL}/v1/models" -H "Authorization: Bearer ${OPENAI_API_KEY:-}" &>/dev/null 2>&1; then
        print_svc "LiteLLM" "running" "(${LITELLM_URL})"
        return 0
    fi
    # 401 means reachable but unauthenticated — still up
    local code
    code=$(curl -so /dev/null -w "%{http_code}" -m 5 "${LITELLM_URL}/v1/models" 2>/dev/null || echo "000")
    if [ "$code" = "401" ] || [ "$code" = "200" ]; then
        print_svc "LiteLLM" "running" "(${LITELLM_URL})"
        return 0
    fi
    print_svc "LiteLLM" "stopped" "(not reachable at ${LITELLM_URL})"
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

    # LiteLLM
    local litellm_code
    litellm_code=$(curl -so /dev/null -w "%{http_code}" -m 5 "${LITELLM_URL}/v1/models" 2>/dev/null || echo "000")
    if [ "$litellm_code" = "200" ] || [ "$litellm_code" = "401" ]; then
        echo -e "  ${GREEN}●${NC} LiteLLM:    ${GREEN}OK${NC}"
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
