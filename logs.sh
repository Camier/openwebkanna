#!/bin/bash

###############################################################################
# OpenWebUI + LiteLLM RAG Logs Viewer
# Shows logs for the surviving compose baseline and selected services.
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

COMPOSE_FILE="config/compose/docker-compose.yml"
MULTIMODAL_RETRIEVAL_API_STATE_DIR="${MULTIMODAL_RETRIEVAL_API_STATE_DIR:-${SCRIPT_DIR}/.codex-state/multimodal_retrieval_api}"
MULTIMODAL_RETRIEVAL_API_PID_FILE="${MULTIMODAL_RETRIEVAL_API_PID_FILE:-${MULTIMODAL_RETRIEVAL_API_STATE_DIR}/multimodal_retrieval_api.pid}"
MULTIMODAL_RETRIEVAL_API_LOG_FILE="${MULTIMODAL_RETRIEVAL_API_LOG_FILE:-${SCRIPT_DIR}/logs/multimodal_retrieval_api.log}"

show_multimodal_retrieval_logs() {
    local lines=${1:-50}
    print_section "Canonical Retrieval Logs (last $lines lines)"
    if [ -f "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" ]; then
        tail -n "$lines" "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}"
        return 0
    fi
    print_info "No canonical retrieval log file found at ${MULTIMODAL_RETRIEVAL_API_LOG_FILE}"
}

follow_multimodal_retrieval_logs() {
    print_section "Canonical Retrieval Logs (live tail)"
    print_info "Press Ctrl+C to stop following"
    echo
    if [ -f "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" ]; then
        tail -n 50 -f "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}"
        return 0
    fi
    print_info "No canonical retrieval log file found at ${MULTIMODAL_RETRIEVAL_API_LOG_FILE}"
}

show_docker_logs() {
    local service=${1:-}
    local lines=${2:-50}

    if [ -n "$service" ]; then
        print_section "Docker Compose Logs: $service (last $lines lines)"
        docker_compose logs --tail="$lines" "$service"
    else
        print_section "Docker Compose Logs (all services, last $lines lines)"
        docker_compose logs --tail="$lines"
    fi
}

follow_docker_logs() {
    local service=${1:-}

    if [ -n "$service" ]; then
        print_section "Docker Compose Logs: $service (live tail)"
        print_info "Press Ctrl+C to stop following"
        echo
        docker_compose logs -f "$service"
    else
        print_section "Docker Compose Logs (all services, live tail)"
        print_info "Press Ctrl+C to stop following"
        echo
        docker_compose logs -f
    fi
}

show_openwebui_logs() {
    local lines=${1:-50}
    print_section "OpenWebUI Logs (last $lines lines)"
    docker_compose logs --tail="$lines" openwebui
}

follow_openwebui_logs() {
    print_section "OpenWebUI Logs (live tail)"
    print_info "Press Ctrl+C to stop following"
    echo
    docker_compose logs -f openwebui
}

show_all_logs() {
    local lines=${1:-50}
    print_header
    show_docker_logs "" "$lines"
}

show_log_summary() {
    print_section "Log Summary"
    if docker info &>/dev/null 2>&1 && [ -f "$COMPOSE_FILE" ]; then
        echo -e "${GREEN}Docker Compose:${NC}"
        docker_compose ps 2>/dev/null | while read -r line; do
            echo "  $line"
        done
    else
        echo -e "${YELLOW}Docker Compose:${NC} Not running or compose file not found"
    fi

    if [ -f "${MULTIMODAL_RETRIEVAL_API_PID_FILE}" ]; then
        local pid
        pid="$(cat "${MULTIMODAL_RETRIEVAL_API_PID_FILE}" 2>/dev/null || true)"
        if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
            echo -e "${GREEN}Canonical Retrieval:${NC} Running (pid=${pid})"
            echo "  log: ${MULTIMODAL_RETRIEVAL_API_LOG_FILE}"
        else
            echo -e "${YELLOW}Canonical Retrieval:${NC} PID file present but process not running"
        fi
    else
        echo -e "${YELLOW}Canonical Retrieval:${NC} No managed PID file"
    fi
}

search_logs() {
    local pattern=$1
    local service=${2:-all}

    print_section "Searching logs for: $pattern"

    if [ "$service" = "docker" ] || [ "$service" = "all" ]; then
        echo -e "\n${MAGENTA}--- Docker Compose Logs ---${NC}"
        if docker info &>/dev/null 2>&1; then
            docker_compose logs 2>/dev/null | grep -i "$pattern" || echo "No matches found in Docker logs"
        else
            echo "Docker is not running"
        fi
    fi

    if [ "$service" = "multimodal_retrieval_api" ] || [ "$service" = "canonical" ] || [ "$service" = "all" ]; then
        echo -e "\n${MAGENTA}--- Canonical Retrieval Logs ---${NC}"
        if [ -f "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" ]; then
            grep -i "$pattern" "${MULTIMODAL_RETRIEVAL_API_LOG_FILE}" || echo "No matches found in canonical retrieval logs"
        else
            echo "Canonical retrieval log file not found"
        fi
    fi

    if [ "$service" = "openwebui" ]; then
        echo -e "\n${MAGENTA}--- OpenWebUI Logs ---${NC}"
        if docker info &>/dev/null 2>&1; then
            docker_compose logs openwebui 2>/dev/null | grep -i "$pattern" || echo "No matches found in OpenWebUI logs"
        else
            echo "Docker is not running"
        fi
    fi
}

show_help() {
    cat <<'EOF'
Usage: ./logs.sh [OPTIONS] [SERVICE]

Options:
  -f, --follow         Follow logs live
  -n, --lines N        Number of lines to show (default: 50)
  -s, --search TERM    Search logs for TERM
  --summary            Show compose log summary
  -h, --help           Show this help

Examples:
  ./logs.sh
  ./logs.sh openwebui
  ./logs.sh multimodal_retrieval_api
  ./logs.sh --follow
  ./logs.sh --follow openwebui
  ./logs.sh --follow multimodal_retrieval_api
  ./logs.sh --search error
EOF
}

main() {
    local follow=false
    local summary=false
    local lines=50
    local search_term=""
    local service=""

    while [[ $# -gt 0 ]]; do
        case $1 in
            -f | --follow)
                follow=true
                shift
                ;;
            -n | --lines)
                lines="$2"
                shift 2
                ;;
            -s | --search)
                search_term="$2"
                shift 2
                ;;
            --summary)
                summary=true
                shift
                ;;
            -h | --help)
                show_help
                exit 0
                ;;
            *)
                service="$1"
                shift
                ;;
        esac
    done

    if ! init_compose_cmd; then
        print_error "Docker Compose is not available"
        exit 1
    fi

    if [ "$summary" = true ]; then
        show_log_summary
        exit 0
    fi

    if [ -n "$search_term" ]; then
        search_logs "$search_term" "${service:-all}"
        exit 0
    fi

    if [ "$follow" = true ]; then
        if [ "$service" = "multimodal_retrieval_api" ] || [ "$service" = "canonical" ]; then
            follow_multimodal_retrieval_logs
            exit 0
        fi
        if [ -n "$service" ]; then
            follow_docker_logs "$service"
        else
            follow_docker_logs
        fi
        exit 0
    fi

    if [ -n "$service" ]; then
        if [ "$service" = "multimodal_retrieval_api" ] || [ "$service" = "canonical" ]; then
            show_multimodal_retrieval_logs "$lines"
            exit 0
        fi
        show_docker_logs "$service" "$lines"
    else
        show_all_logs "$lines"
    fi
}

main "$@"
