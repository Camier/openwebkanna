#!/bin/bash

###############################################################################
# OpenWebUI + vLLM RAG Logs Viewer
# This script displays logs from all services
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
cd "$SCRIPT_DIR"

# Configuration
# shellcheck disable=SC2034
VLLM_PORT=8000
# shellcheck disable=SC2034
OPENWEBUI_PORT=3000
COMPOSE_FILE="docker-compose.yml"
VLLM_LOG_FILE="logs/vllm.log"
CLIPROXYAPI_LOG_FILE="${CLIPROXYAPI_LOG_FILE:-logs/cliproxyapi.log}"

###############################################################################
# Script-specific Helper Functions
###############################################################################

get_file_modified() {
    local file=$1
    if date -r "$file" "+%Y-%m-%d %H:%M:%S" >/dev/null 2>&1; then
        date -r "$file" "+%Y-%m-%d %H:%M:%S"
    else
        echo "unknown"
    fi
}

show_vllm_logs() {
    local lines=${1:-50}

    print_section "vLLM Logs (last $lines lines)"

    if [ ! -f "$VLLM_LOG_FILE" ]; then
        echo -e "${YELLOW}vLLM log file not found: $VLLM_LOG_FILE${NC}"
        echo -e "${YELLOW}vLLM may not be running or logging to a different location${NC}"
        return 1
    fi

    local file_size
    file_size=$(du -h "$VLLM_LOG_FILE" | cut -f1)
    print_info "Log file: $VLLM_LOG_FILE ($file_size)"

    # Show logs with tail
    tail -n "$lines" "$VLLM_LOG_FILE"
}

follow_vllm_logs() {
    print_section "vLLM Logs (live tail)"
    print_info "Press Ctrl+C to stop following"
    echo

    if [ ! -f "$VLLM_LOG_FILE" ]; then
        echo -e "${YELLOW}vLLM log file not found: $VLLM_LOG_FILE${NC}"
        echo -e "${YELLOW}Waiting for log file to be created...${NC}"
        mkdir -p logs
    fi

    tail -f "$VLLM_LOG_FILE" 2>/dev/null
}

show_cliproxyapi_logs() {
    local lines=${1:-50}

    print_section "CLIProxyAPI Logs (last $lines lines)"

    if [ ! -f "$CLIPROXYAPI_LOG_FILE" ]; then
        echo -e "${YELLOW}CLIProxyAPI log file not found: $CLIPROXYAPI_LOG_FILE${NC}"
        echo -e "${YELLOW}CLIProxyAPI may not be running${NC}"
        return 1
    fi

    local total_lines
    total_lines=$(wc -l <"$CLIPROXYAPI_LOG_FILE")
    print_info "Log file: $CLIPROXYAPI_LOG_FILE ($total_lines total lines)"
    tail -n "$lines" "$CLIPROXYAPI_LOG_FILE"
}

follow_cliproxyapi_logs() {
    print_section "CLIProxyAPI Logs (live tail)"
    print_info "Press Ctrl+C to stop following"
    echo

    if [ ! -f "$CLIPROXYAPI_LOG_FILE" ]; then
        echo -e "${YELLOW}CLIProxyAPI log file not found: $CLIPROXYAPI_LOG_FILE${NC}"
        echo -e "${YELLOW}Waiting for log file to be created...${NC}"
        mkdir -p "$(dirname "$CLIPROXYAPI_LOG_FILE")"
    fi

    tail -f "$CLIPROXYAPI_LOG_FILE" 2>/dev/null
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

# shellcheck disable=SC2120
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

    # Show vLLM logs
    show_vllm_logs "$lines"

    # Show CLIProxyAPI logs
    show_cliproxyapi_logs "$lines"

    # Show Docker logs
    show_docker_logs "" "$lines"
}

show_log_summary() {
    print_section "Log Summary"

    # vLLM log info
    if [ -f "$VLLM_LOG_FILE" ]; then
        local vllm_size
        local vllm_modified
        vllm_size=$(du -h "$VLLM_LOG_FILE" | cut -f1)
        vllm_modified=$(get_file_modified "$VLLM_LOG_FILE")
        echo -e "${GREEN}vLLM:${NC} size=$vllm_size, modified=$vllm_modified"
    else
        echo -e "${YELLOW}vLLM:${NC} No log file found"
    fi

    # CLIProxyAPI log info
    if [ -f "$CLIPROXYAPI_LOG_FILE" ]; then
        local cliproxyapi_size cliproxyapi_lines
        cliproxyapi_size=$(du -h "$CLIPROXYAPI_LOG_FILE" | cut -f1)
        cliproxyapi_lines=$(wc -l <"$CLIPROXYAPI_LOG_FILE")
        echo -e "${GREEN}CLIProxyAPI:${NC} $cliproxyapi_lines lines, $cliproxyapi_size"
    else
        echo -e "${YELLOW}CLIProxyAPI:${NC} No log file found"
    fi

    # Docker log info
    if docker info &>/dev/null 2>&1 && [ -f "$COMPOSE_FILE" ]; then
        echo -e "${GREEN}Docker Compose:${NC}"
        docker_compose ps 2>/dev/null | while read -r line; do
            echo "  $line"
        done
    else
        echo -e "${YELLOW}Docker Compose:${NC} Not running or compose file not found"
    fi
}

search_logs() {
    local pattern=$1
    local service=${2:-all}

    print_section "Searching logs for: $pattern"

    if [ "$service" = "vllm" ] || [ "$service" = "all" ]; then
        echo -e "${MAGENTA}--- vLLM Logs ---${NC}"
        if [ -f "$VLLM_LOG_FILE" ]; then
            grep -i "$pattern" "$VLLM_LOG_FILE" || echo "No matches found in vLLM logs"
        else
            echo "vLLM log file not found"
        fi
    fi

    if [ "$service" = "docker" ] || [ "$service" = "all" ]; then
        echo -e "\n${MAGENTA}--- Docker Compose Logs ---${NC}"
        if docker info &>/dev/null 2>&1; then
            docker_compose logs 2>/dev/null | grep -i "$pattern" || echo "No matches found in Docker logs"
        else
            echo "Docker is not running"
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

    if [ "$service" = "cliproxyapi" ] || [ "$service" = "all" ]; then
        echo -e "\n${MAGENTA}--- CLIProxyAPI Logs ---${NC}"
        if [ -f "$CLIPROXYAPI_LOG_FILE" ]; then
            grep -i "$pattern" "$CLIPROXYAPI_LOG_FILE" || echo "No matches found in CLIProxyAPI logs"
        else
            echo "CLIProxyAPI log file not found"
        fi
    fi
}

show_help() {
    cat <<EOF
Usage: $0 [OPTIONS] [SERVICE]

View logs from OpenWebUI + vLLM RAG system services.

Options:
  -f, --follow        Follow logs (live tail)
  -n, --lines NUM     Show last NUM lines (default: 50)
  -s, --search PATTERN Search logs for pattern (case-insensitive)
  --summary           Show log file summary
  -h, --help          Show this help message

Services:
  all                 Show logs from all services (default)
  vllm                Show only vLLM logs
  cliproxyapi         Show only CLIProxyAPI logs
  docker              Show only Docker Compose logs
  openwebui           Show only OpenWebUI logs

Examples:
  $0                      # Show last 50 lines from all logs
  $0 -f                   # Follow all logs
  $0 -n 100 vllm          # Show last 100 lines of vLLM logs
  $0 -n 100 cliproxyapi   # Show last 100 lines of CLIProxyAPI logs
  $0 -f openwebui         # Follow OpenWebUI logs
  $0 -s "error"           # Search all logs for "error"
  $0 -s "error" vllm      # Search vLLM logs for "error"

EOF
}

###############################################################################
# Main Logs Viewer
###############################################################################

main() {
    local follow_mode=false
    local lines=50
    local search_pattern=""
    local service="all"
    local show_summary_only=false

    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -f | --follow)
                follow_mode=true
                shift
                ;;
            -n | --lines)
                lines="$2"
                shift 2
                ;;
            -s | --search)
                search_pattern="$2"
                shift 2
                ;;
            --summary)
                show_summary_only=true
                shift
                ;;
            -h | --help)
                show_help
                exit 0
                ;;
            vllm | cliproxyapi | docker | openwebui | all)
                service="$1"
                shift
                ;;
            *)
                if [ -z "$search_pattern" ]; then
                    search_pattern="$1"
                else
                    echo -e "${RED}Unknown option or too many arguments: $1${NC}"
                    show_help
                    exit 1
                fi
                shift
                ;;
        esac
    done

    # Show summary only if requested
    if [ "$show_summary_only" = true ]; then
        print_header
        show_log_summary
        exit 0
    fi

    # Search mode
    if [ -n "$search_pattern" ]; then
        print_header
        search_logs "$search_pattern" "$service"
        exit 0
    fi

    # Follow mode
    if [ "$follow_mode" = true ]; then
        case $service in
            vllm)
                print_header
                follow_vllm_logs
                ;;
            cliproxyapi)
                print_header
                follow_cliproxyapi_logs
                ;;
            docker)
                print_header
                follow_docker_logs
                ;;
            openwebui)
                print_header
                follow_openwebui_logs
                ;;
            all)
                print_header
                print_info "Following all logs (vLLM and Docker Compose)"
                print_info "Press Ctrl+C to stop"
                echo

                # Follow both in parallel using background processes
                trap "kill 0" EXIT INT TERM

                (
                    while true; do
                        if [ -f "$VLLM_LOG_FILE" ]; then
                            tail -f "$VLLM_LOG_FILE" &
                            wait $!
                        fi
                        sleep 1
                    done
                ) &

                (
                    while true; do
                        if [ -f "$CLIPROXYAPI_LOG_FILE" ]; then
                            tail -f "$CLIPROXYAPI_LOG_FILE" &
                            wait $!
                        fi
                        sleep 1
                    done
                ) &

                (
                    while true; do
                        if docker info &>/dev/null 2>&1; then
                            docker_compose logs -f &
                            wait $!
                        fi
                        sleep 1
                    done
                ) &

                wait
                ;;
        esac
        exit 0
    fi

    # Normal display mode
    print_header

    case $service in
        vllm)
            show_vllm_logs "$lines"
            ;;
        cliproxyapi)
            show_cliproxyapi_logs "$lines"
            ;;
        docker)
            show_docker_logs "" "$lines"
            ;;
        openwebui)
            show_openwebui_logs "$lines"
            ;;
        all)
            show_all_logs "$lines"
            ;;
    esac
}

main "$@"
