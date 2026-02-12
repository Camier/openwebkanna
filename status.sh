#!/bin/bash

###############################################################################
# OpenWebUI + vLLM RAG Status Check Script
# This script checks the status of all services
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

load_env_defaults() {
    local env_file=""
    local line=""
    local key=""
    local value=""

    if [ -f "$SCRIPT_DIR/.env" ]; then
        env_file="$SCRIPT_DIR/.env"
    elif [ -f ".env" ]; then
        env_file=".env"
    fi

    if [ -z "$env_file" ]; then
        return 0
    fi

    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%$'\r'}"
        line="${line#"${line%%[![:space:]]*}"}"
        [ -z "$line" ] && continue
        [[ "$line" = \#* ]] && continue
        [[ "$line" != *=* ]] && continue

        key="${line%%=*}"
        value="${line#*=}"
        key="$(printf "%s" "$key" | tr -d '[:space:]')"

        [[ "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        if [ -n "${!key+x}" ]; then
            continue
        fi

        if [[ "$value" == \"*\" ]] && [[ "$value" == *\" ]]; then
            value="${value:1:${#value}-2}"
        elif [[ "$value" == \'*\' ]] && [[ "$value" == *\' ]]; then
            value="${value:1:${#value}-2}"
        fi

        printf -v "$key" "%s" "$value"
        export "$key"
    done < "$env_file"
}

load_env_defaults
cd "$SCRIPT_DIR"

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Configuration
VLLM_PORT=8000
VLLM_HEALTH_URL="http://localhost:${VLLM_PORT}/health"
OPENWEBUI_PORT="${WEBUI_PORT:-3000}"
COMPOSE_FILE="docker-compose.yml"
VLLM_LOG_FILE="logs/vllm.log"
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_DOCKER_MANAGED="${CLIPROXYAPI_DOCKER_MANAGED:-true}"
CLIPROXYAPI_PID_FILE="${CLIPROXYAPI_PID_FILE:-.cliproxyapi.pid}"
CLIPROXYAPI_LOG_FILE="${CLIPROXYAPI_LOG_FILE:-logs/cliproxyapi.log}"
CLIPROXYAPI_CHECK_SCRIPT="./check-cliproxyapi.sh"

###############################################################################
# Helper Functions
###############################################################################

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     OpenWebUI + vLLM RAG System Status                     ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_section() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────${NC}"
}

print_status() {
    local service=$1
    local status=$2
    local details=$3

    if [ "$status" = "running" ]; then
        echo -e "  ${GREEN}●${NC} ${BOLD}$service${NC}: ${GREEN}Running${NC} $details"
    elif [ "$status" = "stopped" ]; then
        echo -e "  ${RED}○${NC} ${BOLD}$service${NC}: ${RED}Stopped${NC} $details"
    elif [ "$status" = "warning" ]; then
        echo -e "  ${YELLOW}◐${NC} ${BOLD}$service${NC}: ${YELLOW}Warning${NC} $details"
    else
        echo -e "  ${YELLOW}?${NC} ${BOLD}$service${NC}: ${YELLOW}Unknown${NC} $details"
    fi
}

print_info() {
    echo -e "    ${CYAN}ℹ${NC} $1"
}

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

cliproxyapi_quiet_healthy() {
    [ "$CLIPROXYAPI_ENABLED" != "false" ] && [ -x "$CLIPROXYAPI_CHECK_SCRIPT" ] && CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1
}

check_vllm() {
    print_section "vLLM Server"

    # Check if PID file exists
    if [ -f .vllm_pid ]; then
        local pid=$(cat .vllm_pid)
        if kill -0 "$pid" 2>/dev/null; then
            # Check if health endpoint is responding
            if curl -sS -f -m 5 "$VLLM_HEALTH_URL" &> /dev/null; then
                local model_info=$(curl -sS -m 5 "$VLLM_HEALTH_URL" 2>/dev/null || echo "{}")
                print_status "vLLM" "running" "(PID: $pid)"
                print_info "Port: $VLLM_PORT"
                print_info "Health: $VLLM_HEALTH_URL"

                # Get model info if available
                local model=$(echo "$model_info" | grep -o '"model":"[^"]*"' | cut -d'"' -f4 2>/dev/null || echo "N/A")
                if [ "$model" != "N/A" ]; then
                    print_info "Model: $model"
                fi
                return 0
            else
                print_status "vLLM" "warning" "(PID: $pid exists but not responding)"
                print_info "Process exists but health check failing"
                return 1
            fi
        else
            print_status "vLLM" "stopped" "(stale PID file)"
            print_info "Remove .vllm_pid to clear"
            return 1
        fi
    else
        # Check for any vLLM processes
        local pids=$(pgrep -f "vllm.entrypoints.openai.api_server" 2>/dev/null || true)
        if [ -n "$pids" ]; then
            print_status "vLLM" "warning" "(running but no PID file)"
            print_info "PIDs: $pids"
            print_info "Consider creating .vllm_pid"
        else
            # Check if port is in use
            if lsof -i :$VLLM_PORT &> /dev/null 2>&1; then
                print_status "vLLM" "warning" "(port $VLLM_PORT in use by unknown process)"
            else
                print_status "vLLM" "stopped"
            fi
        fi
        return 1
    fi
}

check_docker() {
    print_section "Docker"

    if ! docker info &> /dev/null 2>&1; then
        print_status "Docker Daemon" "stopped"
        print_info "Docker is not running. Start it with: sudo systemctl start docker"
        return 1
    fi

    print_status "Docker Daemon" "running"

    # Get Docker version
    local version=$(docker --version | sed 's/Docker version //')
    print_info "Version: $version"

    # Check disk usage
    local df_output=$(docker system df --format "{{.Size}}" 2>/dev/null | head -1 || echo "N/A")
    print_info "Docker disk usage: $df_output"

    return 0
}

check_docker_compose() {
    print_section "Docker Compose Services"

    if [ ! -f "$COMPOSE_FILE" ]; then
        print_status "Docker Compose" "stopped" "(no compose file found)"
        print_info "Expected file: $COMPOSE_FILE"
        return 1
    fi

    if ! docker info &> /dev/null 2>&1; then
        print_status "Docker Compose" "stopped" "(docker not running)"
        return 1
    fi

    # Get container status
    local compose_ps=$(docker-compose ps 2>&1)
    local exit_code=$?

    if [ $exit_code -ne 0 ]; then
        print_status "Docker Compose" "stopped" "(cannot get status)"
        print_info "Error: $compose_ps"
        return 1
    fi

    # Parse and display each service
    local services_running=0
    local services_total=0

    # Get list of services
    local services=$(docker-compose ps --services 2>/dev/null || echo "")

    if [ -z "$services" ]; then
        print_status "Docker Compose" "stopped" "(no services defined)"
        return 1
    fi

    echo -e "\n  ${BOLD}Services:${NC}"

    for service in $services; do
        services_total=$((services_total + 1))

        # Get container info for this service
        local container_info=$(docker-compose ps -q "$service" 2>/dev/null)

        if [ -n "$container_info" ]; then
            # Check if container is running
            local state=$(docker inspect -f '{{.State.Status}}' "$container_info" 2>/dev/null || echo "unknown")

            if [ "$state" = "running" ]; then
                services_running=$((services_running + 1))

                # Get container details
                local container_name=$(docker inspect -f '{{.Name}}' "$container_info" 2>/dev/null | sed 's|/||')
                local ports=$(docker port "$container_info" 2>/dev/null | head -1 || echo "N/A")
                local uptime=$(docker inspect -f '{{.State.StartedAt}}' "$container_info" 2>/dev/null || echo "N/A")

                # Format uptime
                if [ "$uptime" != "N/A" ]; then
                    uptime=$(date -d "$uptime" "+%Y-%m-%d %H:%M:%S" 2>/dev/null || echo "$uptime")
                fi

                echo -e "    ${GREEN}●${NC} ${BOLD}$service${NC} (${GREEN}running${NC})"
                echo -e "       Container: $container_name"
                echo -e "       Port: ${ports:-N/A}"
                echo -e "       Started: $uptime"

                # Health check if available
                local health=$(docker inspect -f '{{.State.Health.Status}}' "$container_info" 2>/dev/null || echo "")
                if [ -n "$health" ]; then
                    if [ "$health" = "healthy" ]; then
                        echo -e "       Health: ${GREEN}$health${NC}"
                    else
                        echo -e "       Health: ${YELLOW}$health${NC}"
                    fi
                fi
            else
                echo -e "    ${RED}○${NC} ${BOLD}$service${NC} (${RED}$state${NC})"
            fi
        else
            echo -e "    ${YELLOW}?${NC} ${BOLD}$service${NC} (${YELLOW}not created${NC})"
        fi
    done

    echo
    print_info "Services running: $services_running/$services_total"

    if [ $services_running -eq $services_total ] && [ $services_total -gt 0 ]; then
        return 0
    else
        return 1
    fi
}

check_openwebui() {
    print_section "OpenWebUI"

    # Check if container is running via docker-compose
    if docker info &> /dev/null 2>&1 && [ -f "$COMPOSE_FILE" ]; then
        local openwebui_container=$(docker-compose ps -q openwebui 2>/dev/null || echo "")

        if [ -n "$openwebui_container" ]; then
            local state=$(docker inspect -f '{{.State.Status}}' "$openwebui_container" 2>/dev/null || echo "unknown")

            if [ "$state" = "running" ]; then
                # Try to access the web interface
                local openwebui_url="http://localhost:${OPENWEBUI_PORT}"

                if curl -sS -f -m 5 "$openwebui_url" &> /dev/null || \
                   curl -sS -f -m 5 "${openwebui_url}/api/health" &> /dev/null; then
                    print_status "OpenWebUI" "running" "(accessible)"
                    print_info "URL: http://localhost:${OPENWEBUI_PORT}"

                    # Get image version
                    local image=$(docker inspect -f '{{.Config.Image}}' "$openwebui_container" 2>/dev/null || echo "N/A")
                    print_info "Image: $image"
                    return 0
                else
                    print_status "OpenWebUI" "warning" "(container running but not accessible)"
                    print_info "Container is running but web interface is not responding"
                    return 1
                fi
            else
                print_status "OpenWebUI" "stopped" "(container $state)"
                return 1
            fi
        else
            print_status "OpenWebUI" "stopped" "(container not found)"
            return 1
        fi
    else
        print_status "OpenWebUI" "stopped" "(docker-compose not available)"
        return 1
    fi
}

check_cliproxyapi() {
    print_section "CLIProxyAPI"

    if [ "$CLIPROXYAPI_ENABLED" = "false" ]; then
        print_status "CLIProxyAPI" "stopped" "(disabled by CLIPROXYAPI_ENABLED=false)"
        print_info "Set CLIPROXYAPI_ENABLED=true to restore required upstream routing"
        return 1
    fi

    if [ ! -x "$CLIPROXYAPI_CHECK_SCRIPT" ]; then
        print_status "CLIProxyAPI" "stopped" "(check script not found)"
        print_info "Expected script: $CLIPROXYAPI_CHECK_SCRIPT"
        return 1
    fi

    local service_running=false
    local service_pid=""
    if [ -f "$CLIPROXYAPI_PID_FILE" ]; then
        service_pid=$(cat "$CLIPROXYAPI_PID_FILE")
        if [[ "$service_pid" =~ ^[0-9]+$ ]] && kill -0 "$service_pid" 2>/dev/null; then
            service_running=true
        fi
    fi

    if CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1; then
        if [ "$service_running" = true ]; then
            print_status "CLIProxyAPI" "running" "(PID: $service_pid)"
        elif is_true "$CLIPROXYAPI_DOCKER_MANAGED"; then
            print_status "CLIProxyAPI" "running" "(docker-managed)"
            print_info "Compose service: cliproxyapi"
        else
            print_status "CLIProxyAPI" "warning" "(health OK but PID file missing)"
            print_info "Run ./start-cliproxyapi.sh to launch managed mode"
        fi
        print_info "Health command: ./check-cliproxyapi.sh"
        print_info "Log file: $CLIPROXYAPI_LOG_FILE"
        return 0
    fi

    if [ "$service_running" = true ]; then
        print_status "CLIProxyAPI" "warning" "(PID: $service_pid exists but health check failing)"
        print_info "Check CLIProxyAPI logs: $CLIPROXYAPI_LOG_FILE"
    else
        print_status "CLIProxyAPI" "stopped"
    fi

    return 1
}

check_baseline_profile() {
    print_section "OpenWebUI Baseline Profile"

    print_info "RAG_EMBEDDING_MODEL=${RAG_EMBEDDING_MODEL:-<unset>}"
    print_info "RAG_TOP_K=${RAG_TOP_K:-<unset>} | TOP_K=${TOP_K:-<unset>}"
    print_info "CHUNK_SIZE=${CHUNK_SIZE:-<unset>} | CHUNK_OVERLAP=${CHUNK_OVERLAP:-<unset>}"
    print_info "VECTOR_DB=${VECTOR_DB:-<unset>}"

    local web_search_enabled=false
    if is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; then
        web_search_enabled=true
    fi

    if [ "$web_search_enabled" = true ]; then
        print_info "Web search: enabled"
        print_info "SEARXNG_QUERY_URL=${SEARXNG_QUERY_URL:-<unset>}"

        if [ -n "${SEARXNG_QUERY_URL:-}" ]; then
            local probe_url="${SEARXNG_QUERY_URL//\{query\}/status%20probe}"
            local host_probe_url="$probe_url"
            local probe_response
            local container_probe_ok=false

            if echo "$host_probe_url" | grep -q "host.docker.internal"; then
                host_probe_url="${host_probe_url//host.docker.internal/127.0.0.1}"
            fi

            probe_response=$(curl -sS -m 12 "$host_probe_url" 2>/dev/null || true)

            if [ -z "$probe_response" ] && [ "$host_probe_url" != "$probe_url" ]; then
                probe_response=$(curl -sS -m 12 "$probe_url" 2>/dev/null || true)
            fi

            if [ -n "$probe_response" ] && echo "$probe_response" | grep -q "results"; then
                print_status "SearXNG (host-local)" "running" "(reachable on configured URL)"
            else
                print_status "SearXNG (host-local)" "warning" "(configured but probe failed)"
                print_info "Ensure local SearXNG is listening on port 8888"
            fi

            # OpenWebUI runs in Docker; the important path is container -> host SearXNG.
            # Some hosts firewall docker-to-host connections, so check from inside container.
            if docker info >/dev/null 2>&1 && docker ps -q -f "name=^openwebui$" | grep -q .; then
                local container_probe_url="$probe_url"
                local container_probe
                container_probe="$(docker exec openwebui sh -lc "curl -sS -m 8 \"$container_probe_url\" 2>/dev/null || true" 2>/dev/null || true)"

                if [ -n "$container_probe" ] && echo "$container_probe" | grep -q "results"; then
                    container_probe_ok=true
                fi

                if [ "$container_probe_ok" = true ]; then
                    print_status "SearXNG (container path)" "running" "(OpenWebUI can reach host-local SearXNG)"
                else
                    print_status "SearXNG (container path)" "warning" "(OpenWebUI cannot reach host-local SearXNG)"
                    print_info "If host probe works but container probe fails, allow Docker networks to reach port 8888 on the host"
                fi
            fi
        else
            print_status "SearXNG (host-local)" "warning" "(web search enabled but URL is unset)"
        fi
    else
        print_info "Web search: disabled"
    fi
}

check_resources() {
    print_section "System Resources"

    # Check vLLM log file size
    if [ -f "$VLLM_LOG_FILE" ]; then
        local log_size=$(du -h "$VLLM_LOG_FILE" | cut -f1)
        local log_lines=$(wc -l < "$VLLM_LOG_FILE")
        print_info "vLLM log file: $log_size ($log_lines lines)"
    fi

    # Check disk space
    local disk_usage=$(df -h . | tail -1 | awk '{print $5 " used, " $4 " available"}')
    print_info "Disk space: $disk_usage"

    # Check memory (if free command available)
    if command -v free &> /dev/null; then
        local mem_info=$(free -h | awk '/^Mem:/ {print $3 "/" $2 " used"}')
        print_info "Memory: $mem_info"
    fi

    # Check GPU (if nvidia-smi available)
    if command -v nvidia-smi &> /dev/null; then
        local gpu_info=$(nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader 2>/dev/null | head -1 || echo "N/A")
        if [ "$gpu_info" != "N/A" ]; then
            print_info "GPU: $gpu_info"
        fi
    fi
}

print_summary() {
    echo -e "\n${CYAN}${BOLD}════════════════════════════════════════════════════════════${NC}"
    echo -e "${CYAN}${BOLD}                      Summary                                  ${NC}"
    echo -e "${CYAN}${BOLD}════════════════════════════════════════════════════════════${NC}\n"

    local all_good=true
    local cliproxy_ready=false

    # Check vLLM
    if [ -f .vllm_pid ] && kill -0 "$(cat .vllm_pid)" 2>/dev/null && curl -sS -f -m 5 "$VLLM_HEALTH_URL" &> /dev/null; then
        echo -e "  ${GREEN}●${NC} vLLM: ${GREEN}Operational${NC}"
    elif cliproxyapi_quiet_healthy; then
        echo -e "  ${CYAN}◌${NC} vLLM: ${CYAN}Not running (optional)${NC}"
    else
        echo -e "  ${RED}○${NC} vLLM: ${RED}Not operational${NC}"
        all_good=false
    fi

    # Check Docker
    if docker info &> /dev/null 2>&1; then
        echo -e "  ${GREEN}●${NC} Docker: ${GREEN}Running${NC}"
    else
        echo -e "  ${RED}○${NC} Docker: ${RED}Not running${NC}"
        all_good=false
    fi

    # Check OpenWebUI
    if curl -sS -f -m 5 "http://localhost:${OPENWEBUI_PORT}" &> /dev/null 2>&1; then
        echo -e "  ${GREEN}●${NC} OpenWebUI: ${GREEN}Accessible${NC}"
    else
        echo -e "  ${RED}○${NC} OpenWebUI: ${RED}Not accessible${NC}"
        all_good=false
    fi

    # Check CLIProxyAPI
    if [ "$CLIPROXYAPI_ENABLED" = "false" ]; then
        echo -e "  ${RED}○${NC} CLIProxyAPI: ${RED}Disabled (not allowed)${NC}"
        all_good=false
    elif cliproxyapi_quiet_healthy; then
        echo -e "  ${GREEN}●${NC} CLIProxyAPI: ${GREEN}Operational${NC}"
        cliproxy_ready=true
    else
        echo -e "  ${RED}○${NC} CLIProxyAPI: ${RED}Not operational${NC}"
        all_good=false
    fi

    echo

    if [ "$all_good" = true ]; then
        echo -e "  ${GREEN}${BOLD}✓ All systems operational!${NC}"
        echo -e "  ${CYAN}Access OpenWebUI at: http://localhost:${OPENWEBUI_PORT}${NC}"
        if [ "$cliproxy_ready" = true ]; then
            echo -e "  ${CYAN}Primary upstream: CLIProxyAPI (vLLM optional)${NC}"
        fi
    else
        echo -e "  ${YELLOW}${BOLD}⚠ Some services are not running properly${NC}"
        echo -e "  ${YELLOW}Run ./deploy.sh to start all services${NC}"
    fi

    echo
}

###############################################################################
# Main Status Check
###############################################################################

main() {
    print_header

    # Parse command line arguments
    WATCH_MODE=false
    QUIET_MODE=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            -w|--watch)
                WATCH_MODE=true
                shift
                ;;
            -q|--quiet)
                QUIET_MODE=true
                shift
                ;;
            -h|--help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  -w, --watch    Watch mode - refresh status every 5 seconds"
                echo "  -q, --quiet    Quiet mode - show summary only"
                echo "  -h, --help     Show this help message"
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
    done

    if [ "$WATCH_MODE" = true ]; then
        while true; do
            clear
            print_header

            if [ "$QUIET_MODE" = false ]; then
                check_vllm || true
                check_docker || true
                check_docker_compose || true
                check_openwebui || true
                check_cliproxyapi || true
                check_baseline_profile || true
                check_resources || true
            fi

            print_summary

            echo -e "${CYAN}Press Ctrl+C to exit${NC}"
            sleep 5
        done
    else
        if [ "$QUIET_MODE" = false ]; then
            check_vllm || true
            check_docker || true
            check_docker_compose || true
            check_openwebui || true
            check_cliproxyapi || true
            check_baseline_profile || true
            check_resources || true
        fi

        print_summary
    fi
}

main "$@"
