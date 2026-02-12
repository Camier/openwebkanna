#!/bin/bash

###############################################################################
# OpenWebUI + vLLM RAG Deployment Script
# This script deploys and manages the complete RAG system
###############################################################################

set -e  # Exit on error

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
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Configuration
VLLM_PORT=8000
VLLM_HEALTH_URL="http://localhost:${VLLM_PORT}/health"
VLLM_MODEL="${VLLM_MODEL:-meta-llama/Llama-3.1-8B-Instruct}"
OPENWEBUI_PORT="${WEBUI_PORT:-3000}"
COMPOSE_FILE="docker-compose.yml"
VLLM_LOG_FILE="logs/vllm.log"
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_LOG_FILE="${CLIPROXYAPI_LOG_FILE:-logs/cliproxyapi.log}"
CLIPROXYAPI_START_SCRIPT="./start-cliproxyapi.sh"
CLIPROXYAPI_CHECK_SCRIPT="./check-cliproxyapi.sh"
CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"
CLIPROXYAPI_DOCKER_MANAGED="${CLIPROXYAPI_DOCKER_MANAGED:-true}"
CLIPROXYAPI_DOCKER_SERVICE="${CLIPROXYAPI_DOCKER_SERVICE:-cliproxyapi}"
AUTO_SKIP_VLLM_IF_CLIPROXYAPI_HEALTHY="${AUTO_SKIP_VLLM_IF_CLIPROXYAPI_HEALTHY:-true}"
AUTO_START_VLLM_ON_CLIPROXYAPI_FAILURE="${AUTO_START_VLLM_ON_CLIPROXYAPI_FAILURE:-true}"

###############################################################################
# Helper Functions
###############################################################################

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     OpenWebUI + vLLM RAG Deployment Manager                ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}" >&2
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${MAGENTA}ℹ $1${NC}"
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        print_error "$1 is not installed. Please install it first."
        exit 1
    fi
}

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

wait_for_service() {
    local url=$1
    local service_name=$2
    local max_attempts=${3:-60}
    local attempt=1

    print_step "Waiting for $service_name to be ready..."

    while [ $attempt -le $max_attempts ]; do
        if curl -s -f "$url" &> /dev/null; then
            print_success "$service_name is ready!"
            return 0
        fi
        echo -n "."
        sleep 2
        attempt=$((attempt + 1))
    done

    echo
    print_error "$service_name failed to start within expected time."
    return 1
}

###############################################################################
# vLLM Management
###############################################################################

start_vllm() {
    print_step "Checking vLLM status..."

    if curl -s -f "$VLLM_HEALTH_URL" &> /dev/null; then
        print_success "vLLM is already running on port ${VLLM_PORT}"
        return 0
    fi

    print_info "vLLM is not running. Starting vLLM server..."

    # Create logs directory if it doesn't exist
    mkdir -p logs

    # Check if vLLM is installed
    if ! python3 -c "import vllm" 2>/dev/null; then
        print_error "vLLM is not installed. Please install it first:"
        echo "  pip install vllm"
        exit 1
    fi

    # Start vLLM in background
    nohup python3 -m vllm.entrypoints.openai.api_server \
        --model "$VLLM_MODEL" \
        --host 0.0.0.0 \
        --port $VLLM_PORT \
        > "$VLLM_LOG_FILE" 2>&1 &

    local vllm_pid=$!
    echo $vllm_pid > .vllm_pid

    print_info "vLLM started with PID: $vllm_pid"
    print_info "Logs are being written to: $VLLM_LOG_FILE"

    # Wait for vLLM to be ready
    if wait_for_service "$VLLM_HEALTH_URL" "vLLM" 90; then
        print_success "vLLM is healthy and responding"
        return 0
    else
        print_error "vLLM failed to start. Check logs at: $VLLM_LOG_FILE"
        return 1
    fi
}

stop_vllm() {
    if [ -f .vllm_pid ]; then
        local pid=$(cat .vllm_pid)
        if kill -0 "$pid" 2>/dev/null; then
            print_step "Stopping vLLM (PID: $pid)..."
            kill "$pid"
            rm -f .vllm_pid
            print_success "vLLM stopped"
        else
            print_warning "vLLM PID file exists but process is not running"
            rm -f .vllm_pid
        fi
    else
        print_info "vLLM PID file not found. Checking for running vLLM processes..."
        local pids=$(pgrep -f "vllm.entrypoints.openai.api_server" || true)
        if [ -n "$pids" ]; then
            print_step "Stopping vLLM processes: $pids"
            echo "$pids" | xargs kill
            print_success "vLLM stopped"
        else
            print_info "No vLLM processes found"
        fi
    fi
}

###############################################################################
# CLIProxyAPI Management
###############################################################################

start_cliproxyapi() {
    if [ "$CLIPROXYAPI_ENABLED" = "false" ]; then
        print_error "CLIPROXYAPI_ENABLED=false is not supported for this deployment"
        return 1
    fi

    if is_true "$CLIPROXYAPI_DOCKER_MANAGED"; then
        print_step "Starting CLIProxyAPI Docker service..."
        if ! docker-compose up -d "$CLIPROXYAPI_DOCKER_SERVICE"; then
            print_error "Failed to start Docker service: $CLIPROXYAPI_DOCKER_SERVICE"
            return 1
        fi

        if wait_for_service "${CLIPROXYAPI_BASE_URL%/}/" "CLIProxyAPI" 60 && cliproxyapi_is_healthy; then
            print_success "CLIProxyAPI Docker service is healthy"
            return 0
        fi

        print_error "CLIProxyAPI Docker service failed readiness checks"
        return 1
    fi

    if [ ! -x "$CLIPROXYAPI_CHECK_SCRIPT" ]; then
        print_error "CLIProxyAPI check script not found or not executable: $CLIPROXYAPI_CHECK_SCRIPT"
        return 1
    fi

    if CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1; then
        print_success "CLIProxyAPI service is already healthy"
        return 0
    fi

    print_step "Starting CLIProxyAPI service..."

    if [ ! -x "$CLIPROXYAPI_START_SCRIPT" ]; then
        print_error "CLIProxyAPI start script not found or not executable: $CLIPROXYAPI_START_SCRIPT"
        return 1
    fi

    if CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_START_SCRIPT"; then
        if CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1; then
            print_success "CLIProxyAPI service started and healthy"
            return 0
        fi

        print_error "CLIProxyAPI service started but health checks are failing"
        return 1
    fi

    print_error "Failed to start CLIProxyAPI service"
    return 1
}

cliproxyapi_http_is_healthy() {
    local base_url="${CLIPROXYAPI_BASE_URL%/}"
    local health_url="${base_url}/"
    local models_url="${base_url}/v1/models"
    local code=""
    local payload_file=""

    code="$(curl -sS -m 8 -o /dev/null -w "%{http_code}" "$health_url" || true)"
    if [ "$code" != "200" ]; then
        return 1
    fi

    if [ -z "${CLIPROXYAPI_API_KEY:-}" ]; then
        return 0
    fi

    payload_file="$(mktemp)"
    code="$(curl -sS -m 10 -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" -o "$payload_file" -w "%{http_code}" "$models_url" || true)"
    if [ "$code" != "200" ]; then
        rm -f "$payload_file"
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        if ! jq -e '.data and (.data | type == "array") and (.data | length > 0)' "$payload_file" >/dev/null 2>&1; then
            rm -f "$payload_file"
            return 1
        fi
    else
        if ! grep -q '"id"' "$payload_file"; then
            rm -f "$payload_file"
            return 1
        fi
    fi

    rm -f "$payload_file"
    return 0
}

cliproxyapi_is_healthy() {
    if [ "$CLIPROXYAPI_ENABLED" = "false" ]; then
        return 1
    fi

    if is_true "$CLIPROXYAPI_DOCKER_MANAGED"; then
        cliproxyapi_http_is_healthy
        return $?
    fi

    if [ ! -x "$CLIPROXYAPI_CHECK_SCRIPT" ]; then
        return 1
    fi

    CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1
}

###############################################################################
# Docker Compose Management
###############################################################################

start_docker_compose() {
    print_step "Starting Docker Compose services..."

    if [ ! -f "$COMPOSE_FILE" ]; then
        print_error "Docker Compose file not found: $COMPOSE_FILE"
        exit 1
    fi

    # Check if Docker is running
    if ! docker info &> /dev/null; then
        print_error "Docker is not running. Please start Docker first."
        exit 1
    fi

    # Start services
    docker-compose up -d

    if [ $? -eq 0 ]; then
        print_success "Docker Compose services started"
    else
        print_error "Failed to start Docker Compose services"
        exit 1
    fi
}

wait_for_openwebui() {
    local openwebui_url="http://localhost:${OPENWEBUI_PORT}"
    print_step "Waiting for OpenWebUI to be ready..."

    local max_attempts=60
    local attempt=1

    while [ $attempt -le $max_attempts ]; do
        if curl -s -f "$openwebui_url" &> /dev/null || \
           curl -s -f "${openwebui_url}/api/health" &> /dev/null; then
            print_success "OpenWebUI is ready!"
            return 0
        fi
        echo -n "."
        sleep 2
        attempt=$((attempt + 1))
    done

    echo
    print_warning "OpenWebUI may still be starting up. Check logs with: ./logs.sh"
    return 1
}

show_docker_status() {
    print_step "Docker Containers Status:"
    docker-compose ps
}

show_access_info() {
    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}${BOLD}              🚀 Deployment Successful! 🚀                      ${NC}"
    echo -e "${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"

    echo -e "${BOLD}Service URLs:${NC}"
    echo -e "  ${CYAN}OpenWebUI:${NC}     http://localhost:${OPENWEBUI_PORT}"
    echo -e "  ${CYAN}vLLM API:${NC}      http://localhost:${VLLM_PORT}"
    echo -e "  ${CYAN}CLIProxyAPI:${NC}   ${CLIPROXYAPI_BASE_URL}"

    echo -e "\n${BOLD}Quick Commands:${NC}"
    echo -e "  ${YELLOW}View logs:${NC}     ./logs.sh"
    echo -e "  ${YELLOW}Check status:${NC}  ./status.sh"
    echo -e "  ${YELLOW}Stop all:${NC}     ./cleanup.sh"

    echo -e "\n${BOLD}vLLM Model:${NC} ${VLLM_MODEL}"
    echo -e "${BOLD}vLLM Logs:${NC}  ${VLLM_LOG_FILE}"
    echo -e "${BOLD}CLIProxyAPI Logs:${NC} ${CLIPROXYAPI_LOG_FILE}"

    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"
}

show_logs() {
    print_step "Recent logs from all services:"

    echo -e "\n${YELLOW}--- Docker Compose Logs ---${NC}"
    docker-compose logs --tail=20

    if [ -f "$VLLM_LOG_FILE" ]; then
        echo -e "\n${YELLOW}--- vLLM Logs (last 20 lines) ---${NC}"
        tail -20 "$VLLM_LOG_FILE"
    fi

    if [ -f "$CLIPROXYAPI_LOG_FILE" ]; then
        echo -e "\n${YELLOW}--- CLIProxyAPI Logs (last 20 lines) ---${NC}"
        tail -20 "$CLIPROXYAPI_LOG_FILE"
    fi
}

###############################################################################
# Main Deployment Flow
###############################################################################

main() {
    print_header

    # Parse command line arguments
    SKIP_VLLM=false
    SKIP_DOCKER=false
    SKIP_CLIPROXYAPI=false
    SHOW_LOGS=true
    VLLM_SKIP_REASON=""

    while [[ $# -gt 0 ]]; do
        case $1 in
            --skip-vllm)
                SKIP_VLLM=true
                VLLM_SKIP_REASON="--skip-vllm flag used"
                shift
                ;;
            --skip-docker)
                SKIP_DOCKER=true
                shift
                ;;
            --skip-cliproxyapi)
                SKIP_CLIPROXYAPI=true
                shift
                ;;
            --no-logs)
                SHOW_LOGS=false
                shift
                ;;
            -h|--help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --skip-vllm    Skip starting vLLM (use if already running)"
                echo "  --skip-docker   Skip starting Docker Compose (use if already running)"
                echo "  --skip-cliproxyapi Skip starting CLIProxyAPI service (must already be healthy)"
                echo "  --no-logs       Don't show logs after deployment"
                echo "  -h, --help      Show this help message"
                echo ""
                echo "Environment variables:"
                echo "  AUTO_SKIP_VLLM_IF_CLIPROXYAPI_HEALTHY  Auto-skip vLLM start if CLIProxyAPI is already healthy (default: true)"
                echo "  AUTO_START_VLLM_ON_CLIPROXYAPI_FAILURE Start vLLM only when CLIProxyAPI startup fails (default: true)"
                echo "  CLIPROXYAPI_DOCKER_MANAGED             Run CLIProxyAPI from docker-compose service (default: true)"
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
    done

    # Check prerequisites
    print_step "Checking prerequisites..."
    check_command "docker"
    check_command "docker-compose"
    check_command "curl"
    print_success "All prerequisites met"

    # Start CLIProxyAPI first to avoid unnecessary vLLM startups.
    CLIPROXYAPI_READY=false
    if [ "$SKIP_CLIPROXYAPI" = false ]; then
        if start_cliproxyapi; then
            CLIPROXYAPI_READY=true
            echo
        else
            print_warning "CLIProxyAPI startup failed before vLLM fallback decision"
        fi
    else
        print_info "Skipping CLIProxyAPI startup (--skip-cliproxyapi flag used)"
        if [ "$CLIPROXYAPI_ENABLED" = "false" ]; then
            print_error "CLIPROXYAPI_ENABLED=false is not supported for this deployment"
            exit 1
        fi
        if [ ! -x "$CLIPROXYAPI_CHECK_SCRIPT" ] || ! CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1; then
            print_error "CLIProxyAPI must already be healthy when using --skip-cliproxyapi"
            exit 1
        fi
        CLIPROXYAPI_READY=true
        print_success "CLIProxyAPI health check passed"
    fi

    if [ "$SKIP_VLLM" = false ] && [ "$AUTO_SKIP_VLLM_IF_CLIPROXYAPI_HEALTHY" = "true" ] && [ "$CLIPROXYAPI_READY" = true ]; then
        SKIP_VLLM=true
        VLLM_SKIP_REASON="CLIProxyAPI is healthy"
        print_info "CLIProxyAPI is healthy. Skipping vLLM startup."
    fi

    if [ "$CLIPROXYAPI_READY" = false ]; then
        if [ "$AUTO_START_VLLM_ON_CLIPROXYAPI_FAILURE" != "true" ] || [ "$SKIP_VLLM" = true ]; then
            print_error "CLIProxyAPI is not healthy and vLLM fallback is disabled"
            exit 1
        fi

        print_step "Attempting vLLM fallback because CLIProxyAPI is not healthy"
        if ! start_vllm; then
            print_error "Failed to start vLLM fallback"
            exit 1
        fi
        echo

        print_step "Retrying CLIProxyAPI startup after vLLM fallback"
        if ! start_cliproxyapi; then
            print_error "Failed to start CLIProxyAPI service after vLLM fallback"
            exit 1
        fi
        CLIPROXYAPI_READY=true
        echo
    else
        if [ "$SKIP_VLLM" = false ]; then
            print_info "Skipping vLLM startup (CLIProxyAPI is healthy and fallback is unnecessary)"
            SKIP_VLLM=true
            VLLM_SKIP_REASON="CLIProxyAPI healthy after startup"
        fi
    fi

    if [ "$SKIP_VLLM" = true ]; then
        if [ -n "$VLLM_SKIP_REASON" ]; then
            print_info "Skipping vLLM startup (${VLLM_SKIP_REASON})"
        else
            print_info "Skipping vLLM startup"
        fi
    fi

    # Start Docker Compose
    if [ "$SKIP_DOCKER" = false ]; then
        if start_docker_compose; then
            echo
        else
            print_error "Failed to start Docker Compose"
            exit 1
        fi
    else
        print_info "Skipping Docker Compose startup (--skip-docker flag used)"
    fi

    # Wait for services
    if ! wait_for_openwebui; then
        print_error "OpenWebUI failed readiness checks"
        exit 1
    fi

    # Show status
    show_docker_status

    # Show logs if requested
    if [ "$SHOW_LOGS" = true ]; then
        show_logs
    fi

    # Show access information
    show_access_info
}

# Run main function
main "$@"
