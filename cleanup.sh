#!/bin/bash

###############################################################################
# OpenWebUI + LiteLLM RAG Cleanup Script
# This script stops all services, including any optional legacy vLLM fallback, and performs cleanup
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

# Configuration
# shellcheck disable=SC2034  # VLLM_PORT used by lib/init.sh
VLLM_PORT=8000
COMPOSE_FILE="docker-compose.yml"
VLLM_LOG_FILE="logs/vllm.log"
VLLM_PID_FILE="vllm.pid"
VLLM_PID_FILE_LEGACY=".vllm_pid"
OPENWEBUI_IMAGE="${OPENWEBUI_IMAGE:-ghcr.io/open-webui/open-webui:v0.8.10}"
CLIPROXYAPI_LOG_FILE="${CLIPROXYAPI_LOG_FILE:-logs/cliproxyapi.log}"
CLIPROXYAPI_PID_FILE="${CLIPROXYAPI_PID_FILE:-.cliproxyapi.pid}"
CLIPROXYAPI_STOP_SCRIPT="./stop-cliproxyapi.sh"

###############################################################################
# Cleanup Functions
###############################################################################

stop_vllm() {
    print_step "Stopping archived vLLM fallback (repo-managed only)..."

    has_vllm_command() {
        local pid="$1"
        kill -0 "$pid" 2>/dev/null &&
            ps -p "$pid" -o args= 2>/dev/null | grep -q "vllm.entrypoints.openai.api_server"
    }

    local stopped=false
    local pid=""
    local pid_file=""

    if [ -f "$VLLM_PID_FILE" ]; then
        pid_file="$VLLM_PID_FILE"
    elif [ -f "$VLLM_PID_FILE_LEGACY" ]; then
        pid_file="$VLLM_PID_FILE_LEGACY"
    fi

    if [ -z "$pid_file" ]; then
        print_info "No repo-managed archived vLLM PID file found; leaving unrelated host vLLM processes untouched"
        return 0
    fi

    pid=$(cat "$pid_file" 2>/dev/null || true)
    if [[ $pid =~ ^[0-9]+$ ]] && has_vllm_command "$pid"; then
        print_info "Stopping archived vLLM fallback process (PID: $pid)..."
        kill "$pid"
        sleep 2

        if has_vllm_command "$pid"; then
            print_warning "Process still running, sending SIGKILL..."
            kill -9 "$pid"
            sleep 1
        fi

        stopped=true
    elif [ -n "$pid" ]; then
        print_warning "Ignoring stale archived vLLM PID file: $pid_file ($pid)"
    fi

    rm -f "$VLLM_PID_FILE" "$VLLM_PID_FILE_LEGACY"

    if [ "$stopped" = true ]; then
        print_success "Archived vLLM fallback stopped"
    else
        print_info "Archived vLLM fallback was not running"
    fi
}

stop_cliproxyapi() {
    print_step "Stopping CLIProxyAPI..."

    if [ -x "$CLIPROXYAPI_STOP_SCRIPT" ]; then
        if "$CLIPROXYAPI_STOP_SCRIPT"; then
            return 0
        fi
        print_warning "CLIProxyAPI stop script reported an error, attempting fallback stop"
    fi

    local stopped=false

    if [ -f "$CLIPROXYAPI_PID_FILE" ]; then
        local pid
        pid=$(cat "$CLIPROXYAPI_PID_FILE")
        if [[ $pid =~ ^[0-9]+$ ]] && kill -0 "$pid" 2>/dev/null; then
            print_info "Stopping CLIProxyAPI process (PID: $pid)..."
            kill "$pid"
            sleep 1
            if kill -0 "$pid" 2>/dev/null; then
                print_warning "Process still running, sending SIGKILL..."
                kill -9 "$pid"
            fi
            stopped=true
        fi
        rm -f "$CLIPROXYAPI_PID_FILE"
    fi

    if [ "$stopped" = true ]; then
        print_success "CLIProxyAPI stopped"
    else
        print_info "CLIProxyAPI was not running"
    fi
}

stop_docker_compose() {
    print_step "Stopping Docker Compose services..."

    if [ ! -f "$COMPOSE_FILE" ]; then
        print_info "No Docker Compose file found"
        return 0
    fi

    if ! docker info &>/dev/null 2>&1; then
        print_info "Docker is not running"
        return 0
    fi

    if ! init_compose_cmd; then
        print_error "Docker Compose is not available; cannot stop services"
        return 1
    fi

    # Check if any containers are running
    local running
    if ! running=$(docker_compose ps -q 2>/dev/null); then
        print_error "Failed to query Docker Compose services"
        return 1
    fi
    if [ -z "$running" ]; then
        print_info "No Docker Compose containers are running"
        return 0
    fi

    print_info "Stopping Docker Compose..."
    if docker_compose down; then
        print_success "Docker Compose services stopped"
    else
        print_error "Failed to stop Docker Compose services"
        return 1
    fi
}

cleanup_docker() {
    if ! confirm_action "Remove Docker volumes and images? This will delete all data."; then
        return 0
    fi

    print_step "Cleaning up Docker resources..."

    if ! init_compose_cmd; then
        print_error "Docker Compose is not available; cannot remove compose-managed resources"
        return 1
    fi

    # Remove containers
    print_info "Removing containers..."
    if ! docker_compose down -v --remove-orphans; then
        print_error "Failed to remove Docker Compose containers/volumes"
        return 1
    fi

    # Ask about images
    if confirm_action "Remove OpenWebUI Docker image as well?"; then
        print_info "Removing OpenWebUI image..."
        docker rmi "$OPENWEBUI_IMAGE" 2>/dev/null ||
            docker rmi ghcr.io/open-webui/open-webui:v0.8.10 2>/dev/null ||
            docker rmi ghcr.io/open-webui/open-webui:v0.8.7 2>/dev/null ||
            docker rmi ghcr.io/open-webui/open-webui:main 2>/dev/null ||
            docker rmi ghcr.io/open-webui/open-webui:latest 2>/dev/null ||
            print_info "Image not found or already removed"
    fi

    # Clean up dangling images and build cache
    print_info "Cleaning up dangling images and build cache..."
    docker image prune -f 2>/dev/null || true
    docker builder prune -f 2>/dev/null || true

    print_success "Docker cleanup complete"
}

cleanup_logs() {
    if ! confirm_action "Clear all log files?"; then
        return 0
    fi

    print_step "Cleaning up log files..."

    # Clear legacy vLLM fallback log
    if [ -f "$VLLM_LOG_FILE" ]; then
        print_info "Clearing vLLM fallback log file..."
        : >"$VLLM_LOG_FILE"
        print_success "vLLM fallback log cleared"
    fi

    # Clear CLIProxyAPI log
    if [ -f "$CLIPROXYAPI_LOG_FILE" ]; then
        print_info "Clearing CLIProxyAPI log file..."
        : >"$CLIPROXYAPI_LOG_FILE"
        print_success "CLIProxyAPI log cleared"
    fi

    # Clear Docker logs (requires container restart)
    print_info "Note: Docker logs are managed by Docker daemon"
    print_info "To clear Docker logs completely, restart the services"

    print_success "Log cleanup complete"
}

cleanup_cache() {
    if ! confirm_action "Clean up cache files?"; then
        return 0
    fi

    print_step "Cleaning up cache files..."

    # Clean Python cache
    print_info "Removing Python cache files..."
    find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
    find . -type f -name "*.pyc" -delete 2>/dev/null || true

    # Clean temporary files
    print_info "Removing temporary files..."
    find . -type f -name "*.tmp" -delete 2>/dev/null || true
    find . -type f -name ".DS_Store" -delete 2>/dev/null || true

    print_success "Cache cleanup complete"
}

show_cleanup_summary() {
    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}${BOLD}              🧹 Cleanup Complete! 🧹                         ${NC}"
    echo -e "${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"

    echo -e "${BOLD}Services Stopped:${NC}"
    echo -e "  ${CYAN}CLIProxyAPI:${NC}   Stopped"
    echo -e "  ${CYAN}vLLM fallback:${NC} Stopped"
    echo -e "  ${CYAN}Docker Compose:${NC} Stopped"

    echo -e "\n${BOLD}To restart services:${NC}"
    echo -e "  ${YELLOW}./deploy.sh${NC}"

    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"
}

###############################################################################
# Main Cleanup Flow
###############################################################################

main() {
    print_header

    # Parse command line arguments
    CLEAN_DOCKER=false
    CLEAN_LOGS=false
    CLEAN_CACHE=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            -f | --force)
                shift
                ;;
            --docker)
                CLEAN_DOCKER=true
                shift
                ;;
            --logs)
                CLEAN_LOGS=true
                shift
                ;;
            --cache)
                CLEAN_CACHE=true
                shift
                ;;
            --all)
                CLEAN_DOCKER=true
                CLEAN_LOGS=true
                CLEAN_CACHE=true
                shift
                ;;
            -h | --help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Stop and cleanup the OpenWebUI stack, including any optional legacy vLLM fallback."
                echo ""
                echo "Options:"
                echo "  -f, --force    Skip confirmation prompts"
                echo "  --docker       Remove Docker volumes and images"
                echo "  --logs         Clear log files"
                echo "  --cache        Clean cache files"
                echo "  --all          Perform all cleanup operations"
                echo "  -h, --help     Show this help message"
                echo ""
                echo "Examples:"
                echo "  $0                  # Stop all services (default)"
                echo "  $0 -f              # Stop without confirmation"
                echo "  $0 --all           # Stop and clean everything"
                echo "  $0 --docker --logs # Stop and clean Docker and logs"
                exit 0
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
    done

    # Show what will be done
    echo -e "${BOLD}This script will:${NC}"
    echo "  • Stop CLIProxyAPI service"
    echo "  • Stop optional legacy vLLM fallback"
    echo "  • Stop Docker Compose services"

    if [ "$CLEAN_DOCKER" = true ]; then
        echo "  • Remove Docker volumes and images"
    fi

    if [ "$CLEAN_LOGS" = true ]; then
        echo "  • Clear log files"
    fi

    if [ "$CLEAN_CACHE" = true ]; then
        echo "  • Clean cache files"
    fi

    echo

    # Confirm before proceeding
    if ! confirm_action "Proceed with cleanup?"; then
        print_info "Cleanup cancelled"
        exit 0
    fi

    # Stop services
    stop_cliproxyapi
    stop_vllm
    stop_docker_compose

    # Perform additional cleanup if requested
    if [ "$CLEAN_DOCKER" = true ]; then
        cleanup_docker
    fi

    if [ "$CLEAN_LOGS" = true ]; then
        cleanup_logs
    fi

    if [ "$CLEAN_CACHE" = true ]; then
        cleanup_cache
    fi

    # Show summary
    show_cleanup_summary
}

main "$@"
