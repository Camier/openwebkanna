#!/bin/bash

###############################################################################
# OpenWebUI + vLLM RAG Cleanup Script
# This script stops all services and performs cleanup
###############################################################################

set -e

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
COMPOSE_FILE="docker-compose.yml"
VLLM_LOG_FILE="logs/vllm.log"
VLLM_PID_FILE=".vllm_pid"
CLIPROXYAPI_LOG_FILE="${CLIPROXYAPI_LOG_FILE:-logs/cliproxyapi.log}"
CLIPROXYAPI_PID_FILE="${CLIPROXYAPI_PID_FILE:-.cliproxyapi.pid}"
CLIPROXYAPI_STOP_SCRIPT="./stop-cliproxyapi.sh"

###############################################################################
# Helper Functions
###############################################################################

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     OpenWebUI + vLLM RAG Cleanup Manager                  ║"
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

confirm_action() {
    local message=$1
    local default=${2:-n}

    if [ "$FORCE" = true ]; then
        return 0
    fi

    local prompt
    if [ "$default" = "y" ]; then
        prompt="$message [Y/n]"
    else
        prompt="$message [y/N]"
    fi

    while true; do
        echo -en "${YELLOW}$prompt ${NC}"
        read -r response
        response=${response:-$default}

        case $response in
            [Yy]|[Yy][Ee][Ss])
                return 0
                ;;
            [Nn]|[Nn][Oo])
                return 1
                ;;
            *)
                echo -e "${RED}Please respond with yes or no${NC}"
                ;;
        esac
    done
}

###############################################################################
# Cleanup Functions
###############################################################################

stop_vllm() {
    print_step "Stopping vLLM..."

    local stopped=false

    # Check PID file first
    if [ -f "$VLLM_PID_FILE" ]; then
        local pid=$(cat "$VLLM_PID_FILE")
        if kill -0 "$pid" 2>/dev/null; then
            print_info "Stopping vLLM process (PID: $pid)..."
            kill "$pid"
            sleep 2

            # Force kill if still running
            if kill -0 "$pid" 2>/dev/null; then
                print_warning "Process still running, sending SIGKILL..."
                kill -9 "$pid"
                sleep 1
            fi

            stopped=true
        fi
        rm -f "$VLLM_PID_FILE"
    fi

    # Check for any vLLM processes
    local pids=$(pgrep -f "vllm.entrypoints.openai.api_server" 2>/dev/null || true)
    if [ -n "$pids" ]; then
        print_info "Stopping additional vLLM processes: $pids"
        echo "$pids" | xargs kill
        sleep 2

        # Force kill if needed
        pids=$(pgrep -f "vllm.entrypoints.openai.api_server" 2>/dev/null || true)
        if [ -n "$pids" ]; then
            print_warning "Some processes still running, force killing..."
            echo "$pids" | xargs kill -9
        fi

        stopped=true
    fi

    if [ "$stopped" = true ]; then
        print_success "vLLM stopped"
    else
        print_info "vLLM was not running"
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
        local pid=$(cat "$CLIPROXYAPI_PID_FILE")
        if [[ "$pid" =~ ^[0-9]+$ ]] && kill -0 "$pid" 2>/dev/null; then
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

    if ! docker info &> /dev/null 2>&1; then
        print_info "Docker is not running"
        return 0
    fi

    # Check if any containers are running
    local running=$(docker-compose ps -q 2>/dev/null || true)
    if [ -z "$running" ]; then
        print_info "No Docker Compose containers are running"
        return 0
    fi

    print_info "Stopping Docker Compose..."
    docker-compose down

    if [ $? -eq 0 ]; then
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

    # Remove containers
    print_info "Removing containers..."
    docker-compose down -v --remove-orphans 2>/dev/null || true

    # Ask about images
    if confirm_action "Remove OpenWebUI Docker image as well?"; then
        print_info "Removing OpenWebUI image..."
        docker rmi ghcr.io/open-webui/open-webui:main 2>/dev/null || \
        docker rmi ghcr.io/open-webui/open-webui:latest 2>/dev/null || \
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

    # Clear vLLM log
    if [ -f "$VLLM_LOG_FILE" ]; then
        print_info "Clearing vLLM log file..."
        > "$VLLM_LOG_FILE"
        print_success "vLLM log cleared"
    fi

    # Clear CLIProxyAPI log
    if [ -f "$CLIPROXYAPI_LOG_FILE" ]; then
        print_info "Clearing CLIProxyAPI log file..."
        > "$CLIPROXYAPI_LOG_FILE"
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
    echo -e "  ${CYAN}vLLM:${NC}          Stopped"
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
    FORCE=false
    CLEAN_DOCKER=false
    CLEAN_LOGS=false
    CLEAN_CACHE=false
    STOP_ONLY=true

    while [[ $# -gt 0 ]]; do
        case $1 in
            -f|--force)
                FORCE=true
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
                STOP_ONLY=false
                shift
                ;;
            -h|--help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Stop and cleanup OpenWebUI + vLLM RAG system."
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
    echo "  • Stop vLLM server"
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
