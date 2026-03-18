#!/bin/bash

###############################################################################
# OpenWebUI + LiteLLM RAG Cleanup Script
# Stops the current compose baseline and optionally removes logs, cache, and
# compose-managed data.
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

COMPOSE_FILE="config/compose/docker-compose.yml"
OPENWEBUI_IMAGE="${OPENWEBUI_IMAGE:-ghcr.io/open-webui/open-webui:v0.8.10}"

confirm_action() {
    local message=$1
    local default=${2:-n}

    if [ "${FORCE:-false}" = true ]; then
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
            [Yy] | [Yy][Ee][Ss]) return 0 ;;
            [Nn] | [Nn][Oo]) return 1 ;;
            *) echo -e "${RED}Please respond with yes or no${NC}" ;;
        esac
    done
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

    if ! docker_compose ps -q >/dev/null 2>&1; then
        print_error "Failed to query Docker Compose services"
        return 1
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
    if ! confirm_action "Remove Docker volumes and images? This will delete all compose-managed data."; then
        return 0
    fi

    print_step "Cleaning up Docker resources..."

    if ! init_compose_cmd; then
        print_error "Docker Compose is not available; cannot remove compose-managed resources"
        return 1
    fi

    if ! docker_compose down -v --remove-orphans; then
        print_error "Failed to remove Docker Compose containers/volumes"
        return 1
    fi

    if confirm_action "Remove OpenWebUI Docker image as well?"; then
        print_info "Removing OpenWebUI image..."
        docker rmi "$OPENWEBUI_IMAGE" 2>/dev/null || print_info "Image not found or already removed"
    fi

    print_info "Cleaning up dangling images and build cache..."
    docker image prune -f 2>/dev/null || true
    docker builder prune -f 2>/dev/null || true

    print_success "Docker cleanup complete"
}

cleanup_logs() {
    if ! confirm_action "Clear tracked log files under logs/?"; then
        return 0
    fi

    print_step "Cleaning up log files..."
    if [ -d logs ]; then
        find logs -maxdepth 1 -type f -name '*.log' -exec truncate -s 0 {} \;
    fi
    print_info "Docker daemon logs are not truncated by this script"
    print_success "Log cleanup complete"
}

cleanup_cache() {
    if ! confirm_action "Clean up cache files?"; then
        return 0
    fi

    print_step "Cleaning up cache files..."
    find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
    find . -type f -name "*.pyc" -delete 2>/dev/null || true
    find . -type f -name "*.tmp" -delete 2>/dev/null || true
    find . -type f -name ".DS_Store" -delete 2>/dev/null || true
    print_success "Cache cleanup complete"
}

show_cleanup_summary() {
    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}${BOLD}              🧹 Cleanup Complete! 🧹                         ${NC}"
    echo -e "${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"
    echo -e "${BOLD}Services Stopped:${NC}"
    echo -e "  ${CYAN}Docker Compose:${NC} Stopped"
    echo -e "\n${BOLD}To restart services:${NC}"
    echo -e "  ${YELLOW}./deploy.sh${NC}"
    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"
}

main() {
    print_header

    CLEAN_DOCKER=false
    CLEAN_LOGS=false
    CLEAN_CACHE=false
    FORCE=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            -f | --force)
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
                shift
                ;;
            -h | --help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Stop the OpenWebUI compose baseline and optionally clean data, logs, and cache."
                echo ""
                echo "Options:"
                echo "  -f, --force    Skip confirmation prompts"
                echo "  --docker       Remove Docker volumes and images"
                echo "  --logs         Clear log files"
                echo "  --cache        Clean cache files"
                echo "  --all          Perform all cleanup operations"
                echo "  -h, --help     Show this help message"
                exit 0
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
    done

    echo -e "${BOLD}This script will:${NC}"
    echo "  • Stop Docker Compose services"
    [ "$CLEAN_DOCKER" = true ] && echo "  • Remove Docker volumes and images"
    [ "$CLEAN_LOGS" = true ] && echo "  • Clear log files"
    [ "$CLEAN_CACHE" = true ] && echo "  • Clean cache files"
    echo

    if ! confirm_action "Proceed with cleanup?"; then
        print_info "Cleanup cancelled"
        exit 0
    fi

    stop_docker_compose
    [ "$CLEAN_DOCKER" = true ] && cleanup_docker
    [ "$CLEAN_LOGS" = true ] && cleanup_logs
    [ "$CLEAN_CACHE" = true ] && cleanup_cache

    show_cleanup_summary
}

main "$@"
