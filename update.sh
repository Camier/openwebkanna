#!/bin/bash

###############################################################################
# OpenWebUI + vLLM RAG Update Script
# This script updates OpenWebUI and restarts services
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

# Configuration
OPENWEBUI_IMAGE="${OPENWEBUI_IMAGE:-ghcr.io/open-webui/open-webui:v0.8.3}"
OPENWEBUI_ALT_IMAGE="${OPENWEBUI_ALT_IMAGE:-ghcr.io/open-webui/open-webui:v0.8.2}"
COMPOSE_FILE="docker-compose.yml"
BACKUP_DIR="backups"
OPENWEBUI_VOLUME="${OPENWEBUI_VOLUME:-openwebui_data}"
OPENWEBUI_PORT="${WEBUI_PORT:-3000}"
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_DOCKER_MANAGED="${CLIPROXYAPI_DOCKER_MANAGED:-true}"
CLIPROXYAPI_START_SCRIPT="./start-cliproxyapi.sh"
CLIPROXYAPI_RESTART_SCRIPT="./restart-cliproxyapi.sh"
CLIPROXYAPI_CHECK_SCRIPT="./check-cliproxyapi.sh"
COMPOSE_CMD=()

###############################################################################
# Helper Functions
###############################################################################

# Note: print_header, print_step, print_success, print_error, print_warning, print_info
# come from lib/print-utils.sh via lib/init.sh

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
            [Yy] | [Yy][Ee][Ss])
                return 0
                ;;
            [Nn] | [Nn][Oo])
                return 1
                ;;
            *)
                echo -e "${RED}Please respond with yes or no${NC}"
                ;;
        esac
    done
}

###############################################################################
# Backup Functions
###############################################################################

create_backup() {
    print_step "Creating backup..."

    # Create backup directory
    mkdir -p "$BACKUP_DIR"

    local timestamp
    local backup_file
    timestamp=$(date +%Y%m%d_%H%M%S)
    backup_file="$BACKUP_DIR/openwebui_backup_$timestamp.tar.gz"

    print_info "Backup file: $backup_file"

    # Check if Docker volume exists
    local volume_exists=false
    if docker info &>/dev/null 2>&1; then
        if docker volume ls --format '{{.Name}}' | grep -Fxq "${OPENWEBUI_VOLUME}"; then
            volume_exists=true
        fi
    fi

    if [ "$volume_exists" = false ]; then
        print_warning "No Docker volume found to backup"
        return 0
    fi

    # Create temporary container for backup
    print_info "Creating backup of OpenWebUI data..."

    if docker run --rm \
        -v "${OPENWEBUI_VOLUME}:/data" \
        -v "$(pwd)/$BACKUP_DIR:/backup" \
        alpine tar czf "/backup/$(basename "$backup_file")" -C /data . \
        2>/dev/null; then
        local size
        size=$(du -h "$backup_file" | cut -f1)
        print_success "Backup created successfully ($size)"

        # Keep only last 5 backups
        print_info "Cleaning old backups (keeping last 5)..."
        local -a backups=()
        mapfile -d '' -t backups < <(
            find "$BACKUP_DIR" -maxdepth 1 -type f -name 'openwebui_backup_*.tar.gz' -printf '%T@ %p\0' |
                sort -z -nr
        )
        if ((${#backups[@]} > 5)); then
            local backup_entry
            local backup_path
            for backup_entry in "${backups[@]:5}"; do
                backup_path="${backup_entry#* }"
                rm -f -- "$backup_path"
            done
        fi

        return 0
    else
        print_error "Failed to create backup"
        return 1
    fi
}

list_backups() {
    print_step "Available backups:"

    if [ ! -d "$BACKUP_DIR" ] || [ -z "$(ls -A "$BACKUP_DIR" 2>/dev/null)" ]; then
        print_info "No backups found"
        return 0
    fi

    echo
    for backup in "$BACKUP_DIR"/openwebui_backup_*.tar.gz; do
        if [ -f "$backup" ]; then
            local size
            local date
            size=$(du -h "$backup" | cut -f1)
            date=$(basename "$backup" | sed 's/openwebui_backup_//' | sed 's/\.tar\.gz$//' | tr '_' ' ')
            echo -e "  ${CYAN}•${NC} $date (${GREEN}$size${NC})"
        fi
    done
}

restore_backup() {
    local backup_file=$1
    local backup_base=""

    if [ ! -f "$backup_file" ]; then
        print_error "Backup file not found: $backup_file"
        return 1
    fi

    print_step "Restoring from backup..."
    print_info "Backup file: $backup_file"

    if ! confirm_action "This will replace all current data. Continue?"; then
        print_info "Restore cancelled"
        return 0
    fi

    # Stop services
    print_info "Stopping services..."
    docker_compose stop openwebui 2>/dev/null || true

    # Restore from backup
    print_info "Restoring data..."
    backup_base="$(basename "$backup_file")"
    if [ -z "$backup_base" ]; then
        print_error "Invalid backup filename: $backup_file"
        return 1
    fi

    if docker run --rm \
        -v "${OPENWEBUI_VOLUME}:/data" \
        -v "$(pwd)/$BACKUP_DIR:/backup" \
        alpine sh -c 'rm -rf /data/* && tar xzf "$1" -C /data' _ "/backup/$backup_base"; then
        print_success "Backup restored successfully"

        # Start services
        print_info "Starting services..."
        docker_compose start openwebui
        print_success "Services restarted"
        return 0
    else
        print_error "Failed to restore backup"
        return 1
    fi
}

###############################################################################
# Update Functions
###############################################################################

check_docker() {
    print_step "Checking Docker..."

    if ! command -v docker &>/dev/null; then
        print_error "Docker is not installed"
        exit 1
    fi

    if ! docker info &>/dev/null 2>&1; then
        print_error "Docker is not running"
        exit 1
    fi

    if ! init_compose_cmd; then
        exit 1
    fi

    print_success "Docker is running"
}

check_current_version() {
    print_step "Checking current OpenWebUI version..."

    if docker info &>/dev/null 2>&1; then
        local current_container
        current_container=$(docker_compose ps -q openwebui 2>/dev/null || echo "")

        if [ -n "$current_container" ]; then
            local current_image
            local created
            current_image=$(docker inspect -f '{{.Config.Image}}' "$current_container" 2>/dev/null || echo "unknown")
            print_info "Current image: $current_image"

            # Get image creation date
            created=$(docker inspect -f '{{.Created}}' "$current_container" 2>/dev/null || echo "unknown")
            if [ "$created" != "unknown" ]; then
                created=$(date -d "$created" "+%Y-%m-%d %H:%M:%S" 2>/dev/null || echo "$created")
                print_info "Image created: $created"
            fi
        else
            print_info "OpenWebUI container is not running"
        fi
    fi
}

pull_latest_image() {
    print_step "Pulling latest OpenWebUI image..."
    print_info "Image: $OPENWEBUI_IMAGE"

    # Pull the image
    if docker pull "$OPENWEBUI_IMAGE"; then
        print_success "Image pulled successfully"

        # Get new image info
        local new_image_id
        new_image_id=$(docker inspect -f '{{.Id}}' "$OPENWEBUI_IMAGE" 2>/dev/null || echo "unknown")
        print_info "New image ID: ${new_image_id:0:12}"
        return 0
    else
        print_error "Failed to pull image"

        # Try alternative image
        print_info "Trying alternative image: $OPENWEBUI_ALT_IMAGE"
        if docker pull "$OPENWEBUI_ALT_IMAGE"; then
            print_success "Alternative image pulled successfully"
            OPENWEBUI_IMAGE="$OPENWEBUI_ALT_IMAGE"
            return 0
        fi

        return 1
    fi
}

persist_image_selection() {
    local env_file=".env"

    if [ ! -f "$env_file" ]; then
        print_warning ".env not found; using image override only for this update run"
        return 0
    fi

    if grep -q '^OPENWEBUI_IMAGE=' "$env_file"; then
        sed -i "s|^OPENWEBUI_IMAGE=.*|OPENWEBUI_IMAGE=${OPENWEBUI_IMAGE}|" "$env_file"
    else
        printf "\n# Pinned OpenWebUI image for reproducible deployments\nOPENWEBUI_IMAGE=%s\n" "$OPENWEBUI_IMAGE" >>"$env_file"
    fi

    print_success "Persisted OPENWEBUI_IMAGE in .env"
}

restart_services() {
    print_step "Restarting services..."

    if [ ! -f "$COMPOSE_FILE" ]; then
        print_error "Docker Compose file not found: $COMPOSE_FILE"
        return 1
    fi

    # Recreate container with new image
    print_info "Recreating OpenWebUI container..."
    if OPENWEBUI_IMAGE="$OPENWEBUI_IMAGE" docker_compose up -d --force-recreate openwebui; then
        print_success "Services restarted"
        return 0
    else
        print_error "Failed to restart services"
        return 1
    fi
}

wait_for_openwebui() {
    local openwebui_url="http://localhost:${OPENWEBUI_PORT}"
    print_step "Waiting for OpenWebUI to be ready..."

    local max_attempts=60
    local attempt=1

    while [ $attempt -le $max_attempts ]; do
        if curl -s -f "$openwebui_url" &>/dev/null ||
            curl -s -f "${openwebui_url}/api/health" &>/dev/null; then
            print_success "OpenWebUI is ready!"
            return 0
        fi
        echo -n "."
        sleep 2
        attempt=$((attempt + 1))
    done

    echo
    print_warning "OpenWebUI may still be starting. Check logs with: ./logs.sh"
    return 1
}

ensure_cliproxyapi() {
    if [ "$CLIPROXYAPI_ENABLED" = "false" ]; then
        print_error "CLIPROXYAPI_ENABLED=false is not supported for this update flow"
        return 1
    fi

    print_step "Ensuring CLIProxyAPI is healthy..."

    # Keep update flow stable: chat completion checks can fail due to transient upstream
    # throttling (HTTP 429) while the proxy itself is healthy.
    if [ -x "$CLIPROXYAPI_CHECK_SCRIPT" ] &&
        CLIPROXYAPI_ENABLED=true CLIPROXYAPI_CHECK_CHAT_COMPLETION=false "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1; then
        print_success "CLIProxyAPI is already healthy"
        return 0
    fi

    if [ "$CLIPROXYAPI_DOCKER_MANAGED" = "true" ]; then
        print_info "Restarting CLIProxyAPI container (docker-managed)..."
        if ! docker_compose up -d --force-recreate cliproxyapi; then
            print_error "Failed to restart docker-managed CLIProxyAPI container"
            return 1
        fi
    else
        if [ -x "$CLIPROXYAPI_RESTART_SCRIPT" ]; then
            print_info "Restarting CLIProxyAPI service..."
            CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_RESTART_SCRIPT"
        elif [ -x "$CLIPROXYAPI_START_SCRIPT" ]; then
            print_info "Starting CLIProxyAPI service..."
            CLIPROXYAPI_ENABLED=true "$CLIPROXYAPI_START_SCRIPT"
        else
            print_error "CLIProxyAPI management scripts not found"
            return 1
        fi
    fi

    if [ -x "$CLIPROXYAPI_CHECK_SCRIPT" ] &&
        CLIPROXYAPI_ENABLED=true CLIPROXYAPI_CHECK_CHAT_COMPLETION=false "$CLIPROXYAPI_CHECK_SCRIPT" --quiet >/dev/null 2>&1; then
        print_success "CLIProxyAPI is healthy"
        return 0
    fi

    print_error "CLIProxyAPI health check failed after restart"
    return 1
}

show_update_summary() {
    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}${BOLD}              🚀 Update Successful! 🚀                          ${NC}"
    echo -e "${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"

    echo -e "${BOLD}Updated Services:${NC}"
    echo -e "  ${CYAN}OpenWebUI:${NC}    Updated image: ${OPENWEBUI_IMAGE}"
    echo -e "  ${CYAN}vLLM:${NC}         Unchanged (still running)"
    echo -e "  ${CYAN}CLIProxyAPI:${NC}  Checked/restarted as needed"

    echo -e "\n${BOLD}Quick Commands:${NC}"
    echo -e "  ${YELLOW}View logs:${NC}     ./logs.sh"
    echo -e "  ${YELLOW}Check status:${NC}  ./status.sh"

    echo -e "\n${BOLD}Access OpenWebUI:${NC} http://localhost:${OPENWEBUI_PORT}"

    echo -e "\n${GREEN}${BOLD}═══════════════════════════════════════════════════════════${NC}\n"
}

###############################################################################
# Main Update Flow
###############################################################################

main() {
    print_header

    # Parse command line arguments
    FORCE=false
    SKIP_BACKUP=false
    RESTORE_MODE=false
    BACKUP_FILE=""
    LIST_BACKUPS=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            -f | --force)
                FORCE=true
                shift
                ;;
            --no-backup)
                SKIP_BACKUP=true
                shift
                ;;
            --restore)
                if [ -z "${2:-}" ] || [[ $2 == --* ]]; then
                    print_error "Please specify a backup file to restore"
                    echo "Usage: $0 --restore <backup_file>"
                    echo "Use --list-backups to see available backups"
                    exit 1
                fi
                RESTORE_MODE=true
                BACKUP_FILE="$2"
                shift 2
                ;;
            --list-backups)
                LIST_BACKUPS=true
                shift
                ;;
            -h | --help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Update OpenWebUI to the latest version."
                echo ""
                echo "Options:"
                echo "  -f, --force        Skip confirmation prompts"
                echo "  --no-backup        Skip creating backup"
                echo "  --restore FILE     Restore from backup file"
                echo "  --list-backups     List available backups"
                echo "  -h, --help         Show this help message"
                echo ""
                echo "Examples:"
                echo "  $0                  # Update with backup and confirmation"
                echo "  $0 -f              # Update without confirmation"
                echo "  $0 --no-backup     # Update without creating backup"
                echo "  $0 --list-backups  # List available backups"
                echo "  $0 --restore backups/openwebui_backup_20240101_120000.tar.gz"
                exit 0
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
    done

    # List backups mode
    if [ "$LIST_BACKUPS" = true ]; then
        list_backups
        exit 0
    fi

    # Restore mode
    if [ "$RESTORE_MODE" = true ]; then
        if [ -z "$BACKUP_FILE" ]; then
            print_error "Please specify a backup file to restore"
            echo "Usage: $0 --restore <backup_file>"
            echo "Use --list-backups to see available backups"
            exit 1
        fi
        restore_backup "$BACKUP_FILE"
        exit $?
    fi

    # Update mode
    # Show what will be done
    echo -e "${BOLD}This script will:${NC}"
    echo "  • Check current OpenWebUI version"

    if [ "$SKIP_BACKUP" = false ]; then
        echo "  • Create backup of current data"
    else
        echo "  • ${YELLOW}Skip backup${NC} (--no-backup flag used)"
    fi

    echo "  • Pull latest OpenWebUI image"
    echo "  • Restart OpenWebUI container"
    echo "  • Ensure CLIProxyAPI service is healthy"
    echo "  • vLLM will continue running (no interruption)"
    echo

    # Confirm before proceeding
    if ! confirm_action "Proceed with update?"; then
        print_info "Update cancelled"
        exit 0
    fi

    # Check prerequisites
    check_docker

    # Ensure CLIProxyAPI is healthy before restarting OpenWebUI
    if ! ensure_cliproxyapi; then
        print_error "CLIProxyAPI is required and must be healthy before update"
        exit 1
    fi

    # Show current version
    check_current_version

    # Create backup unless skipped
    if [ "$SKIP_BACKUP" = false ]; then
        if ! create_backup; then
            if ! confirm_action "Backup failed. Continue anyway?"; then
                print_error "Update cancelled due to backup failure"
                exit 1
            fi
        fi
    else
        print_warning "Skipping backup (--no-backup flag used)"
    fi

    # Pull latest image
    if ! pull_latest_image; then
        print_error "Failed to pull latest image"
        exit 1
    fi

    # Persist chosen image tag in .env for reproducible subsequent runs
    persist_image_selection

    # Restart services
    if ! restart_services; then
        print_error "Failed to restart services"
        exit 1
    fi

    # Wait for OpenWebUI to be ready
    if ! wait_for_openwebui; then
        print_error "OpenWebUI failed readiness checks"
        exit 1
    fi

    # Ensure CLIProxyAPI service is healthy after OpenWebUI restart
    if ! ensure_cliproxyapi; then
        print_error "Failed to ensure CLIProxyAPI health"
        exit 1
    fi

    # Show summary
    show_update_summary
}

main "$@"
