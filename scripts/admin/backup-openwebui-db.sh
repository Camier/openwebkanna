#!/bin/bash
###############################################################################
# OpenWebUI DB backup helper
#
# Creates a timestamped backup of the OpenWebUI SQLite database from the
# running `openwebui` container into ./backups (gitignored).
#
# Usage:
#   ./scripts/admin/backup-openwebui-db.sh
#   ./scripts/admin/backup-openwebui-db.sh --output backups/webui.db.backup.sqlite
###############################################################################

set -e

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_SERVICE="${OPENWEBUI_SERVICE:-openwebui}"
OPENWEBUI_CONTAINER_NAME="${OPENWEBUI_CONTAINER_NAME:-$OPENWEBUI_SERVICE}"

usage() {
    cat <<'EOF'
Usage:
  ./scripts/admin/backup-openwebui-db.sh
  ./scripts/admin/backup-openwebui-db.sh --output backups/webui.db.backup.sqlite
EOF
}

OUTPUT_PATH=""
while [ $# -gt 0 ]; do
    case "$1" in
        --output)
            OUTPUT_PATH="${2:-}"
            shift 2
            ;;
        -h | --help)
            usage
            exit 0
            ;;
        *)
            print_error "Unknown argument: $1"
            usage
            exit 2
            ;;
    esac
done

print_header "OpenWebUI DB Backup"

if ! docker info >/dev/null 2>&1; then
    print_error "Docker daemon is not running."
    exit 1
fi

if ! init_compose_cmd; then
    print_error "Docker Compose command unavailable"
    exit 1
fi

if ! compose_service_running "$OPENWEBUI_SERVICE"; then
    print_error "OpenWebUI service '$OPENWEBUI_SERVICE' is not running."
    print_error "Start it first: docker compose up -d ${OPENWEBUI_SERVICE}"
    exit 1
fi

CONTAINER_ID="$(resolve_container_id "$OPENWEBUI_SERVICE" "$OPENWEBUI_CONTAINER_NAME" || true)"
if [ -z "$CONTAINER_ID" ]; then
    print_error "Unable to resolve running container for service '$OPENWEBUI_SERVICE'."
    exit 1
fi

mkdir -p backups

if [ -z "$OUTPUT_PATH" ]; then
    ts="$(date +%Y%m%d-%H%M%S)"
    OUTPUT_PATH="backups/webui.db.${ts}.sqlite"
fi

print_step "Copying webui.db from container to: ${OUTPUT_PATH}"
docker cp "${CONTAINER_ID}:/app/backend/data/webui.db" "$OUTPUT_PATH"

print_success "Backup created: ${OUTPUT_PATH}"
