#!/bin/bash
###############################################################################
# OpenWebUI DB backup helper
#
# Creates a timestamped backup of the OpenWebUI SQLite database from the
# running `openwebui` container into ./backups (gitignored).
#
# Usage:
#   ./backup-openwebui-db.sh
#   ./backup-openwebui-db.sh --output backups/webui.db.backup.sqlite
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

usage() {
    cat <<'EOF'
Usage:
  ./backup-openwebui-db.sh
  ./backup-openwebui-db.sh --output backups/webui.db.backup.sqlite
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

if ! docker ps --format '{{.Names}}' | grep -qx 'openwebui'; then
    print_error "Container 'openwebui' is not running."
    print_error "Start it first: docker compose up -d openwebui"
    exit 1
fi

mkdir -p backups

if [ -z "$OUTPUT_PATH" ]; then
    ts="$(date +%Y%m%d-%H%M%S)"
    OUTPUT_PATH="backups/webui.db.${ts}.sqlite"
fi

print_step "Copying webui.db from container to: ${OUTPUT_PATH}"
docker cp "openwebui:/app/backend/data/webui.db" "$OUTPUT_PATH"

print_success "Backup created: ${OUTPUT_PATH}"
