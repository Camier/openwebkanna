#!/bin/bash

###############################################################################
# CLIProxyAPI Server Restart Script
# Stops then starts CLIProxyAPI
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STOP_SCRIPT="${SCRIPT_DIR}/stop-cliproxyapi.sh"
START_SCRIPT="${SCRIPT_DIR}/start-cliproxyapi.sh"
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_RESTART_WAIT_SECONDS="${CLIPROXYAPI_RESTART_WAIT_SECONDS:-1}"

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     CLIProxyAPI Server Restart                            ║"
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

show_help() {
    cat <<EOF_HELP
Usage: $0 [START_OPTIONS]

Restart CLIProxyAPI by stopping then starting it.
All START_OPTIONS are forwarded to start-cliproxyapi.sh.
EOF_HELP
}

main() {
    if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
        show_help
        exit 0
    fi

    print_header

    if ! is_true "$CLIPROXYAPI_ENABLED"; then
        print_error "CLIProxyAPI lifecycle is disabled (CLIPROXYAPI_ENABLED=$CLIPROXYAPI_ENABLED)"
        exit 1
    fi

    if [ ! -x "$STOP_SCRIPT" ] || [ ! -x "$START_SCRIPT" ]; then
        print_error "Missing required scripts: $STOP_SCRIPT or $START_SCRIPT"
        exit 1
    fi

    print_step "Stopping CLIProxyAPI"
    "$STOP_SCRIPT"

    if [ "$CLIPROXYAPI_RESTART_WAIT_SECONDS" -gt 0 ] 2>/dev/null; then
        sleep "$CLIPROXYAPI_RESTART_WAIT_SECONDS"
    fi

    print_step "Starting CLIProxyAPI"
    "$START_SCRIPT" "$@"

    print_success "CLIProxyAPI restarted"
}

main "$@"
