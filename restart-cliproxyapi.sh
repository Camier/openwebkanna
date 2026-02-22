#!/bin/bash

###############################################################################
# CLIProxyAPI Server Restart Script
# Stops then starts CLIProxyAPI
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

STOP_SCRIPT="${SCRIPT_DIR}/stop-cliproxyapi.sh"
START_SCRIPT="${SCRIPT_DIR}/start-cliproxyapi.sh"
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_RESTART_WAIT_SECONDS="${CLIPROXYAPI_RESTART_WAIT_SECONDS:-1}"

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

    print_header "CLIProxyAPI Server Restart"

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
