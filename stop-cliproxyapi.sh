#!/bin/bash

###############################################################################
# CLIProxyAPI Server Stop Script
# Stops CLIProxyAPI using PID file and optional port fallback
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"

# Configuration
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_PID_FILE="${CLIPROXYAPI_PID_FILE:-.cliproxyapi.pid}"
CLIPROXYAPI_PORT="${CLIPROXYAPI_PORT:-8317}"
CLIPROXYAPI_STOP_TIMEOUT="${CLIPROXYAPI_STOP_TIMEOUT:-20}"
CLIPROXYAPI_FORCE_KILL="${CLIPROXYAPI_FORCE_KILL:-true}"
CLIPROXYAPI_CMD="${CLIPROXYAPI_CMD:-./bin/cliproxyapi}"

if [[ $CLIPROXYAPI_CMD == */* ]]; then
    CLIPROXYAPI_CMD="$(resolve_path "$CLIPROXYAPI_CMD")"
fi
CLIPROXYAPI_PID_FILE="$(resolve_path "$CLIPROXYAPI_PID_FILE")"

pid_matches_cli_proxy() {
    local pid="$1"
    local cmd_fragment
    cmd_fragment="$(basename "$CLIPROXYAPI_CMD")"

    if ! kill -0 "$pid" 2>/dev/null; then
        return 1
    fi

    if command_exists ps; then
        local running_cmd
        running_cmd="$(ps -p "$pid" -o command= 2>/dev/null || true)"
        [[ $running_cmd == *"$cmd_fragment"* ]]
        return $?
    fi

    return 0
}

stop_pid() {
    local pid="$1"
    local wait_count=0

    if ! kill -0 "$pid" 2>/dev/null; then
        return 0
    fi

    kill "$pid" 2>/dev/null || true

    while kill -0 "$pid" 2>/dev/null && [ "$wait_count" -lt "$CLIPROXYAPI_STOP_TIMEOUT" ]; do
        sleep 1
        wait_count=$((wait_count + 1))
    done

    if kill -0 "$pid" 2>/dev/null; then
        if is_true "$CLIPROXYAPI_FORCE_KILL"; then
            print_warning "Process still running after ${CLIPROXYAPI_STOP_TIMEOUT}s; sending SIGKILL"
            kill -9 "$pid" 2>/dev/null || true
            sleep 1
        else
            print_error "Process still running and CLIPROXYAPI_FORCE_KILL=false"
            return 1
        fi
    fi

    if kill -0 "$pid" 2>/dev/null; then
        print_error "Failed to stop process $pid"
        return 1
    fi

    return 0
}

main() {
    print_header "CLIProxyAPI Server Stopper"

    if ! is_true "$CLIPROXYAPI_ENABLED"; then
        print_error "CLIProxyAPI lifecycle is disabled (CLIPROXYAPI_ENABLED=$CLIPROXYAPI_ENABLED)"
        return 1
    fi

    print_step "Stopping CLIProxyAPI"

    local stopped=false

    if [ -f "$CLIPROXYAPI_PID_FILE" ]; then
        local pid
        pid="$(cat "$CLIPROXYAPI_PID_FILE")"

        if [[ $pid =~ ^[0-9]+$ ]]; then
            if pid_matches_cli_proxy "$pid"; then
                print_info "Stopping PID from file: $pid"
                if stop_pid "$pid"; then
                    stopped=true
                fi
            else
                print_warning "PID file exists but process does not match CLIProxyAPI: $pid"
            fi
        else
            print_warning "Invalid PID value in $CLIPROXYAPI_PID_FILE"
        fi

        rm -f "$CLIPROXYAPI_PID_FILE"
    fi

    if [ "$stopped" = false ] && command_exists lsof; then
        local port_pid
        port_pid="$(lsof -nP -iTCP:"$CLIPROXYAPI_PORT" -sTCP:LISTEN -t 2>/dev/null | head -n 1 || true)"
        if [ -n "$port_pid" ] && [[ $port_pid =~ ^[0-9]+$ ]]; then
            if pid_matches_cli_proxy "$port_pid"; then
                print_info "Stopping CLIProxyAPI listener on port $CLIPROXYAPI_PORT (PID: $port_pid)"
                if stop_pid "$port_pid"; then
                    stopped=true
                fi
            else
                print_warning "Port $CLIPROXYAPI_PORT is in use by non-CLIProxyAPI process (PID: $port_pid)"
            fi
        fi
    fi

    if [ "$stopped" = true ]; then
        print_success "CLIProxyAPI stopped"
    else
        print_info "CLIProxyAPI was not running"
    fi
}

main "$@"
