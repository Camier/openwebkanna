#!/bin/bash
# Script to restart vLLM server

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

# Configuration
# shellcheck disable=SC2034
PID_FILE="vllm.pid"

print_header "vLLM Server Restart"

# Stop vLLM if running
print_step "Stopping vLLM server..."
if "${SCRIPT_DIR}/stop-vllm.sh"; then
    print_success "vLLM server stopped."
else
    print_warning "vLLM was not running or already stopped."
fi

print_step "Waiting for port to be released..."
sleep 2

# Start vLLM
print_step "Starting vLLM server..."
if "${SCRIPT_DIR}/start-vllm.sh"; then
    print_success "vLLM server restarted successfully!"
else
    print_error "Failed to start vLLM server."
    exit 1
fi
