#!/bin/bash
# Script to restart vLLM server

set -e

# Configuration
PID_FILE="vllm.pid"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=== vLLM Server Restart Script ==="
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Stop vLLM if running
echo "Stopping vLLM server..."
if "${SCRIPT_DIR}/stop-vllm.sh"; then
    echo -e "${GREEN}vLLM server stopped.${NC}"
else
    echo -e "${YELLOW}vLLM was not running or already stopped.${NC}"
fi

echo ""
echo "Waiting for port to be released..."
sleep 2

# Start vLLM
echo ""
echo "Starting vLLM server..."
if "${SCRIPT_DIR}/start-vllm.sh"; then
    echo ""
    echo -e "${GREEN}vLLM server restarted successfully!${NC}"
else
    echo -e "${RED}Failed to start vLLM server.${NC}"
    exit 1
fi
