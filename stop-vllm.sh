#!/bin/bash
# Script to gracefully stop vLLM server

set -e

# Configuration
PID_FILE="vllm.pid"
TIMEOUT=30

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=== vLLM Server Stop Script ==="

# Check if PID file exists
if [ ! -f "${PID_FILE}" ]; then
    echo -e "${YELLOW}No PID file found. vLLM may not be running.${NC}"

    # Try to find vLLM process by port
    PORT=8000
    VLLM_PID=$(lsof -ti:${PORT} 2>/dev/null || true)

    if [ -n "${VLLM_PID}" ]; then
        echo -e "${YELLOW}Found vLLM process (PID: ${VLLM_PID}) listening on port ${PORT}.${NC}"
        read -p "Do you want to stop it? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            kill ${VLLM_PID}
            echo -e "${GREEN}Sent SIGTERM to vLLM process (PID: ${VLLM_PID}).${NC}"
        fi
    else
        echo "No vLLM process found on port ${PORT}."
    fi
    exit 0
fi

# Read PID from file
VLLM_PID=$(cat "${PID_FILE}")

# Check if process is running
if ! ps -p ${VLLM_PID} > /dev/null 2>&1; then
    echo -e "${YELLOW}vLLM process (PID: ${VLLM_PID}) is not running.${NC}"
    rm -f "${PID_FILE}"
    exit 0
fi

echo "Found vLLM process with PID: ${VLLM_PID}"
echo "Attempting graceful shutdown..."

# Send SIGTERM for graceful shutdown
kill -TERM ${VLLM_PID}

# Wait for process to terminate
for i in $(seq 1 ${TIMEOUT}); do
    if ! ps -p ${VLLM_PID} > /dev/null 2>&1; then
        echo -e "${GREEN}vLLM server stopped gracefully.${NC}"
        rm -f "${PID_FILE}"
        exit 0
    fi
    sleep 1
    echo -n "."
done
echo ""

# Process didn't stop, force kill
echo -e "${YELLOW}vLLM did not stop gracefully within ${TIMEOUT} seconds.${NC}"
read -p "Force kill the process? (y/n) " -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]; then
    if ps -p ${VLLM_PID} > /dev/null 2>&1; then
        kill -9 ${VLLM_PID}
        echo -e "${YELLOW}vLLM server force killed.${NC}"
        rm -f "${PID_FILE}"
    else
        echo -e "${GREEN}vLLM server stopped.${NC}"
        rm -f "${PID_FILE}"
    fi
else
    echo -e "${YELLOW}vLLM server is still running. You can stop it manually with: kill ${VLLM_PID}${NC}"
    exit 1
fi
