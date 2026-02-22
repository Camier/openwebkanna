#!/bin/bash
# Script to check if vLLM is running and responding

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

# Configuration
HOST="localhost"
PORT=8000
PID_FILE="vllm.pid"
HEALTH_URL="http://${HOST}:${PORT}/health"
MODELS_URL="http://${HOST}:${PORT}/v1/models"

echo "=== vLLM Server Status Check ==="
echo ""

VLLM_RUNNING=false
VLLM_RESPONDING=false

# Check if PID file exists
if [ -f "${PID_FILE}" ]; then
    VLLM_PID=$(cat "${PID_FILE}")
    echo -e "${BLUE}PID file found: ${VLLM_PID}${NC}"

    # Check if process is running
    if ps -p ${VLLM_PID} >/dev/null 2>&1; then
        echo -e "${GREEN}OK vLLM process is running (PID: ${VLLM_PID})${NC}"
        VLLM_RUNNING=true

        # Show process details
        echo -e "${BLUE}Process details:${NC}"
        ps -p ${VLLM_PID} -o pid,ppid,cmd,etime,pcpu,pmem --no-headers | while read line; do
            echo "  ${line}"
        done
    else
        echo -e "${RED}X vLLM process is NOT running (stale PID file)${NC}"
    fi
else
    echo -e "${YELLOW}No PID file found.${NC}"
fi

echo ""

# Check if port is in use
echo -e "${BLUE}Port check:${NC}"
if lsof -Pi :${PORT} -sTCP:LISTEN -t >/dev/null 2>&1; then
    PORT_PID=$(lsof -ti:${PORT})
    echo -e "${GREEN}OK Port ${PORT} is in use (process: ${PORT_PID})${NC}"

    # Check if it matches our PID file
    if [ -f "${PID_FILE}" ]; then
        VLLM_PID=$(cat "${PID_FILE}")
        if [ "${PORT_PID}" == "${VLLM_PID}" ]; then
            echo -e "${GREEN}OK Port process matches PID file${NC}"
        else
            echo -e "${YELLOW}! Port process (${PORT_PID}) does not match PID file (${VLLM_PID})${NC}"
        fi
    fi
else
    echo -e "${RED}X Port ${PORT} is NOT in use${NC}"
fi

echo ""

# Check if server is responding
echo -e "${BLUE}Server health check:${NC}"

# Check if curl is available
if ! command_exists curl; then
    echo -e "${YELLOW}curl not found. Skipping HTTP health check.${NC}"
else
    # Try health endpoint
    HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" "${HEALTH_URL}" 2>/dev/null || echo "000")

    if [ "${HTTP_CODE}" == "200" ]; then
        echo -e "${GREEN}OK Health endpoint responding (HTTP ${HTTP_CODE})${NC}"
        VLLM_RESPONDING=true

        # Try models endpoint to show loaded model
        MODELS_RESPONSE=$(curl -s "${MODELS_URL}" 2>/dev/null || echo "")
        if [ -n "${MODELS_RESPONSE}" ]; then
            echo -e "${BLUE}Available models:${NC}"
            echo "${MODELS_RESPONSE}" | python3 -m json.tool 2>/dev/null || echo "${MODELS_RESPONSE}"
        fi
    else
        echo -e "${RED}X Health endpoint not responding (HTTP ${HTTP_CODE})${NC}"
        echo "  URL: ${HEALTH_URL}"
    fi
fi

echo ""
echo "=== Summary ==="

if ${VLLM_RUNNING} && ${VLLM_RESPONDING}; then
    echo -e "${GREEN}vLLM server is running and responding correctly!${NC}"
    echo ""
    echo "You can access the API at:"
    echo "  - Base URL: http://${HOST}:${PORT}"
    echo "  - Health: ${HEALTH_URL}"
    echo "  - Models: ${MODELS_URL}"
    exit 0
elif ${VLLM_RUNNING}; then
    echo -e "${YELLOW}vLLM process is running but not responding to HTTP requests.${NC}"
    echo "The server might still be starting up. Check the logs:"
    echo "  tail -f vllm.log"
    exit 1
else
    echo -e "${RED}vLLM server is NOT running.${NC}"
    echo ""
    echo "To start the server, run:"
    echo "  ./start-vllm.sh"
    exit 1
fi
