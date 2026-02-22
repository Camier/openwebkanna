#!/bin/bash
# Script to start vLLM server with proper configuration

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"

# vLLM environment path
VLLM_ENV="${HOME}/.conda/envs/vllm"
VLLM_PYTHON="${VLLM_ENV}/bin/python"

# Configuration (override via environment variables when needed)
MODEL_NAME="${VLLM_MODEL_NAME:-meta-llama/Llama-3.1-8B-Instruct}"
HOST="${VLLM_HOST:-127.0.0.1}"
PORT="${VLLM_PORT:-8000}"
GPU_MEMORY_UTILIZATION="${VLLM_GPU_MEMORY_UTILIZATION:-0.25}"
MAX_MODEL_LEN="${VLLM_MAX_MODEL_LEN:-4096}"
VLLM_EXTRA_ARGS="${VLLM_EXTRA_ARGS:-}"
LOG_FILE="vllm.log"
PID_FILE="vllm.pid"
LEGACY_PID_FILE=".vllm_pid"

echo "=== vLLM Server Startup Script ==="
echo "Model: ${MODEL_NAME}"
echo "Host: ${HOST}"
echo "Port: ${PORT}"
echo "GPU Memory Utilization: ${GPU_MEMORY_UTILIZATION}"
echo "Max Model Len: ${MAX_MODEL_LEN}"
if [ -n "${VLLM_EXTRA_ARGS}" ]; then
    echo "Extra Args: ${VLLM_EXTRA_ARGS}"
fi
if [ "${HOST}" = "0.0.0.0" ] || [ "${HOST}" = "::" ]; then
    echo -e "${YELLOW}Warning: vLLM is configured to listen on all interfaces (${HOST}).${NC}"
fi
echo ""

# Check if vLLM is already running
if [ -f "${PID_FILE}" ] || [ -f "${LEGACY_PID_FILE}" ]; then
    if [ -f "${PID_FILE}" ]; then
        PID=$(cat "${PID_FILE}")
    else
        PID=$(cat "${LEGACY_PID_FILE}")
        # Migrate legacy PID file after process validation succeeds.
    fi

    if ps -p "${PID}" >/dev/null 2>&1; then
        COMMAND=$(ps -p "${PID}" -o args= 2>/dev/null || true)
        if [[ $COMMAND == *"vllm.entrypoints.openai.api_server"* ]]; then
            echo -e "${YELLOW}vLLM is already running with PID ${PID}.${NC}"
            echo "Use 'stop-vllm.sh' to stop it first, or 'check-vllm.sh' to check its status."
            exit 1
        fi
    fi

    echo -e "${YELLOW}Removing stale PID file.${NC}"
    rm -f "${PID_FILE}" "${LEGACY_PID_FILE}"
fi

# Check if port is already in use
if lsof -Pi ":${PORT}" -sTCP:LISTEN -t >/dev/null 2>&1; then
    echo -e "${RED}Error: Port ${PORT} is already in use.${NC}"
    echo "Please check what's using the port and stop it first."
    exit 1
fi

# Check if vLLM environment exists
if [ ! -d "${VLLM_ENV}" ]; then
    echo -e "${RED}Error: vLLM environment not found at ${VLLM_ENV}${NC}"
    exit 1
fi

if [ ! -f "${VLLM_PYTHON}" ]; then
    echo -e "${RED}Error: Python not found in vLLM environment${NC}"
    exit 1
fi

if ! "${VLLM_PYTHON}" -c "import vllm" 2>/dev/null; then
    echo -e "${RED}Error: vLLM module not found in ${VLLM_ENV}${NC}"
    exit 1
fi

echo -e "${GREEN}Starting vLLM server...${NC}"
echo "Using vLLM environment: ${VLLM_ENV}"
echo "Using vLLM version: $("${VLLM_PYTHON}" -c 'import vllm; print(vllm.__version__)')"

declare -a EXTRA_ARGS_ARRAY
if [ -n "${VLLM_EXTRA_ARGS}" ]; then
    # Split optional extra args on whitespace, e.g. '--dtype float16 --enforce-eager'
    read -r -a EXTRA_ARGS_ARRAY <<<"${VLLM_EXTRA_ARGS}"
else
    EXTRA_ARGS_ARRAY=()
fi

# Start vLLM server in background
nohup "${VLLM_PYTHON}" -m vllm.entrypoints.openai.api_server \
    --model "${MODEL_NAME}" \
    --host "${HOST}" \
    --port "${PORT}" \
    --gpu-memory-utilization "${GPU_MEMORY_UTILIZATION}" \
    --max-model-len "${MAX_MODEL_LEN}" \
    "${EXTRA_ARGS_ARRAY[@]}" \
    >"${LOG_FILE}" 2>&1 &

VLLM_PID=$!
echo "${VLLM_PID}" >"${PID_FILE}"
echo "${VLLM_PID}" >"${LEGACY_PID_FILE}"

echo "vLLM started with PID: ${VLLM_PID}"
echo "Logs are being written to: ${LOG_FILE}"
echo ""

# Wait a moment and check if the process is still running
sleep 3
if ps -p "${VLLM_PID}" >/dev/null 2>&1; then
    echo -e "${GREEN}vLLM server is running!${NC}"
    echo ""
    echo "To check the logs:"
    echo "  tail -f ${LOG_FILE}"
    echo ""
    echo "To check if the server is responding:"
    echo "  ./check-vllm.sh"
    echo ""
    echo "To stop the server:"
    echo "  ./stop-vllm.sh"
else
    echo -e "${RED}vLLM server failed to start. Check the log file: ${LOG_FILE}${NC}"
    rm -f "${PID_FILE}"
    exit 1
fi
