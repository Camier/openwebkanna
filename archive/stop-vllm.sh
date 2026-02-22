#!/bin/bash
# Script to gracefully stop vLLM server

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"

# Configuration
PORT="${VLLM_PORT:-8000}"
VLLM_PID_FILE="vllm.pid"
VLLM_PID_FILE_LEGACY=".vllm_pid"
TIMEOUT="${VLLM_STOP_TIMEOUT:-30}"
FORCE_KILL_ON_TIMEOUT="${STOP_VLLM_FORCE_KILL:-true}"

echo "=== vLLM Server Stop Script ==="

has_vllm_command() {
    local pid=$1

    if ! kill -0 "$pid" 2>/dev/null; then
        return 1
    fi

    ps -p "$pid" -o args= 2>/dev/null | grep -q "vllm.entrypoints.openai.api_server"
}

resolve_vllm_pid() {
    local pid_file=$1
    if [ -f "$pid_file" ]; then
        cat "$pid_file"
    fi
}

collect_vllm_pids_from_port() {
    local pids=""
    local candidate=""
    local fallback_pids=""

    if ! command -v lsof >/dev/null 2>&1; then
        fallback_pids="$(pgrep -f "vllm.entrypoints.openai.api_server" 2>/dev/null || true)"
        echo "$fallback_pids"
        return
    fi

    for candidate in $(lsof -tiTCP:"${PORT}" -sTCP:LISTEN 2>/dev/null || true); do
        if has_vllm_command "$candidate"; then
            pids+=" $candidate"
        fi
    done

    echo "$pids"
}

stop_vllm_pid() {
    local pid=$1
    local attempt=1

    if ! has_vllm_command "$pid"; then
        echo -e "${YELLOW}PID ${pid} does not match a running vLLM process.${NC}"
        return 1
    fi

    echo -e "${YELLOW}Attempting graceful shutdown for PID: ${pid}.${NC}"
    kill -TERM "$pid" || true

    while [ "$attempt" -le "${TIMEOUT}" ]; do
        if ! kill -0 "$pid" 2>/dev/null; then
            return 0
        fi
        sleep 1
        echo -n "."
        attempt=$((attempt + 1))
    done
    echo ""

    if [ "${FORCE_KILL_ON_TIMEOUT}" != "true" ]; then
        echo -e "${YELLOW}vLLM process did not exit within ${TIMEOUT} seconds.${NC}"
        return 1
    fi

    echo -e "${YELLOW}vLLM did not stop gracefully. Sending SIGKILL.${NC}"
    if kill -0 "$pid" 2>/dev/null; then
        kill -9 "$pid" || true
        return 0
    fi

    return 1
}

vllm_pid_file=""
if [ -f "${VLLM_PID_FILE}" ]; then
    vllm_pid_file="${VLLM_PID_FILE}"
elif [ -f "${VLLM_PID_FILE_LEGACY}" ]; then
    vllm_pid_file="${VLLM_PID_FILE_LEGACY}"
fi

target_pids=""

if [ -n "${vllm_pid_file}" ]; then
    vllm_pid="$(resolve_vllm_pid "$vllm_pid_file")"
    if has_vllm_command "$vllm_pid"; then
        target_pids="${vllm_pid}"
    else
        echo -e "${YELLOW}Removing stale PID file.${NC}"
        rm -f "${VLLM_PID_FILE}" "${VLLM_PID_FILE_LEGACY}"
    fi
fi

if [ -z "${target_pids}" ]; then
    from_port=$(collect_vllm_pids_from_port)
    if [ -n "${from_port}" ]; then
        echo -e "${YELLOW}Found vLLM process(es) on port ${PORT}: ${from_port}.${NC}"
        target_pids="${from_port}"
    else
        echo -e "${YELLOW}No vLLM process found to stop.${NC}"
        exit 0
    fi
fi

stop_failed=0
for pid in ${target_pids}; do
    if stop_vllm_pid "$pid"; then
        echo -e "${GREEN}vLLM process stopped (PID: ${pid}).${NC}"
    else
        echo -e "${RED}Failed to stop vLLM process (PID: ${pid}).${NC}"
        stop_failed=1
    fi
done

rm -f "${VLLM_PID_FILE}" "${VLLM_PID_FILE_LEGACY}"

if [ "${stop_failed}" -eq 1 ]; then
    exit 1
fi

echo -e "${GREEN}vLLM server stop flow complete.${NC}"
exit 0
