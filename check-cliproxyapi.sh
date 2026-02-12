#!/bin/bash

###############################################################################
# CLIProxyAPI Server Health Check Script
# Checks process, port, and API health for CLIProxyAPI
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

load_env_defaults() {
    local env_file=""
    local line=""
    local key=""
    local value=""

    if [ -f "$SCRIPT_DIR/.env" ]; then
        env_file="$SCRIPT_DIR/.env"
    elif [ -f ".env" ]; then
        env_file=".env"
    fi

    if [ -z "$env_file" ]; then
        return 0
    fi

    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%$'\r'}"
        line="${line#"${line%%[![:space:]]*}"}"
        [ -z "$line" ] && continue
        [[ "$line" = \#* ]] && continue
        [[ "$line" != *=* ]] && continue

        key="${line%%=*}"
        value="${line#*=}"
        key="$(printf "%s" "$key" | tr -d '[:space:]')"

        [[ "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        if [ -n "${!key+x}" ]; then
            continue
        fi

        if [[ "$value" == \"*\" ]] && [[ "$value" == *\" ]]; then
            value="${value:1:${#value}-2}"
        elif [[ "$value" == \'*\' ]] && [[ "$value" == *\' ]]; then
            value="${value:1:${#value}-2}"
        fi

        printf -v "$key" "%s" "$value"
        export "$key"
    done < "$env_file"
}

load_env_defaults

# Configuration
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_DOCKER_MANAGED="${CLIPROXYAPI_DOCKER_MANAGED:-true}"
CLIPROXYAPI_CMD="${CLIPROXYAPI_CMD:-./bin/cliproxyapi}"
CLIPROXYAPI_CONFIG="${CLIPROXYAPI_CONFIG:-./cliproxyapi/config.yaml}"
CLIPROXYAPI_PID_FILE="${CLIPROXYAPI_PID_FILE:-.cliproxyapi.pid}"
CLIPROXYAPI_PORT="${CLIPROXYAPI_PORT:-8317}"
CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
CLIPROXYAPI_HEALTH_PATH="${CLIPROXYAPI_HEALTH_PATH:-/}"
CLIPROXYAPI_MODELS_PATH="${CLIPROXYAPI_MODELS_PATH:-/v1/models}"
CLIPROXYAPI_CHAT_PATH="${CLIPROXYAPI_CHAT_PATH:-/v1/chat/completions}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"
CLIPROXYAPI_CHECK_TIMEOUT="${CLIPROXYAPI_CHECK_TIMEOUT:-8}"
CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY="${CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY:-true}"
CLIPROXYAPI_CHECK_CHAT_COMPLETION="${CLIPROXYAPI_CHECK_CHAT_COMPLETION:-true}"
CLIPROXYAPI_CHAT_MODEL="${CLIPROXYAPI_CHAT_MODEL:-${CLIPROXYAPI_UPSTREAM_MODEL:-}}"
CLIPROXYAPI_CHAT_PROMPT="${CLIPROXYAPI_CHAT_PROMPT:-ping}"
CLIPROXYAPI_CHAT_MAX_TOKENS="${CLIPROXYAPI_CHAT_MAX_TOKENS:-12}"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

QUIET=false

resolve_path() {
    local candidate="$1"
    if [[ "$candidate" = /* ]]; then
        printf "%s" "$candidate"
    else
        printf "%s/%s" "$SCRIPT_DIR" "$candidate"
    fi
}

normalize_base_url() {
    local value="$1"
    printf "%s" "${value%/}"
}

normalize_path() {
    local value="$1"
    if [ -z "$value" ]; then
        printf "/"
        return
    fi
    if [[ "$value" != /* ]]; then
        value="/$value"
    fi
    printf "%s" "$value"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

executable_exists() {
    local executable="$1"
    if [[ "$executable" == */* ]]; then
        [ -x "$executable" ]
        return $?
    fi
    command_exists "$executable"
}

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

print_header() {
    $QUIET && return 0
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     CLIProxyAPI Server Health Check                       ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    $QUIET && return 0
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    $QUIET && return 0
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    $QUIET && return 0
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    $QUIET && return 0
    echo -e "${CYAN}ℹ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}" >&2
}

is_false() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "false" ] || [ "$value" = "0" ] || [ "$value" = "no" ]
}

http_status() {
    local url="$1"
    local header="${2:-}"

    if [ -n "$header" ]; then
        curl -sS -m "$CLIPROXYAPI_CHECK_TIMEOUT" -o /dev/null -w "%{http_code}" -H "$header" "$url" || true
    else
        curl -sS -m "$CLIPROXYAPI_CHECK_TIMEOUT" -o /dev/null -w "%{http_code}" "$url" || true
    fi
}

http_status_with_body() {
    local url="$1"
    local out_file="$2"
    local header="${3:-}"

    if [ -n "$header" ]; then
        curl -sS -m "$CLIPROXYAPI_CHECK_TIMEOUT" -o "$out_file" -w "%{http_code}" -H "$header" "$url" || true
    else
        curl -sS -m "$CLIPROXYAPI_CHECK_TIMEOUT" -o "$out_file" -w "%{http_code}" "$url" || true
    fi
}

extract_models_count() {
    local payload_file="$1"

    if command_exists jq; then
        jq -er '.data | if type == "array" then length else empty end' "$payload_file" 2>/dev/null || return 1
        return 0
    fi

    if ! grep -q '"data"[[:space:]]*:' "$payload_file"; then
        return 1
    fi

    grep -o '"id"[[:space:]]*:[[:space:]]*"[^"]*"' "$payload_file" | wc -l | tr -d '[:space:]'
}

extract_first_model() {
    local payload_file="$1"

    if command_exists jq; then
        jq -r '.data[]?.id // empty' "$payload_file" 2>/dev/null | head -n 1
        return 0
    fi

    grep -o '"id"[[:space:]]*:[[:space:]]*"[^"]*"' "$payload_file" | head -n 1 | cut -d'"' -f4
}

check_process_and_port() {
    local ok=true

    if [ -f "$CLIPROXYAPI_PID_FILE" ]; then
        local pid
        pid="$(cat "$CLIPROXYAPI_PID_FILE")"
        if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
            print_success "PID file process is running (PID: $pid)"
        else
            print_warning "PID file exists but process is not running"
            ok=false
        fi
    else
        print_warning "PID file not found: $CLIPROXYAPI_PID_FILE"
    fi

    local listener_pid
    listener_pid="$(lsof -nP -iTCP:"$CLIPROXYAPI_PORT" -sTCP:LISTEN -t 2>/dev/null | head -n 1 || true)"
    if [ -n "$listener_pid" ]; then
        print_success "Port $CLIPROXYAPI_PORT is listening (PID: $listener_pid)"
    else
        print_warning "Port $CLIPROXYAPI_PORT is not listening"
        ok=false
    fi

    if [ "$ok" = true ]; then
        return 0
    fi
    return 1
}

check_http_health() {
    local health_url="${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_HEALTH_PATH}"
    local code

    print_step "Checking health endpoint"
    code="$(http_status "$health_url")"
    if [ "$code" = "200" ]; then
        print_success "Health endpoint OK ($health_url)"
        return 0
    fi

    print_error "Health endpoint failed ($health_url), status: $code"
    return 1
}

check_models() {
    local header="Authorization: Bearer $CLIPROXYAPI_API_KEY"
    local models_url="${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_MODELS_PATH}"
    local tmp_file
    tmp_file="$(mktemp)"

    print_step "Checking models endpoint"

    local code
    code="$(http_status_with_body "$models_url" "$tmp_file" "$header")"
    if [ "$code" != "200" ]; then
        print_error "Models endpoint failed ($models_url), status: $code"
        rm -f "$tmp_file"
        return 1
    fi

    local models_count
    models_count="$(extract_models_count "$tmp_file" 2>/dev/null || true)"
    if [ -z "$models_count" ]; then
        print_error "Unable to parse models payload"
        rm -f "$tmp_file"
        return 1
    fi

    if is_true "$CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY" && [ "$models_count" -le 0 ] 2>/dev/null; then
        print_error "Models list is empty"
        rm -f "$tmp_file"
        return 1
    fi

    local first_model
    first_model="$(extract_first_model "$tmp_file" || true)"
    if [ -n "$first_model" ]; then
        print_success "Models endpoint OK ($models_count models, first: $first_model)"
    else
        print_success "Models endpoint OK ($models_count models)"
    fi

    if [ -z "$CLIPROXYAPI_CHAT_MODEL" ] && [ -n "$first_model" ]; then
        CLIPROXYAPI_CHAT_MODEL="$first_model"
    fi

    rm -f "$tmp_file"
    return 0
}

check_chat_completion() {
    if ! is_true "$CLIPROXYAPI_CHECK_CHAT_COMPLETION"; then
        print_info "Chat completion check skipped"
        return 0
    fi

    if [ -z "$CLIPROXYAPI_CHAT_MODEL" ]; then
        print_error "Chat completion check enabled but no model is available"
        return 1
    fi

    local chat_url="${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_CHAT_PATH}"
    local header="Authorization: Bearer $CLIPROXYAPI_API_KEY"
    local tmp_file
    tmp_file="$(mktemp)"

    local payload
    if command_exists jq; then
        payload="$(jq -cn \
            --arg model "$CLIPROXYAPI_CHAT_MODEL" \
            --arg prompt "$CLIPROXYAPI_CHAT_PROMPT" \
            --argjson max_tokens "$CLIPROXYAPI_CHAT_MAX_TOKENS" \
            '{model:$model,max_tokens:$max_tokens,messages:[{role:"user",content:$prompt}]}')"
    else
        payload="{\"model\":\"$CLIPROXYAPI_CHAT_MODEL\",\"max_tokens\":$CLIPROXYAPI_CHAT_MAX_TOKENS,\"messages\":[{\"role\":\"user\",\"content\":\"$CLIPROXYAPI_CHAT_PROMPT\"}]}"
    fi

    print_step "Checking chat completion"

    local code
    code="$(curl -sS -m "$CLIPROXYAPI_CHECK_TIMEOUT" \
        -o "$tmp_file" \
        -w "%{http_code}" \
        -H "$header" \
        -H "Content-Type: application/json" \
        -d "$payload" \
        "$chat_url" || true)"

    if [ "$code" != "200" ]; then
        print_error "Chat completion failed ($chat_url), status: $code"
        rm -f "$tmp_file"
        return 1
    fi

    if command_exists jq; then
        if jq -e '.error != null' "$tmp_file" >/dev/null 2>&1; then
            print_error "Chat completion returned error payload"
            rm -f "$tmp_file"
            return 1
        fi
    fi

    print_success "Chat completion OK (model: $CLIPROXYAPI_CHAT_MODEL)"
    rm -f "$tmp_file"
    return 0
}

show_usage() {
    cat <<'EOF_USAGE'
Usage: ./check-cliproxyapi.sh [--quiet] [--help]

Options:
  --quiet   Minimal output (exit code only)
  --help    Show help

Exit codes:
  0 on success
  1 on any failed check
EOF_USAGE
}

parse_args() {
    while [ $# -gt 0 ]; do
        case "$1" in
            --quiet)
                QUIET=true
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                show_usage
                exit 1
                ;;
        esac
        shift
    done
}

main() {
    parse_args "$@"

    if [[ "$CLIPROXYAPI_CMD" == */* ]]; then
        CLIPROXYAPI_CMD="$(resolve_path "$CLIPROXYAPI_CMD")"
    fi
    CLIPROXYAPI_CONFIG="$(resolve_path "$CLIPROXYAPI_CONFIG")"
    CLIPROXYAPI_PID_FILE="$(resolve_path "$CLIPROXYAPI_PID_FILE")"
    CLIPROXYAPI_BASE_URL="$(normalize_base_url "$CLIPROXYAPI_BASE_URL")"
    CLIPROXYAPI_HEALTH_PATH="$(normalize_path "$CLIPROXYAPI_HEALTH_PATH")"
    CLIPROXYAPI_MODELS_PATH="$(normalize_path "$CLIPROXYAPI_MODELS_PATH")"
    CLIPROXYAPI_CHAT_PATH="$(normalize_path "$CLIPROXYAPI_CHAT_PATH")"

    print_header

    if ! is_true "$CLIPROXYAPI_ENABLED"; then
        print_error "CLIProxyAPI lifecycle disabled (CLIPROXYAPI_ENABLED=$CLIPROXYAPI_ENABLED)"
        return 1
    fi

    print_step "Checking files and command"
    if [ ! -f "$CLIPROXYAPI_CONFIG" ]; then
        print_error "Config file not found: $CLIPROXYAPI_CONFIG"
        return 1
    fi

    if is_false "$CLIPROXYAPI_DOCKER_MANAGED"; then
        if ! executable_exists "$CLIPROXYAPI_CMD"; then
            print_error "CLIProxyAPI command not found or not executable: $CLIPROXYAPI_CMD"
            return 1
        fi
    fi

    if is_false "$CLIPROXYAPI_DOCKER_MANAGED"; then
        print_success "Config and command are present"
    else
        print_success "Config is present (docker-managed runtime)"
    fi

    check_process_and_port || true

    check_http_health

    if [ -n "$CLIPROXYAPI_API_KEY" ]; then
        check_models
        check_chat_completion
    else
        print_warning "CLIPROXYAPI_API_KEY is empty; skipping authenticated models/chat checks"
    fi

    print_success "All checks passed"
    return 0
}

main "$@"
