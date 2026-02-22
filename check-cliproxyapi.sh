#!/bin/bash

###############################################################################
# CLIProxyAPI Server Health Check Script
# Checks process, port, and API health for CLIProxyAPI
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source all library modules
source "${SCRIPT_DIR}/lib/init.sh"

# Load environment
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
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-${OPENAI_API_KEY:-}}"
CLIPROXYAPI_CHECK_TIMEOUT="${CLIPROXYAPI_CHECK_TIMEOUT:-8}"
CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY="${CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY:-true}"
CLIPROXYAPI_CHECK_CHAT_COMPLETION="${CLIPROXYAPI_CHECK_CHAT_COMPLETION:-true}"
CLIPROXYAPI_CHAT_MODEL="${CLIPROXYAPI_CHAT_MODEL:-${CLIPROXYAPI_UPSTREAM_MODEL:-}}"
CLIPROXYAPI_CHAT_PROMPT="${CLIPROXYAPI_CHAT_PROMPT:-ping}"
CLIPROXYAPI_CHAT_MAX_TOKENS="${CLIPROXYAPI_CHAT_MAX_TOKENS:-12}"
CLIPROXYAPI_REQUIRE_AUTH_CHECKS="${CLIPROXYAPI_REQUIRE_AUTH_CHECKS:-true}"
CLIPROXYAPI_EXPECT_ALIASES="${CLIPROXYAPI_EXPECT_ALIASES:-}"
CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH="${CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH:-true}"

QUIET=false

###############################################################################
# Script-specific Functions
###############################################################################

executable_exists() {
    local executable="$1"
    if [[ $executable == */* ]]; then
        [ -x "$executable" ]
        return $?
    fi
    command_exists "$executable"
}

is_false() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "false" ] || [ "$value" = "0" ] || [ "$value" = "no" ]
}

print_check_header() {
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

# Additional HTTP helper specific to this script
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

assert_config_api_key_alignment() {
    local config_api_key=""

    if [ -z "${CLIPROXYAPI_API_KEY:-}" ]; then
        return 0
    fi

    config_api_key="$(extract_first_config_api_key "$CLIPROXYAPI_CONFIG" || true)"
    if [ -z "$config_api_key" ]; then
        print_warning "No api-keys entry found in $CLIPROXYAPI_CONFIG"
        if is_true "$CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH"; then
            print_error "Cannot validate authenticated checks without config api-keys"
            return 1
        fi
        return 0
    fi

    if [ "$config_api_key" = "$CLIPROXYAPI_API_KEY" ]; then
        return 0
    fi

    print_error "CLIPROXYAPI_API_KEY does not match first api-keys entry in $CLIPROXYAPI_CONFIG"
    print_info "Set CLIPROXYAPI_API_KEY to config value or update config.yaml api-keys to match OPENAI_API_KEY"
    if is_true "$CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH"; then
        return 1
    fi
    return 0
}

assert_expected_aliases_local() {
    local payload_file="$1"
    local alias=""
    local ids=""
    local expected_aliases=""

    expected_aliases="$(resolved_expected_aliases | sed '/^$/d' || true)"
    if [ -z "$expected_aliases" ]; then
        return 0
    fi

    ids="$(extract_model_ids "$payload_file" | sed '/^$/d' || true)"
    if [ -z "$ids" ]; then
        print_error "Unable to extract model ids for alias assertions"
        return 1
    fi

    for alias in $expected_aliases; do
        if ! printf "%s\n" "$ids" | grep -Fxq "$alias"; then
            print_error "Expected model alias missing from /v1/models: $alias"
            return 1
        fi
    done

    print_success "Expected aliases present"
    return 0
}

check_process_and_port() {
    local ok=true

    # In Docker-managed mode, host port mappings are commonly implemented via
    # iptables/NAT without a userland process listening on the host. Tools like
    # `lsof` can therefore report "not listening" even though the port is
    # reachable. Rely on the HTTP health check instead to avoid false alarms.
    if is_true "$CLIPROXYAPI_DOCKER_MANAGED"; then
        print_info "Docker-managed runtime: skipping PID/port LISTEN checks (health endpoint is authoritative)"
        return 0
    fi

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
    code="$(http_status "$health_url" "$CLIPROXYAPI_CHECK_TIMEOUT")"
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
    code="$(http_status_with_body "$models_url" "$tmp_file" "$CLIPROXYAPI_CHECK_TIMEOUT" "$header")"
    if [ "$code" != "200" ]; then
        print_error "Models endpoint failed ($models_url), status: $code"
        if [ "$code" = "401" ] || [ "$code" = "403" ]; then
            assert_config_api_key_alignment || true
        fi
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

    if ! assert_expected_aliases_local "$tmp_file"; then
        rm -f "$tmp_file"
        return 1
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

Environment:
  CLIPROXYAPI_REQUIRE_AUTH_CHECKS   Fail when CLIPROXYAPI_API_KEY is unset (default: true)
  CLIPROXYAPI_EXPECT_ALIASES        Optional aliases that must exist in /v1/models (auto-derived when empty)
  CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH  Require config api-key to match CLIPROXYAPI_API_KEY (default: true)

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
            -h | --help)
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

    if [[ $CLIPROXYAPI_CMD == */* ]]; then
        CLIPROXYAPI_CMD="$(resolve_path "$CLIPROXYAPI_CMD")"
    fi
    CLIPROXYAPI_CONFIG="$(resolve_path "$CLIPROXYAPI_CONFIG")"
    CLIPROXYAPI_PID_FILE="$(resolve_path "$CLIPROXYAPI_PID_FILE")"
    CLIPROXYAPI_BASE_URL="$(normalize_base_url "$CLIPROXYAPI_BASE_URL")"
    CLIPROXYAPI_HEALTH_PATH="$(normalize_path "$CLIPROXYAPI_HEALTH_PATH")"
    CLIPROXYAPI_MODELS_PATH="$(normalize_path "$CLIPROXYAPI_MODELS_PATH")"
    CLIPROXYAPI_CHAT_PATH="$(normalize_path "$CLIPROXYAPI_CHAT_PATH")"

    print_check_header

    if ! is_true "$CLIPROXYAPI_ENABLED"; then
        print_error "CLIProxyAPI lifecycle disabled (CLIPROXYAPI_ENABLED=$CLIPROXYAPI_ENABLED)"
        return 1
    fi

    print_step "Checking files and command"
    if [ ! -f "$CLIPROXYAPI_CONFIG" ]; then
        print_error "Config file not found: $CLIPROXYAPI_CONFIG"
        return 1
    fi

    if ! assert_config_api_key_alignment; then
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
        if is_true "$CLIPROXYAPI_REQUIRE_AUTH_CHECKS"; then
            print_error "CLIPROXYAPI_API_KEY is empty and CLIPROXYAPI_REQUIRE_AUTH_CHECKS=true"
            return 1
        fi
        print_warning "CLIPROXYAPI_API_KEY is empty; authenticated models/chat checks are disabled"
    fi

    print_success "All checks passed"
    return 0
}

main "$@"
