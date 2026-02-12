#!/bin/bash

###############################################################################
# CLIProxyAPI OAuth Regression Test
# Validates model exposure and chat completion for configured OAuth aliases
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
CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
CLIPROXYAPI_MODELS_PATH="${CLIPROXYAPI_MODELS_PATH:-/v1/models}"
CLIPROXYAPI_CHAT_PATH="${CLIPROXYAPI_CHAT_PATH:-/v1/chat/completions}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"
CLIPROXYAPI_START_SCRIPT="${CLIPROXYAPI_START_SCRIPT:-./start-cliproxyapi.sh}"
CLIPROXYAPI_STOP_SCRIPT="${CLIPROXYAPI_STOP_SCRIPT:-./stop-cliproxyapi.sh}"
CLIPROXYAPI_TEST_ALIASES="${CLIPROXYAPI_TEST_ALIASES:-openai-codex antigravity-oauth qwen-cli kimi-cli}"
CLIPROXYAPI_TEST_PROMPT="${CLIPROXYAPI_TEST_PROMPT:-reply with exactly: ok}"
CLIPROXYAPI_TEST_MAX_TOKENS="${CLIPROXYAPI_TEST_MAX_TOKENS:-12}"
CLIPROXYAPI_REQUEST_TIMEOUT="${CLIPROXYAPI_REQUEST_TIMEOUT:-70}"
CLIPROXYAPI_CHAT_RETRIES="${CLIPROXYAPI_CHAT_RETRIES:-3}"
CLIPROXYAPI_CHAT_RETRY_SLEEP_SECONDS="${CLIPROXYAPI_CHAT_RETRY_SLEEP_SECONDS:-2}"
CLIPROXYAPI_AUTO_START="${CLIPROXYAPI_AUTO_START:-true}"
CLIPROXYAPI_STOP_AFTER_TEST="${CLIPROXYAPI_STOP_AFTER_TEST:-true}"

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

# Test counters
TESTS_TOTAL=0
TESTS_PASSED=0
TESTS_FAILED=0

STARTED_BY_SCRIPT=false
TMP_DIR=""

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

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     CLIProxyAPI OAuth Regression Test                     ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}" >&2
}

print_info() {
    echo -e "${CYAN}ℹ $1${NC}"
}

test_pass() {
    local label="$1"
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    TESTS_PASSED=$((TESTS_PASSED + 1))
    print_success "$label"
}

test_fail() {
    local label="$1"
    local reason="$2"
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    TESTS_FAILED=$((TESTS_FAILED + 1))
    print_error "$label ($reason)"
}

cleanup() {
    if [ -n "$TMP_DIR" ] && [ -d "$TMP_DIR" ]; then
        rm -rf "$TMP_DIR"
    fi

    if [ "$STARTED_BY_SCRIPT" = true ] && is_true "$CLIPROXYAPI_STOP_AFTER_TEST"; then
        if [ -x "$CLIPROXYAPI_STOP_SCRIPT" ]; then
            "$CLIPROXYAPI_STOP_SCRIPT" >/dev/null 2>&1 || true
        fi
    fi
}

get_http_code() {
    local url="$1"
    local output_file="$2"

    if [ -n "$CLIPROXYAPI_API_KEY" ]; then
        curl -sS -m "$CLIPROXYAPI_REQUEST_TIMEOUT" \
            -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
            -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null || true
        return
    fi

    curl -sS -m "$CLIPROXYAPI_REQUEST_TIMEOUT" \
        -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null || true
}

chat_probe() {
    local alias="$1"
    local output_file="$2"
    local payload_file="$3"
    local url="${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_CHAT_PATH}"

    jq -n \
        --arg model "$alias" \
        --arg prompt "$CLIPROXYAPI_TEST_PROMPT" \
        --argjson max_tokens "$CLIPROXYAPI_TEST_MAX_TOKENS" \
        '{
          model: $model,
          messages: [{role: "user", content: $prompt}],
          max_tokens: $max_tokens
        }' > "$payload_file"

    if [ -n "$CLIPROXYAPI_API_KEY" ]; then
        curl -sS -m "$CLIPROXYAPI_REQUEST_TIMEOUT" \
            -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
            -H "Content-Type: application/json" \
            -o "$output_file" -w "%{http_code}" \
            "$url" -d @"$payload_file" 2>/dev/null || true
        return
    fi

    curl -sS -m "$CLIPROXYAPI_REQUEST_TIMEOUT" \
        -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" \
        "$url" -d @"$payload_file" 2>/dev/null || true
}

ensure_ready() {
    local health_url="${CLIPROXYAPI_BASE_URL}/"
    local code
    local health_body="$TMP_DIR/health.json"

    code="$(get_http_code "$health_url" "$health_body")"
    if [ "$code" = "200" ]; then
        print_success "CLIProxyAPI is reachable (${health_url})"
        return 0
    fi

    if ! is_true "$CLIPROXYAPI_AUTO_START"; then
        print_error "CLIProxyAPI is not reachable and CLIPROXYAPI_AUTO_START=false"
        return 1
    fi
    if ! is_true "$CLIPROXYAPI_ENABLED"; then
        print_error "CLIPROXYAPI_ENABLED=false and service is not reachable"
        return 1
    fi
    if [ ! -x "$CLIPROXYAPI_START_SCRIPT" ]; then
        print_error "Start script missing or not executable: $CLIPROXYAPI_START_SCRIPT"
        return 1
    fi

    print_step "Starting CLIProxyAPI for test run"
    "$CLIPROXYAPI_START_SCRIPT"
    STARTED_BY_SCRIPT=true

    code="$(get_http_code "$health_url" "$health_body")"
    if [ "$code" = "200" ]; then
        print_success "CLIProxyAPI became reachable"
        return 0
    fi

    print_error "CLIProxyAPI did not become reachable after start"
    return 1
}

main() {
    CLIPROXYAPI_BASE_URL="$(normalize_base_url "$CLIPROXYAPI_BASE_URL")"
    CLIPROXYAPI_MODELS_PATH="$(normalize_path "$CLIPROXYAPI_MODELS_PATH")"
    CLIPROXYAPI_CHAT_PATH="$(normalize_path "$CLIPROXYAPI_CHAT_PATH")"
    CLIPROXYAPI_START_SCRIPT="$(resolve_path "$CLIPROXYAPI_START_SCRIPT")"
    CLIPROXYAPI_STOP_SCRIPT="$(resolve_path "$CLIPROXYAPI_STOP_SCRIPT")"

    TMP_DIR="$(mktemp -d)"
    trap cleanup EXIT

    print_header
    print_info "Base URL: $CLIPROXYAPI_BASE_URL"
    print_info "Aliases: $CLIPROXYAPI_TEST_ALIASES"
    print_info "Chat retries: $CLIPROXYAPI_CHAT_RETRIES"
    print_info "Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"

    if ! command_exists curl; then
        print_error "curl is required"
        return 1
    fi
    if ! command_exists jq; then
        print_error "jq is required"
        return 1
    fi

    if ! ensure_ready; then
        return 1
    fi

    print_step "Validating /v1/models"
    local models_file="$TMP_DIR/models.json"
    local models_code
    models_code="$(get_http_code "${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_MODELS_PATH}" "$models_file")"
    if [ "$models_code" != "200" ]; then
        test_fail "Models endpoint" "HTTP $models_code"
    else
        test_pass "Models endpoint returned HTTP 200"
    fi

    local alias
    for alias in $CLIPROXYAPI_TEST_ALIASES; do
        if jq -e --arg alias "$alias" '.data[]? | select(.id == $alias)' "$models_file" >/dev/null 2>&1; then
            test_pass "Model listed: $alias"
        else
            test_fail "Model listed: $alias" "missing from /v1/models"
        fi
    done

    print_step "Running chat completion probes"
    for alias in $CLIPROXYAPI_TEST_ALIASES; do
        local out_file="$TMP_DIR/chat_${alias}.json"
        local payload_file="$TMP_DIR/payload_${alias}.json"
        local attempt=1
        local code
        local err=""

        while [ "$attempt" -le "$CLIPROXYAPI_CHAT_RETRIES" ]; do
            code="$(chat_probe "$alias" "$out_file" "$payload_file")"
            if [ "$code" = "200" ]; then
                break
            fi

            err="$(jq -r '.error.message // .message // "unknown error"' "$out_file" 2>/dev/null || true)"
            if [ "$attempt" -lt "$CLIPROXYAPI_CHAT_RETRIES" ]; then
                print_warning "Retrying $alias after HTTP $code (${err}) [attempt ${attempt}/${CLIPROXYAPI_CHAT_RETRIES}]"
                sleep "$CLIPROXYAPI_CHAT_RETRY_SLEEP_SECONDS"
            fi
            attempt=$((attempt + 1))
        done

        if [ "$code" != "200" ]; then
            test_fail "Chat probe: $alias" "HTTP $code (${err})"
            continue
        fi

        test_pass "Chat probe: $alias"
    done

    print_step "Summary"
    print_info "Passed: $TESTS_PASSED"
    print_info "Failed: $TESTS_FAILED"
    print_info "Total:  $TESTS_TOTAL"

    if [ "$TESTS_FAILED" -gt 0 ]; then
        print_error "OAuth regression failed"
        return 1
    fi

    print_success "OAuth regression passed"
    return 0
}

main "$@"
