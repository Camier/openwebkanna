#!/bin/bash

###############################################################################
# OpenWebUI -> CLIProxyAPI Routing Regression Test
# Validates that OpenWebUI exposes and can chat via required CLIProxyAPI aliases
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
WEBUI_PORT="${WEBUI_PORT:-3000}"
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT}}"
OPENWEBUI_HEALTH_PATH="${OPENWEBUI_HEALTH_PATH:-/health}"
OPENWEBUI_MODELS_PATH="${OPENWEBUI_MODELS_PATH:-/api/models}"
OPENWEBUI_CHAT_PATH="${OPENWEBUI_CHAT_PATH:-/api/chat/completions}"
OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"

CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
CLIPROXYAPI_HEALTH_PATH="${CLIPROXYAPI_HEALTH_PATH:-/}"
CLIPROXYAPI_MODELS_PATH="${CLIPROXYAPI_MODELS_PATH:-/v1/models}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"

CLIPROXYAPI_START_SCRIPT="${CLIPROXYAPI_START_SCRIPT:-./start-cliproxyapi.sh}"
CLIPROXYAPI_STOP_SCRIPT="${CLIPROXYAPI_STOP_SCRIPT:-./stop-cliproxyapi.sh}"
CLIPROXYAPI_AUTO_START="${CLIPROXYAPI_AUTO_START:-true}"
CLIPROXYAPI_STOP_AFTER_TEST="${CLIPROXYAPI_STOP_AFTER_TEST:-false}"
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_DOCKER_MANAGED="${CLIPROXYAPI_DOCKER_MANAGED:-true}"
CLIPROXYAPI_DOCKER_SERVICE="${CLIPROXYAPI_DOCKER_SERVICE:-cliproxyapi}"

ROUTING_TEST_ALIASES="${ROUTING_TEST_ALIASES:-openai-codex antigravity-oauth qwen-cli kimi-cli}"
ROUTING_TEST_PROMPT="${ROUTING_TEST_PROMPT:-reply with exactly: routed-ok}"
ROUTING_TEST_MAX_TOKENS="${ROUTING_TEST_MAX_TOKENS:-16}"
ROUTING_REQUEST_TIMEOUT="${ROUTING_REQUEST_TIMEOUT:-90}"
ROUTING_CHAT_RETRIES="${ROUTING_CHAT_RETRIES:-3}"
ROUTING_CHAT_RETRY_SLEEP_SECONDS="${ROUTING_CHAT_RETRY_SLEEP_SECONDS:-2}"
ROUTING_DOCKER_REACHABILITY_CHECK="${ROUTING_DOCKER_REACHABILITY_CHECK:-true}"

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

TESTS_TOTAL=0
TESTS_PASSED=0
TESTS_FAILED=0

TMP_DIR=""
STARTED_CLIPROXYAPI=false
OPENWEBUI_AUTH_TOKEN=""

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

resolve_path() {
    local candidate="$1"
    if [[ "$candidate" = /* ]]; then
        printf "%s" "$candidate"
    else
        printf "%s/%s" "$SCRIPT_DIR" "$candidate"
    fi
}

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  OpenWebUI -> CLIProxyAPI Routing Regression              ║"
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

    if [ "$STARTED_CLIPROXYAPI" = true ] && is_true "$CLIPROXYAPI_STOP_AFTER_TEST"; then
        if [ -x "$CLIPROXYAPI_STOP_SCRIPT" ]; then
            "$CLIPROXYAPI_STOP_SCRIPT" >/dev/null 2>&1 || true
        fi
    fi
}

request_get() {
    local url="$1"
    local output_file="$2"
    local auth_key="$3"

    if [ -n "$auth_key" ]; then
        curl -sS -m "$ROUTING_REQUEST_TIMEOUT" \
            -H "Authorization: Bearer $auth_key" \
            -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null || true
        return
    fi

    curl -sS -m "$ROUTING_REQUEST_TIMEOUT" \
        -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null || true
}

openwebui_chat_probe() {
    local model_id="$1"
    local output_file="$2"
    local payload_file="$3"
    local url="${OPENWEBUI_URL}${OPENWEBUI_CHAT_PATH}"

    jq -n \
        --arg model "$model_id" \
        --arg prompt "$ROUTING_TEST_PROMPT" \
        --argjson max_tokens "$ROUTING_TEST_MAX_TOKENS" \
        '{
          model: $model,
          messages: [{role: "user", content: $prompt}],
          max_tokens: $max_tokens,
          stream: false
        }' > "$payload_file"

    if [ -n "$OPENWEBUI_AUTH_TOKEN" ]; then
        curl -sS -m "$ROUTING_REQUEST_TIMEOUT" \
            -H "Authorization: Bearer $OPENWEBUI_AUTH_TOKEN" \
            -H "Content-Type: application/json" \
            -o "$output_file" -w "%{http_code}" \
            "$url" -d @"$payload_file" 2>/dev/null || true
        return
    fi

    curl -sS -m "$ROUTING_REQUEST_TIMEOUT" \
        -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" \
        "$url" -d @"$payload_file" 2>/dev/null || true
}

openwebui_signin() {
    local output_file="$1"
    local signin_url="${OPENWEBUI_URL}${OPENWEBUI_SIGNIN_PATH}"
    local payload_file="$TMP_DIR/openwebui_signin_payload.json"

    jq -n \
        --arg email "$OPENWEBUI_SIGNIN_EMAIL" \
        --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
        '{email: $email, password: $password}' > "$payload_file"

    curl -sS -m "$ROUTING_REQUEST_TIMEOUT" \
        -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" \
        "$signin_url" -d @"$payload_file" 2>/dev/null || true
}

ensure_openwebui_auth_token() {
    if [ -n "$OPENWEBUI_API_KEY" ]; then
        OPENWEBUI_AUTH_TOKEN="$OPENWEBUI_API_KEY"
        print_info "Using OPENWEBUI_API_KEY for OpenWebUI API calls"
        return 0
    fi

    if ! is_true "$OPENWEBUI_AUTO_AUTH"; then
        print_warning "OPENWEBUI_AUTO_AUTH=false and OPENWEBUI_API_KEY is empty"
        return 0
    fi

    print_step "Acquiring OpenWebUI bearer token via signin"
    local signin_file="$TMP_DIR/openwebui_signin.json"
    local signin_code
    local signin_error=""

    signin_code="$(openwebui_signin "$signin_file")"
    if [ "$signin_code" != "200" ]; then
        signin_error="$(jq -r '.detail // .error.message // .message // "signin failed"' "$signin_file" 2>/dev/null || true)"
        print_warning "OpenWebUI signin failed (HTTP $signin_code: $signin_error)"
        return 0
    fi

    OPENWEBUI_AUTH_TOKEN="$(jq -r '.token // empty' "$signin_file" 2>/dev/null || true)"
    if [ -n "$OPENWEBUI_AUTH_TOKEN" ]; then
        print_success "OpenWebUI bearer token acquired"
        return 0
    fi

    print_warning "OpenWebUI signin response did not include token"
    return 0
}

ensure_cliproxyapi_ready() {
    local health_file="$TMP_DIR/cliproxy_health.json"
    local health_url="${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_HEALTH_PATH}"
    local code

    code="$(request_get "$health_url" "$health_file" "$CLIPROXYAPI_API_KEY")"
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
    if is_true "$CLIPROXYAPI_DOCKER_MANAGED"; then
        if ! command_exists docker; then
            print_error "docker is required for CLIPROXYAPI_DOCKER_MANAGED=true"
            return 1
        fi
        if ! command_exists docker-compose; then
            print_error "docker-compose is required for CLIPROXYAPI_DOCKER_MANAGED=true"
            return 1
        fi

        print_step "Starting CLIProxyAPI Docker service for routing test"
        if ! docker-compose up -d "$CLIPROXYAPI_DOCKER_SERVICE"; then
            print_error "Failed to start docker-compose service: $CLIPROXYAPI_DOCKER_SERVICE"
            return 1
        fi
    else
        if [ ! -x "$CLIPROXYAPI_START_SCRIPT" ]; then
            print_error "Start script missing or not executable: $CLIPROXYAPI_START_SCRIPT"
            return 1
        fi
        print_step "Starting CLIProxyAPI for routing test"
        "$CLIPROXYAPI_START_SCRIPT"
        STARTED_CLIPROXYAPI=true
    fi

    code="$(request_get "$health_url" "$health_file" "$CLIPROXYAPI_API_KEY")"
    if [ "$code" = "200" ]; then
        print_success "CLIProxyAPI became reachable"
        return 0
    fi

    print_error "CLIProxyAPI did not become reachable after start"
    return 1
}

extract_model_ids() {
    local catalog_file="$1"
    jq -r '
        if type == "array" then
            .
        elif (type == "object" and (.data | type == "array")) then
            .data
        else
            []
        end
        | .[]?
        | .id // empty
    ' "$catalog_file" 2>/dev/null || true
}

resolve_openwebui_model_id() {
    local catalog_file="$1"
    local alias="$2"
    local id=""

    while IFS= read -r id; do
        [ -z "$id" ] && continue
        if [ "$id" = "$alias" ]; then
            printf "%s" "$id"
            return 0
        fi
    done < <(extract_model_ids "$catalog_file")

    while IFS= read -r id; do
        [ -z "$id" ] && continue
        case "$id" in
            */"$alias"|*."$alias"|*:"$alias"|*"$alias"*)
                printf "%s" "$id"
                return 0
                ;;
        esac
    done < <(extract_model_ids "$catalog_file")

    return 1
}

build_models_url_from_base() {
    local base_url="$1"
    if [[ "$base_url" == */v1 ]]; then
        printf "%s/models" "$base_url"
        return
    fi
    printf "%s/v1/models" "${base_url%/}"
}

check_openwebui_container_upstream_reachability() {
    local openai_base="${OPENAI_API_BASE_URL:-}"
    local models_url=""
    local code=""

    if ! is_true "$ROUTING_DOCKER_REACHABILITY_CHECK"; then
        return 0
    fi
    if [ -z "$openai_base" ]; then
        return 0
    fi
    if ! command_exists docker; then
        return 0
    fi
    if ! docker ps --format '{{.Names}}' | grep -Fxq "openwebui"; then
        return 0
    fi

    models_url="$(build_models_url_from_base "$openai_base")"
    print_step "Checking OpenWebUI container reachability to upstream"
    print_info "Container URL: $models_url"

    code="$(docker exec openwebui sh -lc "curl -sS -m 8 -o /dev/null -w '%{http_code}' -H 'Authorization: Bearer ${CLIPROXYAPI_API_KEY}' '${models_url}'" 2>/dev/null || true)"
    if [ "$code" = "200" ]; then
        test_pass "OpenWebUI container can reach upstream models URL"
        return 0
    fi

    test_fail "OpenWebUI container upstream reachability" "HTTP ${code:-000} for ${models_url}"
    print_warning "OpenWebUI cannot reach configured OPENAI_API_BASE_URL from inside container"
    return 0
}

main() {
    OPENWEBUI_URL="$(normalize_base_url "$OPENWEBUI_URL")"
    OPENWEBUI_HEALTH_PATH="$(normalize_path "$OPENWEBUI_HEALTH_PATH")"
    OPENWEBUI_MODELS_PATH="$(normalize_path "$OPENWEBUI_MODELS_PATH")"
    OPENWEBUI_CHAT_PATH="$(normalize_path "$OPENWEBUI_CHAT_PATH")"
    OPENWEBUI_SIGNIN_PATH="$(normalize_path "$OPENWEBUI_SIGNIN_PATH")"

    CLIPROXYAPI_BASE_URL="$(normalize_base_url "$CLIPROXYAPI_BASE_URL")"
    CLIPROXYAPI_HEALTH_PATH="$(normalize_path "$CLIPROXYAPI_HEALTH_PATH")"
    CLIPROXYAPI_MODELS_PATH="$(normalize_path "$CLIPROXYAPI_MODELS_PATH")"
    CLIPROXYAPI_START_SCRIPT="$(resolve_path "$CLIPROXYAPI_START_SCRIPT")"
    CLIPROXYAPI_STOP_SCRIPT="$(resolve_path "$CLIPROXYAPI_STOP_SCRIPT")"

    TMP_DIR="$(mktemp -d)"
    trap cleanup EXIT

    print_header
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "CLIProxyAPI URL: $CLIPROXYAPI_BASE_URL"
    print_info "Aliases: $ROUTING_TEST_ALIASES"
    print_info "Chat retries: $ROUTING_CHAT_RETRIES"
    print_info "OpenWebUI auto auth: $OPENWEBUI_AUTO_AUTH"
    print_info "Timestamp (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"

    if ! command_exists curl; then
        print_error "curl is required"
        return 1
    fi
    if ! command_exists jq; then
        print_error "jq is required"
        return 1
    fi

    if ! ensure_cliproxyapi_ready; then
        return 1
    fi

    print_step "Checking OpenWebUI health"
    local openwebui_health_file="$TMP_DIR/openwebui_health.json"
    local openwebui_health_code
    openwebui_health_code="$(request_get "${OPENWEBUI_URL}${OPENWEBUI_HEALTH_PATH}" "$openwebui_health_file" "$OPENWEBUI_API_KEY")"
    if [ "$openwebui_health_code" = "200" ]; then
        test_pass "OpenWebUI health endpoint returned HTTP 200"
    else
        test_fail "OpenWebUI health endpoint" "HTTP $openwebui_health_code (is OpenWebUI running?)"
        print_error "Cannot continue routing checks without OpenWebUI"
        return 1
    fi

    ensure_openwebui_auth_token
    check_openwebui_container_upstream_reachability

    print_step "Fetching model lists from CLIProxyAPI and OpenWebUI"
    local cliproxy_models_file="$TMP_DIR/cliproxy_models.json"
    local openwebui_models_file="$TMP_DIR/openwebui_models.json"
    local cliproxy_models_code
    local openwebui_models_code

    cliproxy_models_code="$(request_get "${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_MODELS_PATH}" "$cliproxy_models_file" "$CLIPROXYAPI_API_KEY")"
    if [ "$cliproxy_models_code" = "200" ]; then
        if [ "$(extract_model_ids "$cliproxy_models_file" | sed '/^$/d' | wc -l | tr -d ' ')" -gt 0 ]; then
            test_pass "CLIProxyAPI models endpoint returned non-empty .data"
        else
            test_fail "CLIProxyAPI models endpoint shape" "expected non-empty model id list"
        fi
    else
        test_fail "CLIProxyAPI models endpoint" "HTTP $cliproxy_models_code"
    fi

    openwebui_models_code="$(request_get "${OPENWEBUI_URL}${OPENWEBUI_MODELS_PATH}" "$openwebui_models_file" "$OPENWEBUI_AUTH_TOKEN")"
    if [ "$openwebui_models_code" = "200" ]; then
        if [ "$(extract_model_ids "$openwebui_models_file" | sed '/^$/d' | wc -l | tr -d ' ')" -gt 0 ]; then
            test_pass "OpenWebUI models endpoint returned non-empty .data"
        else
            test_fail "OpenWebUI models endpoint shape" "expected non-empty model id list"
        fi
    elif [ "$openwebui_models_code" = "401" ]; then
        test_fail "OpenWebUI models endpoint" "HTTP 401 (set OPENWEBUI_API_KEY or provide valid OPENWEBUI_SIGNIN_EMAIL/OPENWEBUI_SIGNIN_PASSWORD)"
    else
        test_fail "OpenWebUI models endpoint" "HTTP $openwebui_models_code"
    fi

    print_step "Checking required aliases in both model catalogs"
    local alias
    local resolved_model_id=""
    for alias in $ROUTING_TEST_ALIASES; do
        if extract_model_ids "$cliproxy_models_file" | grep -Fxq "$alias"; then
            test_pass "CLIProxyAPI model listed: $alias"
        else
            test_fail "CLIProxyAPI model listed: $alias" "missing from /v1/models"
        fi

        resolved_model_id="$(resolve_openwebui_model_id "$openwebui_models_file" "$alias" || true)"
        if [ -n "$resolved_model_id" ]; then
            test_pass "OpenWebUI model listed: $alias"
            if [ "$resolved_model_id" != "$alias" ]; then
                print_info "Resolved OpenWebUI id: $alias -> $resolved_model_id"
            fi
        else
            test_fail "OpenWebUI model listed: $alias" "missing from /api/models"
        fi
    done

    print_step "Running OpenWebUI chat probes per alias"
    for alias in $ROUTING_TEST_ALIASES; do
        resolved_model_id="$(resolve_openwebui_model_id "$openwebui_models_file" "$alias" || true)"
        if [ -z "$resolved_model_id" ]; then
            test_fail "OpenWebUI chat probe: $alias" "no matching OpenWebUI model id"
            continue
        fi

        local out_file="$TMP_DIR/openwebui_chat_${alias}.json"
        local payload_file="$TMP_DIR/openwebui_payload_${alias}.json"
        local attempt=1
        local code
        local err=""

        while [ "$attempt" -le "$ROUTING_CHAT_RETRIES" ]; do
            code="$(openwebui_chat_probe "$resolved_model_id" "$out_file" "$payload_file")"
            if [ "$code" = "200" ]; then
                if jq -e '.choices and (.choices | type == "array") and (.choices | length > 0)' "$out_file" >/dev/null 2>&1; then
                    break
                fi
                err="missing choices array"
            else
                err="$(jq -r '.error.message // .detail // .message // "unknown error"' "$out_file" 2>/dev/null || true)"
            fi

            if [ "$attempt" -lt "$ROUTING_CHAT_RETRIES" ]; then
                print_warning "Retrying OpenWebUI chat for $alias after HTTP $code (${err}) [attempt ${attempt}/${ROUTING_CHAT_RETRIES}]"
                sleep "$ROUTING_CHAT_RETRY_SLEEP_SECONDS"
            fi
            attempt=$((attempt + 1))
        done

        if [ "$code" != "200" ]; then
            test_fail "OpenWebUI chat probe: $alias" "HTTP $code (${err})"
            continue
        fi

        if jq -e '.choices and (.choices | type == "array") and (.choices | length > 0)' "$out_file" >/dev/null 2>&1; then
            test_pass "OpenWebUI chat probe: $alias"
        else
            test_fail "OpenWebUI chat probe: $alias" "response missing choices"
        fi
    done

    print_step "Summary"
    print_info "Passed: $TESTS_PASSED"
    print_info "Failed: $TESTS_FAILED"
    print_info "Total:  $TESTS_TOTAL"

    if [ "$TESTS_FAILED" -gt 0 ]; then
        print_error "OpenWebUI -> CLIProxyAPI routing regression failed"
        return 1
    fi

    print_success "OpenWebUI -> CLIProxyAPI routing regression passed"
    return 0
}

main "$@"
