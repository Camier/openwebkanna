#!/bin/bash

###############################################################################
# OpenWebUI Tool-by-Tool Functional Audit
# Audits OpenWebUI tool endpoints (tools, functions, code execution, web search)
###############################################################################

set -euo pipefail

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
        export "${key?}"
    done < "$env_file"
}

load_env_defaults

# Color definitions (standard repository set)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
CURL_TIMEOUT="${CURL_TIMEOUT:-45}"
VERBOSE="${VERBOSE:-false}"
OUTPUT_DIR="${OUTPUT_DIR:-/tmp/openwebui_tool_audit_$(date +%Y%m%d_%H%M%S)}"
AUDIT_BOOTSTRAP_FUNCTION="${AUDIT_BOOTSTRAP_FUNCTION:-true}"
AUDIT_BOOTSTRAP_FUNCTION_ID="${AUDIT_BOOTSTRAP_FUNCTION_ID:-audit_ping_pipe}"
AUDIT_BOOTSTRAP_FUNCTION_NAME="${AUDIT_BOOTSTRAP_FUNCTION_NAME:-Audit Ping Pipe}"
AUDIT_BOOTSTRAP_FUNCTION_DESCRIPTION="${AUDIT_BOOTSTRAP_FUNCTION_DESCRIPTION:-Bootstrap function created by endpoint audit}"
AUDIT_BOOTSTRAP_CLEANUP="${AUDIT_BOOTSTRAP_CLEANUP:-false}"

# Runtime state
API_TOKEN=""
CURRENT_TEST=""
RESPONSE_CODE=""
RESPONSE_FILE=""
RESPONSE_CURL_ERR_FILE=""
BOOTSTRAP_FUNCTION_CREATED="false"

# Counters
TESTS_TOTAL=0
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_WARNINGS=0
TESTS_SKIPPED=0

declare -a FAILURES=()
declare -a WARNINGS=()

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║       OpenWebUI Tool-by-Tool Functional Audit             ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    echo -e "${BLUE}${BOLD}▶ $1${NC}"
}

print_section() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}" >&2
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

require_dependencies() {
    local missing=0
    for cmd in curl jq; do
        if ! command_exists "$cmd"; then
            print_error "Missing dependency: $cmd"
            missing=1
        fi
    done

    if [ "$missing" -ne 0 ]; then
        exit 1
    fi
}

start_test() {
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    CURRENT_TEST="$1"
    print_step "$CURRENT_TEST"
}

pass_test() {
    TESTS_PASSED=$((TESTS_PASSED + 1))
    print_success "$CURRENT_TEST"
}

fail_test() {
    local reason="$1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILURES+=("$CURRENT_TEST :: $reason")
    print_error "$CURRENT_TEST :: $reason"
}

warn_test() {
    local reason="$1"
    TESTS_WARNINGS=$((TESTS_WARNINGS + 1))
    WARNINGS+=("$CURRENT_TEST :: $reason")
    print_warning "$CURRENT_TEST :: $reason"
}

skip_test() {
    local reason="$1"
    TESTS_SKIPPED=$((TESTS_SKIPPED + 1))
    print_info "Skipped: $CURRENT_TEST :: $reason"
}

uri_encode() {
    local value="$1"
    jq -nr --arg v "$value" '$v|@uri'
}

response_detail() {
    local detail
    detail="$(jq -r '.detail // empty' "$RESPONSE_FILE" 2>/dev/null || true)"
    if [ -n "$detail" ]; then
        printf "%s" "$detail"
    else
        head -c 240 "$RESPONSE_FILE" 2>/dev/null || true
    fi
}

request_api() {
    local method="$1"
    local path="$2"
    local payload_file="${3:-}"
    local label="$4"
    local url="${OPENWEBUI_URL%/}${path}"

    RESPONSE_FILE="$OUTPUT_DIR/${label}.json"
    RESPONSE_CURL_ERR_FILE="$OUTPUT_DIR/${label}.curl.err"

    local -a curl_args=(
        -sS
        -m "$CURL_TIMEOUT"
        -X "$method"
        "$url"
        -H "Accept: application/json"
        -o "$RESPONSE_FILE"
        -w "%{http_code}"
    )

    if [ -n "$API_TOKEN" ]; then
        curl_args+=(-H "Authorization: Bearer $API_TOKEN")
    fi

    if [ -n "$payload_file" ]; then
        curl_args+=(-H "Content-Type: application/json" --data @"$payload_file")
    fi

    RESPONSE_CODE="$(curl "${curl_args[@]}" 2>"$RESPONSE_CURL_ERR_FILE" || true)"
}

sign_in() {
    local signin_output payload_file signin_url signin_code role

    signin_output="$OUTPUT_DIR/auth_signin.json"
    payload_file="$OUTPUT_DIR/auth_signin.payload.json"
    signin_url="${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"

    jq -n \
        --arg email "$OPENWEBUI_SIGNIN_EMAIL" \
        --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
        '{email: $email, password: $password}' > "$payload_file"

    signin_code="$(curl -sS -m "$CURL_TIMEOUT" -X POST "$signin_url" \
        -H "Content-Type: application/json" \
        --data @"$payload_file" \
        -o "$signin_output" -w "%{http_code}" 2>"$OUTPUT_DIR/auth_signin.curl.err" || true)"

    if [ "$signin_code" != "200" ]; then
        print_error "OpenWebUI signin failed (HTTP $signin_code)"
        if [ -s "$signin_output" ]; then
            cat "$signin_output"
        fi
        exit 1
    fi

    API_TOKEN="$(jq -r '.token // empty' "$signin_output")"
    role="$(jq -r '.role // empty' "$signin_output")"

    if [ -z "$API_TOKEN" ]; then
        print_error "Signin succeeded but no token returned"
        exit 1
    fi

    if [ -n "$role" ]; then
        print_info "Authenticated with role: $role"
    else
        print_info "Authenticated successfully"
    fi
}

bootstrap_function_if_registry_empty() {
    local function_count payload_file function_content detail

    if ! is_true "$AUDIT_BOOTSTRAP_FUNCTION"; then
        return 0
    fi

    function_count="$(jq 'length' "$OUTPUT_DIR/functions_list.json" 2>/dev/null || printf '0')"
    if [ "$function_count" -gt 0 ]; then
        return 0
    fi

    start_test "POST /api/v1/functions/create bootstraps minimal function"

    function_content="$(cat <<'PY'
from typing import Optional


class Pipe:
    def pipe(self, body: dict, __user__: Optional[dict] = None) -> str:
        return "audit-ping-ok"
PY
)"

    payload_file="$OUTPUT_DIR/functions_bootstrap.payload.json"
    jq -n \
        --arg id "$AUDIT_BOOTSTRAP_FUNCTION_ID" \
        --arg name "$AUDIT_BOOTSTRAP_FUNCTION_NAME" \
        --arg content "$function_content" \
        --arg description "$AUDIT_BOOTSTRAP_FUNCTION_DESCRIPTION" \
        '{
            id: $id,
            name: $name,
            content: $content,
            meta: {description: $description}
        }' > "$payload_file"

    request_api "POST" "/api/v1/functions/create" "$payload_file" "functions_bootstrap_create"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e --arg id "$AUDIT_BOOTSTRAP_FUNCTION_ID" '.id == $id' "$RESPONSE_FILE" >/dev/null 2>&1; then
        BOOTSTRAP_FUNCTION_CREATED="true"
        pass_test
        return 0
    fi

    detail="$(response_detail)"
    if [ "$RESPONSE_CODE" = "400" ] && [[ "$detail" == *"ID_TAKEN"* || "$detail" == *"ID already taken"* ]]; then
        warn_test "Bootstrap function id already exists"
        return 0
    fi

    fail_test "HTTP $RESPONSE_CODE :: $detail"
}

cleanup_bootstrap_function_if_needed() {
    local function_id_encoded detail

    if ! is_true "$AUDIT_BOOTSTRAP_CLEANUP"; then
        return 0
    fi

    if [ "$BOOTSTRAP_FUNCTION_CREATED" != "true" ]; then
        print_info "Bootstrap cleanup skipped (no function created in this run)"
        return 0
    fi

    function_id_encoded="$(uri_encode "$AUDIT_BOOTSTRAP_FUNCTION_ID")"

    start_test "DELETE /api/v1/functions/id/{id}/delete cleans bootstrap function"
    request_api "DELETE" "/api/v1/functions/id/${function_id_encoded}/delete" "" "functions_bootstrap_delete"

    if [ "$RESPONSE_CODE" = "200" ] && jq -e '. == true' "$RESPONSE_FILE" >/dev/null 2>&1; then
        BOOTSTRAP_FUNCTION_CREATED="false"
        pass_test
        return 0
    fi

    detail="$(response_detail)"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e '. == false' "$RESPONSE_FILE" >/dev/null 2>&1; then
        warn_test "Bootstrap function delete returned false"
        return 0
    fi

    fail_test "HTTP $RESPONSE_CODE :: $detail"
}

test_tools_endpoints() {
    local tools_count first_tool_id first_tool_id_encoded

    print_section "Tools Endpoints"

    start_test "GET /api/v1/tools/list returns tool access list"
    request_api "GET" "/api/v1/tools/list" "" "tools_list"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "array"' "$RESPONSE_FILE" >/dev/null 2>&1; then
        tools_count="$(jq 'length' "$RESPONSE_FILE")"
        print_info "tools/list count: $tools_count"
        pass_test
    else
        fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        return
    fi

    start_test "GET /api/v1/tools/ returns user tool view"
    request_api "GET" "/api/v1/tools/" "" "tools_root"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "array"' "$RESPONSE_FILE" >/dev/null 2>&1; then
        pass_test
    else
        fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
    fi

    first_tool_id="$(jq -r '.[0].id // empty' "$OUTPUT_DIR/tools_list.json")"
    if [ -n "$first_tool_id" ]; then
        first_tool_id_encoded="$(uri_encode "$first_tool_id")"
        start_test "GET /api/v1/tools/id/{id} returns tool details"
        request_api "GET" "/api/v1/tools/id/${first_tool_id_encoded}" "" "tools_first_id"
        if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "object"' "$RESPONSE_FILE" >/dev/null 2>&1; then
            pass_test
        else
            fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        fi
    else
        start_test "GET /api/v1/tools/id/{id} returns tool details"
        skip_test "No tools available in /api/v1/tools/list"
    fi
}

test_functions_endpoints() {
    local functions_count first_function_id first_function_id_encoded

    print_section "Functions Endpoints"

    start_test "GET /api/v1/functions/list returns function registry (admin)"
    request_api "GET" "/api/v1/functions/list" "" "functions_list"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "array"' "$RESPONSE_FILE" >/dev/null 2>&1; then
        functions_count="$(jq 'length' "$RESPONSE_FILE")"
        print_info "functions/list count: $functions_count"
        pass_test
    else
        fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        return
    fi

    bootstrap_function_if_registry_empty

    request_api "GET" "/api/v1/functions/list" "" "functions_list_after_bootstrap"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "array"' "$RESPONSE_FILE" >/dev/null 2>&1; then
        functions_count="$(jq 'length' "$RESPONSE_FILE")"
        cp "$RESPONSE_FILE" "$OUTPUT_DIR/functions_list.json"
        if [ "$BOOTSTRAP_FUNCTION_CREATED" = "true" ]; then
            print_info "functions/list count after bootstrap: $functions_count"
        fi
    else
        fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        return
    fi

    start_test "GET /api/v1/functions/ returns active functions"
    request_api "GET" "/api/v1/functions/" "" "functions_root"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "array"' "$RESPONSE_FILE" >/dev/null 2>&1; then
        pass_test
    else
        fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
    fi

    first_function_id="$(jq -r '.[0].id // empty' "$OUTPUT_DIR/functions_list.json")"
    if [ -n "$first_function_id" ]; then
        first_function_id_encoded="$(uri_encode "$first_function_id")"
        start_test "GET /api/v1/functions/id/{id} returns function details"
        request_api "GET" "/api/v1/functions/id/${first_function_id_encoded}" "" "functions_first_id"
        if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "object"' "$RESPONSE_FILE" >/dev/null 2>&1; then
            pass_test
        else
            fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        fi
    else
        start_test "GET /api/v1/functions/id/{id} returns function details"
        skip_test "No functions available in /api/v1/functions/list"
    fi
}

test_code_execution_endpoints() {
    local code_exec_enabled code_exec_engine payload_file detail

    print_section "Code Execution Endpoints"

    start_test "GET /api/v1/configs/code_execution returns code runtime config"
    request_api "GET" "/api/v1/configs/code_execution" "" "code_execution_config"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "object"' "$RESPONSE_FILE" >/dev/null 2>&1; then
        code_exec_enabled="$(jq -r '.ENABLE_CODE_EXECUTION // false' "$RESPONSE_FILE")"
        code_exec_engine="$(jq -r '.CODE_EXECUTION_ENGINE // empty' "$RESPONSE_FILE")"
        print_info "ENABLE_CODE_EXECUTION=$code_exec_enabled CODE_EXECUTION_ENGINE=$code_exec_engine"
        pass_test
    else
        fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        return
    fi

    payload_file="$OUTPUT_DIR/code_execute.payload.json"
    jq -n --arg code "print(2 + 3)" '{code: $code}' > "$payload_file"

    start_test "POST /api/v1/utils/code/execute executes Python snippet"
    request_api "POST" "/api/v1/utils/code/execute" "$payload_file" "code_execute"
    if [ "$RESPONSE_CODE" = "200" ]; then
        pass_test
        return
    fi

    detail="$(response_detail)"
    if [ "$RESPONSE_CODE" = "400" ] && [[ "$detail" == *"Code execution engine not supported"* ]]; then
        if [ "$code_exec_enabled" = "true" ] && [ "$code_exec_engine" != "jupyter" ]; then
            fail_test "Config mismatch: ENABLE_CODE_EXECUTION=true but engine=$code_exec_engine is unsupported by /utils/code/execute (expects jupyter)"
        else
            warn_test "Engine does not support /utils/code/execute in current config"
        fi
        return
    fi

    fail_test "HTTP $RESPONSE_CODE :: $detail"
}

test_web_search_endpoints() {
    local web_search_enabled web_search_engine searxng_url ssl_verification payload_file detail loaded_count

    print_section "Web Search Endpoints"

    start_test "GET /api/v1/retrieval/config returns web search config"
    request_api "GET" "/api/v1/retrieval/config" "" "retrieval_config"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e '.web and (.web | type == "object")' "$RESPONSE_FILE" >/dev/null 2>&1; then
        web_search_enabled="$(jq -r '.web.ENABLE_WEB_SEARCH // false' "$RESPONSE_FILE")"
        web_search_engine="$(jq -r '.web.WEB_SEARCH_ENGINE // empty' "$RESPONSE_FILE")"
        searxng_url="$(jq -r '.web.SEARXNG_QUERY_URL // empty' "$RESPONSE_FILE")"
        ssl_verification="$(jq -r '.web.ENABLE_WEB_LOADER_SSL_VERIFICATION // empty' "$RESPONSE_FILE")"
        print_info "ENABLE_WEB_SEARCH=$web_search_enabled WEB_SEARCH_ENGINE=$web_search_engine SSL_VERIFY=$ssl_verification"
        if [ "$web_search_engine" = "searxng" ] && [ -n "$searxng_url" ]; then
            print_info "SEARXNG_QUERY_URL=$searxng_url"
        fi
        pass_test
    else
        fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        return
    fi

    payload_file="$OUTPUT_DIR/web_search.payload.json"
    jq -n '{queries: ["OpenWebUI tool endpoint audit"]}' > "$payload_file"

    start_test "POST /api/v1/retrieval/process/web/search returns loaded documents"
    request_api "POST" "/api/v1/retrieval/process/web/search" "$payload_file" "web_search_execute"

    if [ "$web_search_enabled" != "true" ]; then
        if [ "$RESPONSE_CODE" = "403" ]; then
            skip_test "Web search disabled by config"
        else
            warn_test "Web search disabled but endpoint returned HTTP $RESPONSE_CODE"
        fi
        return
    fi

    if [ "$RESPONSE_CODE" = "200" ] && jq -e '.status == true' "$RESPONSE_FILE" >/dev/null 2>&1; then
        loaded_count="$(jq -r '.loaded_count // 0' "$RESPONSE_FILE")"
        if [ "$loaded_count" -ge 1 ]; then
            print_info "loaded_count=$loaded_count"
            pass_test
            return
        fi
        fail_test "HTTP 200 but loaded_count=$loaded_count"
        return
    fi

    detail="$(response_detail)"
    if [[ "$detail" == *"CERTIFICATE_VERIFY_FAILED"* ]] || [[ "$detail" == *"SSL"* ]] || [[ "$detail" == *"certificate"* ]]; then
        fail_test "SSL/CA failure in web search pipeline :: $detail"
    elif [[ "$detail" == *"No results found from web search"* ]]; then
        fail_test "No results returned by configured search engine :: $detail"
    else
        fail_test "HTTP $RESPONSE_CODE :: $detail"
    fi
}

write_summary() {
    local summary_file
    summary_file="$OUTPUT_DIR/summary.txt"

    {
        echo "OpenWebUI Tool-by-Tool Functional Audit"
        echo "Timestamp: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        echo "Target: $OPENWEBUI_URL"
        echo "Output Directory: $OUTPUT_DIR"
        echo
        echo "Total:   $TESTS_TOTAL"
        echo "Passed:  $TESTS_PASSED"
        echo "Failed:  $TESTS_FAILED"
        echo "Warnings:$TESTS_WARNINGS"
        echo "Skipped: $TESTS_SKIPPED"
        echo
        if [ "${#FAILURES[@]}" -gt 0 ]; then
            echo "Failures:"
            for item in "${FAILURES[@]}"; do
                echo "- $item"
            done
            echo
        fi
        if [ "${#WARNINGS[@]}" -gt 0 ]; then
            echo "Warnings:"
            for item in "${WARNINGS[@]}"; do
                echo "- $item"
            done
        fi
    } > "$summary_file"

    print_section "Summary"
    print_info "Output directory: $OUTPUT_DIR"
    print_info "Total=$TESTS_TOTAL Passed=$TESTS_PASSED Failed=$TESTS_FAILED Warnings=$TESTS_WARNINGS Skipped=$TESTS_SKIPPED"
    if [ "${#FAILURES[@]}" -gt 0 ]; then
        for item in "${FAILURES[@]}"; do
            print_error "$item"
        done
    fi
    if [ "${#WARNINGS[@]}" -gt 0 ]; then
        for item in "${WARNINGS[@]}"; do
            print_warning "$item"
        done
    fi
}

main() {
    mkdir -p "$OUTPUT_DIR"

    print_header
    print_info "Target OpenWebUI URL: $OPENWEBUI_URL"

    require_dependencies
    sign_in

    test_tools_endpoints
    test_functions_endpoints
    test_code_execution_endpoints
    test_web_search_endpoints
    cleanup_bootstrap_function_if_needed
    write_summary

    if [ "$TESTS_FAILED" -gt 0 ]; then
        exit 1
    fi
}

main "$@"
