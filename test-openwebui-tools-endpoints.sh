#!/bin/bash

###############################################################################
# OpenWebUI Tool-by-Tool Functional Audit
# Audits OpenWebUI tool endpoints (tools, functions, code execution, web search)
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-}"
OPENWEBUI_TOKEN="${OPENWEBUI_TOKEN:-}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
OPENWEBUI_DOCKER_SERVICE="${OPENWEBUI_DOCKER_SERVICE:-openwebui}"
COMPOSE_FILE="${COMPOSE_FILE:-docker-compose.yml}"
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

read_dot007_key() {
    local wanted="$1"
    local dotfile="$HOME/.007"

    [ -f "$dotfile" ] || return 1

    awk -v k="$wanted" '
        /^[[:space:]]*#/ { next }
        /^[[:space:]]*$/ { next }
        {
            line=$0
            sub(/^[[:space:]]*export[[:space:]]+/, "", line)
            sep = index(line, "=") ? "=" : (index(line, ":") ? ":" : "")
            if (sep == "") next
            split(line, parts, sep)
            key=parts[1]
            sub(/^[[:space:]]+/, "", key); sub(/[[:space:]]+$/, "", key)
            if (key != k) next
            val = substr(line, length(parts[1]) + 2)
            sub(/^[[:space:]]+/, "", val); sub(/[[:space:]]+$/, "", val)
            if ((val ~ /^".*"$/) || (val ~ /^\x27.*\x27$/)) {
                val = substr(val, 2, length(val)-2)
            }
            print val
            exit 0
        }
    ' "$dotfile"
}

docker_compose_ps_q() {
    local service="$1"

    if command_exists docker && docker compose version >/dev/null 2>&1; then
        docker compose -f "$COMPOSE_FILE" ps -q "$service" 2>/dev/null || true
        return 0
    fi

    if command_exists docker-compose; then
        docker-compose -f "$COMPOSE_FILE" ps -q "$service" 2>/dev/null || true
        return 0
    fi

    return 1
}

try_mint_openwebui_admin_jwt() {
    if ! command_exists python3 || ! command_exists docker; then
        return 1
    fi

    if [ -z "${WEBUI_SECRET_KEY:-}" ]; then
        return 1
    fi

    local cid=""
    cid="$(docker_compose_ps_q "$OPENWEBUI_DOCKER_SERVICE" | tr -d '\r' | head -n 1)"
    if [ -z "$cid" ]; then
        cid="$(docker ps -q --filter "name=^/${OPENWEBUI_DOCKER_SERVICE}$" 2>/dev/null | head -n 1 || true)"
    fi
    if [ -z "$cid" ]; then
        return 1
    fi

    local admin_id=""
    admin_id="$(docker exec "$cid" python3 -c "import sqlite3; con=sqlite3.connect('/app/backend/data/webui.db'); cur=con.cursor(); cur.execute(\"SELECT id FROM user WHERE role='admin' ORDER BY created_at ASC LIMIT 1\"); row=cur.fetchone(); print(row[0] if row else '')" 2>/dev/null | tr -d '\r' | head -n 1 || true)"
    if [ -z "$admin_id" ]; then
        return 1
    fi

    local token=""
    token="$(
        WEBUI_SECRET_KEY="$WEBUI_SECRET_KEY" OPENWEBUI_ADMIN_ID="$admin_id" python3 - <<'PY'
import base64
import hashlib
import hmac
import json
import os
import time
import uuid

secret = os.environ.get("WEBUI_SECRET_KEY", "")
user_id = os.environ.get("OPENWEBUI_ADMIN_ID", "")
if not secret or not user_id:
    raise SystemExit(1)

def b64url(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).decode("utf-8").rstrip("=")

header = {"alg": "HS256", "typ": "JWT"}
payload = {
    "id": user_id,
    "exp": int(time.time()) + 60 * 60 * 24 * 30,
    "jti": str(uuid.uuid4()),
}

header_b64 = b64url(json.dumps(header, separators=(",", ":"), sort_keys=True).encode("utf-8"))
payload_b64 = b64url(json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8"))
signing_input = f"{header_b64}.{payload_b64}".encode("utf-8")
sig = hmac.new(secret.encode("utf-8"), signing_input, hashlib.sha256).digest()
token = f"{header_b64}.{payload_b64}.{b64url(sig)}"
print(token)
PY
    )"

    if [ -z "$token" ]; then
        return 1
    fi

    API_TOKEN="$token"
    print_info "Authenticated via locally minted admin JWT (WEBUI_SECRET_KEY)"
    return 0
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
    local dot_token=""
    dot_token="$(read_dot007_key "OPENWEBUI_TOKEN" 2>/dev/null || true)"
    if [ -z "$dot_token" ]; then
        dot_token="$(read_dot007_key "OPENWEBUI_API_KEY" 2>/dev/null || true)"
    fi
    if [ -z "$dot_token" ]; then
        dot_token="$(read_dot007_key "OPENWEBUI_LOCAL_API_KEY" 2>/dev/null || true)"
    fi

    if [ -n "$OPENWEBUI_TOKEN" ]; then
        API_TOKEN="$OPENWEBUI_TOKEN"
        print_info "Using OPENWEBUI_TOKEN for auth"
        return 0
    fi
    if [ -n "$OPENWEBUI_API_KEY" ]; then
        API_TOKEN="$OPENWEBUI_API_KEY"
        print_info "Using OPENWEBUI_API_KEY for auth"
        return 0
    fi
    if [ -n "$dot_token" ]; then
        API_TOKEN="$dot_token"
        print_info "Using ~/.007 OpenWebUI token for auth"
        return 0
    fi

    if is_true "$OPENWEBUI_AUTO_AUTH"; then
        local signin_output payload_file signin_url signin_code role

        signin_output="$OUTPUT_DIR/auth_signin.json"
        payload_file="$OUTPUT_DIR/auth_signin.payload.json"
        signin_url="${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"

        jq -n \
            --arg email "$OPENWEBUI_SIGNIN_EMAIL" \
            --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
            '{email: $email, password: $password}' >"$payload_file"

        signin_code="$(curl -sS -m "$CURL_TIMEOUT" -X POST "$signin_url" \
            -H "Content-Type: application/json" \
            --data @"$payload_file" \
            -o "$signin_output" -w "%{http_code}" 2>"$OUTPUT_DIR/auth_signin.curl.err" || true)"

        if [ "$signin_code" = "200" ]; then
            API_TOKEN="$(jq -r '.token // empty' "$signin_output" 2>/dev/null || true)"
            role="$(jq -r '.role // empty' "$signin_output" 2>/dev/null || true)"
            if [ -n "$API_TOKEN" ]; then
                if [ -n "$role" ]; then
                    print_info "Authenticated with role: $role"
                else
                    print_info "Authenticated successfully"
                fi
                return 0
            fi
            print_warning "Signin succeeded but no token returned"
        else
            local detail=""
            detail="$(jq -r '.detail // .error.message // .message // empty' "$signin_output" 2>/dev/null || true)"
            print_warning "OpenWebUI signin failed (HTTP ${signin_code}${detail:+: $detail}); trying local JWT fallback"
        fi
    else
        print_warning "OPENWEBUI_AUTO_AUTH=false and no token provided; trying local JWT fallback"
    fi

    if try_mint_openwebui_admin_jwt; then
        return 0
    fi

    # If auth is not required, allow proceeding without a token.
    local probe_file="$OUTPUT_DIR/auth_probe_models.json"
    local probe_code=""
    probe_code="$(curl -sS -m "$CURL_TIMEOUT" -o "$probe_file" -w "%{http_code}" "${OPENWEBUI_URL%/}/api/models" 2>"$OUTPUT_DIR/auth_probe_models.curl.err" || true)"
    if [ "$probe_code" = "200" ]; then
        print_warning "Proceeding without auth token (OpenWebUI /api/models is accessible unauthenticated)"
        return 0
    fi

    print_error "Unable to authenticate to OpenWebUI (set OPENWEBUI_TOKEN/OPENWEBUI_API_KEY, configure OPENWEBUI_SIGNIN_EMAIL/OPENWEBUI_SIGNIN_PASSWORD, or ensure WEBUI_SECRET_KEY is available for local JWT minting)"
    exit 1
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

    function_content="$(
        cat <<'PY'
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
        }' >"$payload_file"

    request_api "POST" "/api/v1/functions/create" "$payload_file" "functions_bootstrap_create"
    if [ "$RESPONSE_CODE" = "200" ] && jq -e --arg id "$AUDIT_BOOTSTRAP_FUNCTION_ID" '.id == $id' "$RESPONSE_FILE" >/dev/null 2>&1; then
        BOOTSTRAP_FUNCTION_CREATED="true"
        pass_test
        return 0
    fi

    detail="$(response_detail)"
    if [ "$RESPONSE_CODE" = "400" ] && [[ $detail == *"ID_TAKEN"* || $detail == *"ID already taken"* ]]; then
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
    local tool_id tool_id_encoded idx
    local details_label details_file
    local spec_names expected_specs missing

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

    # Tool-by-tool valves + specs sanity checks.
    idx=0
    while IFS= read -r tool_id; do
        [ -z "$tool_id" ] && continue

        tool_id_encoded="$(uri_encode "$tool_id")"
        details_label="tool_${idx}_details"
        details_file="$OUTPUT_DIR/${details_label}.json"

        start_test "GET /api/v1/tools/id/{id} returns tool details (${tool_id})"
        request_api "GET" "/api/v1/tools/id/${tool_id_encoded}" "" "$details_label"
        if [ "$RESPONSE_CODE" = "200" ] && jq -e 'type == "object"' "$RESPONSE_FILE" >/dev/null 2>&1; then
            pass_test
        else
            fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
            idx=$((idx + 1))
            continue
        fi

        start_test "Tool specs do not expose private'_' helpers (${tool_id})"
        if jq -r '.specs[]?.name // empty' "$details_file" 2>/dev/null | grep -q '^_'; then
            spec_names="$(jq -r '.specs[]?.name // empty' "$details_file" 2>/dev/null | tr '\n' ' ' | sed 's/[[:space:]]*$//')"
            fail_test "Found private spec(s): ${spec_names}"
        else
            pass_test
        fi

        # Enforce expected public API per local tool (guards accidental helper exposure).
        expected_specs=""
        case "$tool_id" in
            llm_council) expected_specs="council" ;;
            web_scraper) expected_specs="search_web scrape_url" ;;
            run_code) expected_specs="run_code" ;;
            tavily) expected_specs="search" ;;
            pubmed) expected_specs="search_pubmed" ;;
            arxiv) expected_specs="search_arxiv" ;;
            github) expected_specs="search_repos" ;;
        esac

        if [ -n "$expected_specs" ]; then
            start_test "Tool specs match expected public functions (${tool_id})"
            missing=""
            for fn in $expected_specs; do
                if ! jq -e --arg fn "$fn" '[.specs[]?.name // empty] | index($fn) != null' "$details_file" >/dev/null 2>&1; then
                    missing="${missing}${missing:+,}${fn}"
                fi
            done
            if [ -n "$missing" ]; then
                spec_names="$(jq -r '.specs[]?.name // empty' "$details_file" 2>/dev/null | tr '\n' ' ' | sed 's/[[:space:]]*$//')"
                fail_test "Missing expected function(s): ${missing} (actual: ${spec_names})"
            else
                # Ensure no unexpected extra functions exist.
                if [ "$(jq -r '.specs | length' "$details_file" 2>/dev/null || echo 0)" -ne "$(printf "%s\n" $expected_specs | wc -l | tr -d ' ')" ]; then
                    spec_names="$(jq -r '.specs[]?.name // empty' "$details_file" 2>/dev/null | tr '\n' ' ' | sed 's/[[:space:]]*$//')"
                    fail_test "Unexpected extra spec(s) present (actual: ${spec_names})"
                else
                    pass_test
                fi
            fi
        fi

        start_test "GET /api/v1/tools/id/{id}/valves/spec returns valves schema (${tool_id})"
        request_api "GET" "/api/v1/tools/id/${tool_id_encoded}/valves/spec" "" "tool_${idx}_valves_spec"
        if [ "$RESPONSE_CODE" = "200" ] && jq -e '(type == "object") or (. == null)' "$RESPONSE_FILE" >/dev/null 2>&1; then
            # For local tools, Valves should exist (schema object, not null).
            if [[ $tool_id != server:* ]] && jq -e '. == null' "$RESPONSE_FILE" >/dev/null 2>&1; then
                fail_test "Valves schema is null (tool likely missing Valves class)"
            else
                pass_test
            fi
        else
            fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        fi

        start_test "GET /api/v1/tools/id/{id}/valves returns valves values (${tool_id})"
        request_api "GET" "/api/v1/tools/id/${tool_id_encoded}/valves" "" "tool_${idx}_valves"
        if [ "$RESPONSE_CODE" = "200" ] && jq -e '(type == "object") or (. == null)' "$RESPONSE_FILE" >/dev/null 2>&1; then
            pass_test
        else
            fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        fi

        start_test "GET /api/v1/tools/id/{id}/valves/user/spec returns user valves schema (${tool_id})"
        request_api "GET" "/api/v1/tools/id/${tool_id_encoded}/valves/user/spec" "" "tool_${idx}_user_valves_spec"
        if [ "$RESPONSE_CODE" = "200" ] && jq -e '(type == "object") or (. == null)' "$RESPONSE_FILE" >/dev/null 2>&1; then
            pass_test
        else
            fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        fi

        start_test "GET /api/v1/tools/id/{id}/valves/user returns user valves values (${tool_id})"
        request_api "GET" "/api/v1/tools/id/${tool_id_encoded}/valves/user" "" "tool_${idx}_user_valves"
        if [ "$RESPONSE_CODE" = "200" ] && jq -e '(type == "object") or (. == null)' "$RESPONSE_FILE" >/dev/null 2>&1; then
            pass_test
        else
            fail_test "HTTP $RESPONSE_CODE :: $(response_detail)"
        fi

        idx=$((idx + 1))
    done < <(jq -r '.[].id // empty' "$OUTPUT_DIR/tools_list.json" 2>/dev/null || true)
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
    jq -n --arg code "print(2 + 3)" '{code: $code}' >"$payload_file"

    start_test "POST /api/v1/utils/code/execute executes Python snippet"
    request_api "POST" "/api/v1/utils/code/execute" "$payload_file" "code_execute"
    if [ "$RESPONSE_CODE" = "200" ]; then
        pass_test
        return
    fi

    detail="$(response_detail)"
    if [ "$RESPONSE_CODE" = "400" ] && [[ $detail == *"Code execution engine not supported"* ]]; then
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
    jq -n '{queries: ["OpenWebUI tool endpoint audit"]}' >"$payload_file"

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
    if [[ $detail == *"CERTIFICATE_VERIFY_FAILED"* ]] || [[ $detail == *"SSL"* ]] || [[ $detail == *"certificate"* ]]; then
        fail_test "SSL/CA failure in web search pipeline :: $detail"
    elif [[ $detail == *"No results found from web search"* ]]; then
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
    } >"$summary_file"

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
