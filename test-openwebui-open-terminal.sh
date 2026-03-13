#!/bin/bash

###############################################################################
# OpenWebUI Open Terminal smoke test
#
# Validates:
# - OpenWebUI health/auth
# - Terminal server config endpoint
# - Proxied Open Terminal endpoints (/health, /openapi.json, /execute)
# - Optional interactive terminal WebSocket auth path
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENWEBUI_HEALTH_PATH="${OPENWEBUI_HEALTH_PATH:-/health}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-}"
OPENWEBUI_TOKEN="${OPENWEBUI_TOKEN:-}"
OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-}"
OPENWEBUI_JWT="${OPENWEBUI_JWT:-}"
OPENWEBUI_SERVICE="${OPENWEBUI_SERVICE:-${OPENWEBUI_DOCKER_SERVICE:-openwebui}}"
OPENWEBUI_CONTAINER_NAME="${OPENWEBUI_CONTAINER_NAME:-${OPENWEBUI_DOCKER_CONTAINER:-$OPENWEBUI_SERVICE}}"

OPEN_TERMINAL_ENABLED="${OPEN_TERMINAL_ENABLED:-false}"
OPEN_TERMINAL_SERVER_ID="${OPEN_TERMINAL_SERVER_ID:-open-terminal}"
OPEN_TERMINAL_SERVER_NAME="${OPEN_TERMINAL_SERVER_NAME:-Open Terminal}"
OPEN_TERMINAL_BACKEND_URL="${OPEN_TERMINAL_BACKEND_URL:-http://open-terminal:8000}"
OPEN_TERMINAL_AUTH_TYPE="${OPEN_TERMINAL_AUTH_TYPE:-bearer}"
OPEN_TERMINAL_PATH="${OPEN_TERMINAL_PATH:-/openapi.json}"
OPEN_TERMINAL_API_KEY="${OPEN_TERMINAL_API_KEY:-}"
OPEN_TERMINAL_ACCESS_PRINCIPAL_TYPE="${OPEN_TERMINAL_ACCESS_PRINCIPAL_TYPE:-user}"
OPEN_TERMINAL_ACCESS_PRINCIPAL_ID="${OPEN_TERMINAL_ACCESS_PRINCIPAL_ID:-*}"
OPEN_TERMINAL_ACCESS_PERMISSION="${OPEN_TERMINAL_ACCESS_PERMISSION:-read}"
OPEN_TERMINAL_SMOKE_APPLY_CONFIG="${OPEN_TERMINAL_SMOKE_APPLY_CONFIG:-true}"
OPEN_TERMINAL_SMOKE_RESTORE_CONFIG="${OPEN_TERMINAL_SMOKE_RESTORE_CONFIG:-true}"
OPEN_TERMINAL_SMOKE_WS="${OPEN_TERMINAL_SMOKE_WS:-true}"
OPEN_TERMINAL_TIMEOUT="${OPEN_TERMINAL_TIMEOUT:-45}"

API_TOKEN=""
JWT_TOKEN=""
TMP_DIR=""
PREV_CONFIG_FILE=""
PENDING_CONFIG_RESTORE=false
SESSION_ID=""

RESPONSE_CODE=""
RESPONSE_FILE=""
RESPONSE_ERR_FILE=""

request_api() {
    local method="$1"
    local path="$2"
    local payload_file="${3:-}"
    local label="$4"
    local url="${OPENWEBUI_URL%/}${path}"

    RESPONSE_FILE="${TMP_DIR}/${label}.json"
    RESPONSE_ERR_FILE="${TMP_DIR}/${label}.curl.err"

    local -a args=(
        -sS
        -m "$OPEN_TERMINAL_TIMEOUT"
        -X "$method"
        "$url"
        -H "Accept: application/json"
        -o "$RESPONSE_FILE"
        -w "%{http_code}"
    )

    if [ -n "$API_TOKEN" ]; then
        args+=(-H "Authorization: Bearer $API_TOKEN")
    fi

    if [ -n "$payload_file" ]; then
        args+=(-H "Content-Type: application/json" --data @"$payload_file")
    fi

    RESPONSE_CODE="$(curl "${args[@]}" 2>"$RESPONSE_ERR_FILE" || true)"
}

request_signin_jwt() {
    local payload_file="${TMP_DIR}/signin_payload.json"
    local out_file="${TMP_DIR}/signin_response.json"
    local err_file="${TMP_DIR}/signin_response.curl.err"
    local code=""

    jq -n \
        --arg email "$OPENWEBUI_SIGNIN_EMAIL" \
        --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
        '{email:$email,password:$password}' >"$payload_file"

    code="$(curl -sS -m "$OPEN_TERMINAL_TIMEOUT" \
        -X POST "${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}" \
        -H "Content-Type: application/json" \
        -H "Accept: application/json" \
        --data @"$payload_file" \
        -o "$out_file" \
        -w "%{http_code}" 2>"$err_file" || true)"

    if [ "$code" != "200" ]; then
        return 1
    fi

    JWT_TOKEN="$(jq -r '.token // empty' "$out_file" 2>/dev/null || true)"
    [ -n "$JWT_TOKEN" ] || return 1
    return 0
}

restore_terminal_config() {
    [ "$PENDING_CONFIG_RESTORE" = true ] || return 0
    [ -n "$PREV_CONFIG_FILE" ] || return 0
    [ -f "$PREV_CONFIG_FILE" ] || return 0
    [ -n "$API_TOKEN" ] || return 0

    local restore_out="${TMP_DIR}/restore_config.json"
    local restore_err="${TMP_DIR}/restore_config.curl.err"
    local restore_code=""
    local endpoint="${OPENWEBUI_URL%/}/api/v1/configs/terminal_servers"

    restore_code="$(curl -sS -m "$OPEN_TERMINAL_TIMEOUT" \
        -X POST "$endpoint" \
        -H "Authorization: Bearer $API_TOKEN" \
        -H "Content-Type: application/json" \
        -H "Accept: application/json" \
        --data @"$PREV_CONFIG_FILE" \
        -o "$restore_out" \
        -w "%{http_code}" 2>"$restore_err" || true)"

    if [ "$restore_code" = "200" ]; then
        print_info "Restored previous terminal server configuration"
    else
        print_warning "Could not restore previous terminal configuration (HTTP ${restore_code})"
    fi
}

cleanup() {
    set +e

    if [ -n "$SESSION_ID" ] && [ -n "$API_TOKEN" ] && [ -n "$TMP_DIR" ]; then
        curl -sS -m "$OPEN_TERMINAL_TIMEOUT" \
            -X DELETE "${OPENWEBUI_URL%/}/api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/api/terminals/${SESSION_ID}" \
            -H "Authorization: Bearer $API_TOKEN" \
            -H "Accept: application/json" \
            -o "${TMP_DIR}/session_delete.json" \
            -w "%{http_code}" >/dev/null 2>"${TMP_DIR}/session_delete.curl.err" || true
    fi

    if is_true "$OPEN_TERMINAL_SMOKE_RESTORE_CONFIG"; then
        restore_terminal_config || true
    fi

    if [ -n "$TMP_DIR" ] && [ -d "$TMP_DIR" ]; then
        rm -rf "$TMP_DIR"
    fi
}

assert_code() {
    local expected_regex="$1"
    local context="$2"
    if ! printf "%s" "$RESPONSE_CODE" | grep -Eq "^(${expected_regex})$"; then
        print_error "${context} failed (HTTP ${RESPONSE_CODE})"
        if [ -s "$RESPONSE_ERR_FILE" ]; then
            print_info "curl error: $(head -c 200 "$RESPONSE_ERR_FILE" 2>/dev/null)"
        fi
        print_info "response: $(head -c 240 "$RESPONSE_FILE" 2>/dev/null)"
        exit 1
    fi
}

require_prereqs() {
    local required=(curl jq docker)
    local cmd
    for cmd in "${required[@]}"; do
        if ! command_exists "$cmd"; then
            print_error "Required command not found: ${cmd}"
            exit 1
        fi
    done
}

authenticate() {
    if [ -n "$OPENWEBUI_TOKEN" ]; then
        API_TOKEN="$OPENWEBUI_TOKEN"
        print_info "Using OPENWEBUI_TOKEN for API auth"
    elif [ -n "$OPENWEBUI_API_KEY" ]; then
        API_TOKEN="$OPENWEBUI_API_KEY"
        print_info "Using OPENWEBUI_API_KEY for API auth"
    fi

    if [ -n "$OPENWEBUI_JWT" ]; then
        JWT_TOKEN="$OPENWEBUI_JWT"
        print_info "Using OPENWEBUI_JWT for websocket auth"
    fi

    if [ -z "$JWT_TOKEN" ] && [ -n "$OPENWEBUI_SIGNIN_EMAIL" ] && [ -n "$OPENWEBUI_SIGNIN_PASSWORD" ]; then
        if request_signin_jwt; then
            print_success "Obtained JWT via signin"
        else
            print_warning "Signin did not return JWT; websocket check may be skipped"
        fi
    fi

    if [ -z "$API_TOKEN" ] && [ -n "$JWT_TOKEN" ]; then
        API_TOKEN="$JWT_TOKEN"
        print_info "Using JWT token for API auth"
    fi

    if [ -z "$API_TOKEN" ]; then
        print_error "No API auth token available. Set OPENWEBUI_API_KEY or OPENWEBUI_TOKEN (or valid signin credentials)."
        exit 1
    fi
}

main() {
    print_header "OpenWebUI Open Terminal Smoke Test"
    require_prereqs

    TMP_DIR="$(mktemp -d /tmp/open_terminal_smoke_XXXXXX)"
    PREV_CONFIG_FILE="${TMP_DIR}/terminal_config_prev.json"
    trap cleanup EXIT

    print_step "Checking OpenWebUI health"
    if ! curl -sS -f -m "$OPEN_TERMINAL_TIMEOUT" "${OPENWEBUI_URL%/}${OPENWEBUI_HEALTH_PATH}" >/dev/null 2>&1; then
        print_error "OpenWebUI health check failed at ${OPENWEBUI_URL%/}${OPENWEBUI_HEALTH_PATH}"
        exit 1
    fi
    print_success "OpenWebUI is healthy"

    authenticate

    if is_true "$OPEN_TERMINAL_ENABLED"; then
        local openwebui_container_id=""

        openwebui_container_id="$(resolve_container_id "$OPENWEBUI_SERVICE" "$OPENWEBUI_CONTAINER_NAME" || true)"
        if [ -z "$openwebui_container_id" ]; then
            print_error "Could not resolve a running OpenWebUI container for service '${OPENWEBUI_SERVICE}'"
            exit 1
        fi

        print_step "Checking Open Terminal reachability from OpenWebUI container"
        if ! docker exec "$openwebui_container_id" sh -lc "curl -sS -m 8 ${OPEN_TERMINAL_BACKEND_URL%/}/health >/dev/null"; then
            print_error "OpenWebUI container cannot reach ${OPEN_TERMINAL_BACKEND_URL%/}/health"
            exit 1
        fi
        print_success "OpenWebUI container can reach Open Terminal backend"
    fi

    print_step "Reading current terminal server config"
    request_api "GET" "/api/v1/configs/terminal_servers" "" "terminal_config_get"
    assert_code "200" "GET /api/v1/configs/terminal_servers"
    cp "$RESPONSE_FILE" "$PREV_CONFIG_FILE"

    if is_true "$OPEN_TERMINAL_SMOKE_APPLY_CONFIG"; then
        print_step "Applying terminal server config for smoke test"
        if [ -z "$OPEN_TERMINAL_API_KEY" ]; then
            print_error "OPEN_TERMINAL_API_KEY is required when OPEN_TERMINAL_SMOKE_APPLY_CONFIG=true"
            exit 1
        fi

        local config_payload="${TMP_DIR}/terminal_config_apply.json"
        jq -n \
            --arg id "$OPEN_TERMINAL_SERVER_ID" \
            --arg name "$OPEN_TERMINAL_SERVER_NAME" \
            --arg url "$OPEN_TERMINAL_BACKEND_URL" \
            --arg path "$OPEN_TERMINAL_PATH" \
            --arg key "$OPEN_TERMINAL_API_KEY" \
            --arg auth_type "$OPEN_TERMINAL_AUTH_TYPE" \
            --arg principal_type "$OPEN_TERMINAL_ACCESS_PRINCIPAL_TYPE" \
            --arg principal_id "$OPEN_TERMINAL_ACCESS_PRINCIPAL_ID" \
            --arg permission "$OPEN_TERMINAL_ACCESS_PERMISSION" \
            '{
                TERMINAL_SERVER_CONNECTIONS: [
                    {
                        id: $id,
                        name: $name,
                        enabled: true,
                        url: $url,
                        path: $path,
                        key: $key,
                        auth_type: $auth_type,
                        config: {
                            access_grants: [
                                {
                                    principal_type: $principal_type,
                                    principal_id: $principal_id,
                                    permission: $permission
                                }
                            ]
                        }
                    }
                ]
            }' >"$config_payload"

        request_api "POST" "/api/v1/configs/terminal_servers" "$config_payload" "terminal_config_set"
        assert_code "200" "POST /api/v1/configs/terminal_servers"
        PENDING_CONFIG_RESTORE=true
        print_success "Terminal server config applied"
    else
        print_info "Skipping config apply (OPEN_TERMINAL_SMOKE_APPLY_CONFIG=false)"
    fi

    print_step "Listing terminal servers via OpenWebUI proxy"
    request_api "GET" "/api/v1/terminals/" "" "terminals_list"
    assert_code "200" "GET /api/v1/terminals/"
    if ! jq -e --arg id "$OPEN_TERMINAL_SERVER_ID" '[.. | objects | select(.id? == $id)] | length > 0' "$RESPONSE_FILE" >/dev/null 2>&1; then
        if ! is_true "$OPEN_TERMINAL_SMOKE_APPLY_CONFIG" && ! is_true "$OPEN_TERMINAL_ENABLED"; then
            print_warning "Terminal server id not registered while Open Terminal is disabled; skipping proxied checks"
            print_info "response: $(head -c 400 "$RESPONSE_FILE" 2>/dev/null)"
            print_success "Open Terminal smoke test passed (disabled/no registered terminal server)"
            exit 0
        fi
        print_error "Terminal server id not found in /api/v1/terminals/: ${OPEN_TERMINAL_SERVER_ID}"
        print_info "response: $(head -c 400 "$RESPONSE_FILE" 2>/dev/null)"
        exit 1
    fi
    print_success "Terminal server id is visible: ${OPEN_TERMINAL_SERVER_ID}"

    print_step "Checking proxied Open Terminal health"
    request_api "GET" "/api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/health" "" "terminal_health"
    assert_code "200" "GET /api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/health"
    print_success "Proxied terminal health endpoint OK"

    print_step "Checking proxied Open Terminal OpenAPI document"
    request_api "GET" "/api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/openapi.json" "" "terminal_openapi"
    assert_code "200" "GET /api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/openapi.json"
    if ! jq -e '.paths and (.paths | type == "object") and ((.paths | keys | length) > 0)' "$RESPONSE_FILE" >/dev/null 2>&1; then
        print_error "OpenAPI payload did not contain endpoints"
        exit 1
    fi
    print_success "OpenAPI proxy payload is valid"

    print_step "Executing command through OpenWebUI terminal proxy"
    local nonce
    nonce="open_terminal_smoke_$(date +%s)"
    local execute_payload="${TMP_DIR}/terminal_execute_payload.json"
    jq -n --arg cmd "echo ${nonce} && pwd && whoami" '{command:$cmd}' >"$execute_payload"
    request_api "POST" "/api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/execute?wait=5" "$execute_payload" "terminal_execute"
    assert_code "200" "POST /api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/execute?wait=5"
    if ! grep -q "$nonce" "$RESPONSE_FILE"; then
        print_error "Execute response did not include expected nonce: ${nonce}"
        print_info "response: $(head -c 400 "$RESPONSE_FILE" 2>/dev/null)"
        exit 1
    fi
    print_success "Execute command returned expected nonce"

    if ! is_true "$OPEN_TERMINAL_SMOKE_WS"; then
        print_info "Skipping websocket check (OPEN_TERMINAL_SMOKE_WS=false)"
        print_success "Open Terminal smoke test passed"
        exit 0
    fi

    if ! command_exists node; then
        print_warning "Node.js not found; skipping websocket check"
        print_success "Open Terminal smoke test passed (without websocket check)"
        exit 0
    fi

    if [ -z "$JWT_TOKEN" ]; then
        print_warning "No JWT token available; skipping websocket check"
        print_success "Open Terminal smoke test passed (without websocket check)"
        exit 0
    fi

    print_step "Creating terminal session"
    request_api "POST" "/api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/api/terminals" "" "terminal_session_create"
    assert_code "200|201" "POST /api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/api/terminals"
    SESSION_ID="$(jq -r '.id // .session_id // empty' "$RESPONSE_FILE" 2>/dev/null || true)"
    if [ -z "$SESSION_ID" ]; then
        print_error "Could not parse terminal session id"
        print_info "response: $(head -c 300 "$RESPONSE_FILE" 2>/dev/null)"
        exit 1
    fi
    print_success "Created terminal session: ${SESSION_ID}"

    print_step "Verifying websocket auth path"
    local ws_base=""
    case "$OPENWEBUI_URL" in
        https://*) ws_base="${OPENWEBUI_URL/https:\/\//wss://}" ;;
        http://*) ws_base="${OPENWEBUI_URL/http:\/\//ws://}" ;;
        *) ws_base="ws://${OPENWEBUI_URL}" ;;
    esac
    ws_base="${ws_base%/}"
    local ws_url="${ws_base}/api/v1/terminals/${OPEN_TERMINAL_SERVER_ID}/api/terminals/${SESSION_ID}"

    WS_URL="$ws_url" OPENWEBUI_JWT="$JWT_TOKEN" node <<'NODE'
const wsUrl = process.env.WS_URL;
const token = process.env.OPENWEBUI_JWT;

if (typeof WebSocket === "undefined") {
  console.error("Global WebSocket is unavailable in this Node runtime");
  process.exit(1);
}

let opened = false;
const ws = new WebSocket(wsUrl);

const fail = (msg) => {
  console.error(msg);
  process.exit(1);
};

const timer = setTimeout(() => fail("WebSocket timeout"), 12000);

ws.onopen = () => {
  opened = true;
  ws.send(JSON.stringify({ type: "auth", token }));
  setTimeout(() => ws.close(), 1000);
};

ws.onerror = () => fail("WebSocket error");

ws.onclose = (ev) => {
  clearTimeout(timer);
  if (!opened) fail("WebSocket never opened");
  if ([4001, 4003, 4004].includes(ev.code)) {
    fail(`WebSocket auth/session failure (code=${ev.code})`);
  }
  process.exit(0);
};
NODE
    print_success "WebSocket proxy/auth check passed"

    print_success "Open Terminal smoke test passed"
}

main "$@"
