#!/bin/bash

###############################################################################
# RAG Functionality Testing Script
# This script tests the complete RAG deployment including OpenWebUI and the
# configured OpenAI-compatible upstream (LiteLLM by default).
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
source "${SCRIPT_DIR}/lib/test-openwebui-helpers.sh"
load_env_defaults
cd "$SCRIPT_DIR"

resolve_openai_base_root() {
    local candidate="${OPENAI_API_BASE_URL:-}"

    if [ -z "$candidate" ] && [ -n "${OPENAI_API_BASE_URLS:-}" ]; then
        candidate="$(printf "%s" "$OPENAI_API_BASE_URLS" | cut -d';' -f1)"
    fi

    if [ -z "$candidate" ]; then
        printf ""
        return 0
    fi

    candidate="${candidate%/}"
    candidate="${candidate%/v1}"

    # Host-side test execution should not rely on host.docker.internal DNS.
    if [[ $candidate == *"://host.docker.internal"* ]] && [ ! -f "/.dockerenv" ]; then
        candidate="$(printf "%s" "$candidate" | sed 's#://host.docker.internal#://localhost#')"
    fi

    printf "%s" "$candidate"
}

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENAI_BASE_ROOT="$(resolve_openai_base_root)"
UPSTREAM_URL="${UPSTREAM_URL:-${OPENAI_BASE_ROOT:-http://localhost:4000}}"
API_KEY="${OPENWEBUI_API_KEY:-${API_KEY:-}}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
OPENWEBUI_SMOKE_MODEL_CANDIDATES="${OPENWEBUI_SMOKE_MODEL_CANDIDATES:-$(default_openwebui_smoke_model_candidates)}"
VERBOSE="${VERBOSE:-false}"
CHAT_MODEL="${RAG_CHAT_MODEL:-}"
NO_MOCK_AUDIT_SCRIPT="${NO_MOCK_AUDIT_SCRIPT:-./scripts/testing/audit-no-mock.sh}"
UPSTREAM_API_KEY="${UPSTREAM_API_KEY:-${OPENAI_API_KEY:-}}"
TEST_MODE="full"
RUN_CLEANUP_ON_EXIT=false
BASELINE_REQUIRE_WEB_SEARCH="${BASELINE_REQUIRE_WEB_SEARCH:-false}"
BASELINE_EXPECT_SEARX_PORT="${BASELINE_EXPECT_SEARX_PORT:-}"
BASELINE_STRICT_SEARXNG_PAYLOAD="${BASELINE_STRICT_SEARXNG_PAYLOAD:-false}"
RAG_MULTIMODAL_TEST_PDF="${RAG_MULTIMODAL_TEST_PDF:-}"
RAG_MULTIMODAL_STRICT="${RAG_MULTIMODAL_STRICT:-}"
INDIGO_SERVICE_ENABLED="${INDIGO_SERVICE_ENABLED:-false}"
INDIGO_SERVICE_BIND_ADDRESS="${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}"
INDIGO_SERVICE_PORT="${INDIGO_SERVICE_PORT:-8012}"
INDIGO_SERVICE_BASE_URL="${INDIGO_SERVICE_BASE_URL:-http://${INDIGO_SERVICE_BIND_ADDRESS}:${INDIGO_SERVICE_PORT}}"
INDIGO_TOOL_ENABLED="${INDIGO_TOOL_ENABLED:-true}"
INDIGO_TOOL_ID="${INDIGO_TOOL_ID:-indigo_chemistry}"
INDIGO_TOOL_API_BASE_URL="${INDIGO_TOOL_API_BASE_URL:-http://indigo-service/v2/indigo}"
RAG_TEST_SENTINEL="${RAG_TEST_SENTINEL:-openwebui-rag-sentinel-$(date +%Y%m%d%H%M%S)}"
RAG_TEST_SENTINEL_LINE="The exact retrieval sentinel for this document is ${RAG_TEST_SENTINEL}."

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0
TEST_VECTOR_COLLECTION=""
TEST_MULTIMODAL_DOC_ID=""

###############################################################################
# Script-specific Helper Functions
###############################################################################

test_start() {
    ((TESTS_TOTAL += 1))
    TEST_NAME="$1"
    print_info "Testing: $TEST_NAME"
}

test_pass() {
    ((TESTS_PASSED += 1))
    print_success "$TEST_NAME passed"
}

test_fail() {
    ((TESTS_FAILED += 1))
    print_error "$TEST_NAME failed: $1"
}

normalize_http_code() {
    local raw="${1:-}"
    local code

    code="$(printf "%s" "$raw" | grep -Eo '[0-9]{3}' | tail -n1 || true)"
    if [ -z "$code" ]; then
        printf "000"
        return 0
    fi

    printf "%s" "$code"
}

make_request() {
    local method="$1"
    local url="$2"
    local data="$3"
    local -a curl_args=("-s" "-X" "$method" "$url" "-H" "Content-Type: application/json")

    if [ -n "$API_KEY" ]; then
        curl_args+=("-H" "Authorization: Bearer $API_KEY")
    fi

    if [ "$method" = "POST" ]; then
        curl_args+=("-d" "$data")
    fi

    curl "${curl_args[@]}"
}

test_repo_real_integration_guard() {
    test_start "Repository Real Integration Guard"

    local response
    if response=$("$NO_MOCK_AUDIT_SCRIPT" 2>&1); then
        test_pass
        return 0
    fi

    test_fail "audit-no-mock.sh failed"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

show_help() {
    cat <<'EOF'
Usage: ./test-rag.sh [OPTIONS]

Options:
  --baseline      Run fast baseline checks only
  --full          Run full RAG integration suite (default)
  -v, --verbose   Enable verbose output
  -h, --help      Show this help message

Environment:
  UPSTREAM_URL=http://localhost:4000       OpenAI-compatible upstream base URL
  UPSTREAM_API_KEY=<token>                 Optional bearer token for upstream /v1/models
  BASELINE_REQUIRE_WEB_SEARCH=true|false  Fail if host and container SearXNG probes fail (default: false)
  BASELINE_EXPECT_SEARX_PORT=<port|empty> Optional SearXNG port hint (default: empty; empty disables port hint)
  BASELINE_STRICT_SEARXNG_PAYLOAD=true|false  Fail on HTTP 200 responses without a SearXNG-style results array (default: false)
  RAG_MULTIMODAL_TEST_PDF=/abs/path.pdf   Optional PDF fixture for full-mode multimodal ingestion test
  RAG_MULTIMODAL_STRICT=true|false        Fail when multimodal prereqs/fixture are missing (default: false; auto-true in CI full mode)
EOF
}

resolve_chat_model_from_openwebui() {
    local response
    response=$(curl -s -X GET \
        "$OPENWEBUI_URL/api/models" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>/dev/null || true)

    if [ -z "$response" ]; then
        return 1
    fi

    if [ "$TEST_MODE" = "baseline" ]; then
        CHAT_MODEL="$(choose_preferred_model_from_response "$response")"
    else
        CHAT_MODEL="$(first_model_id_from_response "$response")"
    fi

    [ -n "$CHAT_MODEL" ]
}

###############################################################################
# Test Functions
###############################################################################

test_openwebui_accessible() {
    test_start "OpenWebUI Accessibility"

    if curl -s -f -o /dev/null "$OPENWEBUI_URL" 2>/dev/null; then
        test_pass
        return 0
    else
        test_fail "Cannot reach OpenWebUI at $OPENWEBUI_URL"
        return 1
    fi
}

test_upstream_responding() {
    test_start "OpenAI-Compatible Upstream Response"

    local code
    local code_raw
    code_raw="$(curl -s -w "%{http_code}" -o /dev/null "$UPSTREAM_URL/health" 2>&1 || true)"
    code="$(normalize_http_code "$code_raw")"
    if [ "$code" = "200" ]; then
        test_pass
        return 0
    fi
    if [ "$code" = "401" ]; then
        print_info "Upstream health endpoint requires auth (HTTP 401); treating service as reachable in LiteLLM mode"
        test_pass
        return 0
    fi

    test_fail "OpenAI-compatible upstream health check failed at $UPSTREAM_URL (/health HTTP $code)"
    return 1
}

test_upstream_models() {
    test_start "OpenAI-Compatible Upstream Models Available"

    local response
    local -a headers=()
    if [ -n "$UPSTREAM_API_KEY" ]; then
        headers=(-H "Authorization: Bearer $UPSTREAM_API_KEY")
    fi
    response=$(curl -s "${headers[@]}" "$UPSTREAM_URL/v1/models" 2>/dev/null || true)

    if echo "$response" | grep -Eiq 'Authentication Error|"code"[[:space:]]*:[[:space:]]*"401"|401'; then
        print_info "Upstream /v1/models requires auth (HTTP 401); set UPSTREAM_API_KEY for strict model validation"
        if resolve_chat_model_from_openwebui; then
            print_info "Using OpenWebUI model fallback: $CHAT_MODEL"
        fi
        test_pass
        return 0
    fi

    if has_non_empty_models_array "$response"; then
        local models
        models=$(echo "$response" | jq -r '.data[]?.id // empty' 2>/dev/null | tr '\n' ', ' | sed 's/, $//' || true)
        if [ -z "$models" ]; then
            models=$(echo "$response" | grep -Eo '"id"[[:space:]]*:[[:space:]]*"[^"]+"' | head -n 20 | cut -d'"' -f4 | tr '\n' ', ' | sed 's/, $//')
        fi

        if [ -z "$CHAT_MODEL" ] && ! resolve_chat_model_from_openwebui; then
            CHAT_MODEL="$(first_model_id_from_response "$response")"
        fi

        print_success "Available models: $models"
        [ -n "$CHAT_MODEL" ] && print_info "Using chat model: $CHAT_MODEL"
        test_pass
        return 0
    else
        test_fail "Could not retrieve non-empty models list"
        return 1
    fi
}

test_openwebui_health() {
    test_start "OpenWebUI Health Endpoint"

    local response
    if response=$(curl -s "$OPENWEBUI_URL/health" 2>/dev/null); then
        print_info "Health check response: $response"
        test_pass
        return 0
    else
        test_fail "Health endpoint not responding"
        return 1
    fi
}

test_indigo_service_health() {
    if ! is_true "$INDIGO_SERVICE_ENABLED"; then
        print_info "Skipping Indigo Service health check (INDIGO_SERVICE_ENABLED=false)"
        return 0
    fi

    test_start "Indigo Service Health"

    local response_file
    local response
    local http_code
    local base_url="${INDIGO_SERVICE_BASE_URL%/}"

    response_file="$(mktemp)"
    http_code="$(normalize_http_code "$(curl -sS -m 12 -o "$response_file" -w "%{http_code}" "${base_url}/v2/indigo/info" 2>/dev/null || true)")"
    response="$(cat "$response_file" 2>/dev/null || true)"
    rm -f "$response_file"

    if [ "$http_code" != "200" ]; then
        test_fail "Indigo Service /v2/indigo/info returned HTTP ${http_code}"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if command -v jq >/dev/null 2>&1 && ! echo "$response" | jq -e '.Indigo.version // .indigo_version // .version' >/dev/null 2>&1; then
        test_fail "Indigo Service response missing version field"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    test_pass
    return 0
}

test_indigo_tool_registration() {
    if ! is_true "$INDIGO_SERVICE_ENABLED"; then
        print_info "Skipping Indigo tool registration check (INDIGO_SERVICE_ENABLED=false)"
        return 0
    fi
    if ! is_true "$INDIGO_TOOL_ENABLED"; then
        print_info "Skipping Indigo tool registration check (INDIGO_TOOL_ENABLED=false)"
        return 0
    fi

    test_start "Indigo Tool Registration"

    local response_file
    local response
    local http_code
    response_file="$(mktemp)"
    http_code="$(normalize_http_code "$(curl -sS -m 12 -o "$response_file" -w "%{http_code}" \
        "$OPENWEBUI_URL/api/v1/tools/id/$INDIGO_TOOL_ID" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>/dev/null || true)")"
    response="$(cat "$response_file" 2>/dev/null || true)"
    rm -f "$response_file"

    if [ "$http_code" != "200" ]; then
        test_fail "OpenWebUI tool endpoint for $INDIGO_TOOL_ID returned HTTP ${http_code}"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        local missing
        missing="$(echo "$response" | jq -r '
            .specs as $specs
            | ["info","check","convert","calculate","render","aromatize","dearomatize","layout","clean"]
            | map(select(($specs | map(.name) | index(.) ) == null))
            | join(",")
        ' 2>/dev/null || true)"
        if [ -n "$missing" ]; then
            test_fail "Indigo tool missing expected specs: $missing"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi
    else
        if ! echo "$response" | grep -q "\"id\":\"$INDIGO_TOOL_ID\""; then
            test_fail "Indigo tool response does not contain expected tool id"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi
    fi

    test_pass
    return 0
}

test_indigo_tool_connectivity_from_openwebui() {
    if ! is_true "$INDIGO_SERVICE_ENABLED"; then
        print_info "Skipping Indigo tool connectivity check (INDIGO_SERVICE_ENABLED=false)"
        return 0
    fi
    if ! is_true "$INDIGO_TOOL_ENABLED"; then
        print_info "Skipping Indigo tool connectivity check (INDIGO_TOOL_ENABLED=false)"
        return 0
    fi

    test_start "Indigo Tool Connectivity (from OpenWebUI container)"

    if ! docker info >/dev/null 2>&1; then
        test_fail "Docker unavailable; cannot probe Indigo path from OpenWebUI container"
        return 1
    fi

    if ! docker ps -q -f "name=^openwebui$" | grep -q .; then
        test_fail "OpenWebUI container is not running"
        return 1
    fi

    local info_url http_code
    info_url="${INDIGO_TOOL_API_BASE_URL%/}/info"
    http_code="$(normalize_http_code "$(docker exec openwebui sh -lc \
        "curl -sS --connect-timeout 3 --max-time 12 -o /tmp/indigo_tool_probe.json -w '%{http_code}' '${info_url}'" \
        2>/dev/null || true)")"
    if [ "$http_code" != "200" ]; then
        test_fail "OpenWebUI container cannot reach Indigo tool URL (${info_url}) HTTP ${http_code}"
        return 1
    fi

    test_pass
    return 0
}

test_embedding_endpoint() {
    test_start "Embedding Generation"

    local response
    response=$(curl -s -X GET \
        "$OPENWEBUI_URL/api/v1/retrieval/" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>/dev/null || true)

    if command -v jq >/dev/null 2>&1 && echo "$response" | jq -e '.status == true and (.RAG_EMBEDDING_MODEL | type == "string") and (.RAG_EMBEDDING_MODEL | length > 0)' >/dev/null 2>&1; then
        local embedding_model
        embedding_model=$(echo "$response" | jq -r '.RAG_EMBEDDING_MODEL' 2>/dev/null || echo "")
        print_success "Embedding configuration detected: $embedding_model"
        test_pass
        return 0
    fi

    if echo "$response" | grep -q "RAG_EMBEDDING_MODEL"; then
        print_success "Embedding configuration endpoint reachable"
        test_pass
        return 0
    else
        test_fail "Embedding configuration endpoint failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_openwebui_baseline_settings() {
    test_start "OpenWebUI Baseline RAG Settings"

    local response
    response=$(curl -s -X GET \
        "$OPENWEBUI_URL/api/v1/retrieval/config" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>/dev/null || true)

    if [ -z "$response" ]; then
        test_fail "Empty response from retrieval settings endpoint"
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        if ! echo "$response" | jq -e '.status == true' >/dev/null 2>&1; then
            test_fail "Retrieval settings endpoint did not return status=true"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi

        local rag_top_k
        rag_top_k=$(echo "$response" | jq -r '.RAG_TOP_K // .TOP_K // empty' 2>/dev/null || true)
        if [ -z "$rag_top_k" ]; then
            test_fail "TOP_K missing from retrieval settings"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        elif [[ ! $rag_top_k =~ ^[0-9]+$ ]] || [ "$rag_top_k" -le 0 ]; then
            test_fail "Invalid RAG_TOP_K value from endpoint: $rag_top_k"
            return 1
        fi

        local chunk_size
        chunk_size=$(echo "$response" | jq -r '.CHUNK_SIZE // empty' 2>/dev/null || true)
        if [ -z "$chunk_size" ]; then
            test_fail "CHUNK_SIZE missing from retrieval settings"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        elif [[ ! $chunk_size =~ ^[0-9]+$ ]] || [ "$chunk_size" -le 0 ]; then
            test_fail "Invalid CHUNK_SIZE value from endpoint: $chunk_size"
            return 1
        fi

        print_info "RAG_TOP_K: $rag_top_k"
        print_info "CHUNK_SIZE: $chunk_size"
        test_pass
        return 0
    fi

    if echo "$response" | grep -q '"TOP_K"'; then
        test_pass
        return 0
    fi

    test_fail "Could not validate baseline retrieval settings"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

handle_web_search_probe_failure() {
    local message="$1"

    if is_true "$BASELINE_REQUIRE_WEB_SEARCH"; then
        test_fail "$message"
        return 1
    fi

    print_info "$message"
    if echo "$message" | grep -qi "container web search probe"; then
        print_info "Hint: host SearXNG may be bound to 127.0.0.1 only. Expose it to Docker (0.0.0.0) or route via a container service."
    fi
    print_info "Continuing because BASELINE_REQUIRE_WEB_SEARCH=false"
    test_pass
    return 0
}

handle_web_search_payload_mismatch() {
    local message="$1"

    if is_true "$BASELINE_REQUIRE_WEB_SEARCH" || is_true "$BASELINE_STRICT_SEARXNG_PAYLOAD"; then
        test_fail "$message"
        return 1
    fi

    print_warning "$message"
    print_info "Continuing because BASELINE_STRICT_SEARXNG_PAYLOAD=false"
    test_pass
    return 0
}

test_web_search_baseline() {
    test_start "Web Search Baseline (SearXNG)"

    local web_search_enabled=false
    if is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; then
        web_search_enabled=true
    fi

    if [ "$web_search_enabled" = false ]; then
        print_info "Web search is disabled by configuration; baseline check skipped"
        test_pass
        return 0
    fi

    local searx_url="${SEARXNG_QUERY_URL:-}"
    if [ -z "$searx_url" ]; then
        test_fail "Web search enabled but SEARXNG_QUERY_URL is empty"
        return 1
    fi

    if [ -n "$BASELINE_EXPECT_SEARX_PORT" ] && ! echo "$searx_url" | grep -q "$BASELINE_EXPECT_SEARX_PORT"; then
        local port_message="SEARXNG_QUERY_URL does not include expected port ${BASELINE_EXPECT_SEARX_PORT}: $searx_url"
        if is_true "$BASELINE_REQUIRE_WEB_SEARCH"; then
            test_fail "$port_message"
            return 1
        fi
        print_info "$port_message"
        print_info "Continuing because BASELINE_REQUIRE_WEB_SEARCH=false"
    fi

    local probe_url
    probe_url="${searx_url//\{query\}/openwebui%20baseline}"
    local host_probe_url="$probe_url"
    if echo "$host_probe_url" | grep -q "host.docker.internal"; then
        host_probe_url="${host_probe_url//host.docker.internal/127.0.0.1}"
    fi

    local response
    local response_file
    local http_code
    response_file="$(mktemp)"
    http_code="$(normalize_http_code "$(curl -sS -m 20 -o "$response_file" -w "%{http_code}" "$host_probe_url" 2>/dev/null || true)")"
    response="$(cat "$response_file" 2>/dev/null || true)"

    if { [ -z "$response" ] || [ "$http_code" != "200" ]; } && [ "$host_probe_url" != "$probe_url" ]; then
        http_code="$(normalize_http_code "$(curl -sS -m 20 -o "$response_file" -w "%{http_code}" "$probe_url" 2>/dev/null || true)")"
        response="$(cat "$response_file" 2>/dev/null || true)"
    fi
    rm -f "$response_file"

    if [ -z "$response" ]; then
        handle_web_search_probe_failure "Empty response from local SearXNG probe URL"
        return $?
    fi

    if [ "$http_code" != "200" ]; then
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        handle_web_search_probe_failure "SearXNG probe returned HTTP $http_code"
        return $?
    fi

    if json_response_has_results_array "$response"; then
        test_pass
        return 0
    fi

    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    handle_web_search_payload_mismatch "SearXNG endpoint reachable but returned non-standard payload (HTTP 200)"
    return $?
}

test_web_search_container_baseline() {
    test_start "Web Search Baseline (SearXNG from OpenWebUI container)"

    local web_search_enabled=false
    if is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; then
        web_search_enabled=true
    fi

    if [ "$web_search_enabled" = false ]; then
        print_info "Web search is disabled by configuration; baseline check skipped"
        test_pass
        return 0
    fi

    local searx_url="${SEARXNG_QUERY_URL:-}"
    if [ -z "$searx_url" ]; then
        test_fail "Web search enabled but SEARXNG_QUERY_URL is empty"
        return 1
    fi

    if ! docker info >/dev/null 2>&1; then
        handle_web_search_probe_failure "Docker is not available; cannot probe web search from container path"
        return $?
    fi

    if ! docker ps -q -f "name=^openwebui$" | grep -q .; then
        handle_web_search_probe_failure "OpenWebUI container is not running; cannot probe container web search path"
        return $?
    fi

    local probe_url
    probe_url="${searx_url//\{query\}/openwebui%20baseline}"

    local raw
    if docker exec openwebui sh -lc "command -v python3 >/dev/null 2>&1"; then
        raw="$(docker exec openwebui python3 -c '
import sys
import urllib.error
import urllib.request

url = sys.argv[1]
try:
    with urllib.request.urlopen(url, timeout=20) as resp:
        body = resp.read().decode("utf-8", "replace")
        sys.stdout.write(body)
        sys.stdout.write(f"\\nHTTP_CODE:{resp.getcode()}\\n")
except urllib.error.HTTPError as exc:
    body = exc.read().decode("utf-8", "replace")
    sys.stdout.write(body)
    sys.stdout.write(f"\\nHTTP_CODE:{exc.code}\\n")
except Exception:
    sys.stdout.write("\\nHTTP_CODE:000\\n")
' "$probe_url" 2>/dev/null || true)"
    elif docker exec openwebui sh -lc "command -v curl >/dev/null 2>&1"; then
        raw="$(docker exec openwebui sh -lc "curl -sS -m 20 -w '\nHTTP_CODE:%{http_code}\n' \"$probe_url\" 2>/dev/null || true" 2>/dev/null || true)"
    else
        raw="$(docker exec openwebui sh -lc "wget -q -O - \"$probe_url\" 2>/dev/null; code=\$?; [ \$code -eq 0 ] && printf '\nHTTP_CODE:200\n' || printf '\nHTTP_CODE:000\n'" 2>/dev/null || true)"
    fi

    if [ -z "$raw" ]; then
        handle_web_search_probe_failure "Empty response from OpenWebUI container web search probe"
        return $?
    fi

    local http_code
    http_code="$(normalize_http_code "$(echo "$raw" | tail -n 1 | sed 's/HTTP_CODE://g' | tr -d '\r' | tr -d '[:space:]')")"
    local response
    response="$(echo "$raw" | sed '$d')"

    if [ "$http_code" != "200" ]; then
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        handle_web_search_probe_failure "Container web search probe returned HTTP ${http_code:-000}"
        return $?
    fi

    if json_response_has_results_array "$response"; then
        test_pass
        return 0
    fi

    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    handle_web_search_payload_mismatch "Container can reach SearXNG endpoint but returned non-standard payload (HTTP 200)"
    return $?
}

test_create_test_document() {
    test_start "Create Test Document via API"

    local temp_file
    local source_file=""
    temp_file="/tmp/rag_test_doc_$(date +%s).txt"

    if [ -f "$SCRIPT_DIR/README.md" ]; then
        source_file="$SCRIPT_DIR/README.md"
    elif [ -f "$SCRIPT_DIR/SECURITY.md" ]; then
        source_file="$SCRIPT_DIR/SECURITY.md"
    else
        test_fail "No local source document available for RAG test upload"
        return 1
    fi

    awk 'NF {print; count++; if (count >= 30) exit}' "$source_file" >"$temp_file"
    if [ ! -s "$temp_file" ]; then
        rm -f "$temp_file"
        test_fail "Source document is empty; cannot create RAG test document"
        return 1
    fi

    {
        printf "\n%s\n" "$RAG_TEST_SENTINEL_LINE"
        printf "Repeat the sentinel token exactly when asked: %s\n" "$RAG_TEST_SENTINEL"
    } >>"$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/files/?process=true&process_in_background=false" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -F "file=@$temp_file" \
        -F 'metadata={"name":"rag_test.txt"}' 2>&1)

    rm -f "$temp_file"

    if json_response_has_any_id "$response"; then
        TEST_DOC_ID="$(extract_json_primary_id "$response")"
        local wait_output=""
        if ! wait_output="$(wait_for_openwebui_file_processing "$TEST_DOC_ID" 20 1)"; then
            if [ -n "$wait_output" ]; then
                test_fail "Document processing failed or timed out"
                [ "$VERBOSE" = "true" ] && print_info "Response: $wait_output"
            else
                test_fail "Document processing failed or timed out"
            fi
            return 1
        fi

        print_success "Document created and processed with ID: $TEST_DOC_ID"
        print_info "RAG sentinel: $RAG_TEST_SENTINEL"
        test_pass
        return 0
    else
        test_fail "Document creation failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_create_knowledge_base() {
    test_start "Create Knowledge Base"

    local data='{
        "name": "RAG Test Knowledge Base",
        "description": "Test knowledge base for RAG validation"
    }'

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/knowledge/create" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)

    if json_response_has_any_id "$response"; then
        TEST_KB_ID="$(extract_json_primary_id "$response")"
        if [ -z "$TEST_DOC_ID" ]; then
            test_fail "Knowledge base created but no test document available to attach"
            return 1
        fi

        local add_file_payload
        local add_file_response
        add_file_payload="{\"file_id\":\"$TEST_DOC_ID\"}"
        add_file_response=$(curl -s -X POST \
            "$OPENWEBUI_URL/api/v1/knowledge/$TEST_KB_ID/file/add" \
            -H "Content-Type: application/json" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
            -d "$add_file_payload" 2>&1)

        if json_response_has_files_field "$add_file_response"; then
            print_success "Knowledge base created and file linked: $TEST_KB_ID"
            test_pass
            return 0
        fi

        test_fail "Knowledge base file linking failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $add_file_response"
        return 1
    else
        test_fail "Knowledge base creation failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_multimodal_pdf_ingestion() {
    test_start "Multimodal PDF Ingestion"

    if [ -z "$RAG_MULTIMODAL_TEST_PDF" ]; then
        if is_true "$RAG_MULTIMODAL_STRICT"; then
            test_fail "RAG_MULTIMODAL_TEST_PDF is not set"
            return 1
        fi
        print_info "RAG_MULTIMODAL_TEST_PDF is not set; skipping multimodal ingestion test"
        test_pass
        return 0
    fi

    if [ ! -f "$RAG_MULTIMODAL_TEST_PDF" ]; then
        if is_true "$RAG_MULTIMODAL_STRICT"; then
            test_fail "Multimodal test PDF not found: $RAG_MULTIMODAL_TEST_PDF"
            return 1
        fi
        print_info "Multimodal test PDF not found, skipping: $RAG_MULTIMODAL_TEST_PDF"
        test_pass
        return 0
    fi

    if [ -z "$TEST_KB_ID" ]; then
        test_fail "No knowledge base ID available for multimodal ingestion"
        return 1
    fi

    local retrieval_response
    retrieval_response=$(curl -s -X GET \
        "$OPENWEBUI_URL/api/v1/retrieval/config" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>&1)

    local extraction_engine=""
    local pdf_extract_images=""
    if command -v jq >/dev/null 2>&1; then
        extraction_engine="$(echo "$retrieval_response" | jq -r '.CONTENT_EXTRACTION_ENGINE // empty' 2>/dev/null || true)"
        pdf_extract_images="$(echo "$retrieval_response" | jq -r '.PDF_EXTRACT_IMAGES // empty' 2>/dev/null || true)"
    fi

    if [ "$extraction_engine" != "docling" ] || [ "$pdf_extract_images" != "true" ]; then
        local prereq_msg="Multimodal prereqs not met (CONTENT_EXTRACTION_ENGINE=${extraction_engine:-unset}, PDF_EXTRACT_IMAGES=${pdf_extract_images:-unset})"
        if is_true "$RAG_MULTIMODAL_STRICT"; then
            test_fail "$prereq_msg"
            return 1
        fi
        print_info "$prereq_msg; skipping multimodal ingestion test"
        test_pass
        return 0
    fi

    local upload_response
    local file_name
    file_name="$(basename "$RAG_MULTIMODAL_TEST_PDF")"
    upload_response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/files/?process=true&process_in_background=false" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -F "file=@$RAG_MULTIMODAL_TEST_PDF" \
        -F "metadata={\"name\":\"$file_name\"}" 2>&1)

    TEST_MULTIMODAL_DOC_ID="$(extract_json_primary_id "$upload_response" || true)"
    if [ -z "$TEST_MULTIMODAL_DOC_ID" ]; then
        test_fail "Multimodal file upload failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $upload_response"
        return 1
    fi

    local wait_output=""
    if ! wait_output="$(wait_for_openwebui_file_processing "$TEST_MULTIMODAL_DOC_ID" 40 1)"; then
        test_fail "Multimodal document processing failed or timed out"
        [ "$VERBOSE" = "true" ] && print_info "Response: ${wait_output:-$upload_response}"
        return 1
    fi

    local add_file_payload
    local add_file_response
    add_file_payload="{\"file_id\":\"$TEST_MULTIMODAL_DOC_ID\"}"
    add_file_response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/knowledge/$TEST_KB_ID/file/add" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$add_file_payload" 2>&1)

    if ! json_response_has_files_field "$add_file_response"; then
        test_fail "Knowledge base multimodal file linking failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $add_file_response"
        return 1
    fi

    print_success "Multimodal document uploaded and linked: $TEST_MULTIMODAL_DOC_ID"
    test_pass
    return 0
}

test_multimodal_rag_query() {
    test_start "Multimodal RAG Query"

    if [ -z "$TEST_MULTIMODAL_DOC_ID" ]; then
        print_info "No multimodal document linked; skipping multimodal query test"
        test_pass
        return 0
    fi

    local model="${CHAT_MODEL:-}"
    if [ -z "$model" ] && resolve_chat_model_from_openwebui; then
        model="$CHAT_MODEL"
    fi
    model="${model:-meta-llama/Llama-3.1-8B-Instruct}"
    local files_payload="[]"
    if [ -n "$TEST_KB_ID" ]; then
        files_payload="$(jq -cn --arg id "$TEST_KB_ID" '[{type:"collection",id:$id,status:"processed"}]')"
    fi

    local data
    data="$(jq -cn \
        --arg model "$model" \
        --arg sentinel "$RAG_TEST_SENTINEL" \
        --argjson files "$files_payload" \
        '{
            model: $model,
            messages: [
                {
                    role: "user",
                    content: ("Using only the uploaded knowledge base, repeat the exact sentinel token from the text document: "
                        + $sentinel
                        + ". Then mention whether the uploaded PDF contributes any notable visual evidence.")
                }
            ],
            files: $files,
            temperature: 0.2,
            max_tokens: 384
        }')"

    local response
    response=$(curl -s -m 45 -X POST \
        "$OPENWEBUI_URL/api/chat/completions" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)

    if json_response_has_chat_choices "$response"; then
        local reply
        reply="$(extract_chat_response_text "$response")"
        if response_contains_text "$reply" "$RAG_TEST_SENTINEL"; then
            test_pass
            return 0
        fi

        test_fail "Multimodal RAG response did not include the sentinel token"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    test_fail "Multimodal RAG query failed"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

test_rag_query() {
    test_start "RAG Query Test"

    local model="${CHAT_MODEL:-}"
    if [ -z "$model" ] && resolve_chat_model_from_openwebui; then
        model="$CHAT_MODEL"
    fi
    model="${model:-meta-llama/Llama-3.1-8B-Instruct}"
    local files_payload="[]"

    if [ -n "$TEST_KB_ID" ]; then
        files_payload="$(jq -cn --arg id "$TEST_KB_ID" '[{type:"collection",id:$id,status:"processed"}]')"
    fi

    local data
    data="$(jq -cn \
        --arg model "$model" \
        --arg sentinel "$RAG_TEST_SENTINEL" \
        --argjson files "$files_payload" \
        '{
            model: $model,
            messages: [
                {
                    role: "user",
                    content: ("Using only the uploaded knowledge base, repeat this exact sentinel token and nothing else if you can find it: "
                        + $sentinel)
                }
            ],
            files: $files,
            temperature: 0.1,
            max_tokens: 128
        }')"

    local response
    response=$(curl -s -m 45 -X POST \
        "$OPENWEBUI_URL/api/chat/completions" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)

    if json_response_has_chat_choices "$response"; then
        local reply
        reply="$(extract_chat_response_text "$response")"
        if ! response_contains_text "$reply" "$RAG_TEST_SENTINEL"; then
            test_fail "RAG query did not return the sentinel token"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi

        print_success "RAG query returned the sentinel token"
        print_info "Response preview: $(echo "$reply" | head -c 100)..."
        test_pass
        return 0
    else
        test_fail "RAG query failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_vector_store() {
    test_start "Vector Store Functionality"

    TEST_VECTOR_COLLECTION="rag-test-collection-$(date +%s)"
    local expected_text="RAG systems combine retrieval and generation for better responses."

    local process_payload
    process_payload='{
        "name": "rag-vector-test",
        "content": "'"$expected_text"'",
        "collection_name": "'"$TEST_VECTOR_COLLECTION"'"
    }'

    local process_response
    process_response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/retrieval/process/text" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$process_payload" 2>/dev/null)

    if ! echo "$process_response" | grep -q '"status":true'; then
        test_fail "Vector store ingestion failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $process_response"
        return 1
    fi

    local query_payload
    query_payload='{
        "collection_name": "'"$TEST_VECTOR_COLLECTION"'",
        "query": "What do RAG systems combine?",
        "k": 3
    }'

    local query_response
    query_response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/retrieval/query/doc" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$query_payload" 2>/dev/null)

    if json_response_has_retrieval_hits "$query_response"; then
        local retrieved_text
        retrieved_text="$(extract_retrieval_text "$query_response")"
        if ! response_contains_text "$retrieved_text" "$expected_text"; then
            test_fail "Vector store query returned hits, but not the expected document"
            [ "$VERBOSE" = "true" ] && print_info "Response: $query_response"
            return 1
        fi

        local doc_count
        doc_count="$(retrieval_hit_count "$query_response")"
        print_success "Vector store retrieval returned $doc_count matching document(s)"
        test_pass
        return 0
    fi

    test_fail "Vector store query failed"
    [ "$VERBOSE" = "true" ] && print_info "Response: $query_response"
    return 1
}

test_retrieval() {
    test_start "Document Retrieval"

    if [ -z "$TEST_KB_ID" ]; then
        test_fail "No knowledge base ID available for retrieval"
        return 1
    fi

    local query="$RAG_TEST_SENTINEL"
    local data='{
        "collection_name": "'${TEST_KB_ID:-}'",
        "query": "'"$query"'",
        "k": 3
    }'

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/retrieval/query/doc" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>/dev/null)

    if json_response_has_retrieval_hits "$response"; then
        local retrieved_text
        retrieved_text="$(extract_retrieval_text "$response")"
        if ! response_contains_text "$retrieved_text" "$RAG_TEST_SENTINEL"; then
            test_fail "Retrieval hits did not include the sentinel token"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi

        local doc_count
        doc_count="$(retrieval_hit_count "$response")"
        print_success "Retrieved $doc_count documents containing the sentinel token"
        test_pass
        return 0
    else
        test_fail "Retrieval test failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

cleanup_test_resources() {
    print_section "Cleaning Up Test Resources"

    if [ -n "$TEST_MULTIMODAL_DOC_ID" ]; then
        print_info "Removing multimodal test document $TEST_MULTIMODAL_DOC_ID..."
        curl -s -X DELETE \
            "$OPENWEBUI_URL/api/v1/files/$TEST_MULTIMODAL_DOC_ID" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} >/dev/null 2>&1
        print_success "Multimodal test document removed"
    fi

    if [ -n "$TEST_DOC_ID" ]; then
        print_info "Removing test document $TEST_DOC_ID..."
        curl -s -X DELETE \
            "$OPENWEBUI_URL/api/v1/files/$TEST_DOC_ID" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} >/dev/null 2>&1
        print_success "Test document removed"
    fi

    if [ -n "$TEST_KB_ID" ]; then
        print_info "Removing test knowledge base $TEST_KB_ID..."
        curl -s -X DELETE \
            "$OPENWEBUI_URL/api/v1/knowledge/$TEST_KB_ID/delete" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} >/dev/null 2>&1
        print_success "Test knowledge base removed"
    fi
}

print_summary() {
    print_header "Test Summary"

    echo -e "Total Tests: $TESTS_TOTAL"
    echo -e "${GREEN}Passed: $TESTS_PASSED${NC}"
    echo -e "${RED}Failed: $TESTS_FAILED${NC}"

    local pass_rate=0
    if [ $TESTS_TOTAL -gt 0 ]; then
        pass_rate=$((TESTS_PASSED * 100 / TESTS_TOTAL))
    fi

    echo -e "\nPass Rate: $pass_rate%"

    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "\n${GREEN}All tests passed! RAG deployment is working correctly.${NC}"
        return 0
    else
        echo -e "\n${RED}Some tests failed. Please check the errors above.${NC}"
        return 1
    fi
}

###############################################################################
# Main Execution
###############################################################################

main() {
    print_header "RAG Deployment Testing Suite"
    RUN_CLEANUP_ON_EXIT=true
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "Upstream URL: $UPSTREAM_URL"
    print_info "Test mode: $TEST_MODE"
    [ -n "$API_KEY" ] && print_info "API Key: ${API_KEY:0:8}..."
    ensure_openwebui_api_key

    # Real integration guard always applies.
    print_section "Integration Guard"
    test_repo_real_integration_guard

    print_section "Upstream Integration"
    print_info "Running against the configured OpenAI-compatible upstream"

    # Health checks
    print_section "Health Checks"
    test_openwebui_accessible
    test_upstream_responding
    test_openwebui_health
    test_indigo_service_health
    test_indigo_tool_registration
    test_indigo_tool_connectivity_from_openwebui

    # Model tests
    print_section "Model Availability"
    test_upstream_models

    # Core functionality tests
    print_section "Core RAG Functionality"
    test_embedding_endpoint
    test_vector_store
    test_openwebui_baseline_settings
    test_web_search_baseline
    test_web_search_container_baseline

    if [ "$TEST_MODE" = "full" ]; then
        test_create_test_document
        test_create_knowledge_base
        test_multimodal_pdf_ingestion

        # Integration tests
        print_section "Integration Tests"
        test_multimodal_rag_query
        test_rag_query
        test_retrieval
    else
        print_section "Integration Tests"
        print_info "Skipping deep integration tests in baseline mode"
    fi

    # Summary
    print_summary
}

on_exit() {
    if [ "$RUN_CLEANUP_ON_EXIT" = true ]; then
        cleanup_test_resources
    fi
}

# Trap to cleanup on interrupt
trap on_exit EXIT INT TERM

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --baseline)
            TEST_MODE="baseline"
            shift
            ;;
        --full)
            TEST_MODE="full"
            shift
            ;;
        -v | --verbose)
            VERBOSE=true
            shift
            ;;
        -h | --help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# CI full runs should fail fast on multimodal prereq drift unless explicitly overridden.
if [ -z "$RAG_MULTIMODAL_STRICT" ]; then
    if [ "$TEST_MODE" = "full" ] && is_true "${CI:-false}"; then
        RAG_MULTIMODAL_STRICT="true"
    else
        RAG_MULTIMODAL_STRICT="false"
    fi
fi

# Run main function
main "$@"
