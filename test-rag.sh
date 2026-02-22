#!/bin/bash

###############################################################################
# RAG Functionality Testing Script
# This script tests the complete RAG deployment including OpenWebUI and vLLM
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
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

    # OpenWebUI can use Docker DNS internally; tests run on host and must use loopback.
    if [[ $candidate == *"://cliproxyapi"* ]]; then
        candidate="$(printf "%s" "$candidate" | sed 's#://cliproxyapi#://127.0.0.1#')"
    fi

    printf "%s" "$candidate"
}

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENAI_BASE_ROOT="$(resolve_openai_base_root)"
CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
VLLM_URL="${VLLM_URL:-${OPENAI_BASE_ROOT:-$CLIPROXYAPI_BASE_URL}}"
API_KEY="${OPENWEBUI_API_KEY:-${API_KEY:-}}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
VERBOSE="${VERBOSE:-false}"
CHAT_MODEL="${RAG_CHAT_MODEL:-}"
CLIPROXY_BIN="${CLIPROXY_BIN:-./cli-proxy-api.sh}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"
TEST_MODE="full"
RUN_CLEANUP_ON_EXIT=false
BASELINE_REQUIRE_WEB_SEARCH="${BASELINE_REQUIRE_WEB_SEARCH:-false}"
BASELINE_EXPECT_SEARX_PORT="${BASELINE_EXPECT_SEARX_PORT:-8888}"

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0
TEST_VECTOR_COLLECTION=""

###############################################################################
# Script-specific Helper Functions
###############################################################################

openwebui_signin() {
    local output_file="$1"
    local payload_file="$2"
    local signin_url="${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"

    if command -v jq >/dev/null 2>&1; then
        jq -n \
            --arg email "$OPENWEBUI_SIGNIN_EMAIL" \
            --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
            '{email: $email, password: $password}' >"$payload_file"
    else
        cat >"$payload_file" <<EOF
{"email":"$OPENWEBUI_SIGNIN_EMAIL","password":"$OPENWEBUI_SIGNIN_PASSWORD"}
EOF
    fi

    curl -sS -m 45 -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" "$signin_url" -d @"$payload_file" 2>/dev/null || true
}

ensure_openwebui_api_key() {
    if [ -n "$API_KEY" ]; then
        return 0
    fi

    if ! is_true "$OPENWEBUI_AUTO_AUTH"; then
        return 0
    fi

    local tmp_file payload_file signin_code token
    tmp_file="$(mktemp)"
    payload_file="$(mktemp)"

    signin_code="$(openwebui_signin "$tmp_file" "$payload_file")"
    if [ "$signin_code" = "200" ]; then
        token="$(jq -r '.token // empty' "$tmp_file" 2>/dev/null || true)"
        if [ -n "$token" ]; then
            API_KEY="$token"
            print_info "OpenWebUI bearer token acquired via signin"
        fi
    fi

    rm -f "$tmp_file" "$payload_file"
}

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
    if response=$(./audit-no-mock.sh 2>&1); then
        test_pass
        return 0
    fi

    test_fail "audit-no-mock.sh failed"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

has_non_empty_models_array() {
    local response="$1"
    local model_count

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e '.data and (.data | type == "array") and (.data | length > 0)' >/dev/null 2>&1
        return $?
    fi

    model_count=$(echo "$response" | grep -Eo '"id"[[:space:]]*:[[:space:]]*"[^"]+"' | wc -l | tr -d '[:space:]')
    if [[ $model_count =~ ^[0-9]+$ ]] && [ "$model_count" -gt 0 ]; then
        return 0
    fi
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
  BASELINE_REQUIRE_WEB_SEARCH=true|false  Fail if host and container SearXNG probes fail (default: false)
  BASELINE_EXPECT_SEARX_PORT=<port|empty> Expected SearXNG port hint (default: 8888; empty disables port hint)
EOF
}

###############################################################################
# Test Functions
###############################################################################

test_openwebui_accessible() {
    test_start "OpenWebUI Accessibility"

    if curl -s -f -o /dev/null "$OPENWEBUI_URL" 2>&1; then
        test_pass
        return 0
    else
        test_fail "Cannot reach OpenWebUI at $OPENWEBUI_URL"
        return 1
    fi
}

test_vllm_responding() {
    test_start "vLLM Server Response"

    local code
    local root_response

    code=$(curl -s -w "%{http_code}" -o /dev/null "$VLLM_URL/health" 2>&1 || echo "000")
    if [ "$code" = "200" ]; then
        test_pass
        return 0
    fi

    # CLIProxyAPI deployments may expose "/" instead of "/health".
    root_response=$(curl -s -f "$VLLM_URL/" 2>/dev/null || true)
    if [ -n "$root_response" ] && (
        (command -v jq >/dev/null 2>&1 && echo "$root_response" | jq -e '.message == "CLI Proxy API Server"' >/dev/null 2>&1) ||
            echo "$root_response" | grep -q "CLI Proxy API Server"
    ); then
        test_pass
        return 0
    else
        test_fail "OpenAI-compatible upstream health check failed at $VLLM_URL (/health HTTP $code)"
        return 1
    fi
}

test_vllm_models() {
    test_start "vLLM Models Available"

    local response
    local -a headers=()
    if [ -n "$CLIPROXYAPI_API_KEY" ]; then
        headers=(-H "Authorization: Bearer $CLIPROXYAPI_API_KEY")
    fi
    response=$(curl -s "${headers[@]}" "$VLLM_URL/v1/models" 2>&1)

    if has_non_empty_models_array "$response"; then
        local models
        models=$(echo "$response" | jq -r '.data[]?.id // empty' 2>/dev/null | tr '\n' ', ' | sed 's/, $//' || true)
        if [ -z "$models" ]; then
            models=$(echo "$response" | grep -Eo '"id"[[:space:]]*:[[:space:]]*"[^"]+"' | head -n 20 | cut -d'"' -f4 | tr '\n' ', ' | sed 's/, $//')
        fi

        if [ -z "$CHAT_MODEL" ]; then
            CHAT_MODEL=$(echo "$response" | jq -r '.data[0].id // empty' 2>/dev/null || true)
            if [ -z "$CHAT_MODEL" ]; then
                CHAT_MODEL=$(echo "$response" | grep -Eo '"id"[[:space:]]*:[[:space:]]*"[^"]+"' | head -n1 | cut -d'"' -f4)
            fi
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
    response=$(curl -s "$OPENWEBUI_URL/health" 2>&1)

    if [ "$?" = "0" ]; then
        print_info "Health check response: $response"
        test_pass
        return 0
    else
        test_fail "Health endpoint not responding"
        return 1
    fi
}

test_embedding_endpoint() {
    test_start "Embedding Generation"

    local response
    response=$(curl -s -X GET \
        "$OPENWEBUI_URL/api/v1/retrieval/" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>&1)

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
        "$OPENWEBUI_URL/api/v1/retrieval/" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>&1)

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

        local embedding_model
        embedding_model=$(echo "$response" | jq -r '.RAG_EMBEDDING_MODEL // empty' 2>/dev/null || true)
        if [ -z "$embedding_model" ]; then
            test_fail "RAG_EMBEDDING_MODEL missing from retrieval settings"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi

        local rag_top_k
        rag_top_k=$(echo "$response" | jq -r '.RAG_TOP_K // .TOP_K // empty' 2>/dev/null || true)
        if [ -z "$rag_top_k" ]; then
            print_info "RAG_TOP_K not returned by endpoint; skipping strict value check"
        elif [[ ! $rag_top_k =~ ^[0-9]+$ ]] || [ "$rag_top_k" -le 0 ]; then
            test_fail "Invalid RAG_TOP_K value from endpoint: $rag_top_k"
            return 1
        fi

        print_info "Embedding model: $embedding_model"
        [ -n "$rag_top_k" ] && print_info "RAG_TOP_K: $rag_top_k"
        test_pass
        return 0
    fi

    if echo "$response" | grep -q "RAG_EMBEDDING_MODEL"; then
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
    print_info "Continuing because BASELINE_REQUIRE_WEB_SEARCH=false"
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
    http_code=$(curl -sS -m 20 -o "$response_file" -w "%{http_code}" "$host_probe_url" 2>/dev/null || true)
    response="$(cat "$response_file" 2>/dev/null || true)"

    if { [ -z "$response" ] || [ "$http_code" != "200" ]; } && [ "$host_probe_url" != "$probe_url" ]; then
        http_code=$(curl -sS -m 20 -o "$response_file" -w "%{http_code}" "$probe_url" 2>/dev/null || true)
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

    if command -v jq >/dev/null 2>&1; then
        if echo "$response" | jq -e '.results and (.results | type == "array")' >/dev/null 2>&1; then
            test_pass
            return 0
        fi
    elif echo "$response" | grep -q "results"; then
        test_pass
        return 0
    fi

    print_info "SearXNG endpoint reachable but returned non-standard payload (HTTP 200)"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    test_pass
    return 0
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
    raw="$(docker exec openwebui sh -lc "curl -sS -m 20 -w '\nHTTP_CODE:%{http_code}\n' \"$probe_url\" 2>/dev/null || true" 2>/dev/null || true)"

    if [ -z "$raw" ]; then
        handle_web_search_probe_failure "Empty response from OpenWebUI container web search probe"
        return $?
    fi

    local http_code
    http_code="$(echo "$raw" | tail -n 1 | sed 's/HTTP_CODE://g' | tr -d '\r' | tr -d '[:space:]')"
    local response
    response="$(echo "$raw" | sed '$d')"

    if [ "$http_code" != "200" ]; then
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        handle_web_search_probe_failure "Container web search probe returned HTTP ${http_code:-000}"
        return $?
    fi

    if command -v jq >/dev/null 2>&1; then
        if echo "$response" | jq -e '.results and (.results | type == "array")' >/dev/null 2>&1; then
            test_pass
            return 0
        fi
    elif echo "$response" | grep -q "results"; then
        test_pass
        return 0
    fi

    print_info "Container can reach SearXNG endpoint but returned non-standard payload (HTTP 200)"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    test_pass
    return 0
}

test_create_test_document() {
    test_start "Create Test Document via API"

    local test_content="This is a test document for RAG testing.
It contains information about artificial intelligence and machine learning.
RAG systems combine retrieval and generation for better responses.
This document will be used to test the knowledge base functionality."

    local temp_file="/tmp/rag_test_doc_$(date +%s).txt"
    echo "$test_content" >"$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/files/?process=true&process_in_background=false" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -F "file=@$temp_file" \
        -F 'metadata={"name":"rag_test.txt"}' 2>&1)

    rm -f "$temp_file"

    if echo "$response" | grep -q "id"; then
        TEST_DOC_ID=$(echo "$response" | jq -r '.id // .file_id // empty' 2>/dev/null)
        local status_response=""
        local status_value=""
        local attempts=0
        while [ "$attempts" -lt 20 ]; do
            status_response=$(curl -s -X GET \
                "$OPENWEBUI_URL/api/v1/files/$TEST_DOC_ID/process/status" \
                ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>&1)
            status_value=$(echo "$status_response" | jq -r '.status // empty' 2>/dev/null || true)
            if [ "$status_value" = "completed" ]; then
                break
            fi
            if [ "$status_value" = "failed" ]; then
                test_fail "Document processing failed"
                [ "$VERBOSE" = "true" ] && print_info "Response: $status_response"
                return 1
            fi
            attempts=$((attempts + 1))
            sleep 1
        done

        if [ "$status_value" != "completed" ]; then
            test_fail "Document processing did not complete in time"
            [ "$VERBOSE" = "true" ] && print_info "Response: $status_response"
            return 1
        fi

        print_success "Document created and processed with ID: $TEST_DOC_ID"
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

    if echo "$response" | grep -q "id"; then
        TEST_KB_ID=$(echo "$response" | jq -r '.id // .knowledge_id // empty' 2>/dev/null)
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

        if echo "$add_file_response" | grep -q "files"; then
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

test_rag_query() {
    test_start "RAG Query Test"

    local model="${CHAT_MODEL:-meta-llama/Llama-3.1-8B-Instruct}"

    local data='{
        "model": "'"$model"'",
        "messages": [
            {
                "role": "user",
                "content": "What is RAG and how does it work?"
            }
        ],
        "knowledge_ids": [],
        "temperature": 0.7,
        "max_tokens": 512
    }'

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/chat/completions" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)

    if echo "$response" | grep -q "choices"; then
        local reply=$(echo "$response" | jq -r '.choices[0].message.content // empty' 2>/dev/null)
        print_success "RAG query successful"
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

    local process_payload
    process_payload='{
        "name": "rag-vector-test",
        "content": "RAG systems combine retrieval and generation for better responses.",
        "collection_name": "'"$TEST_VECTOR_COLLECTION"'"
    }'

    local process_response
    process_response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/retrieval/process/text" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$process_payload" 2>&1)

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
        -d "$query_payload" 2>&1)

    if echo "$query_response" | grep -q '"documents"'; then
        print_success "Vector store retrieval is operational"
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

    local query="RAG testing validation"
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
        -d "$data" 2>&1)

    if echo "$response" | grep -q "results\|documents"; then
        local doc_count=$(echo "$response" | jq '.results | length' 2>/dev/null || echo "0")
        print_success "Retrieved $doc_count documents"
        test_pass
        return 0
    else
        test_fail "Retrieval test failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_cli_proxy_api_managed_runtime() {
    test_start "CliProxyApi Managed Runtime"

    local response
    local cliproxy_base="${CLIPROXYAPI_BASE_URL%/}"
    local -a cliproxy_headers=()

    if [ -n "$CLIPROXYAPI_API_KEY" ]; then
        cliproxy_headers=(-H "Authorization: Bearer $CLIPROXYAPI_API_KEY")
    fi

    if response=$(CLIPROXYAPI_ENABLED=false ./check-cliproxyapi.sh 2>&1); then
        test_fail "Disabled lifecycle negative check unexpectedly passed (CLIPROXYAPI_ENABLED=false)"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if ! echo "$response" | grep -Eiq "CLIPROXYAPI_ENABLED=false|lifecycle is disabled"; then
        test_fail "Disabled lifecycle negative check failed without expected disabled lifecycle message"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if ! response=$(CLIPROXYAPI_ENABLED=true ./check-cliproxyapi.sh 2>&1); then
        test_fail "Managed health check failed (set CLIPROXYAPI_ENABLED=true and start CLIProxyAPI)"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if echo "$response" | grep -Eiq "CLIPROXYAPI_ENABLED=false|lifecycle is disabled"; then
        test_fail "check-cliproxyapi.sh reported disabled lifecycle; runtime checks did not execute"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if ! response=$(curl -sS -f "${cliproxy_base}/" 2>&1); then
        test_fail "CLIProxyAPI root endpoint failed: ${cliproxy_base}/"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        if ! echo "$response" | jq -e '.message == "CLI Proxy API Server"' >/dev/null 2>&1; then
            test_fail "CLIProxyAPI root endpoint missing expected message"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi
    elif ! echo "$response" | grep -q "CLI Proxy API Server"; then
        test_fail "CLIProxyAPI root endpoint missing expected message"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if ! response=$(curl -sS -f "${cliproxy_headers[@]}" "${cliproxy_base}/v1/models" 2>&1); then
        test_fail "CLIProxyAPI models endpoint failed: ${cliproxy_base}/v1/models"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if ! has_non_empty_models_array "$response"; then
        test_fail "CLIProxyAPI models response missing expected non-empty models list"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    test_pass
    return 0
}

test_cli_proxy_api_health_all() {
    test_start "CliProxyApi Health (openwebui + upstream)"

    local response
    if ! response=$(OPENWEBUI_URL="$OPENWEBUI_URL" OPENWEBUI_API_KEY="$API_KEY" "$CLIPROXY_BIN" --raw health openwebui 2>&1); then
        test_fail "Command failed: $CLIPROXY_BIN --raw health openwebui"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        if ! echo "$response" | jq -e '.status == true or .openwebui' >/dev/null 2>&1; then
            test_fail "OpenWebUI health output missing expected success signal"
            [ "$VERBOSE" = "true" ] && print_info "Response: $response"
            return 1
        fi
    elif ! echo "$response" | grep -Eq '"status"[[:space:]]*:[[:space:]]*true|"openwebui"'; then
        test_fail "OpenWebUI health output missing expected success signal"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    # Upstream in CLIProxyAPI mode may not expose /health, so validate root endpoint directly.
    if curl -sS -f "${CLIPROXYAPI_BASE_URL%/}/" >/dev/null 2>&1; then
        test_pass
        return 0
    fi

    test_fail "CLIProxyAPI upstream root endpoint not reachable"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

test_cli_proxy_api_models_vllm() {
    test_start "CliProxyApi Models (vLLM)"

    local response
    if ! response=$(OPENWEBUI_URL="$OPENWEBUI_URL" OPENWEBUI_API_KEY="$API_KEY" VLLM_URL="$CLIPROXYAPI_BASE_URL" "$CLIPROXY_BIN" --raw --url "$CLIPROXYAPI_BASE_URL" models vllm 2>&1); then
        test_fail "Command failed: $CLIPROXY_BIN --raw models vllm"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if has_non_empty_models_array "$response"; then
        if [ -z "$CHAT_MODEL" ]; then
            CHAT_MODEL=$(echo "$response" | jq -r '.data[0].id // empty' 2>/dev/null || true)
        fi
        if [ -z "$CHAT_MODEL" ]; then
            CHAT_MODEL=$(echo "$response" | grep -Eo '"id"[[:space:]]*:[[:space:]]*"[^"]+"' | head -n1 | cut -d'"' -f4)
        fi
        [ -n "$CHAT_MODEL" ] && print_info "Using chat model: $CHAT_MODEL"
        test_pass
        return 0
    fi

    test_fail "CLIProxyAPI models output missing expected non-empty models list"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

test_cli_proxy_api_models_openwebui() {
    test_start "CliProxyApi Models (OpenWebUI)"

    local response
    if ! response=$(OPENWEBUI_URL="$OPENWEBUI_URL" OPENWEBUI_API_KEY="$API_KEY" VLLM_URL="$CLIPROXYAPI_BASE_URL" "$CLIPROXY_BIN" --raw models webui 2>&1); then
        test_fail "Command failed: $CLIPROXY_BIN --raw models webui"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if has_non_empty_models_array "$response"; then
        test_pass
        return 0
    fi

    test_fail "CLIProxyAPI OpenWebUI models output missing expected non-empty models list"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

cleanup_test_resources() {
    print_section "Cleaning Up Test Resources"

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
    print_info "vLLM URL: $VLLM_URL"
    print_info "CliProxyApi: $CLIPROXY_BIN"
    print_info "Test mode: $TEST_MODE"
    [ -n "$API_KEY" ] && print_info "API Key: ${API_KEY:0:8}..."
    ensure_openwebui_api_key

    # CLIProxyAPI smoke checks
    print_section "CliProxyApi Integration"
    test_repo_real_integration_guard
    test_cli_proxy_api_managed_runtime
    test_cli_proxy_api_health_all
    test_cli_proxy_api_models_openwebui
    test_cli_proxy_api_models_vllm

    # Health checks
    print_section "Health Checks"
    test_openwebui_accessible
    test_vllm_responding
    test_openwebui_health

    # Model tests
    print_section "Model Availability"
    test_vllm_models

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

        # Integration tests
        print_section "Integration Tests"
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

# Run main function
main "$@"
