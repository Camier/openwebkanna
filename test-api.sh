#!/bin/bash

###############################################################################
# API Endpoint Testing Script
# This script tests all major API endpoints for the RAG deployment
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
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
VERBOSE="${VERBOSE:-false}"
OUTPUT_DIR="/tmp/openwebui_api_tests_$(date +%Y%m%d_%H%M%S)"
TEST_MODE="full"
RUN_CLEANUP_ON_EXIT=false

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0

# Track created resources
declare -a CREATED_FILES=()
declare -a CREATED_KNOWLEDGE_BASES=()
SELECTED_MODEL=""

###############################################################################
# Helper Functions
###############################################################################

print_endpoint() {
    echo -e "${CYAN}>> $1${NC}"
}

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
    print_endpoint "$TEST_NAME"
    TEST_START_TIME=$(date +%s)
}

test_pass() {
    ((TESTS_PASSED += 1))
    local duration=$(($(date +%s) - TEST_START_TIME))
    print_success "$TEST_NAME (${duration}s)"
}

test_fail() {
    ((TESTS_FAILED += 1))
    local duration=$(($(date +%s) - TEST_START_TIME))
    print_error "$TEST_NAME failed (${duration}s): $1"
}

save_response() {
    local test_name="$1"
    local response="$2"
    mkdir -p "$OUTPUT_DIR"
    echo "$response" >"$OUTPUT_DIR/${test_name}.json"
    if [ "$VERBOSE" = "true" ]; then
        print_info "Response saved to $OUTPUT_DIR/${test_name}.json"
    fi
    return 0
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

make_request() {
    local method="$1"
    local endpoint="$2"
    local data="$3"
    local url="${OPENWEBUI_URL}${endpoint}"
    local -a curl_args=("-sS" "-m" "45" "-X" "$method" "$url" "-H" "Content-Type: application/json")

    if [ -n "$API_KEY" ]; then
        curl_args+=("-H" "Authorization: Bearer $API_KEY")
    fi

    if [ "$method" = "POST" ] || [ "$method" = "PUT" ]; then
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

make_request_with_status() {
    local method="$1"
    local endpoint="$2"
    local data="$3"
    local response_var="$4"
    local status_var="$5"
    local url="${OPENWEBUI_URL}${endpoint}"
    local response_file=""
    local response_body=""
    local status_code=""
    local -a curl_args=("-sS" "-m" "45" "-X" "$method" "$url" "-H" "Content-Type: application/json")

    if [ -n "$API_KEY" ]; then
        curl_args+=("-H" "Authorization: Bearer $API_KEY")
    fi

    if [ "$method" = "POST" ] || [ "$method" = "PUT" ]; then
        curl_args+=("-d" "$data")
    fi

    response_file="$(mktemp)"
    status_code="$(curl "${curl_args[@]}" -o "$response_file" -w "%{http_code}" 2>/dev/null || true)"
    response_body="$(cat "$response_file" 2>/dev/null || true)"
    rm -f "$response_file"

    printf -v "$response_var" "%s" "$response_body"
    printf -v "$status_var" "%s" "$status_code"
}

show_help() {
    cat <<'EOF'
Usage: ./test-api.sh [OPTIONS]

Options:
  --baseline         Run fast baseline API checks
  --full             Run full API regression suite (default)
  -v, --verbose      Enable verbose output
  -u, --url URL      OpenWebUI URL override
  -k, --key KEY      OpenWebUI API key
  -h, --help         Show this help message
EOF
}

###############################################################################
# Test Functions
###############################################################################

# Health Check Tests
test_health_check() {
    test_start "Health Check"

    local response
    local http_code
    make_request_with_status "GET" "/health" "" response http_code
    save_response "health_check" "$response"

    if [ "$http_code" = "200" ] && [ -n "$response" ]; then
        test_pass
        return 0
    else
        test_fail "Health endpoint returned HTTP ${http_code:-000}"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_root_endpoint() {
    test_start "Root Endpoint"

    local response
    response=$(curl -s -I "$OPENWEBUI_URL/" 2>&1 | head -n 1)
    save_response "root_endpoint" "$response"

    if echo "$response" | grep -q "200"; then
        test_pass
        return 0
    else
        test_fail "Root endpoint returned unexpected response"
        return 1
    fi
}

# Authentication Tests
test_api_key_validation() {
    test_start "API Key Validation"

    local response
    local http_code
    make_request_with_status "GET" "/api/auth/status" "" response http_code
    save_response "api_key_validation" "$response"

    if [ "$http_code" = "200" ]; then
        test_pass
        return 0
    elif [ "$http_code" = "401" ] || [ "$http_code" = "403" ]; then
        if is_true "${WEBUI_AUTH:-false}"; then
            test_fail "Authentication is enabled but API key is unauthorized (HTTP $http_code)"
            return 1
        fi
        print_info "API key validation skipped (auth appears disabled; endpoint returned HTTP $http_code)"
        test_pass
        return 0
    else
        test_fail "Auth status endpoint returned HTTP ${http_code:-000}"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

# CliProxyApi Smoke Tests (managed runtime + client wrapper)
test_cli_proxy_api_smoke() {
    local response
    local cliproxy_base="${CLIPROXYAPI_BASE_URL%/}"
    local -a cliproxy_headers=()

    if [ -n "$CLIPROXYAPI_API_KEY" ]; then
        cliproxy_headers=(-H "Authorization: Bearer $CLIPROXYAPI_API_KEY")
    fi

    test_start "CliProxyApi Disabled Lifecycle Negative Check"
    if response=$(CLIPROXYAPI_ENABLED=false ./check-cliproxyapi.sh 2>&1); then
        save_response "cliproxyapi_disabled_lifecycle_check" "$response"
        test_fail "check-cliproxyapi.sh unexpectedly passed with CLIPROXYAPI_ENABLED=false"
    else
        save_response "cliproxyapi_disabled_lifecycle_check" "$response"
        if echo "$response" | grep -Eiq "CLIPROXYAPI_ENABLED=false|lifecycle is disabled"; then
            test_pass
        else
            test_fail "Expected disabled lifecycle failure message was not found"
        fi
    fi

    test_start "CliProxyApi Managed Health Check"
    if response=$(CLIPROXYAPI_ENABLED=true CLIPROXYAPI_CHECK_CHAT_COMPLETION=false ./check-cliproxyapi.sh 2>&1); then
        save_response "cliproxyapi_managed_health" "$response"
        if echo "$response" | grep -Eiq "CLIPROXYAPI_ENABLED=false|lifecycle is disabled"; then
            test_fail "check-cliproxyapi.sh reported disabled lifecycle; real runtime check did not run"
        else
            test_pass
        fi
    else
        save_response "cliproxyapi_managed_health" "$response"
        test_fail "CLIProxyAPI managed health check failed (set CLIPROXYAPI_ENABLED=true and start service)"
    fi

    test_start "CliProxyApi Root Endpoint"
    if response=$(curl -sS -f "${cliproxy_base}/" 2>&1); then
        save_response "cliproxyapi_root" "$response"
        if command -v jq >/dev/null 2>&1; then
            if echo "$response" | jq -e '.message == "CLI Proxy API Server"' >/dev/null 2>&1; then
                test_pass
            else
                test_fail "CLIProxyAPI root response missing expected message"
            fi
        elif echo "$response" | grep -q "CLI Proxy API Server"; then
            test_pass
        else
            test_fail "CLIProxyAPI root response missing expected message"
        fi
    else
        save_response "cliproxyapi_root" "$response"
        test_fail "CLIProxyAPI root endpoint failed: ${cliproxy_base}/"
    fi

    test_start "CliProxyApi Models Endpoint"
    if response=$(curl -sS -f "${cliproxy_headers[@]}" "${cliproxy_base}/v1/models" 2>&1); then
        save_response "cliproxyapi_models" "$response"
        if has_non_empty_models_array "$response"; then
            test_pass
        else
            test_fail "CLIProxyAPI /v1/models output missing expected non-empty models list"
        fi
    else
        save_response "cliproxyapi_models" "$response"
        test_fail "CLIProxyAPI models endpoint failed: ${cliproxy_base}/v1/models"
    fi

    test_start "CliProxyApi Wrapper Help"
    if response=$(./cli-proxy-api.sh --help 2>&1); then
        save_response "cli_proxy_api_help" "$response"
        if echo "$response" | grep -Fq "models [openwebui|webui|vllm|all]" &&
            echo "$response" | grep -Fq "chat --model <id> --message <text>"; then
            test_pass
        else
            test_fail "Help output missing expected command contract"
        fi
    else
        save_response "cli_proxy_api_help" "$response"
        test_fail "Command failed: ./cli-proxy-api.sh --help"
    fi

    test_start "CliProxyApi Wrapper URL Override (models vllm)"
    if response=$(./cli-proxy-api.sh --raw --url "$cliproxy_base" models vllm 2>&1); then
        save_response "cli_proxy_api_models_vllm_override" "$response"
        if has_non_empty_models_array "$response"; then
            test_pass
        else
            test_fail "Wrapper override output missing expected non-empty models list"
        fi
    else
        save_response "cli_proxy_api_models_vllm_override" "$response"
        test_fail "Command failed: ./cli-proxy-api.sh --raw --url \"$cliproxy_base\" models vllm"
    fi

    return 0
}

# Model Tests
test_models_list() {
    test_start "Models List"

    local response
    response=$(make_request "GET" "/api/models" "" "" 2>&1)
    save_response "models_list" "$response"

    if has_non_empty_models_array "$response"; then
        local model_count
        model_count=$(echo "$response" | jq '.data | length' 2>/dev/null || echo "0")
        SELECTED_MODEL=$(echo "$response" | jq -r '.data[0].id // empty' 2>/dev/null || true)
        print_success "Found $model_count models"
        test_pass
        return 0
    else
        test_fail "Could not retrieve models list"
        return 1
    fi
}

test_model_info() {
    test_start "Model Selection"

    local response="{\"selected_model\":\"${SELECTED_MODEL}\"}"
    save_response "model_info" "$response"

    if [ -n "$SELECTED_MODEL" ]; then
        test_pass
        return 0
    else
        test_fail "No model selected from /api/models"
        return 1
    fi
}

# File Upload Tests
test_file_upload_text() {
    test_start "File Upload (Text)"

    local temp_file="/tmp/test_upload_$(date +%s).txt"
    echo "This is a test document for API validation." >"$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/files/?process=true&process_in_background=false" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -F "file=@$temp_file" \
        -F 'metadata={"name":"test_doc.txt"}' 2>&1)
    save_response "file_upload_text" "$response"

    rm -f "$temp_file"

    if echo "$response" | grep -q "id"; then
        local file_id=$(echo "$response" | jq -r '.id // .file_id // empty' 2>/dev/null)
        CREATED_FILES+=("$file_id")
        print_success "File uploaded with ID: $file_id"
        test_pass
        return 0
    else
        test_fail "File upload failed"
        return 1
    fi
}

test_file_upload_pdf() {
    test_start "File Upload (PDF)"

    # Create a minimal PDF
    local temp_file="/tmp/test_upload_$(date +%s).pdf"
    echo "%PDF-1.4
1 0 obj
<<
/Type /Catalog
/Pages 2 0 R
>>
endobj
2 0 obj
<<
/Type /Pages
/Kids [3 0 R]
/Count 1
>>
endobj
3 0 obj
<<
/Type /Page
/Parent 2 0 R
/MediaBox [0 0 612 792]
>>
endobj
xref
0 4
0000000000 65535 f
0000000009 00000 n
0000000058 00000 n
0000000115 00000 n
trailer
<<
/Size 4
/Root 1 0 R
>>
startxref
190
%%EOF" >"$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/files/?process=true&process_in_background=false" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -F "file=@$temp_file" \
        -F 'metadata={"name":"test_doc.pdf"}' 2>&1)
    save_response "file_upload_pdf" "$response"

    rm -f "$temp_file"

    if echo "$response" | grep -q "id"; then
        local file_id=$(echo "$response" | jq -r '.id // .file_id // empty' 2>/dev/null)
        CREATED_FILES+=("$file_id")
        print_success "PDF uploaded with ID: $file_id"
        test_pass
        return 0
    else
        test_fail "PDF upload failed"
        return 1
    fi
}

test_file_list() {
    test_start "File List"

    local response
    response=$(make_request "GET" "/api/v1/files/" "" "" 2>&1)
    save_response "file_list" "$response"

    if echo "$response" | grep -q "id\|filename"; then
        local file_count
        file_count=$(echo "$response" | jq 'length' 2>/dev/null || echo "0")
        print_success "Found $file_count files"
        test_pass
        return 0
    else
        test_fail "Could not retrieve file list"
        return 1
    fi
}

# Knowledge Base Tests
test_knowledge_create() {
    test_start "Knowledge Base Creation"

    local data='{
        "name": "API Test Knowledge Base",
        "description": "Knowledge base created during API testing"
    }'

    local response
    response=$(make_request "POST" "/api/v1/knowledge/create" "$data" "" 2>&1)
    save_response "knowledge_create" "$response"

    if echo "$response" | grep -q "id"; then
        local kb_id=$(echo "$response" | jq -r '.id // .knowledge_id // empty' 2>/dev/null)
        CREATED_KNOWLEDGE_BASES+=("$kb_id")
        print_success "Knowledge base created with ID: $kb_id"
        test_pass
        return 0
    else
        test_fail "Knowledge base creation failed"
        return 1
    fi
}

test_knowledge_list() {
    test_start "Knowledge Base List"

    local response
    response=$(make_request "GET" "/api/v1/knowledge/" "" "" 2>&1)
    save_response "knowledge_list" "$response"

    if echo "$response" | grep -q "id\|name"; then
        local kb_count
        kb_count=$(echo "$response" | jq 'length' 2>/dev/null || echo "0")
        print_success "Found $kb_count knowledge bases"
        test_pass
        return 0
    else
        test_fail "Could not retrieve knowledge base list"
        return 1
    fi
}

test_knowledge_add_file() {
    test_start "Add File to Knowledge Base"

    if [ ${#CREATED_FILES[@]} -eq 0 ] || [ ${#CREATED_KNOWLEDGE_BASES[@]} -eq 0 ]; then
        print_info "Skipping - no file or knowledge base available"
        return 0
    fi

    local file_id="${CREATED_FILES[0]}"
    local kb_id="${CREATED_KNOWLEDGE_BASES[0]}"

    local data="{\"file_id\":\"$file_id\"}"

    local response
    response=$(make_request "POST" "/api/v1/knowledge/$kb_id/file/add" "$data" "" 2>&1)
    save_response "knowledge_add_file" "$response"

    if echo "$response" | grep -q "files"; then
        print_success "File $file_id added to knowledge base $kb_id"
        test_pass
        return 0
    else
        test_fail "Could not add file to knowledge base"
        return 1
    fi
}

# Chat Completion Tests
test_chat_completion_simple() {
    test_start "Chat Completion (Simple)"

    local chat_model="${SELECTED_MODEL:-openai-codex}"
    local data='{
        "model": "'"$chat_model"'",
        "messages": [
            {
                "role": "user",
                "content": "Hello, how are you?"
            }
        ],
        "temperature": 0.7,
        "max_tokens": 100
    }'

    local response
    response=$(make_request "POST" "/api/chat/completions" "$data" "" 2>&1)
    save_response "chat_completion_simple" "$response"

    if echo "$response" | grep -q "choices"; then
        local reply=$(echo "$response" | jq -r '.choices[0].message.content // empty' 2>/dev/null)
        print_success "Got response: $(echo "$reply" | head -c 50)..."
        test_pass
        return 0
    else
        test_fail "Chat completion failed"
        return 1
    fi
}

test_chat_completion_streaming() {
    test_start "Chat Completion (Streaming)"

    local chat_model="${SELECTED_MODEL:-openai-codex}"
    local data='{
        "model": "'"$chat_model"'",
        "messages": [
            {
                "role": "user",
                "content": "Say hello"
            }
        ],
        "stream": true
    }'

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/chat/completions" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)
    save_response "chat_completion_streaming" "$response"

    if echo "$response" | grep -q "data:"; then
        print_success "Streaming response received"
        test_pass
        return 0
    else
        test_fail "Streaming completion failed"
        return 1
    fi
}

# Embedding Tests
test_embedding_single() {
    test_start "Embedding Generation (Single)"

    local response
    response=$(make_request "GET" "/api/v1/retrieval/" "" "" 2>&1)
    save_response "embedding_single" "$response"

    if command -v jq >/dev/null 2>&1 && echo "$response" | jq -e '.status == true and (.RAG_EMBEDDING_MODEL | type == "string") and (.RAG_EMBEDDING_MODEL | length > 0)' >/dev/null 2>&1; then
        local embedding_model
        embedding_model=$(echo "$response" | jq -r '.RAG_EMBEDDING_MODEL' 2>/dev/null || echo "")
        print_success "Embedding configuration detected: $embedding_model"
        test_pass
        return 0
    else
        test_fail "Embedding configuration check failed"
        return 1
    fi
}

test_embedding_batch() {
    test_start "Embedding Generation (Batch)"

    local collection_name="api-test-collection-$(date +%s)"
    local data='{
        "name": "api-test-embed-batch",
        "content": "First text. Second text. Third text.",
        "collection_name": "'"$collection_name"'"
    }'

    local response
    response=$(make_request "POST" "/api/v1/retrieval/process/text" "$data" "" 2>&1)
    save_response "embedding_batch" "$response"

    if echo "$response" | grep -q '"status":true'; then
        print_success "Retrieval embedding pipeline processed text"
        test_pass
        return 0
    else
        test_fail "Retrieval text processing failed"
        return 1
    fi
}

# Vector Database Tests
test_vector_store_status() {
    test_start "Vector Store Status"

    local response
    response=$(make_request "GET" "/api/v1/retrieval/" "" "" 2>&1)
    save_response "vector_store_status" "$response"

    if echo "$response" | grep -q '"status":true'; then
        test_pass
        return 0
    else
        test_fail "Could not get vector store status"
        return 1
    fi
}

test_vector_search() {
    test_start "Vector Search"

    if [ ${#CREATED_KNOWLEDGE_BASES[@]} -eq 0 ]; then
        print_info "Skipping - no knowledge base available"
        return 0
    fi

    local kb_id="${CREATED_KNOWLEDGE_BASES[0]}"
    local data='{
        "collection_name": "'$kb_id'",
        "query": "test document",
        "k": 3
    }'

    local response
    response=$(make_request "POST" "/api/v1/retrieval/query/doc" "$data" "" 2>&1)
    save_response "vector_search" "$response"

    if echo "$response" | grep -q "results\|documents"; then
        test_pass
        return 0
    else
        test_fail "Vector search failed"
        return 1
    fi
}

# Configuration Tests
test_config_get() {
    test_start "Get Configuration"

    local response
    response=$(make_request "GET" "/api/v1/configs/export" "" "" 2>&1)
    save_response "config_get" "$response"

    if echo "$response" | grep -q "{"; then
        test_pass
        return 0
    else
        test_fail "Could not retrieve configuration"
        return 1
    fi
}

# User Management Tests
test_user_profile() {
    test_start "User Profile"

    local response
    response=$(make_request "GET" "/api/v1/users/user/info" "" "" 2>&1)
    save_response "user_profile" "$response"

    if echo "$response" | grep -q "id\|name\|email"; then
        test_pass
        return 0
    else
        print_info "User profile endpoint not available"
        return 0
    fi
}

###############################################################################
# Cleanup Functions
###############################################################################

cleanup_test_resources() {
    print_section "Cleanup Test Resources"

    for kb_id in "${CREATED_KNOWLEDGE_BASES[@]}"; do
        print_info "Removing knowledge base $kb_id..."
        curl -s -X DELETE \
            "$OPENWEBUI_URL/api/v1/knowledge/$kb_id/delete" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} >/dev/null 2>&1
        print_success "Knowledge base $kb_id removed"
    done

    for file_id in "${CREATED_FILES[@]}"; do
        print_info "Removing file $file_id..."
        curl -s -X DELETE \
            "$OPENWEBUI_URL/api/v1/files/$file_id" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} >/dev/null 2>&1
        print_success "File $file_id removed"
    done

    print_info "Test responses saved to: $OUTPUT_DIR"
}

print_summary() {
    echo -e "\n${CYAN}+------------------------------------------------------------+${NC}"
    echo -e "${CYAN}|  API Test Summary                                           |${NC}"
    echo -e "${CYAN}+------------------------------------------------------------+${NC}\n"

    echo -e "Total Tests: $TESTS_TOTAL"
    echo -e "${GREEN}Passed: $TESTS_PASSED${NC}"
    echo -e "${RED}Failed: $TESTS_FAILED${NC}"

    local pass_rate=0
    if [ $TESTS_TOTAL -gt 0 ]; then
        pass_rate=$((TESTS_PASSED * 100 / TESTS_TOTAL))
    fi

    echo -e "\nPass Rate: ${pass_rate}%"

    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "\n${GREEN}All API tests passed!${NC}"
        return 0
    else
        echo -e "\n${RED}Some tests failed. Check responses in $OUTPUT_DIR${NC}"
        return 1
    fi
}

###############################################################################
# Main Execution
###############################################################################

main() {
    echo -e "\n${CYAN}+------------------------------------------------------------+${NC}"
    echo -e "${CYAN}|  OpenWebUI API Testing Suite                                |${NC}"
    echo -e "${CYAN}+------------------------------------------------------------+${NC}\n"

    RUN_CLEANUP_ON_EXIT=true
    ensure_openwebui_api_key
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "vLLM URL: $VLLM_URL"
    print_info "Test mode: $TEST_MODE"
    if [ -n "$API_KEY" ]; then
        print_info "API Key: ${API_KEY:0:8}..."
    fi
    print_info "Output Directory: $OUTPUT_DIR"

    mkdir -p "$OUTPUT_DIR"

    # Health and Auth
    print_section "Health & Authentication"
    test_repo_real_integration_guard
    test_health_check
    test_root_endpoint
    test_api_key_validation

    # CliProxyApi Smoke
    print_section "CliProxyApi Smoke"
    test_cli_proxy_api_smoke

    # Models
    print_section "Models API"
    test_models_list
    test_model_info

    if [ "$TEST_MODE" = "baseline" ]; then
        print_section "Baseline API Core"
        test_chat_completion_simple
        test_embedding_single
        test_vector_store_status
    else
        # Files
        print_section "Files API"
        test_file_upload_text
        test_file_upload_pdf
        test_file_list

        # Knowledge Bases
        print_section "Knowledge Base API"
        test_knowledge_create
        test_knowledge_list
        test_knowledge_add_file

        # Chat
        print_section "Chat Completions API"
        test_chat_completion_simple
        test_chat_completion_streaming

        # Embeddings
        print_section "Embeddings API"
        test_embedding_single
        test_embedding_batch

        # Vector Database
        print_section "Vector Database API"
        test_vector_store_status
        test_vector_search

        # Configuration
        print_section "Configuration API"
        test_config_get

        # User
        print_section "User API"
        test_user_profile
    fi

    # Summary
    print_summary
}

on_exit() {
    if [ "$RUN_CLEANUP_ON_EXIT" = true ]; then
        cleanup_test_resources
    fi
}

# Trap to ensure cleanup
trap on_exit EXIT INT TERM

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v | --verbose)
            VERBOSE=true
            shift
            ;;
        --baseline)
            TEST_MODE="baseline"
            shift
            ;;
        --full)
            TEST_MODE="full"
            shift
            ;;
        -u | --url)
            OPENWEBUI_URL="$2"
            shift 2
            ;;
        -k | --key)
            API_KEY="$2"
            shift 2
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
