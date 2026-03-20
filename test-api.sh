#!/bin/bash

###############################################################################
# API Endpoint Testing Script
# This script tests all major API endpoints for the RAG deployment
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
    if [[ $candidate == *"://host.docker.internal"* ]] && [ ! -f "/.dockerenv" ]; then
        candidate="$(printf "%s" "$candidate" | sed 's#://host.docker.internal#://localhost#')"
    fi
    printf "%s" "$candidate"
}

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
MULTIMODAL_RETRIEVAL_API_URL="${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}"
OPENAI_BASE_ROOT="$(resolve_openai_base_root)"
UPSTREAM_URL="${UPSTREAM_URL:-${OPENAI_BASE_ROOT:-http://localhost:4000}}"
API_KEY="${OPENWEBUI_API_KEY:-${API_KEY:-}}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
OPENWEBUI_TEST_MODEL="${OPENWEBUI_TEST_MODEL:-}"
OPENWEBUI_SMOKE_MODEL_CANDIDATES="${OPENWEBUI_SMOKE_MODEL_CANDIDATES:-$(default_openwebui_smoke_model_candidates)}"
NO_MOCK_AUDIT_SCRIPT="${NO_MOCK_AUDIT_SCRIPT:-./scripts/testing/audit-no-mock.sh}"
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
MODELS_RESPONSE=""

###############################################################################
# Helper Functions
###############################################################################

print_endpoint() {
    echo -e "${CYAN}>> $1${NC}"
}

test_start() {
    ((TESTS_TOTAL += 1))
    TEST_NAME="$1"
    print_endpoint "$TEST_NAME"
    TEST_START_TIME=$(date +%s)
}

test_pass() {
    ((TESTS_PASSED += 1))
    local duration
    duration=$(($(date +%s) - TEST_START_TIME))
    print_success "$TEST_NAME (${duration}s)"
}

test_fail() {
    ((TESTS_FAILED += 1))
    local duration
    duration=$(($(date +%s) - TEST_START_TIME))
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

    curl "${curl_args[@]}" 2>/dev/null
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

make_canonical_request_with_status() {
    local method="$1"
    local endpoint="$2"
    local data="$3"
    local response_var="$4"
    local status_var="$5"
    local url="${MULTIMODAL_RETRIEVAL_API_URL%/}${endpoint}"
    local response_file=""
    local response_body=""
    local status_code=""
    local -a curl_args=("-sS" "-m" "45" "-X" "$method" "$url" "-H" "Content-Type: application/json")

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

    local http_code
    http_code="$(curl -s -I -o /dev/null -w "%{http_code}" "$OPENWEBUI_URL/" 2>/dev/null || true)"
    save_response "root_endpoint" "$http_code"

    if [ "$http_code" = "200" ]; then
        test_pass
        return 0
    else
        test_fail "Root endpoint returned HTTP ${http_code:-000}"
        return 1
    fi
}

test_canonical_health_check() {
    test_start "Canonical Retrieval Health Check"

    local response http_code
    make_canonical_request_with_status "GET" "/health" "" response http_code
    save_response "canonical_health_check" "$response"

    if [ "$http_code" = "200" ] && echo "$response" | jq -e '.status == "ok"' >/dev/null 2>&1; then
        test_pass
        return 0
    fi

    test_fail "Canonical /health returned HTTP ${http_code:-000}"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

test_canonical_ready_check() {
    test_start "Canonical Retrieval Ready Check"

    local response http_code
    make_canonical_request_with_status "GET" "/ready" "" response http_code
    save_response "canonical_ready_check" "$response"

    if [ "$http_code" = "200" ] && echo "$response" | jq -e '.status == "ready"' >/dev/null 2>&1; then
        test_pass
        return 0
    fi

    test_fail "Canonical /ready returned HTTP ${http_code:-000}"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

test_canonical_retrieve() {
    test_start "Canonical Retrieval Query"

    local data response http_code
    data='{
        "query": "sceletium alkaloids multimodal evidence",
        "top_k": 3
    }'
    make_canonical_request_with_status "POST" "/api/v1/retrieve" "$data" response http_code
    save_response "canonical_retrieve" "$response"

    if [ "$http_code" = "200" ]; then
        if echo "$response" | jq -e '
            (.query | type == "string") and
            (.top_k == 3) and
            (.backend | type == "object") and
            (.candidate_hits | type == "array") and
            (.reranked_hits | type == "array") and
            (.evidence_objects | type == "array")
        ' >/dev/null 2>&1; then
            test_pass
            return 0
        fi
    fi

    test_fail "Canonical retrieve returned HTTP ${http_code:-000} or invalid payload"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

# Authentication Tests
test_api_key_validation() {
    test_start "API Key Validation"

    local response
    local http_code
    make_request_with_status "GET" "/api/v1/auths/" "" response http_code
    save_response "api_key_validation" "$response"

    if [ "$http_code" = "200" ]; then
        if echo "$response" | jq -e '
            type == "object" and
            (.id | type == "string") and
            (.email | type == "string") and
            (.token | type == "string")
        ' >/dev/null 2>&1; then
            test_pass
            return 0
        fi
        test_fail "Auth profile returned HTTP 200 but not the expected JSON payload"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    elif [ "$http_code" = "401" ] || [ "$http_code" = "403" ]; then
        if is_true "${WEBUI_AUTH:-false}"; then
            test_fail "Authentication is enabled but API key is unauthorized (HTTP $http_code)"
            return 1
        fi
        print_info "API key validation skipped (auth appears disabled; profile endpoint returned HTTP $http_code)"
        test_pass
        return 0
    else
        test_fail "Auth profile endpoint returned HTTP ${http_code:-000}"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

# Model Tests
test_models_list() {
    test_start "Models List"

    local response
    response=$(make_request "GET" "/api/models" "" "" 2>&1)
    MODELS_RESPONSE="$response"
    save_response "models_list" "$response"

    if has_non_empty_models_array "$response"; then
        local model_count
        local preferred_model
        model_count=$(echo "$response" | jq '.data | length' 2>/dev/null || echo "0")
        SELECTED_MODEL="$(first_model_id_from_response "$response")"
        if [ -n "$OPENWEBUI_TEST_MODEL" ]; then
            if models_response_has_id "$response" "$OPENWEBUI_TEST_MODEL"; then
                SELECTED_MODEL="$OPENWEBUI_TEST_MODEL"
                print_info "Using OPENWEBUI_TEST_MODEL=$OPENWEBUI_TEST_MODEL"
            else
                print_warning "OPENWEBUI_TEST_MODEL not present in /api/models: $OPENWEBUI_TEST_MODEL"
            fi
        elif [ "$TEST_MODE" = "baseline" ]; then
            preferred_model="$(choose_preferred_model_from_response "$response")"
            if [ -n "$preferred_model" ]; then
                SELECTED_MODEL="$preferred_model"
                print_info "Using baseline smoke model: $preferred_model"
            fi
        fi
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

    local temp_file
    temp_file="/tmp/test_upload_$(date +%s).txt"
    echo "This is a test document for API validation." >"$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/v1/files/?process=true&process_in_background=false" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -F "file=@$temp_file" \
        -F 'metadata={"name":"test_doc.txt"}' 2>&1)
    save_response "file_upload_text" "$response"

    rm -f "$temp_file"

    if json_response_has_any_id "$response"; then
        local file_id
        file_id="$(extract_json_primary_id "$response")"
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
    local temp_file
    temp_file="/tmp/test_upload_$(date +%s).pdf"
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

    if json_response_has_any_id "$response"; then
        local file_id
        file_id="$(extract_json_primary_id "$response")"
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

    if json_response_has_any_id "$response"; then
        local kb_id
        kb_id="$(extract_json_primary_id "$response")"
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

    if json_response_has_files_field "$response"; then
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

    choose_chat_fallback_model() {
        local current_model="$1"
        choose_preferred_model_from_response "$MODELS_RESPONSE" "$current_model"
    }

    run_chat_completion_once() {
        local model="$1"
        local payload response_local
        payload='{
            "model": "'"$model"'",
            "messages": [
                {
                    "role": "user",
                    "content": "Hello, how are you?"
                }
            ],
            "temperature": 0.7,
            "max_tokens": 100
        }'
        response_local=$(make_request "POST" "/api/chat/completions" "$payload" "" 2>&1)
        printf "%s" "$response_local"
    }

    local chat_model="${SELECTED_MODEL:-openai-codex}"
    local response fallback_model
    response="$(run_chat_completion_once "$chat_model")"

    if ! json_response_has_chat_choices "$response"; then
        fallback_model="$(choose_chat_fallback_model "$chat_model")"
        if [ -n "$fallback_model" ]; then
            print_warning "Primary model failed ($chat_model). Retrying with fallback: $fallback_model"
            response="$(run_chat_completion_once "$fallback_model")"
            if json_response_has_chat_choices "$response"; then
                chat_model="$fallback_model"
                SELECTED_MODEL="$fallback_model"
            fi
        fi
    fi

    save_response "chat_completion_simple" "$response"

    if json_response_has_chat_choices "$response"; then
        local reply
        reply="$(extract_chat_response_text "$response")"
        print_success "Got response from model $chat_model: $(echo "$reply" | head -c 50)..."
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
    response=$(curl -s -m 45 -X POST \
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
    test_start "Legacy OpenWebUI Retrieval Embedding Configuration"

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
    test_start "Legacy OpenWebUI Retrieval Text Processing"

    local collection_name
    collection_name="api-test-collection-$(date +%s)"
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
    test_start "Legacy OpenWebUI Retrieval Status"

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
    test_start "Legacy OpenWebUI Retrieval Query"

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

    if json_response_has_retrieval_hits "$response"; then
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
    print_info "Canonical Retrieval URL: $MULTIMODAL_RETRIEVAL_API_URL"
    print_info "Upstream URL: $UPSTREAM_URL"
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
    test_canonical_health_check
    test_canonical_ready_check
    test_api_key_validation

    # Models
    print_section "Models API"
    test_models_list
    test_model_info

    if [ "$TEST_MODE" = "baseline" ]; then
        print_section "Baseline API Core"
        test_canonical_retrieve
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
        print_section "Canonical Retrieval API"
        test_canonical_retrieve

        print_section "Legacy OpenWebUI Retrieval API"
        test_embedding_single
        test_embedding_batch

        # Vector Database
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
