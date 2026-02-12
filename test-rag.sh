#!/bin/bash

###############################################################################
# RAG Functionality Testing Script
# This script tests the complete RAG deployment including OpenWebUI and vLLM
###############################################################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:3000}"
VLLM_URL="${VLLM_URL:-http://localhost:8000}"
API_KEY="${OPENWEBUI_API_KEY:-}"
VERBOSE="${VERBOSE:-false}"
CHAT_MODEL="${RAG_CHAT_MODEL:-}"
CLIPROXY_BIN="${CLIPROXY_BIN:-./cli-proxy-api.sh}"
CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0

###############################################################################
# Helper Functions
###############################################################################

print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}\n"
}

print_section() {
    echo -e "\n${YELLOW}>>> $1${NC}\n"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

test_start() {
    ((TESTS_TOTAL++))
    TEST_NAME="$1"
    print_info "Testing: $TEST_NAME"
}

test_pass() {
    ((TESTS_PASSED++))
    print_success "$TEST_NAME passed"
}

test_fail() {
    ((TESTS_FAILED++))
    print_error "$TEST_NAME failed: $1"
}

make_request() {
    local method="$1"
    local url="$2"
    local data="$3"
    local headers="-H 'Content-Type: application/json'"

    if [ -n "$API_KEY" ]; then
        headers="$headers -H 'Authorization: Bearer $API_KEY'"
    fi

    if [ "$method" = "GET" ]; then
        curl -s -X GET "$url" $headers
    elif [ "$method" = "POST" ]; then
        curl -s -X POST "$url" $headers -d "$data"
    elif [ "$method" = "DELETE" ]; then
        curl -s -X DELETE "$url" $headers
    fi
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
    if [[ "$model_count" =~ ^[0-9]+$ ]] && [ "$model_count" -gt 0 ]; then
        return 0
    fi
    return 1
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

    local response
    response=$(curl -s -w "%{http_code}" -o /dev/null "$VLLM_URL/health" 2>&1 || echo "000")

    if [ "$response" = "200" ]; then
        test_pass
        return 0
    else
        test_fail "vLLM health check returned HTTP $response"
        return 1
    fi
}

test_vllm_models() {
    test_start "vLLM Models Available"

    local response
    response=$(curl -s "$VLLM_URL/v1/models" 2>&1)

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

    local data='{
        "model": "text-embedding-ada-002",
        "input": "This is a test text for embedding generation."
    }'

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/embeddings" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)

    if echo "$response" | grep -q "\["; then
        local embedding_length=$(echo "$response" | jq '.[0] | length' 2>/dev/null || echo "0")
        print_success "Embedding generated with length: $embedding_length"
        test_pass
        return 0
    else
        test_fail "Embedding generation failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_create_test_document() {
    test_start "Create Test Document via API"

    local test_content="This is a test document for RAG testing.
It contains information about artificial intelligence and machine learning.
RAG systems combine retrieval and generation for better responses.
This document will be used to test the knowledge base functionality."

    local temp_file="/tmp/rag_test_doc_$(date +%s).txt"
    echo "$test_content" > "$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/files/upload" \
        -H "Content-Type: multipart/form-data" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -F "file=@$temp_file" \
        -F "meta={\"name\":\"rag_test.txt\"}" 2>&1)

    rm -f "$temp_file"

    if echo "$response" | grep -q "id"; then
        TEST_DOC_ID=$(echo "$response" | jq -r '.id // .file_id // empty' 2>/dev/null)
        print_success "Document created with ID: $TEST_DOC_ID"
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
        "$OPENWEBUI_URL/api/knowledge" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)

    if echo "$response" | grep -q "id"; then
        TEST_KB_ID=$(echo "$response" | jq -r '.id // .knowledge_id // empty' 2>/dev/null)
        print_success "Knowledge base created with ID: $TEST_KB_ID"
        test_pass
        return 0
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

    local data='{
        "texts": ["Test text for vector store", "Another test text"],
        "model": "text-embedding-ada-002"
    }'

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/vectordb/embed" \
        -H "Content-Type: application/json" \
        ${API_KEY:+-H "Authorization: Bearer $API_KEY"} \
        -d "$data" 2>&1)

    if echo "$response" | grep -q "embeddings\|vectors"; then
        print_success "Vector store is operational"
        test_pass
        return 0
    else
        test_fail "Vector store test failed"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi
}

test_retrieval() {
    test_start "Document Retrieval"

    local query="RAG testing validation"
    local data='{
        "query": "'"$query"'",
        "top_k": 3,
        "knowledge_id": "'${TEST_KB_ID:-}'"
    }'

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/knowledge/query" \
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
    test_start "CliProxyApi Health (all)"

    local response
    if ! response=$("$CLIPROXY_BIN" --raw health all 2>&1); then
        test_fail "Command failed: $CLIPROXY_BIN --raw health all"
        [ "$VERBOSE" = "true" ] && print_info "Response: $response"
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        if echo "$response" | jq -e '.openwebui and .vllm' >/dev/null 2>&1; then
            test_pass
            return 0
        fi
    elif echo "$response" | grep -q '"openwebui"' && echo "$response" | grep -q '"vllm"'; then
        test_pass
        return 0
    fi

    test_fail "CLIProxyAPI health output missing openwebui/vllm keys"
    [ "$VERBOSE" = "true" ] && print_info "Response: $response"
    return 1
}

test_cli_proxy_api_models_vllm() {
    test_start "CliProxyApi Models (vLLM)"

    local response
    if ! response=$("$CLIPROXY_BIN" --raw models vllm 2>&1); then
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
    if ! response=$("$CLIPROXY_BIN" --raw models webui 2>&1); then
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
            "$OPENWEBUI_URL/api/files/$TEST_DOC_ID" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} > /dev/null 2>&1
        print_success "Test document removed"
    fi

    if [ -n "$TEST_KB_ID" ]; then
        print_info "Removing test knowledge base $TEST_KB_ID..."
        curl -s -X DELETE \
            "$OPENWEBUI_URL/api/knowledge/$TEST_KB_ID" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} > /dev/null 2>&1
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
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "vLLM URL: $VLLM_URL"
    print_info "CliProxyApi: $CLIPROXY_BIN"
    [ -n "$API_KEY" ] && print_info "API Key: ${API_KEY:0:8}..."

    # CLIProxyAPI smoke checks
    print_section "CliProxyApi Integration"
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
    test_create_test_document
    test_create_knowledge_base

    # Integration tests
    print_section "Integration Tests"
    test_rag_query
    test_retrieval

    # Cleanup
    cleanup_test_resources

    # Summary
    print_summary
}

# Trap to cleanup on interrupt
trap cleanup_test_resources EXIT INT TERM

# Run main function
main "$@"
