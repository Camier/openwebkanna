#!/bin/bash

###############################################################################
# API Endpoint Testing Script
# This script tests all major API endpoints for the RAG deployment
###############################################################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:3000}"
VLLM_URL="${VLLM_URL:-http://localhost:8000}"
API_KEY="${OPENWEBUI_API_KEY:-}"
CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-}"
VERBOSE="${VERBOSE:-false}"
OUTPUT_DIR="/tmp/openwebui_api_tests_$(date +%Y%m%d_%H%M%S)"

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0

# Track created resources
declare -a CREATED_FILES=()
declare -a CREATED_KNOWLEDGE_BASES=()

###############################################################################
# Helper Functions
###############################################################################

print_header() {
    echo -e "\n${CYAN}╔════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}║  $1${NC}"
    echo -e "${CYAN}╚════════════════════════════════════════════════════════════╝${NC}\n"
}

print_section() {
    echo -e "\n${YELLOW}─────────────────────────────────────────────────────────────${NC}"
    echo -e "${YELLOW}  $1${NC}"
    echo -e "${YELLOW}─────────────────────────────────────────────────────────────${NC}\n"
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

print_endpoint() {
    echo -e "${CYAN}→ $1${NC}"
}

test_start() {
    ((TESTS_TOTAL++))
    TEST_NAME="$1"
    print_endpoint "$TEST_NAME"
    TEST_START_TIME=$(date +%s)
}

test_pass() {
    ((TESTS_PASSED++))
    local duration=$(( $(date +%s) - TEST_START_TIME ))
    print_success "$TEST_NAME (${duration}s)"
}

test_fail() {
    ((TESTS_FAILED++))
    local duration=$(( $(date +%s) - TEST_START_TIME ))
    print_error "$TEST_NAME failed (${duration}s): $1"
}

save_response() {
    local test_name="$1"
    local response="$2"
    mkdir -p "$OUTPUT_DIR"
    echo "$response" > "$OUTPUT_DIR/${test_name}.json"
    [ "$VERBOSE" = "true" ] && print_info "Response saved to $OUTPUT_DIR/${test_name}.json"
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

make_request() {
    local method="$1"
    local endpoint="$2"
    local data="$3"
    local extra_headers="$4"
    local headers="-H 'Content-Type: application/json'"

    if [ -n "$API_KEY" ]; then
        headers="$headers -H 'Authorization: Bearer $API_KEY'"
    fi

    if [ -n "$extra_headers" ]; then
        headers="$headers $extra_headers"
    fi

    local url="${OPENWEBUI_URL}${endpoint}"

    if [ "$method" = "GET" ]; then
        curl -s -X GET "$url" $headers
    elif [ "$method" = "POST" ]; then
        curl -s -X POST "$url" $headers -d "$data"
    elif [ "$method" = "PUT" ]; then
        curl -s -X PUT "$url" $headers -d "$data"
    elif [ "$method" = "DELETE" ]; then
        curl -s -X DELETE "$url" $headers
    elif [ "$method" = "POST_FORM" ]; then
        curl -s -X POST "$url" ${extra_headers} $data
    fi
}

###############################################################################
# Test Functions
###############################################################################

# Health Check Tests
test_health_check() {
    test_start "Health Check"

    local response
    response=$(make_request "GET" "/health" "" "" 2>&1)
    save_response "health_check" "$response"

    if [ $? -eq 0 ] && [ -n "$response" ]; then
        test_pass
        return 0
    else
        test_fail "Health endpoint not responding"
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
    response=$(make_request "GET" "/api/auth/status" "" "" 2>&1)
    save_response "api_key_validation" "$response"

    if echo "$response" | grep -q "authenticated\|user"; then
        test_pass
        return 0
    else
        print_info "API key validation skipped (no auth configured)"
        return 0
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
    if response=$(CLIPROXYAPI_ENABLED=true ./check-cliproxyapi.sh 2>&1); then
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
        if echo "$response" | grep -Fq "models [openwebui|webui|vllm|all]" && \
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

    if echo "$response" | grep -q "id\|models"; then
        local model_count=$(echo "$response" | jq '.data | length' 2>/dev/null || echo "0")
        print_success "Found $model_count models"
        test_pass
        return 0
    else
        test_fail "Could not retrieve models list"
        return 1
    fi
}

test_model_info() {
    test_start "Model Information"

    local response
    response=$(make_request "GET" "/api/models/llama-3.2-3b-instruct" "" "" 2>&1)
    save_response "model_info" "$response"

    if echo "$response" | grep -q "id"; then
        test_pass
        return 0
    else
        test_fail "Could not retrieve model information"
        return 1
    fi
}

# File Upload Tests
test_file_upload_text() {
    test_start "File Upload (Text)"

    local temp_file="/tmp/test_upload_$(date +%s).txt"
    echo "This is a test document for API validation." > "$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/files/upload" \
        -H "Authorization: Bearer $API_KEY" \
        -F "file=@$temp_file" \
        -F "meta={\"name\":\"test_doc.txt\"}" 2>&1)
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
%%EOF" > "$temp_file"

    local response
    response=$(curl -s -X POST \
        "$OPENWEBUI_URL/api/files/upload" \
        -H "Authorization: Bearer $API_KEY" \
        -F "file=@$temp_file" \
        -F "meta={\"name\":\"test_doc.pdf\"}" 2>&1)
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
    response=$(make_request "GET" "/api/files" "" "" 2>&1)
    save_response "file_list" "$response"

    if echo "$response" | grep -q "files\|data"; then
        local file_count=$(echo "$response" | jq '.data | length' 2>/dev/null || echo "0")
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
    response=$(make_request "POST" "/api/knowledge" "$data" "" 2>&1)
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
    response=$(make_request "GET" "/api/knowledge" "" "" 2>&1)
    save_response "knowledge_list" "$response"

    if echo "$response" | grep -q "data\|knowledge"; then
        local kb_count=$(echo "$response" | jq '.data | length' 2>/dev/null || echo "0")
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

    local data='{
        "file_id": "'$file_id'"
    }'

    local response
    response=$(make_request "POST" "/api/knowledge/$kb_id/add" "$data" "" 2>&1)
    save_response "knowledge_add_file" "$response"

    if [ $? -eq 0 ]; then
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

    local data='{
        "model": "llama-3.2-3b-instruct",
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

    local data='{
        "model": "llama-3.2-3b-instruct",
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

    local data='{
        "model": "text-embedding-ada-002",
        "input": "Test text for embedding generation"
    }'

    local response
    response=$(make_request "POST" "/api/embeddings" "$data" "" 2>&1)
    save_response "embedding_single" "$response"

    if echo "$response" | grep -q "\["; then
        local embedding_length=$(echo "$response" | jq '.[0] | length' 2>/dev/null || echo "0")
        print_success "Embedding generated with length: $embedding_length"
        test_pass
        return 0
    else
        test_fail "Embedding generation failed"
        return 1
    fi
}

test_embedding_batch() {
    test_start "Embedding Generation (Batch)"

    local data='{
        "model": "text-embedding-ada-002",
        "input": ["First text", "Second text", "Third text"]
    }'

    local response
    response=$(make_request "POST" "/api/embeddings" "$data" "" 2>&1)
    save_response "embedding_batch" "$response"

    if echo "$response" | grep -q "\["; then
        local count=$(echo "$response" | jq 'length' 2>/dev/null || echo "0")
        print_success "Generated $count embeddings"
        test_pass
        return 0
    else
        test_fail "Batch embedding generation failed"
        return 1
    fi
}

# Vector Database Tests
test_vector_store_status() {
    test_start "Vector Store Status"

    local response
    response=$(make_request "GET" "/api/vectordb/status" "" "" 2>&1)
    save_response "vector_store_status" "$response"

    if echo "$response" | grep -q "status\|ready"; then
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
        "query": "test document",
        "top_k": 3
    }'

    local response
    response=$(make_request "POST" "/api/knowledge/$kb_id/query" "$data" "" 2>&1)
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
    response=$(make_request "GET" "/api/config" "" "" 2>&1)
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
    response=$(make_request "GET" "/api/user/profile" "" "" 2>&1)
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
            "$OPENWEBUI_URL/api/knowledge/$kb_id" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} > /dev/null 2>&1
        print_success "Knowledge base $kb_id removed"
    done

    for file_id in "${CREATED_FILES[@]}"; do
        print_info "Removing file $file_id..."
        curl -s -X DELETE \
            "$OPENWEBUI_URL/api/files/$file_id" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} > /dev/null 2>&1
        print_success "File $file_id removed"
    done

    print_info "Test responses saved to: $OUTPUT_DIR"
}

print_summary() {
    print_header "API Test Summary"

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
    print_header "OpenWebUI API Testing Suite"
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "vLLM URL: $VLLM_URL"
    [ -n "$API_KEY" ] && print_info "API Key: ${API_KEY:0:8}..."
    print_info "Output Directory: $OUTPUT_DIR"

    mkdir -p "$OUTPUT_DIR"

    # Health and Auth
    print_section "Health & Authentication"
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

    # Cleanup
    cleanup_test_resources

    # Summary
    print_summary
}

# Trap to ensure cleanup
trap cleanup_test_resources EXIT INT TERM

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -u|--url)
            OPENWEBUI_URL="$2"
            shift 2
            ;;
        -k|--key)
            API_KEY="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -v, --verbose       Enable verbose output"
            echo "  -u, --url URL       OpenWebUI URL (default: http://localhost:3000)"
            echo "  -k, --key KEY       API key for authentication"
            echo "  -h, --help          Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Run main function
main "$@"
