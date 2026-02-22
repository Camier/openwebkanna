#!/bin/bash
#
# test-mcp.sh
# Baseline test script for MCP (Model Context Protocol) integration.
#
# Usage:
#   ./test-mcp.sh [--baseline]
#
# Options:
#   --baseline    Run fast baseline checks only (no heavy tests)
#
# Environment:
#   MCPO_BASE_URL         - MCPO endpoint (default: http://127.0.0.1:8000)
#   OPENWEBUI_URL         - OpenWebUI endpoint (default: http://localhost:3000)
#   OPENWEBUI_API_KEY     - Admin API key for OpenWebUI
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/lib/colors.sh"
source "${SCRIPT_DIR}/lib/print-utils.sh"

MCPO_BASE_URL="${MCPO_BASE_URL:-http://127.0.0.1:8000}"
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
BASELINE_MODE=false
PASSED=0
FAILED=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --baseline)
            BASELINE_MODE=true
            shift
            ;;
        --help | -h)
            echo "Usage: $0 [--baseline]"
            echo ""
            echo "Options:"
            echo "  --baseline    Run fast baseline checks only"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

start_test() {
    echo -e "${CYAN}[TEST]${NC} $1"
}

pass_test() {
    echo -e "${GREEN}[PASS]${NC} $1"
    ((PASSED++))
}

fail_test() {
    echo -e "${RED}[FAIL]${NC} $1"
    ((FAILED++))
}

skip_test() {
    echo -e "${YELLOW}[SKIP]${NC} $1"
}

test_mcpo_health() {
    start_test "MCPO health endpoint responds"
    if curl -sS -f "${MCPO_BASE_URL}/docs" >/dev/null 2>&1; then
        pass_test "MCPO docs endpoint accessible"
    else
        fail_test "MCPO not accessible at ${MCPO_BASE_URL}/docs"
        return 1
    fi
}

test_mcpo_openapi_json() {
    start_test "MCPO OpenAPI spec is valid JSON"
    local response
    response=$(curl -sS "${MCPO_BASE_URL}/openapi.json" 2>/dev/null) || {
        fail_test "Failed to fetch OpenAPI spec"
        return 1
    }

    if echo "$response" | jq empty 2>/dev/null; then
        pass_test "OpenAPI spec is valid JSON"
    else
        fail_test "OpenAPI spec is not valid JSON"
        return 1
    fi
}

test_mcp_servers() {
    local servers=("filesystem" "memory" "fetch" "zotero" "time")
    local server

    for server in "${servers[@]}"; do
        start_test "MCP server '${server}' responds"
        local url="${MCPO_BASE_URL}/${server}/openapi.json"
        local response

        if response=$(curl -sS -f "$url" 2>/dev/null); then
            local tool_count
            tool_count=$(echo "$response" | jq '.paths | keys | length' 2>/dev/null || echo "0")
            pass_test "${server}: ${tool_count} tools available"
        else
            fail_test "${server}: endpoint not accessible"
        fi
    done
}

test_openwebui_tool_servers_endpoint() {
    start_test "OpenWebUI tool servers endpoint responds"

    local token
    token="${OPENWEBUI_API_KEY:-}"

    if [[ -z $token ]]; then
        skip_test "No OPENWEBUI_API_KEY set, skipping authenticated endpoint test"
        return 0
    fi

    local response
    if response=$(curl -sS -f \
        -H "Authorization: Bearer ${token}" \
        "${OPENWEBUI_URL}/api/v1/configs/tool_servers" 2>/dev/null); then
        pass_test "Tool servers endpoint accessible"

        # Count configured MCP servers
        local mcp_count
        mcp_count=$(echo "$response" | jq '[.TOOL_SERVER_CONNECTIONS[] | select(.url | contains(":8000"))] | length' 2>/dev/null || echo "0")
        print_info "Configured MCP servers: ${mcp_count}"
    else
        fail_test "Tool servers endpoint not accessible (HTTP error)"
    fi
}

test_mcp_config_syntax() {
    start_test "MCP config files are valid JSON"
    local config_files=(
        "${SCRIPT_DIR}/mcp/config.json"
        "${SCRIPT_DIR}/mcp/config.research.optional.json"
    )
    local file
    local all_valid=true

    for file in "${config_files[@]}"; do
        if [[ -f $file ]]; then
            if jq empty "$file" 2>/dev/null; then
                print_info "$(basename "$file"): valid JSON"
            else
                fail_test "$(basename "$file"): invalid JSON syntax"
                all_valid=false
            fi
        else
            skip_test "$(basename "$file"): file not found"
        fi
    done

    if [[ $all_valid == true ]]; then
        pass_test "All MCP config files are valid JSON"
    fi
}

test_docker_compose_mcpo_service() {
    start_test "MCPO service defined in docker-compose.yml"

    if grep -q "^  mcpo:" "${SCRIPT_DIR}/docker-compose.yml" 2>/dev/null; then
        pass_test "MCPO service defined in docker-compose.yml"
    else
        fail_test "MCPO service not found in docker-compose.yml"
    fi
}

main() {
    print_header "MCP Integration Tests"
    print_info "MCPO URL: ${MCPO_BASE_URL}"
    print_info "OpenWebUI URL: ${OPENWEBUI_URL}"
    if [[ $BASELINE_MODE == true ]]; then
        print_info "Mode: Baseline (fast checks only)"
    fi
    echo ""

    # Always run these fast checks
    test_docker_compose_mcpo_service
    test_mcp_config_syntax

    # Run MCPO endpoint tests if MCPO appears to be running
    if curl -sS -f "${MCPO_BASE_URL}/docs" >/dev/null 2>&1; then
        print_info "MCPO detected, running endpoint tests..."
        test_mcpo_health
        test_mcpo_openapi_json

        if [[ $BASELINE_MODE == false ]]; then
            test_mcp_servers
        else
            skip_test "Individual MCP server tests (use without --baseline for full tests)"
        fi

        test_openwebui_tool_servers_endpoint
    else
        skip_test "MCPO not running at ${MCPO_BASE_URL}, skipping endpoint tests"
        skip_test "Start the stack with: ./deploy.sh --no-logs"
    fi

    echo ""
    print_section "Test Summary"
    echo -e "Passed: ${GREEN}${PASSED}${NC}"
    if [[ $FAILED -gt 0 ]]; then
        echo -e "Failed: ${RED}${FAILED}${NC}"
    else
        echo -e "Failed: ${PASSED}${NC}"
    fi

    if [[ $FAILED -eq 0 ]]; then
        print_success "All MCP tests passed"
        exit 0
    else
        print_error "Some MCP tests failed"
        exit 1
    fi
}

main "$@"
