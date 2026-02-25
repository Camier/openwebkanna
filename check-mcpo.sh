#!/bin/bash

###############################################################################
# check-mcpo.sh - Verify MCPO connectivity and MCP server status
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/print-utils.sh"
source "${SCRIPT_DIR}/lib/docker-helpers.sh"
source "${SCRIPT_DIR}/lib/env-loader.sh"

load_env_defaults
cd "$SCRIPT_DIR"

MCPO_PORT="${MCPO_PORT:-8000}"
MCPO_URL="http://127.0.0.1:${MCPO_PORT}"

MCP_SERVERS=("filesystem" "memory" "fetch" "zotero" "time")

print_header "MCPO Health Check"

check_container() {
    print_step "Checking MCPO container..."

    if ! docker ps --filter "name=mcpo" --filter "status=running" --format "{{.Names}}" | grep -q "mcpo"; then
        print_error "MCPO container is not running"
        return 1
    fi

    local status
    status=$(docker inspect --format='{{.State.Health.Status}}' mcpo 2>/dev/null || echo "unknown")

    if [[ $status == "healthy" ]]; then
        print_success "Container is healthy"
        return 0
    else
        print_warning "Container status: $status"
        return 1
    fi
}

check_docs() {
    print_step "Checking MCPO docs endpoint..."

    if curl -sf "${MCPO_URL}/docs" >/dev/null 2>&1; then
        print_success "Docs endpoint accessible"
        return 0
    else
        print_error "Cannot reach MCPO docs"
        return 1
    fi
}

check_mcp_servers() {
    print_step "Checking MCP server endpoints..."

    local failed=0
    local server title

    for server in "${MCP_SERVERS[@]}"; do
        title=$(curl -sf "${MCPO_URL}/${server}/openapi.json" 2>/dev/null | jq -r '.info.title' 2>/dev/null)

        if [[ -n $title && $title != "null" ]]; then
            echo "  ${GREEN}[OK]${NC} ${server}: ${title}"
        else
            echo "  ${RED}[FAIL]${NC} ${server}: Not accessible"
            ((failed++))
        fi
    done

    return $failed
}

check_openwebui_connectivity() {
    print_step "Checking connectivity from OpenWebUI container..."

    local result
    result=$(docker exec openwebui curl -sf "http://mcpo:8000/docs" 2>&1 | head -c 100)

    if [[ -n $result && $result == *"<!DOCTYPE html>"* ]]; then
        print_success "OpenWebUI can reach MCPO"
        return 0
    else
        print_error "OpenWebUI cannot reach MCPO"
        print_info "Checking network configuration..."

        local openwebui_net mcpo_net
        openwebui_net=$(docker inspect openwebui --format '{{range .NetworkSettings.Networks}}{{.NetworkID}}{{end}}' 2>/dev/null)
        mcpo_net=$(docker inspect mcpo --format '{{range .NetworkSettings.Networks}}{{.NetworkID}}{{end}}' 2>/dev/null)

        if [[ $openwebui_net == *"$mcpo_net"* ]]; then
            print_info "Containers share a network"
        else
            print_warning "Containers may not share a network"
            print_info "OpenWebUI networks: $openwebui_net"
            print_info "MCPO networks: $mcpo_net"
        fi
        return 1
    fi
}

main() {
    local errors=0

    check_container || ((errors++))
    check_docs || ((errors++))
    check_mcp_servers || true
    check_openwebui_connectivity || ((errors++))

    echo ""
    print_section "Summary"

    if [[ $errors -eq 0 ]]; then
        print_success "All checks passed"
        exit 0
    else
        print_error "$errors check(s) failed"
        exit 1
    fi
}

main "$@"
