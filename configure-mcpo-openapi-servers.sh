#!/bin/bash
#
# configure-mcpo-openapi-servers.sh
# Automatically configure MCPO OpenAPI servers in OpenWebUI via admin API.
#
# Usage:
#   ./configure-mcpo-openapi-servers.sh [--dry-run] [--verify]
#
# Environment variables:
#   OPENWEBUI_URL       - OpenWebUI base URL (default: http://localhost:3000)
#   OPENWEBUI_API_KEY   - Admin API key (or will prompt for login)
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/lib/colors.sh"
source "${SCRIPT_DIR}/lib/print-utils.sh"

OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:3000}"
MCPO_URL="${MCPO_URL:-http://mcpo:8000}"
MCPO_HOST_URL="${MCPO_HOST_URL:-http://127.0.0.1:8000}"
DRY_RUN=false
VERIFY_ONLY=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --verify)
            VERIFY_ONLY=true
            shift
            ;;
        --help | -h)
            echo "Usage: $0 [--dry-run] [--verify]"
            echo ""
            echo "Options:"
            echo "  --dry-run    Show what would be configured without making changes"
            echo "  --verify     Only verify existing servers, don't add new ones"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

declare -A MCPO_SERVERS=(
    ["filesystem"]="Filesystem MCP|Read/write access to thesis files"
    ["memory"]="Memory MCP|Session memory and knowledge graph"
    ["fetch"]="Fetch MCP|Web page fetching for research"
    ["zotero"]="Zotero MCP|Local Zotero library (read-only)"
    ["time"]="Time MCP|Date and time utilities"
)

print_header "MCPO OpenAPI Server Configuration"

get_auth_token() {
    if [[ -n ${OPENWEBUI_API_KEY:-} ]]; then
        echo "${OPENWEBUI_API_KEY}"
        return 0
    fi

    if [[ -f "${SCRIPT_DIR}/.env" ]]; then
        local key
        key=$(grep -E "^OPENWEBUI_API_KEY=" "${SCRIPT_DIR}/.env" 2>/dev/null | head -1 | cut -d'=' -f2- | tr -d '"' | tr -d "'")
        if [[ -n $key ]]; then
            echo "$key"
            return 0
        fi
    fi

    print_warning "No OPENWEBUI_API_KEY found. Please set it in environment or .env"
    return 1
}

verify_mcpo() {
    print_step "Verifying MCPO endpoint at ${MCPO_HOST_URL}..."

    if curl -sS -f -o /dev/null "${MCPO_HOST_URL}/docs" 2>/dev/null; then
        print_success "MCPO is accessible"
        return 0
    else
        print_error "MCPO not accessible at ${MCPO_HOST_URL}"
        print_info "Make sure the stack is running: ./deploy.sh --no-logs"
        return 1
    fi
}

get_current_servers() {
    local token="$1"
    local response

    response=$(curl -sS -f \
        -H "Authorization: Bearer ${token}" \
        "${OPENWEBUI_URL}/api/v1/configs/tool_servers" 2>/dev/null) || {
        print_error "Failed to get current tool servers"
        return 1
    }

    echo "$response"
}

build_server_config() {
    local server_name="$1"
    local description="$2"
    local url="${MCPO_URL}/${server_name}"

    jq -n \
        --arg url "$url" \
        --arg path "openapi.json" \
        --arg type "openapi" \
        --arg auth_type "none" \
        --arg name "${server_name^} Server" \
        --arg desc "$description" \
        '{
            "url": $url,
            "path": $path,
            "type": $type,
            "auth_type": $auth_type,
            "key": "",
            "config": {
                "enable": true
            }
        }'
}

verify_server_endpoint() {
    local server_name="$1"
    local url="${MCPO_HOST_URL}/${server_name}/openapi.json"

    print_step "Verifying ${server_name} OpenAPI at ${url}..."

    if curl -sS -f -o /dev/null "${url}" 2>/dev/null; then
        local tool_count
        tool_count=$(curl -sS "${url}" 2>/dev/null | jq '.paths | keys | length' 2>/dev/null || echo "?")
        print_success "${server_name}: ${tool_count} tools available"
        return 0
    else
        print_error "${server_name}: endpoint not accessible"
        return 1
    fi
}

configure_servers() {
    local token="$1"
    local current_config
    # shellcheck disable=SC2034
    local new_servers=()
    local existing_count=0
    local added_count=0

    print_step "Fetching current tool server configuration..."
    current_config=$(get_current_servers "$token")

    if [[ -z $current_config ]]; then
        print_error "Empty response from API"
        return 1
    fi

    local existing_urls
    existing_urls=$(echo "$current_config" | jq -r '.TOOL_SERVER_CONNECTIONS[].url // empty' 2>/dev/null || true)

    local all_servers="[]"

    if [[ -n $existing_urls ]]; then
        while IFS= read -r existing_url; do
            local is_mcpo=false
            for server_name in "${!MCPO_SERVERS[@]}"; do
                if [[ $existing_url == "${MCPO_URL}/${server_name}" ]]; then
                    is_mcpo=true
                    break
                fi
            done

            if [[ $is_mcpo == "false" ]]; then
                all_servers=$(echo "$current_config" | jq \
                    --arg url "$existing_url" \
                    '.TOOL_SERVER_CONNECTIONS | map(select(.url == $url))[0]' 2>/dev/null |
                    jq -s '. + input' "$all_servers" - 2>/dev/null || echo "$all_servers")
            fi
        done <<<"$existing_urls"
    fi

    for server_name in "${!MCPO_SERVERS[@]}"; do
        local description="${MCPO_SERVERS[$server_name]}"
        local server_url="${MCPO_URL}/${server_name}"

        if echo "$existing_urls" | grep -qF "$server_url" 2>/dev/null; then
            print_info "${server_name}: already configured"
            ((existing_count++))
        else
            local server_json
            server_json=$(build_server_config "$server_name" "$description")
            all_servers=$(echo "$all_servers" | jq ". + [$server_json]" 2>/dev/null)
            print_info "${server_name}: will be added"
            ((added_count++))
        fi
    done

    if [[ $DRY_RUN == "true" ]]; then
        print_section "Dry Run - Configuration Preview"
        echo "$all_servers" | jq '.'
        print_info "Would add ${added_count} new servers, ${existing_count} already exist"
        return 0
    fi

    if [[ $added_count -eq 0 ]]; then
        print_success "All MCPO servers already configured (${existing_count} servers)"
        return 0
    fi

    print_step "Applying configuration (${added_count} new servers)..."
    local payload
    payload=$(jq -n --argjson connections "$all_servers" '{"TOOL_SERVER_CONNECTIONS": $connections}')

    local response
    response=$(curl -sS -f \
        -H "Authorization: Bearer ${token}" \
        -H "Content-Type: application/json" \
        -d "$payload" \
        "${OPENWEBUI_URL}/api/v1/configs/tool_servers" 2>/dev/null) || {
        print_error "Failed to update tool servers configuration"
        return 1
    }

    print_success "Configuration updated successfully"
    print_info "Added: ${added_count}, Existing: ${existing_count}"
}

verify_all_servers() {
    print_section "Verifying MCPO Server Endpoints"
    local success_count=0
    local fail_count=0

    for server_name in "${!MCPO_SERVERS[@]}"; do
        if verify_server_endpoint "$server_name"; then
            ((success_count++))
        else
            ((fail_count++))
        fi
    done

    echo ""
    print_section "Verification Summary"
    print_info "Successful: ${success_count}"
    [[ $fail_count -gt 0 ]] && print_error "Failed: ${fail_count}"

    [[ $fail_count -eq 0 ]] && return 0 || return 1
}

main() {
    if ! verify_mcpo; then
        exit 1
    fi

    if ! verify_all_servers; then
        print_warning "Some endpoints not accessible, but continuing..."
    fi

    if [[ $VERIFY_ONLY == "true" ]]; then
        print_success "Verification complete"
        exit 0
    fi

    local token
    if ! token=$(get_auth_token); then
        print_error "Could not obtain authentication token"
        print_info "Set OPENWEBUI_API_KEY or ensure admin session is active"
        exit 1
    fi

    if ! configure_servers "$token"; then
        exit 1
    fi

    print_section "Next Steps"
    echo "1. Refresh OpenWebUI in your browser"
    echo "2. In a chat, click the tools icon to see MCP servers"
    echo "3. Enable desired tools for your thesis workflow"
    echo ""
    print_info "MCPO documentation: http://${MCPO_HOST_URL}/docs"
}

main "$@"
