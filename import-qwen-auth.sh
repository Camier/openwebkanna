#!/bin/bash

###############################################################################
# Qwen OAuth Credentials Importer
# Imports ~/.qwen/oauth_creds.json into CLIProxyAPI auth-dir format
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

load_env_defaults() {
    local env_file=""
    local line=""
    local key=""
    local value=""

    if [ -f "$SCRIPT_DIR/.env" ]; then
        env_file="$SCRIPT_DIR/.env"
    elif [ -f ".env" ]; then
        env_file=".env"
    fi

    if [ -z "$env_file" ]; then
        return 0
    fi

    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%$'\r'}"
        line="${line#"${line%%[![:space:]]*}"}"
        [ -z "$line" ] && continue
        [[ "$line" = \#* ]] && continue
        [[ "$line" != *=* ]] && continue

        key="${line%%=*}"
        value="${line#*=}"
        key="$(printf "%s" "$key" | tr -d '[:space:]')"

        [[ "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        if [ -n "${!key+x}" ]; then
            continue
        fi

        if [[ "$value" == \"*\" ]] && [[ "$value" == *\" ]]; then
            value="${value:1:${#value}-2}"
        elif [[ "$value" == \'*\' ]] && [[ "$value" == *\' ]]; then
            value="${value:1:${#value}-2}"
        fi

        printf -v "$key" "%s" "$value"
        export "$key"
    done < "$env_file"
}

load_env_defaults

# Configuration
CLIPROXYAPI_AUTH_DIR="${CLIPROXYAPI_AUTH_DIR:-./cliproxyapi/auth}"
QWEN_SOURCE_FILE="${CLIPROXYAPI_QWEN_SOURCE_FILE:-$HOME/.qwen/oauth_creds.json}"
QWEN_TARGET_FILE="${CLIPROXYAPI_QWEN_TARGET_FILE:-qwen-cli.json}"
QWEN_EMAIL="${CLIPROXYAPI_QWEN_EMAIL:-qwen-cli}"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

resolve_path() {
    local candidate="$1"
    if [[ "$candidate" = /* ]]; then
        printf "%s" "$candidate"
    else
        printf "%s/%s" "$SCRIPT_DIR" "$candidate"
    fi
}

expand_home_path() {
    local candidate="$1"
    if [ "${candidate:0:2}" = "~/" ]; then
        printf "%s/%s" "$HOME" "${candidate:2}"
        return
    fi
    if [ "$candidate" = "~" ]; then
        printf "%s" "$HOME"
        return
    fi
    printf "%s" "$candidate"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     Qwen OAuth Importer                                   ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}" >&2
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

show_usage() {
    cat <<'EOF_USAGE'
Usage: ./import-qwen-auth.sh [--help]

Environment variables:
  CLIPROXYAPI_AUTH_DIR      Target auth-dir (default: ./cliproxyapi/auth)
  CLIPROXYAPI_QWEN_SOURCE_FILE
                            Source file (default: ~/.qwen/oauth_creds.json)
  CLIPROXYAPI_QWEN_TARGET_FILE
                            Output file name under auth-dir (default: qwen-cli.json)
  CLIPROXYAPI_QWEN_EMAIL    Label/email metadata (default: qwen-cli)
EOF_USAGE
}

main() {
    if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
        show_usage
        return 0
    fi

    print_header

    if ! command_exists jq; then
        print_error "jq is required"
        return 1
    fi
    if ! command_exists date; then
        print_error "date is required"
        return 1
    fi

    CLIPROXYAPI_AUTH_DIR="$(resolve_path "$CLIPROXYAPI_AUTH_DIR")"
    QWEN_SOURCE_FILE="$(expand_home_path "$QWEN_SOURCE_FILE")"
    local target_path="$CLIPROXYAPI_AUTH_DIR/$QWEN_TARGET_FILE"

    if [ ! -f "$QWEN_SOURCE_FILE" ]; then
        print_error "Qwen source file not found: $QWEN_SOURCE_FILE"
        return 1
    fi

    print_step "Reading source credentials"

    local access_token
    local refresh_token
    local resource_url
    local token_type
    local expiry_raw
    access_token="$(jq -r '.access_token // empty' "$QWEN_SOURCE_FILE")"
    refresh_token="$(jq -r '.refresh_token // empty' "$QWEN_SOURCE_FILE")"
    resource_url="$(jq -r '.resource_url // empty' "$QWEN_SOURCE_FILE")"
    token_type="$(jq -r '.token_type // "Bearer"' "$QWEN_SOURCE_FILE")"
    expiry_raw="$(jq -r '.expiry_date // empty' "$QWEN_SOURCE_FILE")"

    if [ -z "$access_token" ] || [ -z "$refresh_token" ] || [ -z "$resource_url" ] || [ -z "$expiry_raw" ]; then
        print_error "Source file is missing one of required fields: access_token, refresh_token, resource_url, expiry_date"
        return 1
    fi

    if ! [[ "$expiry_raw" =~ ^[0-9]+$ ]]; then
        print_error "expiry_date must be an integer epoch value, got: $expiry_raw"
        return 1
    fi

    local expiry_seconds="$expiry_raw"
    if [ "$expiry_raw" -gt 9999999999 ]; then
        expiry_seconds=$((expiry_raw / 1000))
    fi

    local expired_at
    local now_utc
    expired_at="$(date -u -d "@$expiry_seconds" +"%Y-%m-%dT%H:%M:%SZ")"
    now_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

    mkdir -p "$CLIPROXYAPI_AUTH_DIR"

    print_step "Writing CLIProxyAPI auth file"

    jq -n \
        --arg access_token "$access_token" \
        --arg refresh_token "$refresh_token" \
        --arg resource_url "$resource_url" \
        --arg token_type "$token_type" \
        --arg email "$QWEN_EMAIL" \
        --arg expired "$expired_at" \
        --arg last_refresh "$now_utc" \
        '{
          access_token: $access_token,
          refresh_token: $refresh_token,
          resource_url: $resource_url,
          token_type: $token_type,
          email: $email,
          type: "qwen",
          expired: $expired,
          last_refresh: $last_refresh
        }' > "$target_path"

    chmod 600 "$target_path"

    print_success "Imported Qwen credentials: $target_path"
    print_warning "Credentials imported from $QWEN_SOURCE_FILE"
    return 0
}

main "$@"
