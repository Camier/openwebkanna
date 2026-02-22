#!/bin/bash
#
# Configure OpenWebUI to use Ollama Cloud (remote Ollama host: https://ollama.com).
# - Reads the API key from ~/.007 (OLLAMA_CLOUD_API_KEY or OLLAMA_API_KEY)
# - Signs into OpenWebUI using OPENWEBUI_SIGNIN_EMAIL/OPENWEBUI_SIGNIN_PASSWORD (from .env)
# - Updates OpenWebUI Ollama connection via POST /ollama/config/update
#
# Notes:
# - This does NOT write secrets into the repo.
# - OpenWebUI stores Ollama connection config (incl. API key) in its persistent config DB.
#
# Docs:
# - Ollama Cloud API access: https://docs.ollama.com/cloud
# - OpenWebUI Ollama connection: https://docs.openwebui.com/getting-started/quick-start/connect-a-provider/starting-with-ollama/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source all library modules
source "${SCRIPT_DIR}/lib/init.sh"
cd "$SCRIPT_DIR"

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
        [[ $line == \#* ]] && continue
        [[ $line != *=* ]] && continue

        key="${line%%=*}"
        value="${line#*=}"
        key="$(printf "%s" "$key" | tr -d '[:space:]')"

        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        if [ -n "${!key+x}" ]; then
            continue
        fi

        if [[ $value == \"*\" ]] && [[ $value == *\" ]]; then
            value="${value:1:${#value}-2}"
        elif [[ $value == \'*\' ]] && [[ $value == *\' ]]; then
            value="${value:1:${#value}-2}"
        fi

        printf -v "$key" "%s" "$value"
        export "${key?}"
    done <"$env_file"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     Configure Ollama Cloud (OpenWebUI)                     ║"
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
Usage: ./configure-ollama-cloud.sh [--help]

What it does:
  - Reads Ollama Cloud API key from ~/.007 (OLLAMA_CLOUD_API_KEY or OLLAMA_API_KEY)
  - Signs into OpenWebUI (OPENWEBUI_SIGNIN_EMAIL / OPENWEBUI_SIGNIN_PASSWORD)
  - Sets Ollama base URL to https://ollama.com and stores the key for the connection

Env vars (optional):
  OPENWEBUI_URL            Base URL for OpenWebUI (default: http://localhost:${WEBUI_PORT:-3000})
  OPENWEBUI_TOKEN          If set, used instead of signing in (Bearer token / JWT)
  OLLAMA_CLOUD_BASE_URL    Remote Ollama host (default: https://ollama.com)
  OPENWEBUI_SIGNIN_EMAIL   Sign-in email (default: from .env or test@example.com)
  OPENWEBUI_SIGNIN_PASSWORD
                           Sign-in password (default: from .env or admin)
EOF_USAGE
}

read_dot007_key() {
    local wanted="$1"
    local dotfile="$HOME/.007"

    [ -f "$dotfile" ] || return 1

    # Supports lines like:
    #   KEY=value
    #   export KEY=value
    #   KEY: value
    awk -v k="$wanted" '
        /^[[:space:]]*#/ { next }
        /^[[:space:]]*$/ { next }
        {
            line=$0
            sub(/^[[:space:]]*export[[:space:]]+/, "", line)
            # Split on first "=" or ":".
            sep = index(line, "=") ? "=" : (index(line, ":") ? ":" : "")
            if (sep == "") next
            split(line, parts, sep)
            key=parts[1]
            sub(/^[[:space:]]+/, "", key); sub(/[[:space:]]+$/, "", key)
            if (key != k) next
            val = substr(line, length(parts[1]) + 2)
            sub(/^[[:space:]]+/, "", val); sub(/[[:space:]]+$/, "", val)
            # Strip surrounding quotes.
            if ((val ~ /^".*"$/) || (val ~ /^\x27.*\x27$/)) {
                val = substr(val, 2, length(val)-2)
            }
            print val
            exit 0
        }
    ' "$dotfile"
}

main() {
    if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
        show_usage
        return 0
    fi

    load_env_defaults

    # Colors
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    CYAN='\033[0;36m'
    NC='\033[0m'
    BOLD='\033[1m'

    print_header

    if ! command_exists curl; then
        print_error "curl is required"
        return 1
    fi
    if ! command_exists jq; then
        print_error "jq is required"
        return 1
    fi

    local openwebui_url="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
    local ollama_cloud_base_url="${OLLAMA_CLOUD_BASE_URL:-${OLLAMA_BASE_URL:-https://ollama.com}}"
    local signin_email="${OPENWEBUI_SIGNIN_EMAIL:-test@example.com}"
    local signin_password="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"

    local api_key=""
    api_key="$(read_dot007_key "OLLAMA_CLOUD_API_KEY" || true)"
    if [ -z "$api_key" ]; then
        api_key="$(read_dot007_key "OLLAMA_API_KEY" || true)"
    fi
    if [ -z "$api_key" ]; then
        print_error "Missing Ollama Cloud API key. Expected OLLAMA_CLOUD_API_KEY or OLLAMA_API_KEY in ~/.007"
        return 1
    fi

    local token="${OPENWEBUI_TOKEN:-}"
    local signin_payload=""
    local signin_resp=""
    local update_payload=""
    umask 077
    trap 'rm -f "$signin_payload" "$signin_resp" "$update_payload"' EXIT

    if [ -z "$token" ]; then
        print_step "Signing into OpenWebUI"
        signin_payload="$(mktemp)"
        signin_resp="$(mktemp)"

        jq -n \
            --arg email "$signin_email" \
            --arg password "$signin_password" \
            '{email:$email, password:$password}' >"$signin_payload"

        local signin_code
        signin_code="$(
            curl -sS \
                -o "$signin_resp" \
                -w '%{http_code}' \
                -X POST \
                "${openwebui_url%/}/api/v1/auths/signin" \
                -H 'Content-Type: application/json' \
                --data-binary "@$signin_payload"
        )"

        if [ "$signin_code" != "200" ]; then
            print_error "OpenWebUI signin failed (HTTP $signin_code)"
            return 1
        fi

        token="$(jq -r '.token // empty' "$signin_resp")"
        if [ -z "$token" ]; then
            print_error "OpenWebUI signin did not return a token"
            return 1
        fi
    else
        print_step "Using OPENWEBUI_TOKEN (skipping signin)"
    fi

    print_step "Updating Ollama connection (Cloud)"
    update_payload="$(mktemp)"

    jq -n \
        --arg url "${ollama_cloud_base_url%/}" \
        --arg key "$api_key" \
        '{
            ENABLE_OLLAMA_API: true,
            OLLAMA_BASE_URLS: [$url],
            OLLAMA_API_CONFIGS: {
                "0": {
                    enable: true,
                    key: $key,
                    connection_type: "cloud"
                }
            }
        }' >"$update_payload"

    local update_code
    update_code="$(
        curl -sS \
            -o /dev/null \
            -w '%{http_code}' \
            -X POST \
            "${openwebui_url%/}/ollama/config/update" \
            -H "Authorization: Bearer $token" \
            -H 'Content-Type: application/json' \
            --data-binary "@$update_payload"
    )"

    if [ "$update_code" != "200" ]; then
        print_error "Failed to update /ollama/config/update (HTTP $update_code)"
        return 1
    fi

    print_step "Verifying model list"
    local count
    count="$(
        curl -sS \
            -H "Authorization: Bearer $token" \
            "${openwebui_url%/}/ollama/api/tags" | jq -r '.models | length'
    )"

    if ! [[ $count =~ ^[0-9]+$ ]]; then
        print_warning "Could not parse model count from /ollama/api/tags"
    else
        print_success "Ollama Cloud connected. Models visible: $count"
    fi

    print_success "Done"
}

main "$@"
