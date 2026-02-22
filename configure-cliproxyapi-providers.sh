#!/bin/bash
#
# Generate a local (gitignored) CLIProxyAPI config from ~/.007 provider keys.
#
# Why:
# - cliproxyapi/config.yaml is tracked and must never contain real secrets.
# - This script writes cliproxyapi/config.local.yaml with real provider keys
#   sourced from ~/.007 and local API keys sourced from .env.
#
# What it configures (by default):
# - Z.ai (GLM-5) via OpenAI-compatible provider
# - MiniMax (MiniMax-M2.5) via OpenAI-compatible provider
# - OAuth aliases for Codex + Qwen if you have oauth creds under cliproxyapi/auth/
#
# Safety:
# - Never prints secret values.
# - Writes the output file with umask 077 (mode 600).
#

set -euo pipefail

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

# Color definitions (standard repository set)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     Configure CLIProxyAPI Providers (local config)         ║"
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

yaml_quote() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    printf '"%s"' "$value"
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
            sep = index(line, "=") ? "=" : (index(line, ":") ? ":" : "")
            if (sep == "") next
            split(line, parts, sep)
            key=parts[1]
            sub(/^[[:space:]]+/, "", key); sub(/[[:space:]]+$/, "", key)
            if (key != k) next
            val = substr(line, length(parts[1]) + 2)
            sub(/^[[:space:]]+/, "", val); sub(/[[:space:]]+$/, "", val)
            if ((val ~ /^".*"$/) || (val ~ /^\x27.*\x27$/)) {
                val = substr(val, 2, length(val)-2)
            }
            print val
            exit 0
        }
    ' "$dotfile"
}

show_usage() {
    cat <<'EOF_USAGE'
Usage: ./configure-cliproxyapi-providers.sh [--help]

Writes:
  cliproxyapi/config.local.yaml

Inputs:
  - ~/.007: ZAI_API_KEY, MINIMAX_API_KEY (or MINIMAX_CODING_PLAN_API_KEY)
  - .env:  CLIPROXYAPI_API_KEY (fallback OPENAI_API_KEY)

Optional env overrides:
  CLIPROXYAPI_LOCAL_CONFIG   Output path (default: ./cliproxyapi/config.local.yaml)
  ZAI_BASE_URL               Default: https://api.z.ai/api/coding/paas/v4
  MINIMAX_BASE_URL           Default: https://api.minimax.io/v1
EOF_USAGE
}

main() {
    if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
        show_usage
        return 0
    fi

    load_env_defaults
    print_header

    local out="${CLIPROXYAPI_LOCAL_CONFIG:-./cliproxyapi/config.local.yaml}"
    local local_key="${CLIPROXYAPI_API_KEY:-${OPENAI_API_KEY:-}}"

    if [ -z "$local_key" ]; then
        print_error "Missing local API key. Set CLIPROXYAPI_API_KEY (or OPENAI_API_KEY) in .env"
        return 1
    fi

    local zai_key minimax_key
    zai_key="$(read_dot007_key "ZAI_API_KEY" || true)"
    if [ -z "$zai_key" ]; then
        zai_key="$(read_dot007_key "ZAI_DEVPACK_API_KEY" || true)"
    fi
    minimax_key="$(read_dot007_key "MINIMAX_API_KEY" || true)"
    if [ -z "$minimax_key" ]; then
        minimax_key="$(read_dot007_key "MINIMAX_CODING_PLAN_API_KEY" || true)"
    fi

    if [ -z "$zai_key" ]; then
        print_error "Missing Z.ai key. Expected ZAI_API_KEY (or ZAI_DEVPACK_API_KEY) in ~/.007"
        return 1
    fi
    if [ -z "$minimax_key" ]; then
        print_error "Missing MiniMax key. Expected MINIMAX_API_KEY (or MINIMAX_CODING_PLAN_API_KEY) in ~/.007"
        return 1
    fi

    # Z.ai offers multiple OpenAI-compatible entrypoints. GLM Coding Plan requires
    # the dedicated coding API base. Using the general API often yields
    # 429 "Insufficient balance or no resource package" even with a valid key.
    local zai_base="${ZAI_BASE_URL:-https://api.z.ai/api/coding/paas/v4}"
    local minimax_base="${MINIMAX_BASE_URL:-https://api.minimax.io/v1}"

    print_step "Writing local config"
    print_success "Output: $out"
    print_success "Z.ai key: present (from ~/.007)"
    print_success "MiniMax key: present (from ~/.007)"
    print_success "Local API key: present (from .env)"

    umask 077
    mkdir -p "$(dirname "$out")"

    cat >"$out" <<EOF_CONFIG
host: ""
port: 8317

tls:
  enable: false
  cert: ""
  key: ""

remote-management:
  allow-remote: false
  secret-key: ""
  disable-control-panel: false

auth-dir: "/CLIProxyAPI/cliproxyapi/auth"

api-keys:
  - $(yaml_quote "$local_key")

debug: false
request-log: false
usage-statistics-enabled: false
logging-to-file: false
logs-max-total-size-mb: 0
error-logs-max-files: 10
proxy-url: ""
force-model-prefix: false
request-retry: 3
max-retry-interval: 30

quota-exceeded:
  switch-project: true
  switch-preview-model: true

routing:
  strategy: "round-robin"

openai-compatibility:
  - name: "zai"
    base-url: $(yaml_quote "$zai_base")
    api-key-entries:
      - api-key: $(yaml_quote "$zai_key")
    models:
      - name: "glm-5"
        alias: "glm-5"

  - name: "minimax"
    base-url: $(yaml_quote "$minimax_base")
    api-key-entries:
      - api-key: $(yaml_quote "$minimax_key")
    models:
      # MiniMax Coding Plan model per docs: MiniMax-M2.5.
      #
      # Preferred model ID for OpenWebUI callers:
      # - minimax/chat-elite
      # Legacy compatibility:
      # - minmax
      # - minimax/chat-thinking
      # - minimax/chat-quality
      - name: "MiniMax-M2.5"
        alias: "minmax"
      - name: "MiniMax-M2.5"
        alias: "minimax/chat-elite"
      - name: "MiniMax-M2.5"
        alias: "minimax/chat-thinking"
      - name: "MiniMax-M2.5"
        alias: "minimax/chat-quality"

# OAuth-backed providers (if you have credentials under cliproxyapi/auth/).
oauth-model-alias:
  codex:
    - name: "gpt-5.3-codex"
      alias: "gpt-5.3-codex"
  qwen:
    - name: "qwen3-coder-flash"
      alias: "qwen3-coder-flash"

ws-auth: false
nonstream-keepalive-interval: 0
EOF_CONFIG

    print_step "Next step"
    echo "Set this in .env for docker-compose variable substitution:"
    echo "  CLIPROXYAPI_CONFIG=$out"
    print_success "Done"
}

main "$@"
