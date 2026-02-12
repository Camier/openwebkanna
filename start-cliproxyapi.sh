#!/bin/bash

###############################################################################
# CLIProxyAPI Server Start Script
# Starts CLIProxyAPI as a managed background process
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
CLIPROXYAPI_ENABLED="${CLIPROXYAPI_ENABLED:-true}"
CLIPROXYAPI_CMD="${CLIPROXYAPI_CMD:-./bin/cliproxyapi}"
CLIPROXYAPI_CONFIG="${CLIPROXYAPI_CONFIG:-./cliproxyapi/config.yaml}"
CLIPROXYAPI_CONFIG_FLAG="${CLIPROXYAPI_CONFIG_FLAG:---config}"
CLIPROXYAPI_SETUP_SCRIPT="${CLIPROXYAPI_SETUP_SCRIPT:-./setup-cliproxyapi.sh}"
CLIPROXYAPI_INSTALL_ON_START="${CLIPROXYAPI_INSTALL_ON_START:-true}"
CLIPROXYAPI_VERSION="${CLIPROXYAPI_VERSION:-latest}"
CLIPROXYAPI_REBUILD="${CLIPROXYAPI_REBUILD:-false}"
CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
CLIPROXYAPI_BIND_HOST="${CLIPROXYAPI_BIND_HOST:-0.0.0.0}"
CLIPROXYAPI_HEALTH_PATH="${CLIPROXYAPI_HEALTH_PATH:-/}"
CLIPROXYAPI_HEALTH_RETRIES="${CLIPROXYAPI_HEALTH_RETRIES:-45}"
CLIPROXYAPI_STARTUP_WAIT_SECONDS="${CLIPROXYAPI_STARTUP_WAIT_SECONDS:-2}"
CLIPROXYAPI_PORT="${CLIPROXYAPI_PORT:-8317}"
CLIPROXYAPI_EXTRA_ARGS="${CLIPROXYAPI_EXTRA_ARGS:-}"
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-layra-cliproxyapi-key}"
CLIPROXYAPI_AUTH_DIR="${CLIPROXYAPI_AUTH_DIR:-./cliproxyapi/auth}"
CLIPROXYAPI_AUTH_DIR_PATH=""
CLIPROXYAPI_CONFIG_AUTOGEN="${CLIPROXYAPI_CONFIG_AUTOGEN:-true}"
CLIPROXYAPI_PROVIDER_MODE="${CLIPROXYAPI_PROVIDER_MODE:-claude-api-key}"
CLIPROXYAPI_UPSTREAM_NAME="${CLIPROXYAPI_UPSTREAM_NAME:-upstream}"
CLIPROXYAPI_UPSTREAM_BASE_URL="${CLIPROXYAPI_UPSTREAM_BASE_URL:-}"
CLIPROXYAPI_UPSTREAM_API_KEY="${CLIPROXYAPI_UPSTREAM_API_KEY:-}"
CLIPROXYAPI_UPSTREAM_MODEL="${CLIPROXYAPI_UPSTREAM_MODEL:-}"
CLIPROXYAPI_UPSTREAM_MODEL_ALIAS="${CLIPROXYAPI_UPSTREAM_MODEL_ALIAS:-}"
CLIPROXYAPI_OAUTH_ANTIGRAVITY_MODEL="${CLIPROXYAPI_OAUTH_ANTIGRAVITY_MODEL:-}"
CLIPROXYAPI_OAUTH_ANTIGRAVITY_ALIAS="${CLIPROXYAPI_OAUTH_ANTIGRAVITY_ALIAS:-}"
CLIPROXYAPI_OAUTH_CODEX_MODEL="${CLIPROXYAPI_OAUTH_CODEX_MODEL:-}"
CLIPROXYAPI_OAUTH_CODEX_ALIAS="${CLIPROXYAPI_OAUTH_CODEX_ALIAS:-}"
CLIPROXYAPI_OAUTH_QWEN_MODEL="${CLIPROXYAPI_OAUTH_QWEN_MODEL:-}"
CLIPROXYAPI_OAUTH_QWEN_ALIAS="${CLIPROXYAPI_OAUTH_QWEN_ALIAS:-}"
CLIPROXYAPI_OAUTH_KIMI_MODEL="${CLIPROXYAPI_OAUTH_KIMI_MODEL:-}"
CLIPROXYAPI_OAUTH_KIMI_ALIAS="${CLIPROXYAPI_OAUTH_KIMI_ALIAS:-}"
CLIPROXYAPI_DISABLE_MANAGEMENT_PANEL="${CLIPROXYAPI_DISABLE_MANAGEMENT_PANEL:-true}"

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

normalize_base_url() {
    local value="$1"
    printf "%s" "${value%/}"
}

normalize_path() {
    local value="$1"
    if [ -z "$value" ]; then
        printf "/"
        return
    fi
    if [[ "$value" != /* ]]; then
        value="/$value"
    fi
    printf "%s" "$value"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

executable_exists() {
    local executable="$1"
    if [[ "$executable" == */* ]]; then
        [ -x "$executable" ]
        return $?
    fi
    command_exists "$executable"
}

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

yaml_quote() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    printf "\"%s\"" "$value"
}

if [[ "$CLIPROXYAPI_CMD" == */* ]]; then
    CLIPROXYAPI_CMD="$(resolve_path "$CLIPROXYAPI_CMD")"
fi
CLIPROXYAPI_CONFIG="$(resolve_path "$CLIPROXYAPI_CONFIG")"
if [[ "$CLIPROXYAPI_SETUP_SCRIPT" == */* ]]; then
    CLIPROXYAPI_SETUP_SCRIPT="$(resolve_path "$CLIPROXYAPI_SETUP_SCRIPT")"
fi
CLIPROXYAPI_LOG_FILE="$(resolve_path "${CLIPROXYAPI_LOG_FILE:-logs/cliproxyapi.log}")"
CLIPROXYAPI_PID_FILE="$(resolve_path "${CLIPROXYAPI_PID_FILE:-.cliproxyapi.pid}")"
if [[ "$CLIPROXYAPI_AUTH_DIR" = /* ]]; then
    CLIPROXYAPI_AUTH_DIR_PATH="$CLIPROXYAPI_AUTH_DIR"
else
    CLIPROXYAPI_AUTH_DIR_PATH="$(resolve_path "$CLIPROXYAPI_AUTH_DIR")"
fi
CLIPROXYAPI_BASE_URL="$(normalize_base_url "$CLIPROXYAPI_BASE_URL")"
CLIPROXYAPI_HEALTH_PATH="$(normalize_path "$CLIPROXYAPI_HEALTH_PATH")"

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     CLIProxyAPI Server Starter                            ║"
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

print_info() {
    echo -e "${CYAN}ℹ $1${NC}"
}

prepare_provider_defaults() {
    CLIPROXYAPI_PROVIDER_MODE="$(printf "%s" "$CLIPROXYAPI_PROVIDER_MODE" | tr '[:upper:]' '[:lower:]')"

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "claude" ]; then
        CLIPROXYAPI_PROVIDER_MODE="claude-api-key"
    fi
    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "openai" ]; then
        CLIPROXYAPI_PROVIDER_MODE="openai-compatibility"
    fi

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "oauth" ]; then
        if [ -z "$CLIPROXYAPI_OAUTH_ANTIGRAVITY_MODEL" ]; then
            CLIPROXYAPI_OAUTH_ANTIGRAVITY_MODEL="gemini-3-pro-high"
        fi
        if [ -z "$CLIPROXYAPI_OAUTH_ANTIGRAVITY_ALIAS" ]; then
            CLIPROXYAPI_OAUTH_ANTIGRAVITY_ALIAS="antigravity-oauth"
        fi
        if [ -z "$CLIPROXYAPI_OAUTH_CODEX_MODEL" ]; then
            CLIPROXYAPI_OAUTH_CODEX_MODEL="gpt-5-codex"
        fi
        if [ -z "$CLIPROXYAPI_OAUTH_CODEX_ALIAS" ]; then
            CLIPROXYAPI_OAUTH_CODEX_ALIAS="openai-codex"
        fi
        if [ -z "$CLIPROXYAPI_OAUTH_QWEN_MODEL" ]; then
            CLIPROXYAPI_OAUTH_QWEN_MODEL="qwen3-coder-plus"
        fi
        if [ -z "$CLIPROXYAPI_OAUTH_QWEN_ALIAS" ]; then
            CLIPROXYAPI_OAUTH_QWEN_ALIAS="qwen-cli"
        fi
        if [ -z "$CLIPROXYAPI_OAUTH_KIMI_MODEL" ]; then
            CLIPROXYAPI_OAUTH_KIMI_MODEL="kimi-k2.5"
        fi
        if [ -z "$CLIPROXYAPI_OAUTH_KIMI_ALIAS" ]; then
            CLIPROXYAPI_OAUTH_KIMI_ALIAS="kimi-cli"
        fi
    fi

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "claude-api-key" ]; then
        if [ -z "$CLIPROXYAPI_UPSTREAM_BASE_URL" ] && [ -n "${ANTHROPIC_BASE_URL:-}" ]; then
            CLIPROXYAPI_UPSTREAM_BASE_URL="$ANTHROPIC_BASE_URL"
        fi
        if [ -z "$CLIPROXYAPI_UPSTREAM_API_KEY" ] && [ -n "${ANTHROPIC_AUTH_TOKEN:-}" ]; then
            CLIPROXYAPI_UPSTREAM_API_KEY="$ANTHROPIC_AUTH_TOKEN"
        fi
        if [ -z "$CLIPROXYAPI_UPSTREAM_MODEL" ] && [ -n "${ANTHROPIC_DEFAULT_SONNET_MODEL:-}" ]; then
            CLIPROXYAPI_UPSTREAM_MODEL="$ANTHROPIC_DEFAULT_SONNET_MODEL"
        fi
        if [ -z "$CLIPROXYAPI_UPSTREAM_MODEL" ]; then
            CLIPROXYAPI_UPSTREAM_MODEL="glm-4.7"
        fi
    fi

    if [ -z "$CLIPROXYAPI_UPSTREAM_MODEL_ALIAS" ]; then
        CLIPROXYAPI_UPSTREAM_MODEL_ALIAS="$CLIPROXYAPI_UPSTREAM_MODEL"
    fi
}

validate_provider_config() {
    case "$CLIPROXYAPI_PROVIDER_MODE" in
        claude-api-key|openai-compatibility|oauth)
            ;;
        *)
            print_error "Unsupported CLIPROXYAPI_PROVIDER_MODE: $CLIPROXYAPI_PROVIDER_MODE"
            print_info "Supported values: claude-api-key, openai-compatibility, oauth"
            return 1
            ;;
    esac

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "oauth" ]; then
        return 0
    fi

    if [ -z "$CLIPROXYAPI_UPSTREAM_BASE_URL" ]; then
        print_error "CLIPROXYAPI_UPSTREAM_BASE_URL is required"
        return 1
    fi
    if [ -z "$CLIPROXYAPI_UPSTREAM_API_KEY" ]; then
        print_error "CLIPROXYAPI_UPSTREAM_API_KEY is required"
        print_info "Or export ANTHROPIC_AUTH_TOKEN when using claude-api-key mode"
        return 1
    fi
    if [ -z "$CLIPROXYAPI_UPSTREAM_MODEL" ]; then
        print_error "CLIPROXYAPI_UPSTREAM_MODEL is required"
        return 1
    fi
    return 0
}

autogenerate_config() {
    mkdir -p "$(dirname "$CLIPROXYAPI_CONFIG")"
    mkdir -p "$CLIPROXYAPI_AUTH_DIR_PATH"

    local disable_panel
    if is_true "$CLIPROXYAPI_DISABLE_MANAGEMENT_PANEL"; then
        disable_panel="true"
    else
        disable_panel="false"
    fi

    print_step "Generating CLIProxyAPI config"
    print_info "Provider mode: $CLIPROXYAPI_PROVIDER_MODE"
    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "oauth" ]; then
        print_info "Antigravity alias: $CLIPROXYAPI_OAUTH_ANTIGRAVITY_ALIAS -> $CLIPROXYAPI_OAUTH_ANTIGRAVITY_MODEL"
        print_info "Codex alias: $CLIPROXYAPI_OAUTH_CODEX_ALIAS -> $CLIPROXYAPI_OAUTH_CODEX_MODEL"
        print_info "Qwen alias: $CLIPROXYAPI_OAUTH_QWEN_ALIAS -> $CLIPROXYAPI_OAUTH_QWEN_MODEL"
        print_info "Kimi alias: $CLIPROXYAPI_OAUTH_KIMI_ALIAS -> $CLIPROXYAPI_OAUTH_KIMI_MODEL"
    else
        print_info "Upstream base URL: $CLIPROXYAPI_UPSTREAM_BASE_URL"
        print_info "Model alias: $CLIPROXYAPI_UPSTREAM_MODEL_ALIAS"
    fi

cat > "$CLIPROXYAPI_CONFIG" <<EOF_CONFIG
host: $(yaml_quote "$CLIPROXYAPI_BIND_HOST")
port: $CLIPROXYAPI_PORT

tls:
  enable: false
  cert: ""
  key: ""

remote-management:
  allow-remote: false
  secret-key: ""
  disable-control-panel: $disable_panel

auth-dir: $(yaml_quote "$CLIPROXYAPI_AUTH_DIR")

api-keys:
  - $(yaml_quote "$CLIPROXYAPI_API_KEY")

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
EOF_CONFIG

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "claude-api-key" ]; then
        cat >> "$CLIPROXYAPI_CONFIG" <<EOF_CONFIG

claude-api-key:
  - api-key: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_API_KEY")
    base-url: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_BASE_URL")
    models:
      - name: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_MODEL")
        alias: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_MODEL_ALIAS")
EOF_CONFIG
    elif [ "$CLIPROXYAPI_PROVIDER_MODE" = "openai-compatibility" ]; then
        cat >> "$CLIPROXYAPI_CONFIG" <<EOF_CONFIG

openai-compatibility:
  - name: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_NAME")
    base-url: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_BASE_URL")
    api-key-entries:
      - api-key: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_API_KEY")
    models:
      - name: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_MODEL")
        alias: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_MODEL_ALIAS")
EOF_CONFIG
    else
        cat >> "$CLIPROXYAPI_CONFIG" <<EOF_CONFIG

oauth-model-alias:
  antigravity:
    - name: $(yaml_quote "$CLIPROXYAPI_OAUTH_ANTIGRAVITY_MODEL")
      alias: $(yaml_quote "$CLIPROXYAPI_OAUTH_ANTIGRAVITY_ALIAS")
  codex:
    - name: $(yaml_quote "$CLIPROXYAPI_OAUTH_CODEX_MODEL")
      alias: $(yaml_quote "$CLIPROXYAPI_OAUTH_CODEX_ALIAS")
  qwen:
    - name: $(yaml_quote "$CLIPROXYAPI_OAUTH_QWEN_MODEL")
      alias: $(yaml_quote "$CLIPROXYAPI_OAUTH_QWEN_ALIAS")
  kimi:
    - name: $(yaml_quote "$CLIPROXYAPI_OAUTH_KIMI_MODEL")
      alias: $(yaml_quote "$CLIPROXYAPI_OAUTH_KIMI_ALIAS")
EOF_CONFIG
    fi

    cat >> "$CLIPROXYAPI_CONFIG" <<'EOF_CONFIG'

ws-auth: false
nonstream-keepalive-interval: 0
EOF_CONFIG

    print_success "Config written: $CLIPROXYAPI_CONFIG"
}

install_official_binary() {
    if [[ "$CLIPROXYAPI_CMD" != */* ]]; then
        print_error "CLIPROXYAPI_CMD must be a filesystem path to install binary: $CLIPROXYAPI_CMD"
        return 1
    fi

    if [ ! -x "$CLIPROXYAPI_SETUP_SCRIPT" ]; then
        print_error "Setup script missing or not executable: $CLIPROXYAPI_SETUP_SCRIPT"
        return 1
    fi

    print_step "Installing official CLIProxyAPI binary"
    if ! CLIPROXYAPI_INSTALL_PATH="$CLIPROXYAPI_CMD" \
        CLIPROXYAPI_VERSION="$CLIPROXYAPI_VERSION" \
        CLIPROXYAPI_FORCE_REINSTALL="$CLIPROXYAPI_REBUILD" \
        "$CLIPROXYAPI_SETUP_SCRIPT"; then
        return 1
    fi

    if ! executable_exists "$CLIPROXYAPI_CMD"; then
        print_error "Installed binary is not executable: $CLIPROXYAPI_CMD"
        return 1
    fi

    print_success "CLIProxyAPI binary ready"
}

ensure_binary() {
    if executable_exists "$CLIPROXYAPI_CMD"; then
        return 0
    fi

    if ! is_true "$CLIPROXYAPI_INSTALL_ON_START"; then
        print_error "CLIProxyAPI binary missing: $CLIPROXYAPI_CMD"
        print_info "Set CLIPROXYAPI_INSTALL_ON_START=true or run ./setup-cliproxyapi.sh"
        return 1
    fi

    install_official_binary
}

wait_for_health() {
    local health_url="${CLIPROXYAPI_BASE_URL}${CLIPROXYAPI_HEALTH_PATH}"
    local attempt=1

    print_step "Waiting for health endpoint"
    print_info "URL: $health_url"

    while [ "$attempt" -le "$CLIPROXYAPI_HEALTH_RETRIES" ]; do
        if curl -sS -m 5 -o /dev/null -w "%{http_code}" "$health_url" | grep -q '^200$'; then
            print_success "CLIProxyAPI is healthy"
            return 0
        fi
        sleep 1
        attempt=$((attempt + 1))
    done

    print_error "CLIProxyAPI health check timed out"
    return 1
}

show_usage() {
    cat <<'EOF_USAGE'
Usage: ./start-cliproxyapi.sh [--help]

Environment variables:
  CLIPROXYAPI_ENABLED                Enable/disable lifecycle (must be true)
  CLIPROXYAPI_CMD                    Binary path (default: ./bin/cliproxyapi)
  CLIPROXYAPI_CONFIG                 Config path (default: ./cliproxyapi/config.yaml)
  CLIPROXYAPI_SETUP_SCRIPT           Installer script path (default: ./setup-cliproxyapi.sh)
  CLIPROXYAPI_INSTALL_ON_START       Install binary if missing (default: true)
  CLIPROXYAPI_VERSION                Release tag or latest (default: latest)
  CLIPROXYAPI_REBUILD                Force reinstall on start (default: false)
  CLIPROXYAPI_CONFIG_AUTOGEN         Generate config from env (default: true)
  CLIPROXYAPI_PROVIDER_MODE          claude-api-key, openai-compatibility, or oauth
  CLIPROXYAPI_BIND_HOST              Listener host written to config (default: 0.0.0.0)
  CLIPROXYAPI_UPSTREAM_BASE_URL      Required upstream base URL
  CLIPROXYAPI_UPSTREAM_API_KEY       Required upstream key (or ANTHROPIC_AUTH_TOKEN fallback)
  CLIPROXYAPI_UPSTREAM_MODEL         Required upstream model id
  CLIPROXYAPI_UPSTREAM_MODEL_ALIAS   Model alias exposed to clients
  CLIPROXYAPI_UPSTREAM_NAME          Provider name for openai-compatibility mode
  CLIPROXYAPI_OAUTH_*                OAuth alias mapping env vars for antigravity/codex/qwen/kimi
  CLIPROXYAPI_API_KEY                Local key accepted by CLIProxyAPI
EOF_USAGE
}

main() {
    if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
        show_usage
        return 0
    fi

    print_header

    print_step "Checking prerequisites"
    if ! command_exists curl; then
        print_error "curl is required"
        return 1
    fi
    if ! command_exists lsof; then
        print_error "lsof is required"
        return 1
    fi

    if ! is_true "$CLIPROXYAPI_ENABLED"; then
        print_error "CLIProxyAPI lifecycle disabled (CLIPROXYAPI_ENABLED=$CLIPROXYAPI_ENABLED)"
        return 1
    fi

    ensure_binary

    prepare_provider_defaults

    if is_true "$CLIPROXYAPI_CONFIG_AUTOGEN"; then
        if ! validate_provider_config; then
            return 1
        fi
        autogenerate_config
    fi

    if [ ! -f "$CLIPROXYAPI_CONFIG" ]; then
        print_error "Config file not found: $CLIPROXYAPI_CONFIG"
        return 1
    fi

    if [ -f "$CLIPROXYAPI_PID_FILE" ]; then
        local existing_pid
        existing_pid="$(cat "$CLIPROXYAPI_PID_FILE")"
        if [ -n "$existing_pid" ] && kill -0 "$existing_pid" 2>/dev/null; then
            print_success "CLIProxyAPI already running (PID: $existing_pid)"
            wait_for_health
            return 0
        fi
        print_warning "Removing stale PID file: $CLIPROXYAPI_PID_FILE"
        rm -f "$CLIPROXYAPI_PID_FILE"
    fi

    local listener
    listener="$(lsof -nP -iTCP:"$CLIPROXYAPI_PORT" -sTCP:LISTEN 2>/dev/null | tail -n +2 || true)"
    if [ -n "$listener" ]; then
        print_error "Port $CLIPROXYAPI_PORT is already in use"
        echo "$listener"
        return 1
    fi

    mkdir -p "$(dirname "$CLIPROXYAPI_LOG_FILE")"

    local -a command=("$CLIPROXYAPI_CMD" "$CLIPROXYAPI_CONFIG_FLAG" "$CLIPROXYAPI_CONFIG")
    if [ -n "$CLIPROXYAPI_EXTRA_ARGS" ]; then
        # shellcheck disable=SC2206
        local -a extra_args=( $CLIPROXYAPI_EXTRA_ARGS )
        command+=("${extra_args[@]}")
    fi

    print_step "Starting CLIProxyAPI"
    print_info "Command: ${command[*]}"
    print_info "Log file: $CLIPROXYAPI_LOG_FILE"

    nohup "${command[@]}" >"$CLIPROXYAPI_LOG_FILE" 2>&1 &
    local pid=$!
    echo "$pid" >"$CLIPROXYAPI_PID_FILE"

    sleep "$CLIPROXYAPI_STARTUP_WAIT_SECONDS"
    if ! kill -0 "$pid" 2>/dev/null; then
        print_error "CLIProxyAPI exited immediately"
        rm -f "$CLIPROXYAPI_PID_FILE"
        tail -n 40 "$CLIPROXYAPI_LOG_FILE" 2>/dev/null || true
        return 1
    fi

    if ! wait_for_health; then
        tail -n 40 "$CLIPROXYAPI_LOG_FILE" 2>/dev/null || true
        return 1
    fi

    print_success "CLIProxyAPI started (PID: $pid)"
    return 0
}

main "$@"
