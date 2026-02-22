#!/bin/bash

###############################################################################
# CLIProxyAPI Server Start Script
# Starts CLIProxyAPI as a managed background process
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

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
CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-${OPENAI_API_KEY:-}}"
CLIPROXYAPI_ALLOW_INSECURE_DEFAULT_KEY="${CLIPROXYAPI_ALLOW_INSECURE_DEFAULT_KEY:-false}"
CLIPROXYAPI_AUTH_DIR="${CLIPROXYAPI_AUTH_DIR:-./cliproxyapi/auth}"
CLIPROXYAPI_AUTH_DIR_PATH=""
CLIPROXYAPI_CONFIG_AUTOGEN="${CLIPROXYAPI_CONFIG_AUTOGEN:-true}"
CLIPROXYAPI_CONFIG_AUTOGEN_OVERWRITE="${CLIPROXYAPI_CONFIG_AUTOGEN_OVERWRITE:-false}"
CLIPROXYAPI_PROVIDER_MODE="${CLIPROXYAPI_PROVIDER_MODE:-claude-api-key}"
CLIPROXYAPI_UPSTREAM_NAME="${CLIPROXYAPI_UPSTREAM_NAME:-upstream}"
CLIPROXYAPI_UPSTREAM_BASE_URL="${CLIPROXYAPI_UPSTREAM_BASE_URL:-}"
CLIPROXYAPI_UPSTREAM_API_KEY="${CLIPROXYAPI_UPSTREAM_API_KEY:-}"
CLIPROXYAPI_UPSTREAM_MODEL="${CLIPROXYAPI_UPSTREAM_MODEL:-}"
CLIPROXYAPI_UPSTREAM_MODEL_ALIAS="${CLIPROXYAPI_UPSTREAM_MODEL_ALIAS:-}"
CLIPROXYAPI_OAUTH_CODEX_MODEL="${CLIPROXYAPI_OAUTH_CODEX_MODEL:-}"
CLIPROXYAPI_OAUTH_CODEX_ALIAS="${CLIPROXYAPI_OAUTH_CODEX_ALIAS:-}"
CLIPROXYAPI_OAUTH_QWEN_MODEL="${CLIPROXYAPI_OAUTH_QWEN_MODEL:-}"
CLIPROXYAPI_OAUTH_QWEN_ALIAS="${CLIPROXYAPI_OAUTH_QWEN_ALIAS:-}"
CLIPROXYAPI_OAUTH_KIMI_MODEL="${CLIPROXYAPI_OAUTH_KIMI_MODEL:-}"
CLIPROXYAPI_OAUTH_KIMI_ALIAS="${CLIPROXYAPI_OAUTH_KIMI_ALIAS:-}"
CLIPROXYAPI_DISABLE_MANAGEMENT_PANEL="${CLIPROXYAPI_DISABLE_MANAGEMENT_PANEL:-true}"
CLIPROXYAPI_VERIFY_MODELS_ON_START="${CLIPROXYAPI_VERIFY_MODELS_ON_START:-true}"
CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY="${CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY:-true}"
CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH="${CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH:-true}"
CLIPROXYAPI_SYNC_CONFIG_API_KEY="${CLIPROXYAPI_SYNC_CONFIG_API_KEY:-true}"

is_weak_local_key() {
    local value="$1"
    [ -z "$value" ] && return 0

    case "$value" in
        layra-cliproxyapi-key | change-me-cliproxyapi-key | replace-with-strong-local-api-key | change-me | changeme | default | test | test-key)
            return 0
            ;;
    esac

    return 1
}

extract_models_count() {
    local payload_file="$1"

    if command_exists jq; then
        jq -er '.data | if type == "array" then length else empty end' "$payload_file" 2>/dev/null || return 1
        return 0
    fi

    if ! grep -q '"data"[[:space:]]*:' "$payload_file"; then
        return 1
    fi

    grep -o '"id"[[:space:]]*:[[:space:]]*"[^"]*"' "$payload_file" | wc -l | tr -d '[:space:]'
}

yaml_quote() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    printf '"%s"' "$value"
}

executable_exists() {
    local executable="$1"
    if [[ $executable == */* ]]; then
        [ -x "$executable" ]
        return $?
    fi
    command_exists "$executable"
}

if [[ $CLIPROXYAPI_CMD == */* ]]; then
    CLIPROXYAPI_CMD="$(resolve_path "$CLIPROXYAPI_CMD")"
fi
CLIPROXYAPI_CONFIG="$(resolve_path "$CLIPROXYAPI_CONFIG")"
if [[ $CLIPROXYAPI_SETUP_SCRIPT == */* ]]; then
    CLIPROXYAPI_SETUP_SCRIPT="$(resolve_path "$CLIPROXYAPI_SETUP_SCRIPT")"
fi
CLIPROXYAPI_LOG_FILE="$(resolve_path "${CLIPROXYAPI_LOG_FILE:-logs/cliproxyapi.log}")"
CLIPROXYAPI_PID_FILE="$(resolve_path "${CLIPROXYAPI_PID_FILE:-.cliproxyapi.pid}")"
if [[ $CLIPROXYAPI_AUTH_DIR == /* ]]; then
    CLIPROXYAPI_AUTH_DIR_PATH="$CLIPROXYAPI_AUTH_DIR"
else
    CLIPROXYAPI_AUTH_DIR_PATH="$(resolve_path "$CLIPROXYAPI_AUTH_DIR")"
fi
CLIPROXYAPI_BASE_URL="$(normalize_base_url "$CLIPROXYAPI_BASE_URL")"
CLIPROXYAPI_HEALTH_PATH="$(normalize_path "$CLIPROXYAPI_HEALTH_PATH")"

prepare_provider_defaults() {
    CLIPROXYAPI_PROVIDER_MODE="$(printf "%s" "$CLIPROXYAPI_PROVIDER_MODE" | tr '[:upper:]' '[:lower:]')"

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "claude" ]; then
        CLIPROXYAPI_PROVIDER_MODE="claude-api-key"
    fi
    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "openai" ]; then
        CLIPROXYAPI_PROVIDER_MODE="openai-compatibility"
    fi

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "oauth" ]; then
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
    fi

    if [ -z "$CLIPROXYAPI_UPSTREAM_MODEL_ALIAS" ]; then
        CLIPROXYAPI_UPSTREAM_MODEL_ALIAS="$CLIPROXYAPI_UPSTREAM_MODEL"
    fi
}

validate_provider_config() {
    case "$CLIPROXYAPI_PROVIDER_MODE" in
        claude-api-key | openai-compatibility | oauth) ;;
        *)
            print_error "Unsupported CLIPROXYAPI_PROVIDER_MODE: $CLIPROXYAPI_PROVIDER_MODE"
            print_info "Supported values: claude-api-key, openai-compatibility, oauth"
            return 1
            ;;
    esac

    if [ "$CLIPROXYAPI_PROVIDER_MODE" = "oauth" ]; then
        if [ ! -d "$CLIPROXYAPI_AUTH_DIR_PATH" ]; then
            print_error "OAuth mode requires auth directory: $CLIPROXYAPI_AUTH_DIR_PATH"
            print_info "Run ./configure-cliproxyapi-oauth.sh to populate credentials"
            return 1
        fi
        if ! find "$CLIPROXYAPI_AUTH_DIR_PATH" -maxdepth 2 -type f -name '*.json' | grep -q .; then
            print_error "OAuth mode requires credential JSON files in $CLIPROXYAPI_AUTH_DIR_PATH"
            print_info "Run ./configure-cliproxyapi-oauth.sh to create/update provider credentials"
            return 1
        fi
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

validate_local_api_key() {
    local require_for_autogen="${1:-false}"

    if [ -z "$CLIPROXYAPI_API_KEY" ]; then
        if is_true "$require_for_autogen" || is_true "$CLIPROXYAPI_VERIFY_MODELS_ON_START"; then
            print_error "CLIPROXYAPI_API_KEY is required for config autogen and startup verification"
            return 1
        fi
        return 0
    fi

    if is_weak_local_key "$CLIPROXYAPI_API_KEY" && ! is_true "$CLIPROXYAPI_ALLOW_INSECURE_DEFAULT_KEY"; then
        print_error "CLIPROXYAPI_API_KEY uses a weak default value"
        print_info "Set a unique local key or set CLIPROXYAPI_ALLOW_INSECURE_DEFAULT_KEY=true (not recommended)"
        return 1
    fi

    return 0
}

autogenerate_config() {
    if [ -f "$CLIPROXYAPI_CONFIG" ] && ! is_true "$CLIPROXYAPI_CONFIG_AUTOGEN_OVERWRITE"; then
        print_error "Config already exists and overwrite is disabled: $CLIPROXYAPI_CONFIG"
        print_info "Set CLIPROXYAPI_CONFIG_AUTOGEN_OVERWRITE=true to regenerate"
        return 1
    fi

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
        print_info "Codex alias: $CLIPROXYAPI_OAUTH_CODEX_ALIAS -> $CLIPROXYAPI_OAUTH_CODEX_MODEL"
        print_info "Qwen alias: $CLIPROXYAPI_OAUTH_QWEN_ALIAS -> $CLIPROXYAPI_OAUTH_QWEN_MODEL"
        print_info "Kimi alias: $CLIPROXYAPI_OAUTH_KIMI_ALIAS -> $CLIPROXYAPI_OAUTH_KIMI_MODEL"
    else
        print_info "Upstream base URL: $CLIPROXYAPI_UPSTREAM_BASE_URL"
        print_info "Model alias: $CLIPROXYAPI_UPSTREAM_MODEL_ALIAS"
    fi

    cat >"$CLIPROXYAPI_CONFIG" <<EOF_CONFIG
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
        cat >>"$CLIPROXYAPI_CONFIG" <<EOF_CONFIG

claude-api-key:
  - api-key: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_API_KEY")
    base-url: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_BASE_URL")
    models:
      - name: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_MODEL")
        alias: $(yaml_quote "$CLIPROXYAPI_UPSTREAM_MODEL_ALIAS")
EOF_CONFIG
    elif [ "$CLIPROXYAPI_PROVIDER_MODE" = "openai-compatibility" ]; then
        cat >>"$CLIPROXYAPI_CONFIG" <<EOF_CONFIG

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
        cat >>"$CLIPROXYAPI_CONFIG" <<EOF_CONFIG

oauth-model-alias:
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

    cat >>"$CLIPROXYAPI_CONFIG" <<'EOF_CONFIG'

ws-auth: false
nonstream-keepalive-interval: 0
EOF_CONFIG

    print_success "Config written: $CLIPROXYAPI_CONFIG"
}

install_official_binary() {
    if [[ $CLIPROXYAPI_CMD != */* ]]; then
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

verify_models_on_start() {
    local models_url="${CLIPROXYAPI_BASE_URL}/v1/models"
    local attempt=1
    local tmp_file=""
    local code=""
    local models_count=""

    if ! is_true "$CLIPROXYAPI_VERIFY_MODELS_ON_START"; then
        return 0
    fi

    if [ -z "$CLIPROXYAPI_API_KEY" ]; then
        print_error "CLIPROXYAPI_API_KEY is required when CLIPROXYAPI_VERIFY_MODELS_ON_START=true"
        return 1
    fi

    print_step "Verifying authenticated models endpoint"
    print_info "URL: $models_url"

    while [ "$attempt" -le "$CLIPROXYAPI_HEALTH_RETRIES" ]; do
        tmp_file="$(mktemp)"
        code="$(curl -sS -m 8 \
            -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
            -o "$tmp_file" -w "%{http_code}" \
            "$models_url" || true)"

        if [ "$code" = "200" ]; then
            models_count="$(extract_models_count "$tmp_file" 2>/dev/null || true)"
            rm -f "$tmp_file"

            if [ -z "$models_count" ]; then
                print_warning "Unable to parse models payload (attempt ${attempt}/${CLIPROXYAPI_HEALTH_RETRIES})"
            elif is_true "$CLIPROXYAPI_EXPECT_MODELS_NON_EMPTY" && [ "$models_count" -le 0 ] 2>/dev/null; then
                print_warning "Models list is empty (attempt ${attempt}/${CLIPROXYAPI_HEALTH_RETRIES})"
            else
                print_success "Authenticated models endpoint is ready (${models_count:-unknown} models)"
                return 0
            fi
        else
            rm -f "$tmp_file"
        fi

        sleep 1
        attempt=$((attempt + 1))
    done

    print_error "Authenticated models endpoint did not become ready: $models_url"
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
  CLIPROXYAPI_CONFIG_AUTOGEN_OVERWRITE Overwrite existing config when autogen is enabled
  CLIPROXYAPI_PROVIDER_MODE          claude-api-key, openai-compatibility, or oauth
  CLIPROXYAPI_BIND_HOST              Listener host written to config (default: 0.0.0.0)
  CLIPROXYAPI_UPSTREAM_BASE_URL      Required upstream base URL
  CLIPROXYAPI_UPSTREAM_API_KEY       Required upstream key (or ANTHROPIC_AUTH_TOKEN fallback)
  CLIPROXYAPI_UPSTREAM_MODEL         Required upstream model id
  CLIPROXYAPI_UPSTREAM_MODEL_ALIAS   Model alias exposed to clients
  CLIPROXYAPI_UPSTREAM_NAME          Provider name for openai-compatibility mode
  CLIPROXYAPI_OAUTH_*                OAuth alias mapping env vars for codex/qwen/kimi
  CLIPROXYAPI_API_KEY                Local key accepted by CLIProxyAPI
  CLIPROXYAPI_VERIFY_MODELS_ON_START Validate authenticated /v1/models after health checks
  CLIPROXYAPI_ALLOW_INSECURE_DEFAULT_KEY Allow known weak default keys (not recommended)
  CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH Require config api-key to match CLIPROXYAPI_API_KEY
  CLIPROXYAPI_SYNC_CONFIG_API_KEY    Sync first config api-key entry with CLIPROXYAPI_API_KEY
EOF_USAGE
}

main() {
    local autogen_needed=false

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

    if is_true "$CLIPROXYAPI_CONFIG_AUTOGEN"; then
        if [ -f "$CLIPROXYAPI_CONFIG" ] && ! is_true "$CLIPROXYAPI_CONFIG_AUTOGEN_OVERWRITE"; then
            print_info "Config autogen skipped (existing config preserved): $CLIPROXYAPI_CONFIG"
        else
            autogen_needed=true
        fi
    fi

    if ! validate_local_api_key "$autogen_needed"; then
        return 1
    fi

    ensure_binary

    prepare_provider_defaults

    if [ "$autogen_needed" = true ]; then
        if ! validate_provider_config; then
            return 1
        fi
        if ! autogenerate_config; then
            return 1
        fi
    fi

    if [ ! -f "$CLIPROXYAPI_CONFIG" ]; then
        print_error "Config file not found: $CLIPROXYAPI_CONFIG"
        return 1
    fi

    if ! ensure_cliproxyapi_config_api_key_alignment; then
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
        local -a extra_args=($CLIPROXYAPI_EXTRA_ARGS)
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

    if ! verify_models_on_start; then
        tail -n 40 "$CLIPROXYAPI_LOG_FILE" 2>/dev/null || true
        return 1
    fi

    print_success "CLIProxyAPI started (PID: $pid)"
    return 0
}

main "$@"
