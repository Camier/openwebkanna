#!/bin/bash

###############################################################################
# CLIProxyAPI OAuth Setup Script
# Configures auth credentials for codex, qwen, and kimi
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

# Configuration
CLIPROXYAPI_CMD="${CLIPROXYAPI_CMD:-./bin/cliproxyapi}"
CLIPROXYAPI_CONFIG="${CLIPROXYAPI_CONFIG:-./cliproxyapi/config.yaml}"
CLIPROXYAPI_CONFIG_FLAG="${CLIPROXYAPI_CONFIG_FLAG:---config}"
CLIPROXYAPI_OAUTH_NO_BROWSER="${CLIPROXYAPI_OAUTH_NO_BROWSER:-false}"
CLIPROXYAPI_QWEN_AUTH_MODE="${CLIPROXYAPI_QWEN_AUTH_MODE:-auto}"
CLIPROXYAPI_CODEX_CALLBACK_PORT="${CLIPROXYAPI_CODEX_CALLBACK_PORT:-1455}"
CLIPROXYAPI_QWEN_CALLBACK_PORT="${CLIPROXYAPI_QWEN_CALLBACK_PORT:-0}"
CLIPROXYAPI_KIMI_CALLBACK_PORT="${CLIPROXYAPI_KIMI_CALLBACK_PORT:-0}"
CLIPROXYAPI_QWEN_SOURCE_FILE="${CLIPROXYAPI_QWEN_SOURCE_FILE:-$HOME/.qwen/oauth_creds.json}"
QWEN_IMPORT_SCRIPT="${QWEN_IMPORT_SCRIPT:-./import-qwen-auth.sh}"

expand_home_path() {
    local candidate="$1"
    # shellcheck disable=SC2088
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

executable_exists() {
    local executable="$1"
    if [[ $executable == */* ]]; then
        [ -x "$executable" ]
        return $?
    fi
    command_exists "$executable"
}

show_usage() {
    cat <<'EOF_USAGE'
Usage: ./configure-cliproxyapi-oauth.sh [all|codex|qwen|kimi] [--help]

Examples:
  ./configure-cliproxyapi-oauth.sh
  ./configure-cliproxyapi-oauth.sh all
  ./configure-cliproxyapi-oauth.sh qwen codex

Behavior:
  - qwen: controlled by CLIPROXYAPI_QWEN_AUTH_MODE
      auto  -> import if source file exists, else native -qwen-login OAuth
      import -> imports ~/.qwen/oauth_creds.json using import-qwen-auth.sh
      oauth  -> runs native -qwen-login OAuth flow
  - codex/kimi: runs CLIProxyAPI OAuth login flows
EOF_USAGE
}

run_oauth_login() {
    local provider="$1"
    local login_flag="$2"
    local callback_port="$3"

    local -a command=("$CLIPROXYAPI_CMD" "$CLIPROXYAPI_CONFIG_FLAG" "$CLIPROXYAPI_CONFIG" "$login_flag")
    if [[ $callback_port =~ ^[0-9]+$ ]] && [ "$callback_port" -gt 0 ]; then
        command+=("-oauth-callback-port" "$callback_port")
    fi
    if is_true "$CLIPROXYAPI_OAUTH_NO_BROWSER"; then
        command+=("-no-browser")
    fi

    print_step "Running $provider OAuth login"
    echo "Command: ${command[*]}"
    "${command[@]}"
    print_success "$provider OAuth login completed"
}

run_qwen_setup() {
    local auth_mode="$CLIPROXYAPI_QWEN_AUTH_MODE"
    auth_mode="$(printf "%s" "$auth_mode" | tr '[:upper:]' '[:lower:]')"
    local qwen_source_file
    qwen_source_file="$(expand_home_path "$CLIPROXYAPI_QWEN_SOURCE_FILE")"

    case "$auth_mode" in
        auto)
            if [ -f "$qwen_source_file" ]; then
                auth_mode="import"
            else
                auth_mode="oauth"
            fi
            ;;
        import | oauth) ;;
        *)
            print_error "Unsupported CLIPROXYAPI_QWEN_AUTH_MODE: $CLIPROXYAPI_QWEN_AUTH_MODE"
            print_error "Allowed values: auto, import, oauth"
            return 1
            ;;
    esac

    if [ "$auth_mode" = "import" ]; then
        if [ ! -x "$QWEN_IMPORT_SCRIPT" ]; then
            print_error "Qwen import script missing or not executable: $QWEN_IMPORT_SCRIPT"
            return 1
        fi
        print_step "Importing Qwen credentials from $qwen_source_file"
        "$QWEN_IMPORT_SCRIPT"
        print_success "Qwen credentials imported"
        return 0
    fi

    run_oauth_login "qwen" "-qwen-login" "$CLIPROXYAPI_QWEN_CALLBACK_PORT"
    return 0
}

main() {
    if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
        show_usage
        return 0
    fi

    print_header

    CLIPROXYAPI_CMD="$(resolve_path "$CLIPROXYAPI_CMD")"
    CLIPROXYAPI_CONFIG="$(resolve_path "$CLIPROXYAPI_CONFIG")"
    QWEN_IMPORT_SCRIPT="$(resolve_path "$QWEN_IMPORT_SCRIPT")"

    if ! executable_exists "$CLIPROXYAPI_CMD"; then
        print_error "CLIProxyAPI binary missing or not executable: $CLIPROXYAPI_CMD"
        return 1
    fi

    if [ ! -f "$CLIPROXYAPI_CONFIG" ]; then
        print_error "Config not found: $CLIPROXYAPI_CONFIG"
        print_error "Run ./start-cliproxyapi.sh once to generate config.yaml in oauth mode"
        return 1
    fi

    local -a providers=()
    if [ $# -eq 0 ] || [ "${1:-}" = "all" ]; then
        providers=("qwen" "codex" "kimi")
    else
        while [ $# -gt 0 ]; do
            case "$1" in
                qwen | codex | kimi)
                    providers+=("$1")
                    ;;
                *)
                    print_error "Unsupported provider: $1"
                    show_usage
                    return 1
                    ;;
            esac
            shift
        done
    fi

    local provider
    for provider in "${providers[@]}"; do
        case "$provider" in
            qwen)
                run_qwen_setup
                ;;
            codex)
                run_oauth_login "codex" "-codex-login" "$CLIPROXYAPI_CODEX_CALLBACK_PORT"
                ;;
            kimi)
                run_oauth_login "kimi" "-kimi-login" "$CLIPROXYAPI_KIMI_CALLBACK_PORT"
                ;;
        esac
    done

    print_success "OAuth setup workflow finished"
    return 0
}

main "$@"
