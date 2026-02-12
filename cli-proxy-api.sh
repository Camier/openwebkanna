#!/bin/bash

###############################################################################
# CLI Proxy API Wrapper
# Unified shell CLI for OpenWebUI and vLLM endpoints
###############################################################################

set -e

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

# Defaults
DEFAULT_OPENWEBUI_URL="http://localhost:3000"
DEFAULT_VLLM_URL="http://localhost:8000"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Runtime configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-}"
VLLM_URL="${VLLM_URL:-}"
OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-}"
URL_OVERRIDE=""
VERBOSE="false"
RAW="false"

###############################################################################
# Helper Functions
###############################################################################

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "+------------------------------------------------------------+"
    echo "|                     CLI Proxy API                          |"
    echo "+------------------------------------------------------------+"
    echo -e "${NC}"
}

print_step() {
    echo -e "${BLUE}${BOLD}[STEP] $1${NC}" >&2
}

print_success() {
    echo -e "${GREEN}[OK] $1${NC}" >&2
}

print_error() {
    echo -e "${RED}[ERROR] $1${NC}" >&2
}

print_warning() {
    echo -e "${YELLOW}[WARN] $1${NC}" >&2
}

print_info() {
    echo -e "${CYAN}[INFO] $1${NC}" >&2
}

usage() {
    cat <<'EOF'
Usage:
  cli-proxy-api.sh [global options] <command> [args]

Global Options:
  --url <base-url>           Override service base URL for the invoked command
  --api-key <token>          OpenWebUI API key (Bearer token)
  --verbose                  Print request details to stderr
  --raw                      Print raw JSON response without pretty formatting
  --help, -h                 Show this help message

Commands:
  health [openwebui|webui|vllm|all]
  models [openwebui|webui|vllm|all]
  chat --model <id> --message <text> [--service openwebui|webui|vllm] [--temperature <float>] [--max-tokens <int>]
  files list
  files upload --path <file>
  knowledge list

Environment:
  - If ".env" exists, it is loaded.
  - OPENWEBUI_URL defaults to http://localhost:3000
  - VLLM_URL defaults to http://localhost:8000
  - OPENWEBUI_API_KEY can be provided via environment variable.

Examples:
  ./cli-proxy-api.sh health all
  ./cli-proxy-api.sh models openwebui
  ./cli-proxy-api.sh models vllm
  ./cli-proxy-api.sh --api-key "$OPENWEBUI_API_KEY" chat --model llama-3.1-8b --message "Summarize this paper."
  ./cli-proxy-api.sh chat --service vllm --model meta-llama/Llama-3.1-8B-Instruct --message "Hello" --temperature 0.2 --max-tokens 200
  ./cli-proxy-api.sh files list
  ./cli-proxy-api.sh files upload --path ./data/pdfs/paper.pdf
  ./cli-proxy-api.sh knowledge list
EOF
}

load_env_file() {
    local env_file=""
    local file_openwebui_url=""
    local file_vllm_url=""
    local file_openai_api_base_url=""
    local file_openai_api_base_urls=""
    local file_webui_port=""
    local file_openwebui_api_key=""

    if [ -f ".env" ]; then
        env_file=".env"
    elif [ -f "${SCRIPT_DIR}/.env" ]; then
        env_file="${SCRIPT_DIR}/.env"
    fi

    if [ -n "$env_file" ]; then
        file_openwebui_url="$(read_env_value "$env_file" "OPENWEBUI_URL")"
        file_vllm_url="$(read_env_value "$env_file" "VLLM_URL")"
        file_openai_api_base_url="$(read_env_value "$env_file" "OPENAI_API_BASE_URL")"
        file_openai_api_base_urls="$(read_env_value "$env_file" "OPENAI_API_BASE_URLS")"
        file_webui_port="$(read_env_value "$env_file" "WEBUI_PORT")"
        file_openwebui_api_key="$(read_env_value "$env_file" "OPENWEBUI_API_KEY")"
    fi

    if [ -z "$file_vllm_url" ] && [ -n "$file_openai_api_base_url" ]; then
        file_vllm_url="$file_openai_api_base_url"
    fi

    if [ -z "$file_vllm_url" ] && [ -n "$file_openai_api_base_urls" ]; then
        # OPENAI_API_BASE_URLS can contain semicolon-separated URLs.
        file_vllm_url="$(printf '%s' "$file_openai_api_base_urls" | cut -d';' -f1)"
    fi

    # Convert OpenAI-compatible /v1 base URL to server root base.
    if [ -n "$file_vllm_url" ]; then
        file_vllm_url="$(printf '%s' "$file_vllm_url" | sed 's#/v1/*$##')"
    fi

    if [ -n "$file_openwebui_url" ] && [ -z "$OPENWEBUI_URL" ]; then
        OPENWEBUI_URL="$file_openwebui_url"
    fi

    if [ -z "$OPENWEBUI_URL" ] && [ -n "$file_webui_port" ]; then
        OPENWEBUI_URL="http://localhost:${file_webui_port}"
    fi

    if [ -n "$file_vllm_url" ] && [ -z "$VLLM_URL" ]; then
        VLLM_URL="$file_vllm_url"
    fi
    if [ -n "$file_openwebui_api_key" ] && [ -z "$OPENWEBUI_API_KEY" ]; then
        OPENWEBUI_API_KEY="$file_openwebui_api_key"
    fi

    OPENWEBUI_URL="${OPENWEBUI_URL:-$DEFAULT_OPENWEBUI_URL}"
    VLLM_URL="${VLLM_URL:-$DEFAULT_VLLM_URL}"
    OPENWEBUI_URL="$(normalize_base_url "$OPENWEBUI_URL")"
    VLLM_URL="$(normalize_base_url "$VLLM_URL")"
    OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-}"
}

read_env_value() {
    local env_file="$1"
    local key="$2"
    local value=""

    if [ ! -f "$env_file" ]; then
        return 0
    fi

    value="$(
        awk -v key="$key" '
            {
                line = $0
                sub(/\r$/, "", line)
            }
            /^[[:space:]]*$/ { next }
            /^[[:space:]]*#/ { next }
            {
                sub(/^[[:space:]]*export[[:space:]]+/, "", line)
                eq_index = index(line, "=")
                if (eq_index == 0) {
                    next
                }

                lhs = substr(line, 1, eq_index - 1)
                rhs = substr(line, eq_index + 1)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", lhs)

                if (lhs != key) {
                    next
                }

                sub(/^[[:space:]]+/, "", rhs)
                sub(/[[:space:]]+$/, "", rhs)
                print rhs
            }
        ' "$env_file" | tail -n 1
    )"

    if [ -z "$value" ]; then
        return 0
    fi

    value="$(printf '%s' "$value" | sed 's/[[:space:]]#.*$//')"
    value="$(printf '%s' "$value" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

    if [[ "$value" == \"*\" && "$value" == *\" ]]; then
        value="${value:1:${#value}-2}"
    elif [[ "$value" == \'*\' && "$value" == *\' ]]; then
        value="${value:1:${#value}-2}"
    fi

    printf '%s' "$value"
}

running_in_container() {
    if [ -f "/.dockerenv" ]; then
        return 0
    fi

    grep -qaE "(docker|containerd|kubepods|podman|lxc)" /proc/1/cgroup 2>/dev/null
}

host_is_resolvable() {
    local host="$1"
    if command -v getent >/dev/null 2>&1; then
        getent hosts "$host" >/dev/null 2>&1
        return $?
    fi
    if command -v nslookup >/dev/null 2>&1; then
        nslookup "$host" >/dev/null 2>&1
        return $?
    fi
    return 1
}

normalize_base_url() {
    local url="$1"
    if [ -z "$url" ]; then
        printf '%s' "$url"
        return 0
    fi

    if [[ "$url" == *"host.docker.internal"* ]]; then
        if ! running_in_container; then
            url="$(printf '%s' "$url" | sed 's#://host\.docker\.internal#://localhost#g')"
        elif ! host_is_resolvable "host.docker.internal"; then
            print_warning "host.docker.internal is not resolvable in container context: $url"
        fi
    fi

    printf '%s' "$url"
}

normalize_service_name() {
    local service="$1"
    service="$(printf '%s' "$service" | tr '[:upper:]' '[:lower:]')"

    case "$service" in
        webui|openwebui)
            printf 'openwebui'
            ;;
        vllm|all)
            printf '%s' "$service"
            ;;
        *)
            printf '%s' "$service"
            ;;
    esac
}

json_escape() {
    printf '%s' "$1" | sed \
        -e 's/\\/\\\\/g' \
        -e 's/"/\\"/g' \
        -e ':a;N;$!ba;s/\n/\\n/g' \
        -e 's/\r/\\r/g' \
        -e 's/\t/\\t/g'
}

emit_json() {
    local payload="$1"

    if [ "$RAW" = "true" ]; then
        printf '%s\n' "$payload"
        return 0
    fi

    if command -v jq >/dev/null 2>&1; then
        if echo "$payload" | jq -e . >/dev/null 2>&1; then
            echo "$payload" | jq .
            return 0
        fi
    fi

    printf '%s\n' "$payload"
}

resolve_service_url() {
    local service="$1"

    if [ -n "$URL_OVERRIDE" ]; then
        echo "$URL_OVERRIDE"
        return 0
    fi

    if [ "$service" = "openwebui" ]; then
        echo "$OPENWEBUI_URL"
    else
        echo "$VLLM_URL"
    fi
}

is_float() {
    local value="$1"
    [[ "$value" =~ ^-?[0-9]+([.][0-9]+)?$ ]]
}

is_int() {
    local value="$1"
    [[ "$value" =~ ^-?[0-9]+$ ]]
}

run_curl() {
    local method="$1"
    local url="$2"
    shift 2

    local -a cmd
    local response=""
    local err_file=""

    cmd=("curl" "--fail" "--show-error" "--silent" "-X" "$method" "$url")
    cmd+=("$@")

    if [ "$VERBOSE" = "true" ]; then
        print_step "$method $url"
    fi

    err_file="$(mktemp)"
    response=$("${cmd[@]}" 2>"$err_file")
    local curl_exit="$?"

    if [ "$curl_exit" -eq 0 ]; then
        rm -f "$err_file"
        printf '%s' "$response"
        return 0
    else
        local err_output=""

        err_output="$(cat "$err_file")"
        rm -f "$err_file"

        case "$curl_exit" in
            6)
                print_error "Unable to resolve host for request: $url"
                ;;
            7)
                print_error "Unable to connect to endpoint: $url"
                ;;
            22)
                print_error "Request returned non-success HTTP status: $method $url"
                ;;
            *)
                print_error "Request failed (${curl_exit}) for $method $url"
                ;;
        esac

        if [ -n "$err_output" ]; then
            print_error "$err_output"
        fi

        return 1
    fi
}

request_openwebui_json() {
    local method="$1"
    local endpoint="$2"
    local payload="${3:-}"
    local base_url="$4"
    local url="${base_url%/}${endpoint}"
    local -a headers

    headers=("-H" "Content-Type: application/json")
    if [ -n "$OPENWEBUI_API_KEY" ]; then
        headers+=("-H" "Authorization: Bearer $OPENWEBUI_API_KEY")
    fi

    if [ -n "$payload" ]; then
        run_curl "$method" "$url" "${headers[@]}" "-d" "$payload"
    else
        run_curl "$method" "$url" "${headers[@]}"
    fi
}

request_openwebui_form_upload() {
    local path="$1"
    local base_url="$2"
    local url="${base_url%/}/api/files/upload"
    local filename
    local escaped_filename
    local meta_json
    local -a headers

    filename="$(basename "$path")"
    escaped_filename="$(json_escape "$filename")"
    meta_json="{\"name\":\"$escaped_filename\"}"

    headers=()
    if [ -n "$OPENWEBUI_API_KEY" ]; then
        headers+=("-H" "Authorization: Bearer $OPENWEBUI_API_KEY")
    fi

    run_curl "POST" "$url" "${headers[@]}" "-F" "file=@${path}" "-F" "meta=${meta_json}"
}

request_vllm_json() {
    local method="$1"
    local endpoint="$2"
    local payload="${3:-}"
    local base_url="$4"
    local url="${base_url%/}${endpoint}"
    local -a headers

    headers=("-H" "Content-Type: application/json")

    if [ -n "$payload" ]; then
        run_curl "$method" "$url" "${headers[@]}" "-d" "$payload"
    else
        run_curl "$method" "$url" "${headers[@]}"
    fi
}

require_arg() {
    local value="$1"
    local message="$2"
    if [ -z "$value" ]; then
        print_error "$message"
        exit 1
    fi
}

###############################################################################
# Command Handlers
###############################################################################

handle_health() {
    local target="${1:-all}"
    local openwebui_base
    local vllm_base
    local openwebui_resp
    local vllm_resp

    target="$(normalize_service_name "$target")"

    if [ "$target" != "openwebui" ] && [ "$target" != "vllm" ] && [ "$target" != "all" ]; then
        print_error "Invalid health target '$target'. Use one of: openwebui|webui|vllm|all."
        usage
        exit 1
    fi

    openwebui_base="$(resolve_service_url "openwebui")"
    vllm_base="$(resolve_service_url "vllm")"

    if [ "$target" = "openwebui" ]; then
        openwebui_resp="$(request_openwebui_json "GET" "/health" "" "$openwebui_base")" || return 1
        emit_json "$openwebui_resp"
        return 0
    fi

    if [ "$target" = "vllm" ]; then
        vllm_resp="$(request_vllm_json "GET" "/health" "" "$vllm_base")" || return 1
        emit_json "$vllm_resp"
        return 0
    fi

    openwebui_resp="$(request_openwebui_json "GET" "/health" "" "$openwebui_base")" || return 1
    vllm_resp="$(request_vllm_json "GET" "/health" "" "$vllm_base")" || return 1
    emit_json "{\"openwebui\":${openwebui_resp},\"vllm\":${vllm_resp}}"
}

handle_models() {
    local target="${1:-openwebui}"
    local openwebui_base
    local vllm_base
    local openwebui_resp
    local vllm_resp

    target="$(normalize_service_name "$target")"

    if [ "$target" != "openwebui" ] && [ "$target" != "vllm" ] && [ "$target" != "all" ]; then
        print_error "Invalid models target '$target'. Use one of: openwebui|webui|vllm|all."
        usage
        exit 1
    fi

    openwebui_base="$(resolve_service_url "openwebui")"
    vllm_base="$(resolve_service_url "vllm")"

    if [ "$target" = "openwebui" ]; then
        openwebui_resp="$(request_openwebui_json "GET" "/api/models" "" "$openwebui_base")" || return 1
        emit_json "$openwebui_resp"
        return 0
    fi

    if [ "$target" = "vllm" ]; then
        vllm_resp="$(request_vllm_json "GET" "/v1/models" "" "$vllm_base")" || return 1
        emit_json "$vllm_resp"
        return 0
    fi

    openwebui_resp="$(request_openwebui_json "GET" "/api/models" "" "$openwebui_base")" || return 1
    vllm_resp="$(request_vllm_json "GET" "/v1/models" "" "$vllm_base")" || return 1
    emit_json "{\"openwebui\":${openwebui_resp},\"vllm\":${vllm_resp}}"
}

handle_chat() {
    local service="openwebui"
    local model=""
    local message=""
    local temperature=""
    local max_tokens=""
    local base_url
    local payload=""
    local escaped_model=""
    local escaped_message=""

    while [ $# -gt 0 ]; do
        case "$1" in
            --service)
                shift
                require_arg "${1:-}" "Missing value for --service"
                service="$1"
                ;;
            --model)
                shift
                require_arg "${1:-}" "Missing value for --model"
                model="$1"
                ;;
            --message)
                shift
                require_arg "${1:-}" "Missing value for --message"
                message="$1"
                ;;
            --temperature)
                shift
                require_arg "${1:-}" "Missing value for --temperature"
                temperature="$1"
                ;;
            --max-tokens)
                shift
                require_arg "${1:-}" "Missing value for --max-tokens"
                max_tokens="$1"
                ;;
            --help|-h)
                usage
                exit 0
                ;;
            *)
                print_error "Unknown chat option: $1"
                usage
                exit 1
                ;;
        esac
        shift
    done

    require_arg "$model" "chat requires --model <id>"
    require_arg "$message" "chat requires --message <text>"

    service="$(normalize_service_name "$service")"
    if [ "$service" = "all" ]; then
        print_error "chat does not support '--service all'. Use: openwebui|webui|vllm."
        exit 1
    fi

    if [ "$service" != "openwebui" ] && [ "$service" != "vllm" ]; then
        print_error "Invalid chat service '$service'. Use: openwebui|webui|vllm."
        exit 1
    fi

    if [ -n "$temperature" ] && ! is_float "$temperature"; then
        print_error "--temperature must be a number"
        exit 1
    fi

    if [ -n "$max_tokens" ] && ! is_int "$max_tokens"; then
        print_error "--max-tokens must be an integer"
        exit 1
    fi

    escaped_model="$(json_escape "$model")"
    escaped_message="$(json_escape "$message")"
    payload="{\"model\":\"${escaped_model}\",\"messages\":[{\"role\":\"user\",\"content\":\"${escaped_message}\"}]"

    if [ -n "$temperature" ]; then
        payload="${payload},\"temperature\":${temperature}"
    fi
    if [ -n "$max_tokens" ]; then
        payload="${payload},\"max_tokens\":${max_tokens}"
    fi
    payload="${payload}}"

    base_url="$(resolve_service_url "$service")"
    if [ "$service" = "openwebui" ]; then
        payload="$(request_openwebui_json "POST" "/api/chat/completions" "$payload" "$base_url")" || return 1
        emit_json "$payload"
    else
        payload="$(request_vllm_json "POST" "/v1/chat/completions" "$payload" "$base_url")" || return 1
        emit_json "$payload"
    fi
}

handle_files() {
    local subcommand="${1:-}"
    shift || true
    local base_url

    base_url="$(resolve_service_url "openwebui")"

    case "$subcommand" in
        list)
            local files_response=""
            files_response="$(request_openwebui_json "GET" "/api/v1/files/" "" "$base_url")" || return 1
            emit_json "$files_response"
            ;;
        upload)
            local path=""
            while [ $# -gt 0 ]; do
                case "$1" in
                    --path)
                        shift
                        require_arg "${1:-}" "Missing value for --path"
                        path="$1"
                        ;;
                    --help|-h)
                        usage
                        exit 0
                        ;;
                    *)
                        print_error "Unknown files upload option: $1"
                        usage
                        exit 1
                        ;;
                esac
                shift
            done

            require_arg "$path" "files upload requires --path <file>"
            if [ ! -f "$path" ]; then
                print_error "File does not exist: $path"
                exit 1
            fi

            local upload_response=""
            upload_response="$(request_openwebui_form_upload "$path" "$base_url")" || return 1
            emit_json "$upload_response"
            ;;
        *)
            print_error "Invalid files subcommand. Use 'files list' or 'files upload --path <file>'."
            usage
            exit 1
            ;;
    esac
}

handle_knowledge() {
    local subcommand="${1:-}"
    local base_url

    base_url="$(resolve_service_url "openwebui")"

    case "$subcommand" in
        list)
            local knowledge_response=""
            knowledge_response="$(request_openwebui_json "GET" "/api/v1/knowledge/" "" "$base_url")" || return 1
            emit_json "$knowledge_response"
            ;;
        *)
            print_error "Invalid knowledge subcommand. Use 'knowledge list'."
            usage
            exit 1
            ;;
    esac
}

###############################################################################
# Main
###############################################################################

load_env_file

# Parse global options
while [ $# -gt 0 ]; do
    case "$1" in
        --url)
            shift
            require_arg "${1:-}" "Missing value for --url"
            URL_OVERRIDE="$1"
            ;;
        --api-key)
            shift
            require_arg "${1:-}" "Missing value for --api-key"
            OPENWEBUI_API_KEY="$1"
            ;;
        --verbose)
            VERBOSE="true"
            ;;
        --raw)
            RAW="true"
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            print_error "Unknown global option: $1"
            usage
            exit 1
            ;;
        *)
            break
            ;;
    esac
    shift
done

COMMAND="${1:-}"
if [ -z "$COMMAND" ]; then
    usage
    exit 1
fi
shift

case "$COMMAND" in
    health)
        handle_health "${1:-all}"
        ;;
    models)
        handle_models "${1:-openwebui}"
        ;;
    chat)
        handle_chat "$@"
        ;;
    files)
        handle_files "$@"
        ;;
    knowledge)
        handle_knowledge "$@"
        ;;
    --help|-h|help)
        usage
        ;;
    *)
        print_error "Unknown command: $COMMAND"
        usage
        exit 1
        ;;
esac
