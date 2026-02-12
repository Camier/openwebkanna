#!/bin/bash

###############################################################################
# OpenWebUI Knowledge Base Importer (PDFs)
# Upload PDFs to OpenWebUI, wait for processing, and attach them to a Knowledge Base.
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

is_uuid() {
    local value="${1:-}"
    echo "$value" | rg -q '^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$'
}

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

# Colors
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
    echo "║  OpenWebUI PDF Import → Knowledge Base                     ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    echo -e "${BLUE}${BOLD}▶ $1${NC}"
}

print_info() {
    echo -e "  ${CYAN}ℹ${NC} $1"
}

print_success() {
    echo -e "  ${GREEN}✓${NC} $1"
}

print_warning() {
    echo -e "  ${YELLOW}⚠${NC} $1"
}

print_error() {
    echo -e "  ${RED}✗${NC} $1" >&2
}

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        print_error "Missing dependency: $cmd"
        exit 1
    fi
}

OPENWEBUI_URL_DEFAULT="http://localhost:${WEBUI_PORT:-3000}"
OPENWEBUI_URL="${OPENWEBUI_URL:-$OPENWEBUI_URL_DEFAULT}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"

API_KEY="${OPENWEBUI_API_KEY:-${API_KEY:-}}"

PDF_DIR="data/pdfs"
KB_NAME="Academic Papers"
KB_DESCRIPTION="Imported PDFs"
GLOB_PATTERN="*.pdf"
LIMIT_COUNT="${LIMIT_COUNT:-0}"
PROCESS_TIMEOUT_SECONDS="${PROCESS_TIMEOUT_SECONDS:-900}"
PROCESS_POLL_SECONDS="${PROCESS_POLL_SECONDS:-3}"

SHOW_HELP=false

show_help() {
    cat <<'EOF'
Usage: ./import-pdfs-to-kb.sh [OPTIONS]

Options:
  --dir PATH              Directory containing PDFs (default: data/pdfs)
  --glob PATTERN          Filename pattern (default: *.pdf)
  --limit N               Process at most N PDFs (0 = no limit; default: 0)
  --kb-name NAME          Knowledge Base name (default: Academic Papers)
  --kb-description TEXT   Knowledge Base description (default: Imported PDFs)
  --url URL               OpenWebUI URL override (default: http://localhost:${WEBUI_PORT})
  --timeout SECONDS       Max seconds to wait for per-file processing (default: 900)
  --poll SECONDS          Poll interval seconds (default: 3)
  -h, --help              Show help

Auth:
  - If OPENWEBUI_API_KEY (or API_KEY) is set, it's used directly.
  - Otherwise, the script signs in via /api/v1/auths/signin using:
    OPENWEBUI_SIGNIN_EMAIL / OPENWEBUI_SIGNIN_PASSWORD
EOF
}

openwebui_signin() {
    local payload_file="$1"
    local output_file="$2"
    local signin_url="${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"

    curl -sS -m 30 \
        -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" \
        "$signin_url" -d @"$payload_file" 2>/dev/null || true
}

ensure_api_key() {
    if [ -n "$API_KEY" ]; then
        return 0
    fi

    if [ "$OPENWEBUI_AUTO_AUTH" = "false" ]; then
        print_error "No API key present and OPENWEBUI_AUTO_AUTH=false; cannot authenticate"
        exit 1
    fi

    local tmp_out payload_file signin_code token
    tmp_out="$(mktemp)"
    payload_file="$(mktemp)"

    cat >"$payload_file" <<EOF
{"email":"$OPENWEBUI_SIGNIN_EMAIL","password":"$OPENWEBUI_SIGNIN_PASSWORD"}
EOF

    signin_code="$(openwebui_signin "$payload_file" "$tmp_out")"
    if [ "$signin_code" != "200" ]; then
        print_error "Signin failed (HTTP $signin_code) to ${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"
        rm -f "$tmp_out" "$payload_file"
        exit 1
    fi

    token="$(cat "$tmp_out" | jq -r '.token // empty' 2>/dev/null || true)"
    rm -f "$tmp_out" "$payload_file"

    if [ -z "$token" ]; then
        print_error "Signin succeeded but no token returned"
        exit 1
    fi

    API_KEY="$token"
    export API_KEY
}

api_get() {
    local path="$1"
    curl -sS -m 45 \
        -H "Authorization: Bearer $API_KEY" \
        "${OPENWEBUI_URL%/}${path}"
}

api_post_json() {
    local path="$1"
    local json="$2"
    curl -sS -m 45 \
        -H "Authorization: Bearer $API_KEY" \
        -H "Content-Type: application/json" \
        -X POST \
        -d "$json" \
        "${OPENWEBUI_URL%/}${path}"
}

api_post_json_with_status() {
    local path="$1"
    local json="$2"
    local out_file="$3"
    curl -sS -m 45 \
        -H "Authorization: Bearer $API_KEY" \
        -H "Content-Type: application/json" \
        -X POST \
        -d "$json" \
        -o "$out_file" -w "%{http_code}" \
        "${OPENWEBUI_URL%/}${path}" 2>/dev/null || true
}

find_or_create_kb() {
    local kb_list kb_id
    kb_list="$(api_get "/api/v1/knowledge/" 2>/dev/null || true)"

    if [ -n "$kb_list" ]; then
        kb_id="$(echo "$kb_list" | jq -r --arg name "$KB_NAME" '.[] | select(.name == $name) | .id' 2>/dev/null | head -n 1 || true)"
        if [ -n "$kb_id" ] && [ "$kb_id" != "null" ]; then
            print_success "Using existing Knowledge Base: $KB_NAME ($kb_id)" >&2
            echo "$kb_id"
            return 0
        fi
    fi

    local create_payload response kb_create_out kb_create_code
    create_payload="$(jq -cn --arg name "$KB_NAME" --arg desc "$KB_DESCRIPTION" '{name:$name, description:$desc}')"
    kb_create_out="$(mktemp)"
    kb_create_code="$(api_post_json_with_status "/api/v1/knowledge/create" "$create_payload" "$kb_create_out")"
    response="$(cat "$kb_create_out" 2>/dev/null || true)"
    rm -f "$kb_create_out"

    if [ -z "$kb_create_code" ] || ! echo "$kb_create_code" | rg -q '^[0-9]{3}$' || [ "$kb_create_code" -lt 200 ] || [ "$kb_create_code" -ge 300 ]; then
        print_error "Knowledge Base create failed (HTTP ${kb_create_code:-unknown})"
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)" >&2
        exit 1
    fi

    kb_id="$(echo "$response" | jq -r '.id // .knowledge_id // empty' 2>/dev/null || true)"

    if [ -z "$kb_id" ]; then
        print_error "Failed to create Knowledge Base"
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)" >&2
        exit 1
    fi

    print_success "Created Knowledge Base: $KB_NAME ($kb_id)" >&2
    echo "$kb_id"
}

upload_pdf() {
    local pdf_path="$1"
    local filename
    filename="$(basename "$pdf_path")"

    curl -sS -m 300 -X POST \
        "${OPENWEBUI_URL%/}/api/v1/files/?process=true&process_in_background=false" \
        -H "Authorization: Bearer $API_KEY" \
        -F "file=@$pdf_path" \
        -F "metadata={\"name\":\"$filename\"}" 2>/dev/null || true
}

file_processing_status() {
    local file_id="$1"
    api_get "/api/v1/files/$file_id/process/status" 2>/dev/null || true
}

wait_for_processing() {
    local file_id="$1"
    local deadline=$((SECONDS + PROCESS_TIMEOUT_SECONDS))

    while [ $SECONDS -lt $deadline ]; do
        local status_resp status_val
        status_resp="$(file_processing_status "$file_id")"
        status_val="$(echo "$status_resp" | jq -r '.status // empty' 2>/dev/null || true)"

        if [ "$status_val" = "completed" ]; then
            return 0
        fi
        if [ "$status_val" = "failed" ]; then
            return 1
        fi

        sleep "$PROCESS_POLL_SECONDS"
    done

    return 1
}

attach_file_to_kb() {
    local kb_id="$1"
    local file_id="$2"
    local payload attach_out attach_code response ok

    payload="$(jq -cn --arg file_id "$file_id" '{file_id:$file_id}')"
    attach_out="$(mktemp)"
    attach_code="$(api_post_json_with_status "/api/v1/knowledge/$kb_id/file/add" "$payload" "$attach_out")"
    response="$(cat "$attach_out" 2>/dev/null || true)"
    rm -f "$attach_out"

    if [ -z "$attach_code" ] || ! echo "$attach_code" | rg -q '^[0-9]{3}$' || [ "$attach_code" -lt 200 ] || [ "$attach_code" -ge 300 ]; then
        print_error "Attach API failed (HTTP ${attach_code:-unknown})" >&2
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)" >&2
        return 1
    fi

    ok="$(echo "$response" | jq -r --arg file_id "$file_id" '
        if (.files? | type) == "array" then
            ([.files[]? | .id? // empty] | index($file_id)) // empty
        else
            empty
        end
    ' 2>/dev/null || true)"

    if [ -n "$ok" ] && [ "$ok" != "null" ]; then
        return 0
    fi

    # Fallback: some versions may return a different object; accept presence of files array.
    if echo "$response" | jq -e '(.files? | type) == "array"' >/dev/null 2>&1; then
        return 0
    fi

    print_error "Attach succeeded but response does not look like a KB file list" >&2
    [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)" >&2
    return 1
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --dir)
                PDF_DIR="$2"
                shift 2
                ;;
            --glob)
                GLOB_PATTERN="$2"
                shift 2
                ;;
            --kb-name)
                KB_NAME="$2"
                shift 2
                ;;
            --kb-description)
                KB_DESCRIPTION="$2"
                shift 2
                ;;
            --limit)
                LIMIT_COUNT="$2"
                shift 2
                ;;
            --url)
                OPENWEBUI_URL="$2"
                shift 2
                ;;
            --timeout)
                PROCESS_TIMEOUT_SECONDS="$2"
                shift 2
                ;;
            --poll)
                PROCESS_POLL_SECONDS="$2"
                shift 2
                ;;
            -h|--help)
                SHOW_HELP=true
                shift 1
                ;;
            *)
                print_error "Unknown argument: $1"
                SHOW_HELP=true
                shift 1
                ;;
        esac
    done
}

main() {
    parse_args "$@"
    if [ "$SHOW_HELP" = true ]; then
        show_help
        exit 0
    fi

    require_cmd curl
    require_cmd jq
    require_cmd rg

    print_header
    print_step "Configuration"
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "PDF directory: $PDF_DIR (pattern: $GLOB_PATTERN)"
    print_info "Knowledge Base: $KB_NAME"
    print_info "Limit: ${LIMIT_COUNT}"
    print_info "Processing timeout (per file): ${PROCESS_TIMEOUT_SECONDS}s"
    print_info "Poll interval: ${PROCESS_POLL_SECONDS}s"

    if ! curl -sS -m 10 -f -o /dev/null "${OPENWEBUI_URL%/}/health" 2>/dev/null; then
        print_error "OpenWebUI health check failed at ${OPENWEBUI_URL%/}/health"
        exit 1
    fi
    print_success "OpenWebUI health check OK"

    ensure_api_key

    if [ ! -d "$PDF_DIR" ]; then
        print_error "Directory not found: $PDF_DIR"
        exit 1
    fi

    local kb_id
    print_step "Knowledge Base"
    kb_id="$(find_or_create_kb)"
    if [ -z "$kb_id" ]; then
        print_error "Knowledge Base ID is empty"
        exit 1
    fi
    if ! is_uuid "$kb_id"; then
        print_error "Knowledge Base ID is not a UUID: $kb_id"
        exit 1
    fi

    shopt -s nullglob
    local pdfs=("$PDF_DIR"/$GLOB_PATTERN)
    shopt -u nullglob

    if [ ${#pdfs[@]} -eq 0 ]; then
        print_warning "No PDFs found in $PDF_DIR matching $GLOB_PATTERN"
        exit 0
    fi

    if ! echo "${LIMIT_COUNT}" | rg -q '^[0-9]+$'; then
        print_error "--limit must be an integer >= 0 (got: ${LIMIT_COUNT})"
        exit 1
    fi

    print_step "Upload + Attach"
    local total=0
    local uploaded=0
    local processed=0
    local attached=0
    local failed=0

    for pdf_path in "${pdfs[@]}"; do
        if [ "${LIMIT_COUNT}" -gt 0 ] && [ "${total}" -ge "${LIMIT_COUNT}" ]; then
            break
        fi
        total=$((total + 1))
        local filename response file_id
        filename="$(basename "$pdf_path")"

        print_info "Uploading: $filename"
        response="$(upload_pdf "$pdf_path")"

        file_id="$(echo "$response" | jq -r '.id // .file_id // empty' 2>/dev/null || true)"
        if [ -z "$file_id" ]; then
            failed=$((failed + 1))
            print_error "Upload failed for: $filename"
            continue
        fi

        uploaded=$((uploaded + 1))
        print_success "Uploaded: $filename (file_id=$file_id)"

        print_info "Waiting for processing: $filename"
        if wait_for_processing "$file_id"; then
            processed=$((processed + 1))
            print_success "Processed: $filename"
        else
            failed=$((failed + 1))
            print_error "Processing failed or timed out: $filename"
            continue
        fi

        print_info "Attaching to Knowledge Base: $filename"
        if attach_file_to_kb "$kb_id" "$file_id"; then
            attached=$((attached + 1))
            print_success "Attached: $filename"
        else
            failed=$((failed + 1))
            print_error "Attach failed: $filename"
            continue
        fi
    done

    print_step "Summary"
    print_info "Total PDFs seen: $total"
    print_info "Uploaded: $uploaded"
    print_info "Processed: $processed"
    print_info "Attached: $attached"
    print_info "Failed: $failed"

    if [ "$failed" -gt 0 ]; then
        exit 1
    fi
    print_success "Import complete"
}

main "$@"
