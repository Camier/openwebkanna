#!/bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_URL_DEFAULT="http://localhost:${WEBUI_PORT:-3000}"
OPENWEBUI_URL="${OPENWEBUI_URL:-$OPENWEBUI_URL_DEFAULT}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
API_KEY="${OPENWEBUI_API_KEY:-${API_KEY:-}}"

KB_NAME="${KB_NAME:-Molecular SMILES High Confidence}"
KEEP_ID="${KEEP_ID:-}"
DRY_RUN=true

show_help() {
    cat <<'EOF'
Usage: ./dedupe-smiles-kb.sh [OPTIONS]

Options:
  --kb-name NAME       Knowledge Base name to dedupe
  --keep-id UUID       Keep this KB id, delete other same-name ids
  --url URL            OpenWebUI base URL
  --apply              Apply deletions (default: dry-run)
  --dry-run            Show what would be deleted
  -h, --help           Show help
EOF
}

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        print_error "Missing dependency: $cmd"
        exit 1
    fi
}

is_uuid() {
    local value="${1:-}"
    echo "$value" | rg -q '^[0-9a-fA-F-]{36}$'
}

slugify() {
    echo "$1" | tr '[:upper:]' '[:lower:]' | tr -cs 'a-z0-9' '-'
}

cached_keep_id() {
    local slug path
    slug="$(slugify "$KB_NAME")"
    path="$SCRIPT_DIR/.state/${slug}.kb_id"
    if [ -f "$path" ]; then
        cat "$path"
    fi
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
        print_error "No API key present and OPENWEBUI_AUTO_AUTH=false"
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
        print_error "Signin failed (HTTP $signin_code)"
        rm -f "$tmp_out" "$payload_file"
        exit 1
    fi

    token="$(jq -r '.token // empty' "$tmp_out" 2>/dev/null || true)"
    rm -f "$tmp_out" "$payload_file"
    if [ -z "$token" ]; then
        print_error "Signin succeeded but no token returned"
        exit 1
    fi

    API_KEY="$token"
    export API_KEY
}

list_kb_json() {
    curl -sS -m 30 \
        -H "Authorization: Bearer $API_KEY" \
        "${OPENWEBUI_URL%/}/api/v1/knowledge/" 2>/dev/null || true
}

delete_kb() {
    local kb_id="$1"
    local code

    code="$(curl -sS -m 30 -o /dev/null -w "%{http_code}" \
        -H "Authorization: Bearer $API_KEY" \
        -X DELETE "${OPENWEBUI_URL%/}/api/v1/knowledge/$kb_id/delete" 2>/dev/null || true)"

    if [ "$code" = "200" ] || [ "$code" = "204" ]; then
        return 0
    fi

    code="$(curl -sS -m 30 -o /dev/null -w "%{http_code}" \
        -H "Authorization: Bearer $API_KEY" \
        -H "Content-Type: application/json" \
        -X POST "${OPENWEBUI_URL%/}/api/v1/knowledge/$kb_id/delete" 2>/dev/null || true)"

    [ "$code" = "200" ] || [ "$code" = "204" ]
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --kb-name)
                KB_NAME="$2"
                shift 2
                ;;
            --keep-id)
                KEEP_ID="$2"
                shift 2
                ;;
            --url)
                OPENWEBUI_URL="$2"
                shift 2
                ;;
            --apply)
                DRY_RUN=false
                shift 1
                ;;
            --dry-run)
                DRY_RUN=true
                shift 1
                ;;
            -h | --help)
                show_help
                exit 0
                ;;
            *)
                print_error "Unknown argument: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

main() {
    parse_args "$@"

    require_cmd curl
    require_cmd jq
    require_cmd rg

    print_header "OpenWebUI KB Dedupe"
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "Target KB name: $KB_NAME"

    ensure_api_key

    if [ -z "$KEEP_ID" ]; then
        KEEP_ID="$(cached_keep_id)"
    fi

    local raw
    raw="$(list_kb_json)"
    if [ -z "$raw" ]; then
        print_error "Unable to load knowledge list"
        exit 1
    fi

    local ids
    ids="$(echo "$raw" | jq -r --arg name "$KB_NAME" '
        if type == "array" then
            .[]
        elif (.items? | type) == "array" then
            .items[]
        elif (.data? | type) == "array" then
            .data[]
        else
            empty
        end
        | select(((.name // "") | ascii_downcase) == ($name | ascii_downcase))
        | .id
    ' 2>/dev/null || true)"

    if [ -z "$ids" ]; then
        print_warning "No KB entries found for name: $KB_NAME"
        exit 0
    fi

    local count
    count="$(echo "$ids" | sed '/^$/d' | wc -l | tr -d ' ')"
    print_info "Found $count matching KB entries"

    if [ -z "$KEEP_ID" ]; then
        KEEP_ID="$(echo "$ids" | head -n 1)"
        print_warning "No keep-id supplied; defaulting to first id: $KEEP_ID"
    fi

    if ! is_uuid "$KEEP_ID"; then
        print_error "Invalid keep-id: $KEEP_ID"
        exit 1
    fi

    local deleted=0 skipped=0 failed=0
    local kb_id
    while IFS= read -r kb_id; do
        [ -z "$kb_id" ] && continue
        if [ "$kb_id" = "$KEEP_ID" ]; then
            print_success "Keeping $kb_id"
            skipped=$((skipped + 1))
            continue
        fi

        if [ "$DRY_RUN" = true ]; then
            print_warning "Would delete duplicate KB: $kb_id"
            continue
        fi

        if delete_kb "$kb_id"; then
            print_success "Deleted duplicate KB: $kb_id"
            deleted=$((deleted + 1))
        else
            print_error "Failed to delete duplicate KB: $kb_id"
            failed=$((failed + 1))
        fi
    done <<<"$ids"

    print_section "Summary"
    print_info "Keep ID: $KEEP_ID"
    print_info "Deleted: $deleted"
    print_info "Kept: $skipped"
    print_info "Failed: $failed"
    if [ "$DRY_RUN" = true ]; then
        print_info "Mode: dry-run"
    else
        print_info "Mode: apply"
    fi

    if [ "$failed" -gt 0 ]; then
        exit 1
    fi
}

main "$@"
