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

SOURCE_ROOT="${SOURCE_ROOT:-/LAB/@thesis/openwebui/prod_max}"
FILE_NAME="${FILE_NAME:-molecules_high_confidence.jsonl}"
KB_NAME="${KB_NAME:-Molecular SMILES High Confidence}"
KB_DESCRIPTION="${KB_DESCRIPTION:-High-confidence molecular extractions for retrieval}"
KB_ID_OVERRIDE="${KB_ID_OVERRIDE:-}"
LIMIT_COUNT="${LIMIT_COUNT:-0}"
PAPER_FILTER="${PAPER_FILTER:-}"
PROCESS_TIMEOUT_SECONDS="${PROCESS_TIMEOUT_SECONDS:-900}"
PROCESS_POLL_SECONDS="${PROCESS_POLL_SECONDS:-3}"

SHOW_HELP=false

show_help() {
    cat <<'EOF'
Usage: ./import-smiles-to-kb.sh [OPTIONS]

Options:
  --source-root PATH      Root directory containing per-paper folders (default: /LAB/@thesis/openwebui/prod_max)
  --file-name NAME        Artifact filename in each paper dir (default: molecules_high_confidence.jsonl)
  --paper-filter TEXT     Keep only paper directories containing TEXT
  --limit N               Process at most N files (0 = no limit)
  --kb-name NAME          Knowledge Base name
  --kb-description TEXT   Knowledge Base description
  --url URL               OpenWebUI URL override
  --timeout SECONDS       Processing timeout per file (default: 900)
  --poll SECONDS          Poll interval seconds (default: 3)
  -h, --help              Show help

Auth:
  - Uses OPENWEBUI_API_KEY if set.
  - Else signs in using OPENWEBUI_SIGNIN_EMAIL / OPENWEBUI_SIGNIN_PASSWORD.
EOF
}

is_uuid() {
    local value="${1:-}"
    echo "$value" | rg -q '^[0-9a-fA-F-]{36}$'
}

kb_cache_file() {
    local slug
    slug="$(echo "$KB_NAME" | tr '[:upper:]' '[:lower:]' | tr -cs 'a-z0-9' '-')"
    echo "$SCRIPT_DIR/.state/${slug}.kb_id"
}

cached_kb_id() {
    local path
    path="$(kb_cache_file)"
    if [ -f "$path" ]; then
        cat "$path"
    fi
}

store_kb_id_cache() {
    local kb_id="$1"
    local path
    path="$(kb_cache_file)"
    mkdir -p "$(dirname "$path")"
    printf "%s" "$kb_id" >"$path"
}

kb_id_exists() {
    local kb_id="$1"
    local code
    code="$(curl -sS -m 20 -H "Authorization: Bearer $API_KEY" -o /dev/null -w "%{http_code}" "${OPENWEBUI_URL%/}/api/v1/knowledge/$kb_id" 2>/dev/null || true)"
    [ "$code" = "200" ]
}

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        print_error "Missing dependency: $cmd"
        exit 1
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

api_get() {
    local path="$1"
    curl -sS -m 45 -H "Authorization: Bearer $API_KEY" "${OPENWEBUI_URL%/}${path}"
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

    if [ -n "$KB_ID_OVERRIDE" ] && is_uuid "$KB_ID_OVERRIDE" && kb_id_exists "$KB_ID_OVERRIDE"; then
        print_success "Using override Knowledge Base ID: $KB_ID_OVERRIDE" >&2
        store_kb_id_cache "$KB_ID_OVERRIDE"
        echo "$KB_ID_OVERRIDE"
        return 0
    fi

    kb_id="$(cached_kb_id)"
    if [ -n "$kb_id" ] && is_uuid "$kb_id" && kb_id_exists "$kb_id"; then
        print_success "Using cached Knowledge Base: $KB_NAME ($kb_id)" >&2
        echo "$kb_id"
        return 0
    fi

    kb_list="$(api_get "/api/v1/knowledge/" 2>/dev/null || true)"
    if [ -n "$kb_list" ]; then
        kb_id="$(echo "$kb_list" | jq -r --arg name "$KB_NAME" '
            if type == "array" then
                .[]
            elif (.data? | type) == "array" then
                .data[]
            elif (.items? | type) == "array" then
                .items[]
            else
                empty
            end
            | select(((.name // "") | ascii_downcase) == ($name | ascii_downcase))
            | .id
        ' 2>/dev/null | head -n 1 || true)"
        if [ -n "$kb_id" ] && [ "$kb_id" != "null" ]; then
            print_success "Using existing Knowledge Base: $KB_NAME ($kb_id)" >&2
            store_kb_id_cache "$kb_id"
            echo "$kb_id"
            return 0
        fi
    fi

    local out code resp payload
    payload="$(jq -cn --arg name "$KB_NAME" --arg desc "$KB_DESCRIPTION" '{name:$name, description:$desc}')"
    out="$(mktemp)"
    code="$(api_post_json_with_status "/api/v1/knowledge/create" "$payload" "$out")"
    resp="$(cat "$out" 2>/dev/null || true)"
    rm -f "$out"

    if [ -z "$code" ] || ! echo "$code" | rg -q '^[0-9]{3}$' || [ "$code" -lt 200 ] || [ "$code" -ge 300 ]; then
        print_error "Knowledge Base create failed (HTTP ${code:-unknown})"
        exit 1
    fi

    kb_id="$(echo "$resp" | jq -r '.id // .knowledge_id // empty' 2>/dev/null || true)"
    if [ -z "$kb_id" ]; then
        print_error "Failed to parse Knowledge Base ID"
        exit 1
    fi

    print_success "Created Knowledge Base: $KB_NAME ($kb_id)" >&2
    store_kb_id_cache "$kb_id"
    echo "$kb_id"
}

upload_file() {
    local file_path="$1"
    local filename
    filename="$(basename "$file_path")"

    curl -sS -m 300 -X POST \
        "${OPENWEBUI_URL%/}/api/v1/files/?process=true&process_in_background=false" \
        -H "Authorization: Bearer $API_KEY" \
        -F "file=@$file_path" \
        -F "metadata={\"name\":\"$filename\"}" 2>/dev/null || true
}

wait_for_processing() {
    local file_id="$1"
    local deadline=$((SECONDS + PROCESS_TIMEOUT_SECONDS))
    while [ $SECONDS -lt $deadline ]; do
        local resp status
        resp="$(api_get "/api/v1/files/$file_id/process/status" 2>/dev/null || true)"
        status="$(echo "$resp" | jq -r '.status // empty' 2>/dev/null || true)"
        if [ "$status" = "completed" ]; then
            return 0
        fi
        if [ "$status" = "failed" ]; then
            return 1
        fi
        sleep "$PROCESS_POLL_SECONDS"
    done
    return 1
}

attach_file_to_kb() {
    local kb_id="$1"
    local file_id="$2"
    local payload out code response detail
    payload="$(jq -cn --arg file_id "$file_id" '{file_id:$file_id}')"
    out="$(mktemp)"
    code="$(api_post_json_with_status "/api/v1/knowledge/$kb_id/file/add" "$payload" "$out")"
    response="$(cat "$out" 2>/dev/null || true)"
    rm -f "$out"
    if [ "$code" = "400" ]; then
        detail="$(echo "$response" | jq -r '.detail // empty' 2>/dev/null || true)"
        if echo "$detail" | rg -qi 'duplicate content'; then
            print_warning "Duplicate content skipped for file_id=$file_id"
            return 0
        fi
    fi
    if [ -z "$code" ] || ! echo "$code" | rg -q '^[0-9]{3}$' || [ "$code" -lt 200 ] || [ "$code" -ge 300 ]; then
        if [ -n "$response" ]; then
            print_error "Attach API failed (HTTP ${code:-unknown}): $(echo "$response" | head -c 240)"
        fi
        return 1
    fi
    return 0
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --source-root)
                SOURCE_ROOT="$2"
                shift 2
                ;;
            --file-name)
                FILE_NAME="$2"
                shift 2
                ;;
            --paper-filter)
                PAPER_FILTER="$2"
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
            --kb-id)
                KB_ID_OVERRIDE="$2"
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
            -h | --help)
                SHOW_HELP=true
                shift
                ;;
            *)
                print_error "Unknown argument: $1"
                SHOW_HELP=true
                shift
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

    print_header "OpenWebUI SMILES Artifact Import"
    print_info "OpenWebUI URL: $OPENWEBUI_URL"
    print_info "Source root: $SOURCE_ROOT"
    print_info "Artifact file: $FILE_NAME"
    print_info "Knowledge Base: $KB_NAME"

    if ! curl -sS -m 10 -f -o /dev/null "${OPENWEBUI_URL%/}/health" 2>/dev/null; then
        print_error "OpenWebUI health check failed at ${OPENWEBUI_URL%/}/health"
        exit 1
    fi

    ensure_api_key
    local kb_id
    kb_id="$(find_or_create_kb)"
    if [ -z "$kb_id" ] || ! is_uuid "$kb_id"; then
        print_error "Invalid KB ID: $kb_id"
        exit 1
    fi

    local -a files=()
    while IFS= read -r f; do
        files+=("$f")
    done < <(find "$SOURCE_ROOT" -mindepth 2 -maxdepth 2 -type f -name "$FILE_NAME" | sort)

    if [ -n "$PAPER_FILTER" ]; then
        local -a filtered=()
        local file
        for file in "${files[@]}"; do
            if echo "$file" | rg -qi "$PAPER_FILTER"; then
                filtered+=("$file")
            fi
        done
        files=("${filtered[@]}")
    fi

    if [ ${#files[@]} -eq 0 ]; then
        print_warning "No matching artifacts found"
        exit 0
    fi

    local total=0 uploaded=0 processed=0 attached=0 failed=0
    local artifact response file_id
    for artifact in "${files[@]}"; do
        if [ "$LIMIT_COUNT" -gt 0 ] && [ "$total" -ge "$LIMIT_COUNT" ]; then
            break
        fi
        total=$((total + 1))
        print_step "[$total] Uploading $(basename "$artifact")"
        response="$(upload_file "$artifact")"
        file_id="$(echo "$response" | jq -r '.id // .file_id // empty' 2>/dev/null || true)"
        if [ -z "$file_id" ]; then
            failed=$((failed + 1))
            print_error "Upload failed: $artifact"
            continue
        fi
        uploaded=$((uploaded + 1))

        if wait_for_processing "$file_id"; then
            processed=$((processed + 1))
        else
            failed=$((failed + 1))
            print_error "Processing failed: $artifact"
            continue
        fi

        if attach_file_to_kb "$kb_id" "$file_id"; then
            attached=$((attached + 1))
            print_success "Attached: $(basename "$artifact")"
        else
            failed=$((failed + 1))
            print_error "Attach failed: $artifact"
        fi
    done

    print_section "Summary"
    print_info "Artifacts seen: $total"
    print_info "Uploaded: $uploaded"
    print_info "Processed: $processed"
    print_info "Attached: $attached"
    print_info "Failed: $failed"
    print_info "KB: $KB_NAME ($kb_id)"

    if [ "$failed" -gt 0 ]; then
        exit 1
    fi
    print_success "SMILES artifact import complete"
}

main "$@"
