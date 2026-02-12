#!/bin/bash

###############################################################################
# OpenWebUI Documents/RAG Tuner
# Applies a reproducible set of admin "Documents" settings via OpenWebUI APIs.
#
# Notes:
# - This changes persistent config stored in OpenWebUI's DB (not just env vars).
# - Avoid changing embedding model on an existing vector store unless you plan to
#   reset/re-ingest documents, otherwise retrieval quality can degrade or break.
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
    echo "║  OpenWebUI Documents Settings Tuner                        ║"
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
API_KEY="${OPENWEBUI_API_KEY:-${API_KEY:-}}"

APPLY_CHANGES=false
SNAPSHOT_DIR="${SNAPSHOT_DIR:-logs}"
RESTORE_FILE=""

# Tunables (defaults chosen for academic PDFs)
TUNE_TOP_K="${TUNE_TOP_K:-8}"
TUNE_CHUNK_SIZE="${TUNE_CHUNK_SIZE:-1800}"
TUNE_CHUNK_OVERLAP="${TUNE_CHUNK_OVERLAP:-200}"
TUNE_CHUNK_MIN_SIZE_TARGET="${TUNE_CHUNK_MIN_SIZE_TARGET:-900}"
TUNE_ALLOWED_EXTENSIONS="${TUNE_ALLOWED_EXTENSIONS:-pdf,txt,md}"
TUNE_EMBEDDING_BATCH_SIZE="${TUNE_EMBEDDING_BATCH_SIZE:-4}"

SHOW_HELP=false

show_help() {
    cat <<'EOF'
Usage: ./tune-openwebui-documents.sh [OPTIONS]

Options:
  --url URL                 OpenWebUI base URL (default: http://localhost:${WEBUI_PORT})
  --apply                   Apply changes (default: print + snapshot only)
  --restore FILE            Restore retrieval+embedding config from a snapshot JSON file
  --snapshot-dir DIR        Snapshot output directory (default: logs)

Tuning knobs (env overrides):
  TUNE_TOP_K
  TUNE_CHUNK_SIZE
  TUNE_CHUNK_OVERLAP
  TUNE_CHUNK_MIN_SIZE_TARGET
  TUNE_ALLOWED_EXTENSIONS   Comma-separated, e.g. "pdf,txt,md"
  TUNE_EMBEDDING_BATCH_SIZE

Auth:
  - Set OPENWEBUI_API_KEY (recommended) or API_KEY to an admin Bearer token.
  - This script intentionally does not attempt to sign in automatically.
    Keep auth manual-only to avoid hidden credential assumptions.

Notes:
  - This modifies OpenWebUI persistent config in the DB.
  - Avoid changing embedding model on an existing corpus unless you plan to re-ingest.
EOF
}

ensure_api_key() {
    if [ -n "$API_KEY" ]; then
        return 0
    fi

    print_error "Missing auth token. Export OPENWEBUI_API_KEY (or API_KEY) with an admin Bearer token."
    print_info "Admin UI path: ${OPENWEBUI_URL%/}/admin/settings/documents"
    print_info "Tip: If you already have a Bearer token, run:"
    print_info "  OPENWEBUI_API_KEY='<token>' ./tune-openwebui-documents.sh --apply"
    exit 1
}

api_get() {
    local path="$1"
    curl -sS -m 45 \
        -H "Authorization: Bearer $API_KEY" \
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

snapshot_configs() {
    local ts out_file retrieval embedding
    ts="$(date -u +%Y%m%dT%H%M%SZ)"
    mkdir -p "$SNAPSHOT_DIR"
    out_file="${SNAPSHOT_DIR%/}/openwebui-documents-snapshot-${ts}.json"

    retrieval="$(api_get "/api/v1/retrieval/config" 2>/dev/null || true)"
    embedding="$(api_get "/api/v1/retrieval/embedding" 2>/dev/null || true)"

    jq -cn --argjson retrieval "$retrieval" --argjson embedding "$embedding" \
        '{retrieval:$retrieval, embedding:$embedding}' >"$out_file"

    print_success "Snapshot written: $out_file"
    echo "$out_file"
}

apply_tuning() {
    local allowed_json retrieval_payload embedding_payload out_file code response

    allowed_json="$(echo "$TUNE_ALLOWED_EXTENSIONS" | tr ',' '\n' | sed 's/^ *//;s/ *$//' | rg -v '^$' | jq -R . | jq -s .)"

    retrieval_payload="$(jq -cn \
        --argjson allowed "$allowed_json" \
        --argjson top_k "$TUNE_TOP_K" \
        --argjson chunk_size "$TUNE_CHUNK_SIZE" \
        --argjson chunk_overlap "$TUNE_CHUNK_OVERLAP" \
        --argjson chunk_min "$TUNE_CHUNK_MIN_SIZE_TARGET" \
        '{TOP_K:$top_k,CHUNK_SIZE:$chunk_size,CHUNK_OVERLAP:$chunk_overlap,CHUNK_MIN_SIZE_TARGET:$chunk_min,ALLOWED_FILE_EXTENSIONS:$allowed}')"

    embedding_payload="$(jq -cn \
        --arg engine "" \
        --arg model "sentence-transformers/all-MiniLM-L6-v2" \
        --argjson batch "$TUNE_EMBEDDING_BATCH_SIZE" \
        --argjson enable_async true \
        '{RAG_EMBEDDING_ENGINE:$engine,RAG_EMBEDDING_MODEL:$model,RAG_EMBEDDING_BATCH_SIZE:$batch,ENABLE_ASYNC_EMBEDDING:$enable_async}')"

    print_step "Applying retrieval/documents tuning"
    out_file="$(mktemp)"
    code="$(api_post_json_with_status "/api/v1/retrieval/config/update" "$retrieval_payload" "$out_file")"
    response="$(cat "$out_file" 2>/dev/null || true)"
    rm -f "$out_file"

    if [ -z "$code" ] || [ "$code" -lt 200 ] || [ "$code" -ge 300 ]; then
        print_error "Retrieval config update failed (HTTP ${code:-unknown})"
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)"
        exit 1
    fi
    print_success "Retrieval config updated"

    print_step "Applying embedding throughput tuning"
    out_file="$(mktemp)"
    code="$(api_post_json_with_status "/api/v1/retrieval/embedding/update" "$embedding_payload" "$out_file")"
    response="$(cat "$out_file" 2>/dev/null || true)"
    rm -f "$out_file"

    if [ -z "$code" ] || [ "$code" -lt 200 ] || [ "$code" -ge 300 ]; then
        print_error "Embedding config update failed (HTTP ${code:-unknown})"
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)"
        exit 1
    fi
    print_success "Embedding config updated"
}

restore_from_snapshot() {
    local file="$1"
    local retrieval embedding out_file code response

    if [ ! -f "$file" ]; then
        print_error "Snapshot file not found: $file"
        exit 1
    fi

    retrieval="$(jq -c '.retrieval' "$file" 2>/dev/null || true)"
    embedding="$(jq -c '.embedding' "$file" 2>/dev/null || true)"

    if [ -z "$retrieval" ] || [ -z "$embedding" ]; then
        print_error "Snapshot missing retrieval or embedding sections"
        exit 1
    fi

    # Restore retrieval config via /config/update. The endpoint expects only the fields
    # it knows; sending the whole object is fine but includes extra keys, so we filter.
    retrieval="$(echo "$retrieval" | jq -c '{
        RAG_TEMPLATE,
        TOP_K,
        BYPASS_EMBEDDING_AND_RETRIEVAL,
        RAG_FULL_CONTEXT,
        ENABLE_RAG_HYBRID_SEARCH,
        ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS,
        TOP_K_RERANKER,
        RELEVANCE_THRESHOLD,
        HYBRID_BM25_WEIGHT,
        CONTENT_EXTRACTION_ENGINE,
        PDF_EXTRACT_IMAGES,
        TEXT_SPLITTER,
        ENABLE_MARKDOWN_HEADER_TEXT_SPLITTER,
        CHUNK_SIZE,
        CHUNK_MIN_SIZE_TARGET,
        CHUNK_OVERLAP,
        FILE_MAX_SIZE,
        FILE_MAX_COUNT,
        FILE_IMAGE_COMPRESSION_WIDTH,
        FILE_IMAGE_COMPRESSION_HEIGHT,
        ALLOWED_FILE_EXTENSIONS
    }')"

    embedding="$(echo "$embedding" | jq -c '{
        RAG_EMBEDDING_ENGINE,
        RAG_EMBEDDING_MODEL,
        RAG_EMBEDDING_BATCH_SIZE,
        ENABLE_ASYNC_EMBEDDING,
        openai_config,
        ollama_config,
        azure_openai_config
    }')"

    print_step "Restoring retrieval config"
    out_file="$(mktemp)"
    code="$(api_post_json_with_status "/api/v1/retrieval/config/update" "$retrieval" "$out_file")"
    response="$(cat "$out_file" 2>/dev/null || true)"
    rm -f "$out_file"
    if [ -z "$code" ] || [ "$code" -lt 200 ] || [ "$code" -ge 300 ]; then
        print_error "Retrieval restore failed (HTTP ${code:-unknown})"
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)"
        exit 1
    fi
    print_success "Retrieval restored"

    print_step "Restoring embedding config"
    out_file="$(mktemp)"
    code="$(api_post_json_with_status "/api/v1/retrieval/embedding/update" "$embedding" "$out_file")"
    response="$(cat "$out_file" 2>/dev/null || true)"
    rm -f "$out_file"
    if [ -z "$code" ] || [ "$code" -lt 200 ] || [ "$code" -ge 300 ]; then
        print_error "Embedding restore failed (HTTP ${code:-unknown})"
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 200)"
        exit 1
    fi
    print_success "Embedding restored"
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --url)
                OPENWEBUI_URL="$2"
                shift 2
                ;;
            --apply)
                APPLY_CHANGES=true
                shift 1
                ;;
            --snapshot-dir)
                SNAPSHOT_DIR="$2"
                shift 2
                ;;
            --restore)
                RESTORE_FILE="$2"
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
    print_info "Apply changes: $APPLY_CHANGES"
    print_info "Snapshot dir: $SNAPSHOT_DIR"

    if ! curl -sS -m 10 -f -o /dev/null "${OPENWEBUI_URL%/}/health" 2>/dev/null; then
        print_error "OpenWebUI health check failed at ${OPENWEBUI_URL%/}/health"
        exit 1
    fi
    print_success "OpenWebUI health check OK"

    ensure_api_key

    print_step "Snapshot (current)"
    local snapshot_file=""
    snapshot_file="$(snapshot_configs)"

    if [ -n "$RESTORE_FILE" ]; then
        print_step "Restore"
        restore_from_snapshot "$RESTORE_FILE"
    elif [ "$APPLY_CHANGES" = true ]; then
        print_step "Apply"
        apply_tuning
    else
        print_warning "Dry-run mode: no changes applied. Re-run with --apply to apply tuning."
    fi

    print_step "Verify (current)"
    api_get "/api/v1/retrieval/config" | jq '{TOP_K,CHUNK_SIZE,CHUNK_OVERLAP,CHUNK_MIN_SIZE_TARGET,ALLOWED_FILE_EXTENSIONS,ENABLE_MARKDOWN_HEADER_TEXT_SPLITTER,RAG_FULL_CONTEXT}'
    api_get "/api/v1/retrieval/embedding" | jq '{RAG_EMBEDDING_ENGINE,RAG_EMBEDDING_MODEL,RAG_EMBEDDING_BATCH_SIZE,ENABLE_ASYNC_EMBEDDING}'

    print_step "Notes"
    print_info "Snapshot saved at: $snapshot_file"
    print_info "Admin UI path: ${OPENWEBUI_URL%/}/admin/settings/documents"
}

main "$@"
