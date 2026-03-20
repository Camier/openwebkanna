#!/bin/bash
#
# Sync OpenWebUI persistent retrieval config with local .env values.
#
# Why:
# - OpenWebUI stores retrieval settings in webui.db.
# - If .env changes (e.g. RAG_TOP_K or CONTENT_EXTRACTION_ENGINE) but webui.db
#   is not updated, multimodal and SMILES ingestion can silently drift.
#
# What it does:
# - Authenticates to OpenWebUI admin API (token env, signin, or local JWT mint).
# - Reads current /api/v1/retrieval/config.
# - Updates retrieval keys used by ingestion/retrieval and writes back via
#   /api/v1/retrieval/config/update.
#

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
cd "$SCRIPT_DIR"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_SERVICE="${OPENWEBUI_SERVICE:-${OPENWEBUI_DOCKER_SERVICE:-openwebui}}"
OPENWEBUI_CONTAINER_NAME="${OPENWEBUI_CONTAINER_NAME:-${OPENWEBUI_DOCKER_CONTAINER:-$OPENWEBUI_SERVICE}}"
OPENWEBUI_TOKEN_INPUT="${OPENWEBUI_TOKEN:-${OPENWEBUI_API_KEY:-}}"
CURL_TIMEOUT="${CURL_TIMEOUT:-45}"

TOP_K_RAW="${RAG_TOP_K:-}"
CHUNK_SIZE_RAW="${CHUNK_SIZE:-}"
CHUNK_OVERLAP_RAW="${CHUNK_OVERLAP:-}"
TOP_K_RERANKER_RAW="${RAG_TOP_K_RERANKER:-}"
RAG_RELEVANCE_THRESHOLD_RAW="${RAG_RELEVANCE_THRESHOLD:-}"
RAG_SYSTEM_CONTEXT_RAW="${RAG_SYSTEM_CONTEXT:-}"
ENABLE_RAG_HYBRID_SEARCH_RAW="${ENABLE_RAG_HYBRID_SEARCH:-}"
ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS_RAW="${ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS:-}"
RAG_HYBRID_BM25_WEIGHT_RAW="${RAG_HYBRID_BM25_WEIGHT:-}"
RETRIEVAL_FUSION_STRATEGY="${RETRIEVAL_FUSION_STRATEGY:-}"
RETRIEVAL_FUSION_RRF_K_RAW="${RETRIEVAL_FUSION_RRF_K:-}"
RETRIEVAL_FUSION_TOP_K_RAW="${RETRIEVAL_FUSION_TOP_K:-}"
RETRIEVAL_CHANNEL_WEIGHTS="${RETRIEVAL_CHANNEL_WEIGHTS:-}"
CONTENT_EXTRACTION_ENGINE="${CONTENT_EXTRACTION_ENGINE:-docling}"
PDF_EXTRACT_IMAGES_RAW="${PDF_EXTRACT_IMAGES:-true}"

API_TOKEN=""

to_bool_json_or_null() {
    local value="$1"
    if [ -z "$value" ]; then
        printf "null"
        return 0
    fi
    if is_true "$value"; then
        printf "true"
    else
        printf "false"
    fi
}

to_nonneg_int_or_null() {
    local value="$1"
    if [[ $value =~ ^[0-9]+$ ]]; then
        printf "%s" "$value"
    else
        printf "null"
    fi
}

to_float_or_null() {
    local value="$1"
    if [[ $value =~ ^[0-9]+([.][0-9]+)?$ ]] || [[ $value =~ ^[0-9]*[.][0-9]+$ ]]; then
        printf "%s" "$value"
    else
        printf "null"
    fi
}

resolve_token() {
    API_TOKEN="$(
        resolve_openwebui_api_token \
            "$OPENWEBUI_TOKEN_INPUT" \
            "$OPENWEBUI_URL" \
            "$OPENWEBUI_SIGNIN_EMAIL" \
            "$OPENWEBUI_SIGNIN_PASSWORD" \
            "$CURL_TIMEOUT" \
            "$OPENWEBUI_SERVICE" \
            "$OPENWEBUI_CONTAINER_NAME" \
            "/api/v1/auths/signin" || true
    )"
    [ -n "$API_TOKEN" ]
}

main() {
    print_header "Sync OpenWebUI Retrieval Config"

    if ! command_exists jq; then
        print_error "jq is required"
        exit 1
    fi

    if ! resolve_token; then
        print_error "Unable to authenticate to OpenWebUI (set OPENWEBUI_API_KEY/OPENWEBUI_TOKEN or ensure WEBUI_SECRET_KEY is available)"
        exit 1
    fi

    local cfg_file cfg_code payload_file payload response_file response_code
    local top_k chunk_size chunk_overlap top_k_reranker
    local relevance_threshold bm25_weight
    local hybrid_search hybrid_search_enriched
    local fusion_k fusion_top_k pdf_extract_images
    cfg_file="$(mktemp)"
    payload_file="$(mktemp)"
    response_file="$(mktemp)"

    cfg_code="$(curl -sS -m "$CURL_TIMEOUT" \
        -H "Authorization: Bearer $API_TOKEN" \
        -H "Accept: application/json" \
        -o "$cfg_file" -w "%{http_code}" \
        "${OPENWEBUI_URL%/}/api/v1/retrieval/config" || true)"

    if [ "$cfg_code" != "200" ]; then
        print_error "Failed to read retrieval config (HTTP $cfg_code)"
        rm -f "$cfg_file" "$payload_file" "$response_file"
        exit 1
    fi

    top_k="$(to_nonneg_int_or_null "$TOP_K_RAW")"
    chunk_size="$(to_nonneg_int_or_null "$CHUNK_SIZE_RAW")"
    chunk_overlap="$(to_nonneg_int_or_null "$CHUNK_OVERLAP_RAW")"
    top_k_reranker="$(to_nonneg_int_or_null "$TOP_K_RERANKER_RAW")"
    relevance_threshold="$(to_float_or_null "$RAG_RELEVANCE_THRESHOLD_RAW")"
    bm25_weight="$(to_float_or_null "$RAG_HYBRID_BM25_WEIGHT_RAW")"
    hybrid_search="$(to_bool_json_or_null "$ENABLE_RAG_HYBRID_SEARCH_RAW")"
    hybrid_search_enriched="$(to_bool_json_or_null "$ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS_RAW")"
    pdf_extract_images="$(to_bool_json_or_null "$PDF_EXTRACT_IMAGES_RAW")"
    fusion_k="$(to_nonneg_int_or_null "$RETRIEVAL_FUSION_RRF_K_RAW")"
    fusion_top_k="$(to_nonneg_int_or_null "$RETRIEVAL_FUSION_TOP_K_RAW")"

    payload="$(jq -c \
        --argjson top_k "$top_k" \
        --argjson chunk_size "$chunk_size" \
        --argjson chunk_overlap "$chunk_overlap" \
        --argjson top_k_reranker "$top_k_reranker" \
        --argjson relevance_threshold "$relevance_threshold" \
        --argjson rag_system_context_json "$(to_bool_json_or_null "$RAG_SYSTEM_CONTEXT_RAW")" \
        --argjson enable_hybrid_search "$hybrid_search" \
        --argjson enable_hybrid_search_enriched_texts "$hybrid_search_enriched" \
        --argjson bm25_weight "$bm25_weight" \
        --arg fusion_strategy "$RETRIEVAL_FUSION_STRATEGY" \
        --argjson fusion_k "$fusion_k" \
        --argjson fusion_top_k "$fusion_top_k" \
        --argjson pdf_extract_images "$pdf_extract_images" \
        --arg extraction_engine "$CONTENT_EXTRACTION_ENGINE" \
        --arg channel_weights "$RETRIEVAL_CHANNEL_WEIGHTS" \
        '
        if $top_k != null then
          .TOP_K = $top_k | .RAG_TOP_K = $top_k
        else
          .
        end
        | if $chunk_size != null then .CHUNK_SIZE = $chunk_size else . end
        | if $chunk_overlap != null then .CHUNK_OVERLAP = $chunk_overlap else . end
        | if $top_k_reranker != null then .RAG_TOP_K_RERANKER = $top_k_reranker else . end
        | if $relevance_threshold != null then .RAG_RELEVANCE_THRESHOLD = $relevance_threshold else . end
        | if $rag_system_context_json != null then .RAG_SYSTEM_CONTEXT = $rag_system_context_json else . end
        | if $enable_hybrid_search != null then .ENABLE_RAG_HYBRID_SEARCH = $enable_hybrid_search else . end
        | if $enable_hybrid_search_enriched_texts != null then .ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS = $enable_hybrid_search_enriched_texts else . end
        | if $bm25_weight != null then .RAG_HYBRID_BM25_WEIGHT = $bm25_weight else . end
        | if $fusion_strategy != "" then .RETRIEVAL_FUSION_STRATEGY = $fusion_strategy else . end
        | if $fusion_k != null then .RETRIEVAL_FUSION_RRF_K = $fusion_k else . end
        | if $fusion_top_k != null then .RETRIEVAL_FUSION_TOP_K = $fusion_top_k else . end
        | if $channel_weights != "" then .RETRIEVAL_CHANNEL_WEIGHTS = $channel_weights else . end
        | if $extraction_engine != "" then .CONTENT_EXTRACTION_ENGINE = $extraction_engine else . end
        | if $pdf_extract_images != null then .PDF_EXTRACT_IMAGES = $pdf_extract_images else . end
        ' "$cfg_file")"

    if [ -z "$payload" ]; then
        print_error "Failed to build retrieval config payload"
        rm -f "$cfg_file" "$payload_file" "$response_file"
        exit 1
    fi

    printf "%s" "$payload" >"$payload_file"

    local before after
    before="$(jq -S -c '
        {
          TOP_K: .TOP_K,
          RAG_TOP_K: .RAG_TOP_K,
          CHUNK_SIZE: .CHUNK_SIZE,
          CHUNK_OVERLAP: .CHUNK_OVERLAP,
          RAG_TOP_K_RERANKER: .RAG_TOP_K_RERANKER,
          RAG_RELEVANCE_THRESHOLD: .RAG_RELEVANCE_THRESHOLD,
          RAG_SYSTEM_CONTEXT: .RAG_SYSTEM_CONTEXT,
          ENABLE_RAG_HYBRID_SEARCH: .ENABLE_RAG_HYBRID_SEARCH,
          ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS: .ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS,
          RAG_HYBRID_BM25_WEIGHT: .RAG_HYBRID_BM25_WEIGHT,
          RETRIEVAL_FUSION_STRATEGY: .RETRIEVAL_FUSION_STRATEGY,
          RETRIEVAL_FUSION_RRF_K: .RETRIEVAL_FUSION_RRF_K,
          RETRIEVAL_FUSION_TOP_K: .RETRIEVAL_FUSION_TOP_K,
          RETRIEVAL_CHANNEL_WEIGHTS: .RETRIEVAL_CHANNEL_WEIGHTS,
          CONTENT_EXTRACTION_ENGINE: .CONTENT_EXTRACTION_ENGINE,
          PDF_EXTRACT_IMAGES: .PDF_EXTRACT_IMAGES
        }' "$cfg_file")"
    after="$(jq -S -c '
        {
          TOP_K: .TOP_K,
          RAG_TOP_K: .RAG_TOP_K,
          CHUNK_SIZE: .CHUNK_SIZE,
          CHUNK_OVERLAP: .CHUNK_OVERLAP,
          RAG_TOP_K_RERANKER: .RAG_TOP_K_RERANKER,
          RAG_RELEVANCE_THRESHOLD: .RAG_RELEVANCE_THRESHOLD,
          RAG_SYSTEM_CONTEXT: .RAG_SYSTEM_CONTEXT,
          ENABLE_RAG_HYBRID_SEARCH: .ENABLE_RAG_HYBRID_SEARCH,
          ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS: .ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS,
          RAG_HYBRID_BM25_WEIGHT: .RAG_HYBRID_BM25_WEIGHT,
          RETRIEVAL_FUSION_STRATEGY: .RETRIEVAL_FUSION_STRATEGY,
          RETRIEVAL_FUSION_RRF_K: .RETRIEVAL_FUSION_RRF_K,
          RETRIEVAL_FUSION_TOP_K: .RETRIEVAL_FUSION_TOP_K,
          RETRIEVAL_CHANNEL_WEIGHTS: .RETRIEVAL_CHANNEL_WEIGHTS,
          CONTENT_EXTRACTION_ENGINE: .CONTENT_EXTRACTION_ENGINE,
          PDF_EXTRACT_IMAGES: .PDF_EXTRACT_IMAGES
        }' "$payload_file")"

    if [ "$before" = "$after" ]; then
        print_info "Retrieval config already in sync"
        rm -f "$cfg_file" "$payload_file" "$response_file"
        exit 0
    fi

    print_step "Updating OpenWebUI retrieval config"
    response_code="$(curl -sS -m "$CURL_TIMEOUT" \
        -H "Authorization: Bearer $API_TOKEN" \
        -H "Content-Type: application/json" \
        -X POST --data @"$payload_file" \
        -o "$response_file" -w "%{http_code}" \
        "${OPENWEBUI_URL%/}/api/v1/retrieval/config/update" || true)"

    if [ -z "$response_code" ] || [ "$response_code" -lt 200 ] || [ "$response_code" -ge 300 ]; then
        print_error "Failed to update retrieval config (HTTP ${response_code:-unknown})"
        print_info "Response: $(head -c 240 "$response_file" 2>/dev/null || true)"
        rm -f "$cfg_file" "$payload_file" "$response_file"
        exit 1
    fi

    print_success "OpenWebUI retrieval config synced"
    if [ -n "$TOP_K_RAW" ]; then
        print_info "TOP_K=$TOP_K_RAW"
    fi
    if [ -n "$CHUNK_SIZE_RAW" ]; then
        print_info "CHUNK_SIZE=$CHUNK_SIZE_RAW"
    fi
    if [ -n "$CHUNK_OVERLAP_RAW" ]; then
        print_info "CHUNK_OVERLAP=$CHUNK_OVERLAP_RAW"
    fi
    if [ -n "$CONTENT_EXTRACTION_ENGINE" ]; then
        print_info "CONTENT_EXTRACTION_ENGINE=$CONTENT_EXTRACTION_ENGINE"
    fi
    if [ -n "$PDF_EXTRACT_IMAGES_RAW" ]; then
        print_info "PDF_EXTRACT_IMAGES=$PDF_EXTRACT_IMAGES_RAW"
    fi

    rm -f "$cfg_file" "$payload_file" "$response_file"
}

main "$@"
