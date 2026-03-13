#!/bin/bash
#
# Sync OpenWebUI persistent retrieval web-search config with local .env values.
#
# Why:
# - OpenWebUI stores retrieval/web settings in webui.db.
# - If .env changes (e.g. SEARXNG_QUERY_URL) but DB still has stale values,
#   /api/v1/retrieval/process/web/search can timeout or fail.
#
# What it does:
# - Authenticates to OpenWebUI admin API (token env, signin, or local JWT mint).
# - Reads current /api/v1/retrieval/config.
# - Updates only web-search related fields and writes them back via
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

WEB_ENABLED=false
WEB_ENGINE="${WEB_SEARCH_ENGINE:-${WEBSEARCH_ENGINE:-searxng}}"
SEARX_URL="${SEARXNG_QUERY_URL:-}"
SEARX_LANGUAGE="${SEARXNG_LANGUAGE:-}"
RESULT_COUNT_RAW="${WEB_SEARCH_RESULT_COUNT:-}"
CONCURRENT_RAW="${WEB_SEARCH_CONCURRENT_REQUESTS:-}"
SSL_VERIFY_RAW="${ENABLE_WEB_LOADER_SSL_VERIFICATION:-}"

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

to_nonneg_int_or_sentinel() {
    local value="$1"
    if [[ $value =~ ^[0-9]+$ ]]; then
        printf "%s" "$value"
    else
        printf -- "-1"
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
    print_header "Sync OpenWebUI Web Search Config"

    if is_true "${ENABLE_WEB_SEARCH:-false}" || is_true "${ENABLE_WEBSEARCH:-false}"; then
        WEB_ENABLED=true
    fi

    if [ "$WEB_ENABLED" = true ] && [ "$WEB_ENGINE" = "searxng" ] && [ -z "$SEARX_URL" ]; then
        print_error "Web search is enabled with searxng engine but SEARXNG_QUERY_URL is empty"
        exit 1
    fi

    if ! resolve_token; then
        print_error "Unable to authenticate to OpenWebUI (set OPENWEBUI_API_KEY/OPENWEBUI_TOKEN or ensure WEBUI_SECRET_KEY is available)"
        exit 1
    fi

    local cfg_file cfg_code payload_file payload response_file response_code
    local result_count concurrent ssl_verify_json
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

    result_count="$(to_nonneg_int_or_sentinel "$RESULT_COUNT_RAW")"
    concurrent="$(to_nonneg_int_or_sentinel "$CONCURRENT_RAW")"
    ssl_verify_json="$(to_bool_json_or_null "$SSL_VERIFY_RAW")"

    payload="$(jq -c \
        --arg engine "$WEB_ENGINE" \
        --arg searx_url "$SEARX_URL" \
        --arg searx_language "$SEARX_LANGUAGE" \
        --argjson web_enabled "$(if [ "$WEB_ENABLED" = true ]; then printf "true"; else printf "false"; fi)" \
        --argjson result_count "$result_count" \
        --argjson concurrent "$concurrent" \
        --argjson ssl_verify "$ssl_verify_json" \
        '
        .web.ENABLE_WEB_SEARCH = $web_enabled
        | .web.WEB_SEARCH_ENGINE = $engine
        | (if $searx_url != "" then .web.SEARXNG_QUERY_URL = $searx_url else . end)
        | (if $searx_language != "" then .web.SEARXNG_LANGUAGE = $searx_language else . end)
        | (if $result_count >= 0 then .web.WEB_SEARCH_RESULT_COUNT = $result_count else . end)
        | (if $concurrent >= 0 then .web.WEB_SEARCH_CONCURRENT_REQUESTS = $concurrent else . end)
        | (if $ssl_verify == null then . else .web.ENABLE_WEB_LOADER_SSL_VERIFICATION = $ssl_verify end)
        ' "$cfg_file")"

    if [ -z "$payload" ]; then
        print_error "Failed to build retrieval web config payload"
        rm -f "$cfg_file" "$payload_file" "$response_file"
        exit 1
    fi

    printf "%s" "$payload" >"$payload_file"

    local before after
    before="$(jq -c '.web | {ENABLE_WEB_SEARCH,WEB_SEARCH_ENGINE,SEARXNG_QUERY_URL,SEARXNG_LANGUAGE,WEB_SEARCH_RESULT_COUNT,WEB_SEARCH_CONCURRENT_REQUESTS,ENABLE_WEB_LOADER_SSL_VERIFICATION}' "$cfg_file")"
    after="$(jq -c '.web | {ENABLE_WEB_SEARCH,WEB_SEARCH_ENGINE,SEARXNG_QUERY_URL,SEARXNG_LANGUAGE,WEB_SEARCH_RESULT_COUNT,WEB_SEARCH_CONCURRENT_REQUESTS,ENABLE_WEB_LOADER_SSL_VERIFICATION}' "$payload_file")"

    if [ "$before" = "$after" ]; then
        print_info "Retrieval web-search config already in sync"
        rm -f "$cfg_file" "$payload_file" "$response_file"
        exit 0
    fi

    print_step "Updating OpenWebUI retrieval web-search config"
    response_code="$(curl -sS -m "$CURL_TIMEOUT" \
        -H "Authorization: Bearer $API_TOKEN" \
        -H "Content-Type: application/json" \
        -X POST --data @"$payload_file" \
        -o "$response_file" -w "%{http_code}" \
        "${OPENWEBUI_URL%/}/api/v1/retrieval/config/update" || true)"

    if [ -z "$response_code" ] || [ "$response_code" -lt 200 ] || [ "$response_code" -ge 300 ]; then
        print_error "Failed to update retrieval config (HTTP ${response_code:-unknown})"
        print_info "Response: $(head -c 200 "$response_file" 2>/dev/null || true)"
        rm -f "$cfg_file" "$payload_file" "$response_file"
        exit 1
    fi

    local updated_url updated_engine
    updated_url="$(jq -r '.web.SEARXNG_QUERY_URL // empty' "$response_file" 2>/dev/null || true)"
    updated_engine="$(jq -r '.web.WEB_SEARCH_ENGINE // empty' "$response_file" 2>/dev/null || true)"

    print_success "Retrieval web-search config synced"
    [ -n "$updated_engine" ] && print_info "WEB_SEARCH_ENGINE=$updated_engine"
    [ -n "$updated_url" ] && print_info "SEARXNG_QUERY_URL=$updated_url"

    rm -f "$cfg_file" "$payload_file" "$response_file"
}

main "$@"
