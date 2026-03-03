#!/bin/bash

openwebui_signin() {
    local output_file="$1"
    local payload_file="$2"
    local signin_url="${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"

    if command -v jq >/dev/null 2>&1; then
        jq -n \
            --arg email "$OPENWEBUI_SIGNIN_EMAIL" \
            --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
            '{email: $email, password: $password}' >"$payload_file"
    else
        cat >"$payload_file" <<EOF
{"email":"$OPENWEBUI_SIGNIN_EMAIL","password":"$OPENWEBUI_SIGNIN_PASSWORD"}
EOF
    fi

    curl -sS -m 45 -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" "$signin_url" -d @"$payload_file" 2>/dev/null || true
}

ensure_openwebui_api_key() {
    if [ -n "$API_KEY" ]; then
        return 0
    fi

    if ! is_true "$OPENWEBUI_AUTO_AUTH"; then
        return 0
    fi

    local tmp_file payload_file signin_code token
    tmp_file="$(mktemp)"
    payload_file="$(mktemp)"

    signin_code="$(openwebui_signin "$tmp_file" "$payload_file")"
    if [ "$signin_code" = "200" ]; then
        token="$(jq -r '.token // empty' "$tmp_file" 2>/dev/null || true)"
        if [ -n "$token" ]; then
            API_KEY="$token"
            print_info "OpenWebUI bearer token acquired via signin"
        fi
    fi

    rm -f "$tmp_file" "$payload_file"
}

has_non_empty_models_array() {
    local response="$1"
    local model_count

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e '.data and (.data | type == "array") and (.data | length > 0)' >/dev/null 2>&1
        return $?
    fi

    model_count="$(echo "$response" | grep -Eo '"id"[[:space:]]*:[[:space:]]*"[^"]+"' | wc -l | tr -d '[:space:]')"
    if [[ $model_count =~ ^[0-9]+$ ]] && [ "$model_count" -gt 0 ]; then
        return 0
    fi
    return 1
}

wait_for_openwebui_file_processing() {
    local file_id="$1"
    local max_attempts="${2:-20}"
    local sleep_seconds="${3:-1}"
    local attempts=0
    local status_response=""
    local status_value=""

    if [ -z "$file_id" ]; then
        return 2
    fi

    while [ "$attempts" -lt "$max_attempts" ]; do
        status_response="$(curl -s -X GET \
            "$OPENWEBUI_URL/api/v1/files/$file_id/process/status" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>&1)"
        status_value="$(echo "$status_response" | jq -r '.status // empty' 2>/dev/null || true)"

        if [ "$status_value" = "completed" ]; then
            return 0
        fi

        if [ "$status_value" = "failed" ]; then
            printf "%s" "$status_response"
            return 1
        fi

        attempts=$((attempts + 1))
        sleep "$sleep_seconds"
    done

    printf "%s" "$status_response"
    return 3
}
