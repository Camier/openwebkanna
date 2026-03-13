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

default_openwebui_smoke_model_candidates() {
    printf "%s" "gpt-oss:20b,kimi-k2-thinking,gemini-3-flash-preview,ministral-3:3b,gemma3:4b"
}

mint_local_admin_jwt() {
    local container_name="${OPENWEBUI_CONTAINER:-openwebui}"
    if ! command -v docker >/dev/null 2>&1 || ! command -v python3 >/dev/null 2>&1; then
        return 1
    fi
    if [ -z "${WEBUI_SECRET_KEY:-}" ]; then
        return 1
    fi
    if ! docker ps --format '{{.Names}}' | grep -qx "$container_name"; then
        return 1
    fi

    local admin_id=""
    admin_id="$(docker exec "$container_name" python3 -c "import sqlite3; con=sqlite3.connect('/app/backend/data/webui.db'); cur=con.cursor(); cur.execute(\"SELECT id FROM user WHERE role='admin' ORDER BY created_at ASC LIMIT 1\"); row=cur.fetchone(); print(row[0] if row else '')" 2>/dev/null | tr -d '\r' | head -n 1 || true)"
    if [ -z "$admin_id" ]; then
        return 1
    fi

    API_KEY="$(
        WEBUI_SECRET_KEY="$WEBUI_SECRET_KEY" OPENWEBUI_ADMIN_ID="$admin_id" python3 - <<'PY'
import base64
import hashlib
import hmac
import json
import os
import time
import uuid

secret = os.environ.get("WEBUI_SECRET_KEY", "")
user_id = os.environ.get("OPENWEBUI_ADMIN_ID", "")
if not secret or not user_id:
    raise SystemExit(1)

def b64url(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).decode("utf-8").rstrip("=")

header = {"alg": "HS256", "typ": "JWT"}
payload = {
    "id": user_id,
    "exp": int(time.time()) + (60 * 60 * 24 * 30),
    "jti": str(uuid.uuid4()),
}

header_b64 = b64url(json.dumps(header, separators=(",", ":"), sort_keys=True).encode("utf-8"))
payload_b64 = b64url(json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8"))
signing_input = f"{header_b64}.{payload_b64}".encode("utf-8")
sig = hmac.new(secret.encode("utf-8"), signing_input, hashlib.sha256).digest()
print(f"{header_b64}.{payload_b64}.{b64url(sig)}")
PY
    )"

    [ -n "$API_KEY" ]
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

    if [ -n "$API_KEY" ]; then
        return 0
    fi

    if mint_local_admin_jwt; then
        print_info "OpenWebUI bearer token acquired via local JWT mint"
        return 0
    fi
}

json_response_has_any_id() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e '
            (.id // .file_id // .knowledge_id // empty)
            | select(type == "string" and length > 0)
        ' >/dev/null 2>&1
        return $?
    fi

    echo "$response" | grep -Eq '"(id|file_id|knowledge_id)"[[:space:]]*:[[:space:]]*"[^"]+"'
}

extract_json_primary_id() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -r '
            .id // .file_id // .knowledge_id // empty
        ' 2>/dev/null
        return 0
    fi

    echo "$response" | grep -Eo '"(id|file_id|knowledge_id)"[[:space:]]*:[[:space:]]*"[^"]+"' | head -n1 | cut -d'"' -f4
}

json_response_has_files_field() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e '
            .files
            | select(type == "array" or type == "object")
        ' >/dev/null 2>&1
        return $?
    fi

    echo "$response" | grep -q '"files"'
}

json_response_has_chat_choices() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e '
            .choices and (.choices | type == "array") and (.choices | length > 0)
        ' >/dev/null 2>&1
        return $?
    fi

    echo "$response" | grep -q '"choices"'
}

extract_chat_response_text() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -r '
            .choices[0].message.content
            // .choices[0].delta.content
            // empty
        ' 2>/dev/null
        return 0
    fi

    printf ""
}

json_response_has_results_array() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e '
            .results and (.results | type == "array")
        ' >/dev/null 2>&1
        return $?
    fi

    echo "$response" | grep -q '"results"'
}

json_response_has_retrieval_hits() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e '
            (
                (.results? | if type == "array" then length else 0 end) > 0
            ) or (
                (.documents? | [.. | strings] | length) > 0
            )
        ' >/dev/null 2>&1
        return $?
    fi

    echo "$response" | grep -q '"results"\|"documents"'
}

extract_retrieval_text() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -r '
            [
                (
                    .results[]? |
                        if type == "object" then
                            (.content // .text // .document // .body // empty)
                        elif type == "string" then
                            .
                        else
                            empty
                        end
                ),
                (
                    .documents[]? |
                        if type == "array" then
                            .[]
                        elif type == "string" then
                            .
                        else
                            empty
                        end
                )
            ]
            | .[]
            | select(type == "string" and length > 0)
        ' 2>/dev/null
        return 0
    fi

    printf "%s" "$response"
}

retrieval_hit_count() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -r '
            if (.results? | type) == "array" then
                (.results | length)
            elif (.documents? | type) == "array" then
                ([.documents[]? | if type == "array" then .[] else . end | select(type == "string" and length > 0)] | length)
            else
                0
            end
        ' 2>/dev/null
        return 0
    fi

    printf "0"
}

response_contains_text() {
    local response="$1"
    local needle="$2"

    if [ -z "$response" ] || [ -z "$needle" ]; then
        return 1
    fi

    printf "%s" "$response" | grep -F -- "$needle" >/dev/null 2>&1
}

extract_openwebui_file_status() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -r '.status // empty' 2>/dev/null
        return 0
    fi

    echo "$response" | grep -Eo '"status"[[:space:]]*:[[:space:]]*"[^"]+"' | head -n1 | cut -d'"' -f4
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

models_response_has_id() {
    local response="$1"
    local model_id="$2"

    if [ -z "$response" ] || [ -z "$model_id" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -e --arg model "$model_id" 'any(.data[]?; .id == $model)' >/dev/null 2>&1
        return $?
    fi

    echo "$response" | grep -F "\"id\":\"$model_id\"" >/dev/null 2>&1
}

first_model_id_from_response() {
    local response="$1"

    if [ -z "$response" ]; then
        return 1
    fi

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -r '
            .data[]?.id
            | select(type == "string" and length > 0)
            | select(. != "arena-model")
        ' 2>/dev/null | head -n 1
        return 0
    fi

    echo "$response" | grep -Eo '"id"[[:space:]]*:[[:space:]]*"[^"]+"' | head -n1 | cut -d'"' -f4
}

choose_preferred_model_from_response() {
    local response="$1"
    local exclude_model="${2:-}"
    local candidates candidate

    if [ -z "$response" ]; then
        return 1
    fi

    candidates="$(printf "%s" "${OPENWEBUI_SMOKE_MODEL_CANDIDATES:-}" | tr ',' '\n')"
    while IFS= read -r candidate; do
        candidate="$(printf "%s" "$candidate" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        if [ -z "$candidate" ] || [ "$candidate" = "$exclude_model" ]; then
            continue
        fi
        if models_response_has_id "$response" "$candidate"; then
            printf "%s" "$candidate"
            return 0
        fi
    done <<EOF
$candidates
EOF

    if command -v jq >/dev/null 2>&1; then
        echo "$response" | jq -r --arg current "$exclude_model" '
            .data[]?.id
            | select(type == "string" and length > 0)
            | select(. != $current)
            | select((contains("*") | not))
            | select((startswith("ollama-cloud/") | not))
            | select(. != "arena-model")
        ' 2>/dev/null | head -n 1
        return 0
    fi

    first_model_id_from_response "$response"
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
        status_response="$(curl -sS -m 20 -X GET \
            "$OPENWEBUI_URL/api/v1/files/$file_id/process/status" \
            ${API_KEY:+-H "Authorization: Bearer $API_KEY"} 2>/dev/null || true)"
        status_value="$(extract_openwebui_file_status "$status_response" || true)"

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
