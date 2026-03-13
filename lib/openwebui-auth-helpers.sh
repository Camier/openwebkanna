#!/usr/bin/env bash

###############################################################################
# OpenWebUI auth helper library
# Provides reusable helpers for acquiring admin API tokens and resolving the
# local admin user ID from the running OpenWebUI container.
###############################################################################

openwebui_write_signin_payload() {
    local output_file="$1"
    local email="$2"
    local password="$3"

    if command_exists jq; then
        jq -n \
            --arg email "$email" \
            --arg password "$password" \
            '{email: $email, password: $password}' >"$output_file"
        return 0
    fi

    cat >"$output_file" <<EOF
{"email":"$email","password":"$password"}
EOF
}

openwebui_signin_token() {
    local base_url="$1"
    local email="$2"
    local password="$3"
    local timeout="${4:-45}"
    local signin_path="${5:-/api/v1/auths/signin}"
    local payload_file=""
    local output_file=""
    local status_code=""
    local token=""

    payload_file="$(mktemp)"
    output_file="$(mktemp)"
    openwebui_write_signin_payload "$payload_file" "$email" "$password"

    status_code="$(curl -sS -m "$timeout" \
        -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" \
        "${base_url%/}${signin_path}" -d @"$payload_file" 2>/dev/null || true)"

    if [ "$status_code" = "200" ]; then
        token="$(jq -r '.token // empty' "$output_file" 2>/dev/null || true)"
    fi

    rm -f "$payload_file" "$output_file"

    [ -n "$token" ] || return 1
    printf "%s" "$token"
}

resolve_openwebui_admin_id() {
    local service="${1:-openwebui}"
    local container_name="${2:-$service}"
    local container_id=""
    local container_python=""
    local admin_id=""

    if ! command_exists docker; then
        return 1
    fi

    init_compose_cmd >/dev/null
    container_id="$(resolve_container_id "$service" "$container_name" || true)"
    if [ -z "$container_id" ]; then
        return 1
    fi

    if ! container_python="$(resolve_container_python "$container_id")"; then
        return 1
    fi

    admin_id="$(docker exec "$container_id" "$container_python" -c "import sqlite3; con=sqlite3.connect('/app/backend/data/webui.db'); cur=con.cursor(); cur.execute(\"SELECT id FROM user WHERE role='admin' ORDER BY created_at ASC LIMIT 1\"); row=cur.fetchone(); print(row[0] if row else '')" 2>/dev/null | tr -d '\r' | head -n 1 || true)"
    [ -n "$admin_id" ] || return 1
    printf "%s" "$admin_id"
}

mint_openwebui_admin_jwt() {
    local service="${1:-openwebui}"
    local container_name="${2:-$service}"
    local host_python=""
    local admin_id=""
    local token=""

    if [ -z "${WEBUI_SECRET_KEY:-}" ]; then
        return 1
    fi

    if ! host_python="$(resolve_host_python)"; then
        return 1
    fi

    if ! admin_id="$(resolve_openwebui_admin_id "$service" "$container_name")"; then
        return 1
    fi

    token="$(
        WEBUI_SECRET_KEY="$WEBUI_SECRET_KEY" OPENWEBUI_ADMIN_ID="$admin_id" "$host_python" - <<'PY'
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

    [ -n "$token" ] || return 1
    printf "%s" "$token"
}

resolve_openwebui_api_token() {
    local token_input="${1:-}"
    local base_url="$2"
    local email="$3"
    local password="$4"
    local timeout="${5:-45}"
    local service="${6:-openwebui}"
    local container_name="${7:-$service}"
    local signin_path="${8:-/api/v1/auths/signin}"
    local token=""

    if [ -n "$token_input" ]; then
        printf "%s" "$token_input"
        return 0
    fi

    token="$(openwebui_signin_token "$base_url" "$email" "$password" "$timeout" "$signin_path" || true)"
    if [ -n "$token" ]; then
        printf "%s" "$token"
        return 0
    fi

    token="$(mint_openwebui_admin_jwt "$service" "$container_name" || true)"
    [ -n "$token" ] || return 1
    printf "%s" "$token"
}
