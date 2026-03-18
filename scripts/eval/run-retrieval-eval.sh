#!/bin/bash

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_URL_DEFAULT="http://localhost:${WEBUI_PORT:-3000}"
OPENWEBUI_URL="${OPENWEBUI_URL:-$OPENWEBUI_URL_DEFAULT}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
OPENWEBUI_SERVICE="${OPENWEBUI_SERVICE:-${OPENWEBUI_DOCKER_SERVICE:-openwebui}}"
OPENWEBUI_CONTAINER_NAME="${OPENWEBUI_CONTAINER_NAME:-${OPENWEBUI_DOCKER_CONTAINER:-$OPENWEBUI_SERVICE}}"

API_TOKEN_INPUT="${OPENWEBUI_TOKEN:-${OPENWEBUI_API_KEY:-}}"
API_TOKEN="$API_TOKEN_INPUT"

if [ -z "$API_TOKEN" ] && [ "$OPENWEBUI_AUTO_AUTH" != "false" ]; then
    API_TOKEN="$(
        resolve_openwebui_api_token \
            "$API_TOKEN_INPUT" \
            "$OPENWEBUI_URL" \
            "$OPENWEBUI_SIGNIN_EMAIL" \
            "$OPENWEBUI_SIGNIN_PASSWORD" \
            "30" \
            "$OPENWEBUI_SERVICE" \
            "$OPENWEBUI_CONTAINER_NAME" \
            "$OPENWEBUI_SIGNIN_PATH" || true
    )"
fi

if [ -z "$API_TOKEN" ]; then
    print_error "Unable to acquire OpenWebUI token for retrieval evaluation"
    exit 1
fi

export OPENWEBUI_URL
export OPENWEBUI_EVAL_TOKEN="$API_TOKEN"

exec python3 "${SELF_DIR}/run_retrieval_eval.py" "$@"
