#!/bin/bash

###############################################################################
# http-helpers.sh - HTTP request and model extraction utilities
#
# Provides functions for HTTP status checks, service waiting, and
# JSON model payload parsing (for OpenAI-compatible /v1/models responses).
#
# Functions:
#   wait_for_service()      - Wait for HTTP endpoint to become available
#   json_models_non_empty() - Check if /v1/models response has models
#   extract_model_ids()     - Extract model IDs from /v1/models response
#   assert_expected_aliases() - Verify expected aliases exist in models
#
# Usage:
#   source /path/to/lib/http-helpers.sh
#   wait_for_service "http://localhost:8080/health" "MyService" 30
###############################################################################

###############################################################################
# wait_for_service
# Wait for an HTTP endpoint to become available.
#
# Arguments:
#   $1 - URL to check
#   $2 - Service name (for display)
#   $3 - Maximum attempts (optional, default: 60)
#
# Returns:
#   0 if service became ready, 1 if timeout
#
# Example:
#   wait_for_service "http://localhost:8080/health" "API Server" 30
###############################################################################
wait_for_service() {
    local url="$1"
    local service_name="$2"
    local max_attempts="${3:-60}"
    local attempt=1

    print_step "Waiting for $service_name to be ready..."

    while [ $attempt -le $max_attempts ]; do
        if curl -s -f "$url" &>/dev/null; then
            print_success "$service_name is ready!"
            return 0
        fi
        echo -n "."
        sleep 2
        attempt=$((attempt + 1))
    done

    echo
    print_error "$service_name failed to start within expected time."
    return 1
}

###############################################################################
# json_models_non_empty
# Check if a JSON file contains a non-empty OpenAI-compatible models array.
#
# Supports both jq and python3 for parsing, with grep fallback.
#
# Arguments:
#   $1 - Path to JSON file containing /v1/models response
#
# Returns:
#   0 if data array exists and has at least one element, 1 otherwise
#
# Example:
#   curl -s "$MODELS_URL" -o /tmp/models.json
#   if json_models_non_empty /tmp/models.json; then ...
###############################################################################
json_models_non_empty() {
    local payload_file="$1"

    if command -v jq >/dev/null 2>&1; then
        jq -e '.data and (.data | type == "array") and (.data | length > 0)' "$payload_file" >/dev/null 2>&1
        return $?
    fi

    if command -v python3 >/dev/null 2>&1; then
        python3 - "$payload_file" <<'PY'
import json
import sys

path = sys.argv[1]
try:
    with open(path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
except Exception:
    raise SystemExit(1)

data = payload.get("data")
if isinstance(data, list) and len(data) > 0:
    raise SystemExit(0)
raise SystemExit(1)
PY
        return $?
    fi

    # Fallback: grep-based check (less reliable)
    if ! grep -q '"data"[[:space:]]*:[[:space:]]*\[' "$payload_file"; then
        return 1
    fi
    grep -q '"id"[[:space:]]*:[[:space:]]*"[^"]\+"' "$payload_file"
}

###############################################################################
# extract_model_ids
# Extract model IDs from a JSON file containing /v1/models response.
#
# Arguments:
#   $1 - Path to JSON file
#
# Output:
#   Model IDs, one per line (stdout)
#
# Example:
#   ids=$(extract_model_ids /tmp/models.json)
###############################################################################
extract_model_ids() {
    local payload_file="$1"

    if command -v jq >/dev/null 2>&1; then
        jq -r '.data[]?.id // empty' "$payload_file" 2>/dev/null || true
        return 0
    fi

    # Fallback: grep-based extraction
    grep -o '"id"[[:space:]]*:[[:space:]]*"[^"]*"' "$payload_file" | cut -d'"' -f4
}

###############################################################################
# assert_expected_aliases
# Verify that all expected model aliases exist in a models response.
#
# Arguments:
#   $1 - Path to JSON file containing /v1/models response
#
# Required globals (set before calling):
#   Uses resolved_expected_aliases() from cliproxyapi-helpers.sh
#
# Returns:
#   0 if all expected aliases are present, 1 if any are missing
#
# Example:
#   CLIPROXYAPI_EXPECT_ALIASES="model-a model-b"
#   if assert_expected_aliases /tmp/models.json; then ...
###############################################################################
assert_expected_aliases() {
    local payload_file="$1"
    local alias=""
    local ids=""
    local expected_aliases=""

    # Get expected aliases (from cliproxyapi-helpers.sh if available)
    if command -v resolved_expected_aliases >/dev/null 2>&1; then
        expected_aliases="$(resolved_expected_aliases | sed '/^$/d' || true)"
    elif [ -n "${CLIPROXYAPI_EXPECT_ALIASES:-}" ]; then
        expected_aliases="$CLIPROXYAPI_EXPECT_ALIASES"
    else
        return 0
    fi

    if [ -z "$expected_aliases" ]; then
        return 0
    fi

    ids="$(extract_model_ids "$payload_file" | sed '/^$/d' || true)"
    if [ -z "$ids" ]; then
        return 1
    fi

    for alias in $expected_aliases; do
        if ! printf "%s\n" "$ids" | grep -Fxq "$alias"; then
            return 1
        fi
    done

    return 0
}

###############################################################################
# http_status
# Get HTTP status code for a URL.
#
# Arguments:
#   $1 - URL to check
#   $2 - Timeout in seconds (optional, default: 10)
#
# Output:
#   HTTP status code (stdout), or "000" on failure
#
# Example:
#   status=$(http_status "http://localhost:8080/health")
#   [ "$status" = "200" ] && echo "OK"
###############################################################################
http_status() {
    local url="$1"
    local timeout="${2:-10}"

    curl -sS -m "$timeout" -o /dev/null -w "%{http_code}" "$url" 2>/dev/null || echo "000"
}

###############################################################################
# http_status_with_body
# Get HTTP status code and save response body to a file.
#
# Arguments:
#   $1 - URL to check
#   $2 - Output file for response body
#   $3 - Timeout in seconds (optional, default: 15)
#   $4 - Optional Authorization header value (e.g., "Bearer token")
#
# Output:
#   HTTP status code (stdout), or "000" on failure
#   Response body written to output file
#
# Example:
#   status=$(http_status_with_body "$URL" /tmp/response.json 10 "Bearer $TOKEN")
###############################################################################
http_status_with_body() {
    local url="$1"
    local output_file="$2"
    local timeout="${3:-15}"
    local auth_header="${4:-}"
    local code=""

    if [ -n "$auth_header" ]; then
        code="$(curl -sS -m "$timeout" -H "Authorization: $auth_header" -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null)" || code="000"
    else
        code="$(curl -sS -m "$timeout" -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null)" || code="000"
    fi

    echo "$code"
}
