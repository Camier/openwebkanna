#!/bin/bash

default_openwebui_smoke_model_candidates() {
    printf "%s" "gpt-oss:20b,kimi-k2-thinking,gemini-3-flash-preview,ministral-3:3b,gemma3:4b"
}

ensure_openwebui_api_key() {
    if [ -n "$API_KEY" ]; then
        return 0
    fi

    if ! is_true "$OPENWEBUI_AUTO_AUTH"; then
        return 0
    fi

    local service="${OPENWEBUI_SERVICE:-openwebui}"
    local container_name="${OPENWEBUI_CONTAINER:-${OPENWEBUI_CONTAINER_NAME:-$service}}"
    API_KEY="$(
        resolve_openwebui_api_token \
            "" \
            "$OPENWEBUI_URL" \
            "$OPENWEBUI_SIGNIN_EMAIL" \
            "$OPENWEBUI_SIGNIN_PASSWORD" \
            "45" \
            "$service" \
            "$container_name" \
            "$OPENWEBUI_SIGNIN_PATH" || true
    )"

    if [ -n "$API_KEY" ]; then
        print_info "OpenWebUI bearer token acquired"
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
