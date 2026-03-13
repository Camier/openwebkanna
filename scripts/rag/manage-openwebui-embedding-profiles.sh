#!/bin/bash

###############################################################################
# OpenWebUI Embedding Profiles Manager
# - Profile registry for embedding models
# - Safe runtime switching
# - Optional model prewarm in OpenWebUI container
# - KB-to-profile bindings and drift diagnostics
###############################################################################

set -e

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
OPENWEBUI_CONTAINER_NAME="${OPENWEBUI_CONTAINER_NAME:-${OPENWEBUI_CONTAINER:-${OPENWEBUI_DOCKER_CONTAINER:-$OPENWEBUI_SERVICE}}}"
POSTGRES_CONTAINER="${POSTGRES_CONTAINER:-openwebui_postgres}"
POSTGRES_DB="${POSTGRES_DB:-openwebui}"
POSTGRES_USER="${POSTGRES_USER:-openwebui}"

API_KEY="${OPENWEBUI_API_KEY:-${API_KEY:-}}"

PROFILES_FILE="${EMBEDDING_PROFILES_FILE:-config/embeddings/profiles.json}"
BINDINGS_FILE="${EMBEDDING_BINDINGS_FILE:-config/embeddings/kb-bindings.json}"

SUBCOMMAND="${1:-}"
KB_ID_FILTER=""

show_help() {
    cat <<'EOF'
Usage: ./scripts/rag/manage-openwebui-embedding-profiles.sh <command> [options]

Commands:
  list
      Show local profiles and current runtime embedding config.

  runtime
      Show current runtime embedding config as JSON.

  lanes
      Show configured KB lanes (kb_key -> kb/profile/model).

  prewarm --profile ID
      Download/cache profile model inside OpenWebUI container.

  apply --profile ID [--prewarm]
      Apply profile to OpenWebUI runtime embedding config.
      Optionally prewarm first.

  bindings
      Show KB bindings from local policy file (profile + optional lane metadata).

  bind-kb --kb-id UUID --profile ID [--kb-name NAME] [--lane KEY] [--replace-lane] [--use-case TEXT]
      Upsert KB binding in config/embeddings/kb-bindings.json.
      `--lane` is the operator-facing lane key.
      `--kb-key` remains supported as a compatibility alias.
      `--replace-lane` transfers an existing lane binding to the target KB.
      `--use-case` is optional legacy metadata.

  use-kb (--kb-id UUID | --lane KEY) [--prewarm] [--profile ID]
      One-shot KB activation:
      - Resolve profile from binding (or explicit --profile)
      - Apply runtime embedding profile
      - Run targeted drift diagnosis

  diagnose [--kb-id UUID]
      Compare:
      - Runtime query model
      - Indexed model per KB (from pgvector metadata)
      - Optional expected model from local bindings

Environment overrides:
  OPENWEBUI_URL
  OPENWEBUI_API_KEY (or API_KEY)
  OPENWEBUI_SIGNIN_EMAIL / OPENWEBUI_SIGNIN_PASSWORD
  OPENWEBUI_CONTAINER
  POSTGRES_CONTAINER / POSTGRES_DB / POSTGRES_USER
  EMBEDDING_PROFILES_FILE
  EMBEDDING_BINDINGS_FILE
EOF
}

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        print_error "Missing dependency: $cmd"
        exit 1
    fi
}

ensure_json_file() {
    local file="$1"
    local fallback="$2"

    if [ -f "$file" ]; then
        return 0
    fi

    mkdir -p "$(dirname "$file")"
    printf "%s\n" "$fallback" >"$file"
}

ensure_files() {
    if [ ! -f "$PROFILES_FILE" ]; then
        print_error "Profiles file not found: $PROFILES_FILE"
        exit 1
    fi
    ensure_json_file "$BINDINGS_FILE" '{"schema_version":1,"bindings":[]}'
}

ensure_api_key() {
    if [ -n "$API_KEY" ]; then
        return 0
    fi

    if ! is_true "$OPENWEBUI_AUTO_AUTH"; then
        print_error "No API key and OPENWEBUI_AUTO_AUTH=false"
        exit 1
    fi

    API_KEY="$(
        resolve_openwebui_api_token \
            "" \
            "$OPENWEBUI_URL" \
            "$OPENWEBUI_SIGNIN_EMAIL" \
            "$OPENWEBUI_SIGNIN_PASSWORD" \
            "30" \
            "$OPENWEBUI_SERVICE" \
            "$OPENWEBUI_CONTAINER_NAME" \
            "$OPENWEBUI_SIGNIN_PATH" || true
    )"

    if [ -z "$API_KEY" ]; then
        print_error "Unable to authenticate to OpenWebUI at ${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"
        exit 1
    fi
}

try_api_key() {
    if [ -n "$API_KEY" ]; then
        return 0
    fi
    if ! is_true "$OPENWEBUI_AUTO_AUTH"; then
        return 1
    fi
    if ensure_api_key 2>/dev/null; then
        return 0
    fi
    return 1
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

profile_exists() {
    local profile_id="$1"
    jq -e --arg id "$profile_id" '.profiles[]? | select(.id == $id)' "$PROFILES_FILE" >/dev/null 2>&1
}

get_profile_json() {
    local profile_id="$1"
    jq -c --arg id "$profile_id" '.profiles[]? | select(.id == $id)' "$PROFILES_FILE"
}

get_profile_model() {
    local profile_id="$1"
    jq -r --arg id "$profile_id" '.profiles[]? | select(.id == $id) | .model // empty' "$PROFILES_FILE"
}

set_active_profile() {
    local profile_id="$1"
    local tmp_file
    tmp_file="$(mktemp)"
    jq --arg id "$profile_id" '.active_profile = $id' "$PROFILES_FILE" >"$tmp_file"
    mv "$tmp_file" "$PROFILES_FILE"
}

print_profile_table() {
    jq -r '
        .profiles[] |
        [
            .id,
            (.model // ""),
            (.batch_size // 0 | tostring),
            (.enable_async // true | tostring),
            (.description // "")
        ] | @tsv
    ' "$PROFILES_FILE" | while IFS=$'\t' read -r id model batch async desc; do
        printf "  %-24s | %-42s | batch=%-2s | async=%-5s | %s\n" "$id" "$model" "$batch" "$async" "$desc"
    done
}

cmd_list() {
    local active_profile active_model runtime_json runtime_model
    active_profile="$(jq -r '.active_profile // empty' "$PROFILES_FILE")"
    active_model=""
    if [ -n "$active_profile" ]; then
        active_model="$(get_profile_model "$active_profile")"
    fi

    print_header "Embedding Profiles Registry"
    print_info "Profiles file: $PROFILES_FILE"
    print_info "Bindings file: $BINDINGS_FILE"
    if [ -n "$active_profile" ]; then
        print_info "Local active profile: $active_profile"
    else
        print_warning "No active_profile set in $PROFILES_FILE"
    fi
    if [ -n "$active_model" ]; then
        print_info "Local active model: $active_model"
    fi
    print_section "Available Profiles"
    print_profile_table

    if try_api_key; then
        runtime_json="$(api_get "/api/v1/retrieval/embedding" 2>/dev/null || true)"
        runtime_model="$(echo "$runtime_json" | jq -r '.RAG_EMBEDDING_MODEL // empty' 2>/dev/null || true)"
        if [ -n "$runtime_model" ]; then
            print_section "Runtime Embedding"
            print_info "OpenWebUI URL: $OPENWEBUI_URL"
            print_info "Runtime model: $runtime_model"
            if [ -n "$active_model" ] && [ "$runtime_model" != "$active_model" ]; then
                print_warning "Runtime model differs from local active profile model"
            fi
        else
            print_warning "Could not read runtime embedding config"
        fi
    else
        print_warning "No auth available; skipping runtime embedding query"
    fi
}

cmd_runtime() {
    local runtime_json runtime_model
    ensure_api_key
    runtime_json="$(api_get "/api/v1/retrieval/embedding" 2>/dev/null || true)"
    runtime_model="$(echo "$runtime_json" | jq -r '.RAG_EMBEDDING_MODEL // empty' 2>/dev/null || true)"
    if [ -z "$runtime_model" ]; then
        print_error "Could not read runtime embedding config"
        exit 1
    fi
    echo "$runtime_json" | jq '{RAG_EMBEDDING_ENGINE,RAG_EMBEDDING_MODEL,RAG_EMBEDDING_BATCH_SIZE,ENABLE_ASYNC_EMBEDDING}'
}

cmd_prewarm() {
    local profile_id="$1"
    local model_id
    model_id="$(get_profile_model "$profile_id")"
    if [ -z "$model_id" ]; then
        print_error "Profile model is empty for profile: $profile_id"
        exit 1
    fi

    local openwebui_container_id=""
    local openwebui_python=""

    init_compose_cmd >/dev/null
    openwebui_container_id="$(resolve_container_id "$OPENWEBUI_SERVICE" "$OPENWEBUI_CONTAINER_NAME" || true)"
    if [ -z "$openwebui_container_id" ]; then
        print_error "OpenWebUI container is not running"
        exit 1
    fi
    if ! openwebui_python="$(resolve_container_python "$openwebui_container_id")"; then
        print_error "No python interpreter found inside OpenWebUI container"
        exit 1
    fi

    print_step "Prewarming model in container '$OPENWEBUI_CONTAINER_NAME': $model_id"
    docker exec -e MODEL_ID="$model_id" "$openwebui_container_id" "$openwebui_python" - <<'PY'
import os
from sentence_transformers import SentenceTransformer

model_id = os.environ["MODEL_ID"]
SentenceTransformer(model_id)
print(f"model_cached={model_id}")
PY
    print_success "Prewarm completed"
}

cmd_apply() {
    local profile_id="$1"
    local do_prewarm="$2"
    local profile_json payload out_file code response

    ensure_api_key
    profile_json="$(get_profile_json "$profile_id")"
    if [ -z "$profile_json" ]; then
        print_error "Unknown profile: $profile_id"
        exit 1
    fi

    if is_true "$do_prewarm"; then
        cmd_prewarm "$profile_id"
    fi

    payload="$(echo "$profile_json" | jq -c '{
        RAG_EMBEDDING_ENGINE: (.engine // ""),
        RAG_EMBEDDING_MODEL: .model,
        RAG_EMBEDDING_BATCH_SIZE: (.batch_size // 4),
        ENABLE_ASYNC_EMBEDDING: (.enable_async // true)
    }')"

    print_step "Applying profile '$profile_id' to OpenWebUI runtime"
    out_file="$(mktemp)"
    code="$(api_post_json_with_status "/api/v1/retrieval/embedding/update" "$payload" "$out_file")"
    response="$(cat "$out_file" 2>/dev/null || true)"
    rm -f "$out_file"

    if [ -z "$code" ] || [ "$code" -lt 200 ] || [ "$code" -ge 300 ]; then
        print_error "Embedding update failed (HTTP ${code:-unknown})"
        [ -n "$response" ] && print_info "Response: $(echo "$response" | head -c 240)"
        exit 1
    fi

    set_active_profile "$profile_id"
    print_success "Applied profile: $profile_id"
    echo "$response" | jq '{RAG_EMBEDDING_ENGINE,RAG_EMBEDDING_MODEL,RAG_EMBEDDING_BATCH_SIZE,ENABLE_ASYNC_EMBEDDING}'
}

cmd_bindings() {
    print_header "KB Embedding Bindings"
    if ! jq -e '.bindings | type == "array"' "$BINDINGS_FILE" >/dev/null 2>&1; then
        print_error "Invalid bindings file format: $BINDINGS_FILE"
        exit 1
    fi

    local total
    total="$(jq -r '.bindings | length' "$BINDINGS_FILE")"
    print_info "Bindings count: $total"
    if [ "$total" = "0" ]; then
        print_warning "No KB bindings configured yet"
        exit 0
    fi

    jq -r '.bindings[] | [(.kb_id // ""),(.kb_name // ""),(.kb_key // ""),(.profile_id // ""),(.use_case // ""),(.updated_at // "")] | @tsv' "$BINDINGS_FILE" |
        while IFS=$'\t' read -r kb_id kb_name kb_key profile_id use_case updated_at; do
            local model_id
            model_id="$(get_profile_model "$profile_id")"
            printf "  kb=%-36s | lane=%-14s | profile=%-24s | model=%-42s | name=%s\n" "$kb_id" "${kb_key:-<none>}" "$profile_id" "${model_id:-<unknown>}" "${kb_name:-<unset>}"
            if [ -n "$use_case" ]; then
                print_info "    legacy_use_case: $use_case"
            fi
            if [ -n "$updated_at" ]; then
                print_info "    updated_at: $updated_at"
            fi
        done
}

cmd_lanes() {
    print_header "KB Lanes"
    if ! jq -e '.bindings | type == "array"' "$BINDINGS_FILE" >/dev/null 2>&1; then
        print_error "Invalid bindings file format: $BINDINGS_FILE"
        exit 1
    fi

    local total keyed
    total="$(jq -r '.bindings | length' "$BINDINGS_FILE")"
    keyed="$(jq -r '[.bindings[]? | select((.kb_key // "") != "")] | length' "$BINDINGS_FILE")"
    print_info "Bindings total: $total"
    print_info "Bindings with lane key: $keyed"
    if [ "$keyed" = "0" ]; then
        print_warning "No lane key set. Add one with: bind-kb --kb-id <UUID> --profile <ID> --lane <KEY>"
        exit 0
    fi

    jq -r '
      [.bindings[]? | select((.kb_key // "") != "")]
      | sort_by(.kb_key)
      | .[]
      | [(.kb_key // ""),(.kb_id // ""),(.profile_id // ""),(.kb_name // ""),(.use_case // "")]
      | @tsv
    ' "$BINDINGS_FILE" | while IFS=$'\t' read -r kb_key kb_id profile_id kb_name use_case; do
        local model_id
        model_id="$(get_profile_model "$profile_id")"
        printf "  lane=%-14s | kb=%-36s | profile=%-24s | model=%-42s\n" "$kb_key" "$kb_id" "$profile_id" "${model_id:-<unknown>}"
        print_info "    name: ${kb_name:-<unset>}"
        if [ -n "$use_case" ]; then
            print_info "    legacy_use_case: $use_case"
        fi
    done
}

cmd_bind_kb() {
    local kb_id="$1"
    local profile_id="$2"
    local kb_name="$3"
    local kb_key="$4"
    local use_case="$5"
    local replace_lane="$6"
    local updated_at tmp_file existing_kb_for_key

    if ! echo "$kb_id" | rg -q '^[0-9a-fA-F-]{36}$'; then
        print_error "Invalid KB ID format: $kb_id"
        exit 1
    fi

    if [ -n "$kb_key" ] && ! echo "$kb_key" | rg -q '^[a-z0-9][a-z0-9_-]{1,31}$'; then
        print_error "Invalid lane key: $kb_key (expected: lowercase letters/digits with _ or -, length 2-32)"
        exit 1
    fi

    if [ -n "$kb_key" ]; then
        existing_kb_for_key="$(lookup_binding_kb_id_by_key "$kb_key")"
        if [ -n "$existing_kb_for_key" ] && [ "$existing_kb_for_key" != "$kb_id" ]; then
            if [ "$replace_lane" = "true" ]; then
                print_warning "Replacing lane binding: $kb_key ($existing_kb_for_key -> $kb_id)"
            else
                print_error "Lane key already used by another KB: $kb_key -> $existing_kb_for_key"
                print_info "Re-run with --replace-lane to transfer the lane to the new KB"
                exit 1
            fi
        fi
    fi

    if ! profile_exists "$profile_id"; then
        print_error "Unknown profile: $profile_id"
        exit 1
    fi

    updated_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    tmp_file="$(mktemp)"

    jq --arg kb_id "$kb_id" \
        --arg kb_name "$kb_name" \
        --arg kb_key "$kb_key" \
        --arg profile_id "$profile_id" \
        --arg use_case "$use_case" \
        --arg updated_at "$updated_at" \
        --arg replace_lane "$replace_lane" \
        '
       .bindings = (
         (.bindings // [])
         | map(
             select(.kb_id != $kb_id)
           )
         | map(
             if ($replace_lane == "true" and (.kb_key // "") == $kb_key and ($kb_key | length) > 0)
             then . + { kb_key: null }
             else .
             end
           )
         + [{
             kb_id: $kb_id,
             kb_name: (if ($kb_name|length) > 0 then $kb_name else null end),
             kb_key: (if ($kb_key|length) > 0 then $kb_key else null end),
             profile_id: $profile_id,
             use_case: (if ($use_case|length) > 0 then $use_case else null end),
             updated_at: $updated_at
           }]
         | sort_by(.kb_id)
       )
       ' "$BINDINGS_FILE" >"$tmp_file"

    mv "$tmp_file" "$BINDINGS_FILE"
    print_success "Binding updated: kb=$kb_id -> profile=$profile_id"
    if [ -n "$kb_key" ]; then
        print_info "Lane: $kb_key"
    fi
}

get_runtime_model() {
    local runtime_json runtime_model
    runtime_json="$(api_get "/api/v1/retrieval/embedding" 2>/dev/null || true)"
    runtime_model="$(echo "$runtime_json" | jq -r '.RAG_EMBEDDING_MODEL // empty' 2>/dev/null || true)"
    printf '%s\n' "$runtime_model"
}

extract_db_models_tsv() {
    docker exec "$POSTGRES_CONTAINER" psql -U "$POSTGRES_USER" -d "$POSTGRES_DB" -At -F $'\t' -c "
        SELECT collection_name,
               split_part(split_part(vmetadata->>'embedding_config', '''model'': ''', 2), '''', 1) AS model,
               COUNT(*) AS chunk_count
        FROM document_chunk
        WHERE vmetadata ? 'embedding_config'
        GROUP BY collection_name, model;
    " 2>/dev/null || true
}

lookup_binding_profile_id() {
    local kb_id="$1"
    jq -r --arg kb "$kb_id" '.bindings[]? | select(.kb_id == $kb) | .profile_id // empty' "$BINDINGS_FILE" | head -n 1
}

lookup_binding_kb_key() {
    local kb_id="$1"
    jq -r --arg kb "$kb_id" '.bindings[]? | select(.kb_id == $kb) | .kb_key // empty' "$BINDINGS_FILE" | head -n 1
}

lookup_binding_kb_id_by_key() {
    local kb_key="$1"
    jq -r --arg key "$kb_key" '.bindings[]? | select((.kb_key // "") == $key) | .kb_id // empty' "$BINDINGS_FILE" | head -n 1
}

lookup_binding_profile_id_by_key() {
    local kb_key="$1"
    jq -r --arg key "$kb_key" '.bindings[]? | select((.kb_key // "") == $key) | .profile_id // empty' "$BINDINGS_FILE" | head -n 1
}

lookup_kb_name() {
    local kb_id="$1"
    local kb_list
    kb_list="$(api_get "/api/v1/knowledge/" 2>/dev/null || true)"
    if [ -z "$kb_list" ]; then
        printf '\n'
        return 0
    fi

    echo "$kb_list" | jq -r --arg kb "$kb_id" '
      if type == "array" then
        .[]
      elif type == "object" and (.items | type == "array") then
        .items[]
      else
        empty
      end
      | select(.id == $kb)
      | (.name // "")
    ' | head -n 1
}

get_kb_index_summary() {
    local kb_id="$1"
    local db_rows
    db_rows="$(extract_db_models_tsv)"

    printf "%s\n" "$db_rows" | awk -F'\t' -v id="$kb_id" '
      $1==id && NF>=3 {
        model = $2
        count = $3 + 0
        total += count
        if (count > max_count) {
          max_count = count
          max_model = model
        }
      }
      END {
        if (total > 0) {
          printf "%s\t%d\n", max_model, total
        }
      }
    '
}

cmd_use_kb() {
    local kb_id="$1"
    local kb_key="$2"
    local do_prewarm="$3"
    local explicit_profile_id="$4"
    local target_kb_id target_kb_name
    local resolved_profile_id binding_profile_id
    local resolved_model runtime_model kb_name
    local db_summary db_model chunk_count

    ensure_api_key

    if [ -n "$kb_id" ] && [ -n "$kb_key" ]; then
        print_error "Use either --kb-id or --lane, not both"
        exit 1
    fi

    if [ -z "$kb_id" ] && [ -z "$kb_key" ]; then
        print_error "Missing target KB. Provide --kb-id or --lane"
        exit 1
    fi

    if [ -n "$kb_key" ]; then
        if ! echo "$kb_key" | rg -q '^[a-z0-9][a-z0-9_-]{1,31}$'; then
            print_error "Invalid lane key: $kb_key"
            exit 1
        fi
        kb_id="$(lookup_binding_kb_id_by_key "$kb_key")"
        if [ -z "$kb_id" ]; then
            print_error "No binding found for lane key: $kb_key"
            print_info "Check available lanes with:"
            print_info "  ./scripts/rag/manage-openwebui-embedding-profiles.sh lanes"
            exit 1
        fi
    fi

    if ! echo "$kb_id" | rg -q '^[0-9a-fA-F-]{36}$'; then
        print_error "Invalid KB ID format: $kb_id"
        exit 1
    fi

    kb_name="$(lookup_kb_name "$kb_id")"
    if [ -z "$kb_name" ]; then
        print_error "KB not found in /api/v1/knowledge/: $kb_id"
        exit 1
    fi
    target_kb_id="$kb_id"
    target_kb_name="$kb_name"

    binding_profile_id="$(lookup_binding_profile_id "$target_kb_id")"
    if [ -n "$kb_key" ]; then
        binding_profile_id="$(lookup_binding_profile_id_by_key "$kb_key")"
    fi

    if [ -n "$explicit_profile_id" ]; then
        profile_exists "$explicit_profile_id" || {
            print_error "Unknown profile: $explicit_profile_id"
            exit 1
        }
        resolved_profile_id="$explicit_profile_id"
        if [ -n "$binding_profile_id" ] && [ "$binding_profile_id" != "$explicit_profile_id" ]; then
            print_warning "Explicit profile overrides binding ($binding_profile_id -> $explicit_profile_id)"
        fi
    else
        resolved_profile_id="$binding_profile_id"
        if [ -z "$resolved_profile_id" ]; then
            print_error "No binding found for KB $target_kb_id"
            print_info "Create one with:"
            print_info "  ./scripts/rag/manage-openwebui-embedding-profiles.sh bind-kb --kb-id $target_kb_id --profile <PROFILE_ID>"
            exit 1
        fi
    fi

    resolved_model="$(get_profile_model "$resolved_profile_id")"
    if [ -z "$resolved_model" ]; then
        print_error "Profile model is empty for profile: $resolved_profile_id"
        exit 1
    fi

    print_header "KB Activation"
    if [ -n "$kb_key" ]; then
        print_info "Lane: $kb_key"
    fi
    print_info "KB ID: $target_kb_id"
    print_info "KB name: $target_kb_name"
    print_info "Resolved profile: $resolved_profile_id"
    print_info "Resolved model: $resolved_model"

    cmd_apply "$resolved_profile_id" "$do_prewarm"

    KB_ID_FILTER="$target_kb_id"
    cmd_diagnose

    runtime_model="$(get_runtime_model)"
    db_summary="$(get_kb_index_summary "$target_kb_id")"
    db_model="$(echo "$db_summary" | awk -F'\t' '{print $1}' | head -n 1)"
    chunk_count="$(echo "$db_summary" | awk -F'\t' '{print $2}' | head -n 1)"

    if [ "$runtime_model" != "$resolved_model" ]; then
        print_error "Runtime model does not match resolved model after apply"
        print_info "Runtime: ${runtime_model:-<empty>}"
        print_info "Resolved: $resolved_model"
        exit 1
    fi

    if [ -n "$db_model" ] && [ "$db_model" != "$runtime_model" ]; then
        print_error "KB indexed model mismatches runtime model"
        print_info "Indexed: $db_model"
        print_info "Runtime: $runtime_model"
        exit 1
    fi

    if [ -z "$chunk_count" ] || [ "$chunk_count" = "0" ]; then
        print_warning "KB has no indexed chunks yet; runtime is aligned but retrieval may return no results"
    fi

    print_success "KB ready: $target_kb_name ($target_kb_id) on profile $resolved_profile_id"
}

cmd_diagnose() {
    ensure_api_key

    local runtime_model kb_list db_rows
    local kb_map_file db_file agg_file filter_file
    runtime_model="$(get_runtime_model)"

    if [ -z "$runtime_model" ]; then
        print_error "Could not read runtime embedding model"
        exit 1
    fi

    kb_list="$(api_get "/api/v1/knowledge/" 2>/dev/null || true)"
    if [ -z "$kb_list" ]; then
        print_error "Could not read knowledge base list"
        exit 1
    fi

    db_rows="$(extract_db_models_tsv)"
    kb_map_file="$(mktemp)"
    db_file="$(mktemp)"
    agg_file="$(mktemp)"
    filter_file="$(mktemp)"

    trap 'rm -f "$kb_map_file" "$db_file" "$agg_file" "$filter_file"' RETURN

    echo "$kb_list" | jq -r '
      if type == "array" then
        .[]
      elif type == "object" and (.items | type == "array") then
        .items[]
      else
        empty
      end
      | [.id, (.name // "")] | @tsv
    ' >"$kb_map_file"
    printf "%s\n" "$db_rows" >"$db_file"

    awk -F'\t' '
      NF >= 3 {
        key = $1
        model = $2
        count = $3 + 0
        total[key] += count
        if (count > max_count[key]) {
          max_count[key] = count
          max_model[key] = model
        }
      }
      END {
        for (k in total) {
          printf "%s\t%s\t%d\n", k, max_model[k], total[k]
        }
      }
    ' "$db_file" | sort >"$agg_file"

    if [ -n "$KB_ID_FILTER" ]; then
        rg -n "^${KB_ID_FILTER}\t" "$kb_map_file" >/dev/null 2>&1 || {
            print_error "KB ID not found in /api/v1/knowledge/: $KB_ID_FILTER"
            exit 1
        }
        rg "^${KB_ID_FILTER}\t" "$kb_map_file" >"$filter_file"
    else
        cp "$kb_map_file" "$filter_file"
    fi

    print_header "Embedding Drift Diagnosis"
    print_info "Runtime model: $runtime_model"
    print_info "Knowledge endpoint: ${OPENWEBUI_URL%/}/api/v1/knowledge/"
    print_section "KB Compatibility"

    while IFS=$'\t' read -r kb_id_row kb_name_row; do
        [ -n "$kb_id_row" ] || continue

        local db_model chunk_count
        local binding_profile_id binding_model binding_kb_key
        local status runtime_status binding_status

        db_model="$(awk -F'\t' -v id="$kb_id_row" '$1==id {print $2}' "$agg_file" | head -n 1)"
        chunk_count="$(awk -F'\t' -v id="$kb_id_row" '$1==id {print $3}' "$agg_file" | head -n 1)"
        binding_profile_id="$(lookup_binding_profile_id "$kb_id_row")"
        binding_kb_key="$(lookup_binding_kb_key "$kb_id_row")"
        binding_model=""
        if [ -n "$binding_profile_id" ]; then
            binding_model="$(get_profile_model "$binding_profile_id")"
        fi

        if [ -z "$chunk_count" ]; then
            status="NO_CHUNKS"
            runtime_status="-"
        elif [ "$db_model" = "$runtime_model" ]; then
            runtime_status="MATCH"
            status="OK"
        else
            runtime_status="MISMATCH"
            status="DRIFT"
        fi

        binding_status="-"
        if [ -n "$binding_model" ] && [ -n "$db_model" ]; then
            if [ "$binding_model" = "$db_model" ]; then
                binding_status="MATCH"
            else
                binding_status="MISMATCH"
            fi
        fi

        printf "  %-8s | kb=%s | chunks=%s\n" "$status" "$kb_id_row" "${chunk_count:-0}"
        print_info "    name: ${kb_name_row:-<unnamed>}"
        if [ -n "$binding_kb_key" ]; then
            print_info "    lane: $binding_kb_key"
        fi
        print_info "    indexed_model: ${db_model:-<none>}"
        print_info "    runtime_model: $runtime_model ($runtime_status)"
        if [ -n "$binding_profile_id" ]; then
            print_info "    binding: ${binding_profile_id} -> ${binding_model:-<unknown>} ($binding_status)"
        else
            print_info "    binding: <none>"
        fi
    done <"$filter_file"
}

main() {
    require_cmd jq
    require_cmd curl
    require_cmd rg
    ensure_files

    if [ -z "$SUBCOMMAND" ]; then
        show_help
        exit 1
    fi

    shift 1
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --url)
                OPENWEBUI_URL="$2"
                shift 2
                ;;
            *)
                break
                ;;
        esac
    done

    case "$SUBCOMMAND" in
        list)
            cmd_list
            ;;
        runtime)
            cmd_runtime
            ;;
        lanes)
            cmd_lanes
            ;;
        prewarm)
            local profile_id=""
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    --profile)
                        profile_id="$2"
                        shift 2
                        ;;
                    *)
                        print_error "Unknown argument: $1"
                        exit 1
                        ;;
                esac
            done
            [ -n "$profile_id" ] || {
                print_error "Missing required argument: --profile"
                exit 1
            }
            profile_exists "$profile_id" || {
                print_error "Unknown profile: $profile_id"
                exit 1
            }
            cmd_prewarm "$profile_id"
            ;;
        apply)
            local profile_id=""
            local do_prewarm="false"
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    --profile)
                        profile_id="$2"
                        shift 2
                        ;;
                    --prewarm)
                        do_prewarm="true"
                        shift 1
                        ;;
                    *)
                        print_error "Unknown argument: $1"
                        exit 1
                        ;;
                esac
            done
            [ -n "$profile_id" ] || {
                print_error "Missing required argument: --profile"
                exit 1
            }
            profile_exists "$profile_id" || {
                print_error "Unknown profile: $profile_id"
                exit 1
            }
            cmd_apply "$profile_id" "$do_prewarm"
            ;;
        bindings)
            cmd_bindings
            ;;
        bind-kb)
            local kb_id=""
            local profile_id=""
            local kb_name=""
            local kb_key=""
            local use_case=""
            local replace_lane="false"
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    --kb-id)
                        kb_id="$2"
                        shift 2
                        ;;
                    --profile)
                        profile_id="$2"
                        shift 2
                        ;;
                    --kb-name)
                        kb_name="$2"
                        shift 2
                        ;;
                    --lane)
                        kb_key="$2"
                        shift 2
                        ;;
                    --kb-key)
                        kb_key="$2"
                        shift 2
                        ;;
                    --use-case)
                        use_case="$2"
                        shift 2
                        ;;
                    --replace-lane)
                        replace_lane="true"
                        shift 1
                        ;;
                    *)
                        print_error "Unknown argument: $1"
                        exit 1
                        ;;
                esac
            done
            [ -n "$kb_id" ] || {
                print_error "Missing required argument: --kb-id"
                exit 1
            }
            [ -n "$profile_id" ] || {
                print_error "Missing required argument: --profile"
                exit 1
            }
            cmd_bind_kb "$kb_id" "$profile_id" "$kb_name" "$kb_key" "$use_case" "$replace_lane"
            ;;
        use-kb)
            local kb_id=""
            local kb_key=""
            local do_prewarm="false"
            local profile_id=""
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    --kb-id)
                        kb_id="$2"
                        shift 2
                        ;;
                    --lane)
                        kb_key="$2"
                        shift 2
                        ;;
                    --kb-key)
                        kb_key="$2"
                        shift 2
                        ;;
                    --prewarm)
                        do_prewarm="true"
                        shift 1
                        ;;
                    --profile)
                        profile_id="$2"
                        shift 2
                        ;;
                    *)
                        print_error "Unknown argument: $1"
                        exit 1
                        ;;
                esac
            done
            if [ -z "$kb_id" ] && [ -z "$kb_key" ]; then
                print_error "Missing target KB: provide --kb-id or --lane"
                exit 1
            fi
            cmd_use_kb "$kb_id" "$kb_key" "$do_prewarm" "$profile_id"
            ;;
        diagnose)
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    --kb-id)
                        KB_ID_FILTER="$2"
                        shift 2
                        ;;
                    *)
                        print_error "Unknown argument: $1"
                        exit 1
                        ;;
                esac
            done
            cmd_diagnose
            ;;
        help | -h | --help)
            show_help
            ;;
        *)
            print_error "Unknown command: $SUBCOMMAND"
            show_help
            exit 1
            ;;
    esac
}

main "$@"
