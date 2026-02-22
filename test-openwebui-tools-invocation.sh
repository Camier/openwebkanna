#!/bin/bash
###############################################################################
# OpenWebUI Tool Invocation Integration Test
#
# Goal:
# - Prove that OpenWebUI tools *actually execute* when triggered through
#   POST /api/chat/completions (vs. only validating CRUD/specs/valves endpoints).
#
# Strategy:
# - Use tool_ids to restrict available tools during chat, then force a call.
# - Avoid external internet by using a local ephemeral HTTP server that returns a
#   random nonce not present in prompts. The nonce must appear in the response,
#   proving real tool execution.
#
# Tools covered:
# - web_scraper: scrape_url against local HTML endpoint (nonce in page body)
# - llm_council: council against local OpenAI-compatible endpoint via valves override
# - run_code: reads nonce from a file created inside the Jupyter container
#
# Safety:
# - Never prints secrets (tokens / API keys).
# - Restores llm_council valves after test.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source all library modules
source "${SCRIPT_DIR}/lib/init.sh"
cd "$SCRIPT_DIR"

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
        [[ $line == \#* ]] && continue
        [[ $line != *=* ]] && continue

        key="${line%%=*}"
        value="${line#*=}"
        key="$(printf "%s" "$key" | tr -d '[:space:]')"

        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        if [ -n "${!key+x}" ]; then
            continue
        fi

        if [[ $value == \"*\" ]] && [[ $value == *\" ]]; then
            value="${value:1:${#value}-2}"
        elif [[ $value == \'*\' ]] && [[ $value == *\' ]]; then
            value="${value:1:${#value}-2}"
        fi

        printf -v "$key" "%s" "$value"
        export "${key?}"
    done <"$env_file"
}

# Color definitions (standard repository set)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

# Configuration
OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENWEBUI_HEALTH_PATH="${OPENWEBUI_HEALTH_PATH:-/health}"
OPENWEBUI_CHAT_PATH="${OPENWEBUI_CHAT_PATH:-/api/chat/completions}"
OPENWEBUI_MODELS_PATH="${OPENWEBUI_MODELS_PATH:-/api/models}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_CHATS_NEW_PATH="${OPENWEBUI_CHATS_NEW_PATH:-/api/v1/chats/new}"
OPENWEBUI_CHAT_BY_ID_PREFIX="${OPENWEBUI_CHAT_BY_ID_PREFIX:-/api/v1/chats}"
OPENWEBUI_TASKS_BY_CHAT_PREFIX="${OPENWEBUI_TASKS_BY_CHAT_PREFIX:-/api/tasks/chat}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-}"
OPENWEBUI_TOKEN="${OPENWEBUI_TOKEN:-}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
OPENWEBUI_DOCKER_SERVICE="${OPENWEBUI_DOCKER_SERVICE:-openwebui}"
COMPOSE_FILE="${COMPOSE_FILE:-docker-compose.yml}"

TOOL_INVOCATION_MODEL="${TOOL_INVOCATION_MODEL:-}"
TOOL_INVOCATION_MAX_TOKENS="${TOOL_INVOCATION_MAX_TOKENS:-512}"
TOOL_INVOCATION_REQUEST_TIMEOUT="${TOOL_INVOCATION_REQUEST_TIMEOUT:-120}"
TOOL_INVOCATION_CHAT_RETRIES="${TOOL_INVOCATION_CHAT_RETRIES:-3}"
TOOL_INVOCATION_CHAT_RETRY_SLEEP_SECONDS="${TOOL_INVOCATION_CHAT_RETRY_SLEEP_SECONDS:-2}"
TOOL_INVOCATION_TREAT_429_AS_WARN="${TOOL_INVOCATION_TREAT_429_AS_WARN:-true}"
TOOL_INVOCATION_FUNCTION_CALLING="${TOOL_INVOCATION_FUNCTION_CALLING:-default}"
TOOL_INVOCATION_FORCE_TOOL_CALL="${TOOL_INVOCATION_FORCE_TOOL_CALL:-true}"
# NOTE:
# - OpenWebUI executes upstream "native" tool_calls only via the streaming handler.
#   If you set function_calling=native, this test forces stream=true.
# - For deterministic local regression tests, function_calling=default + stream=false
#   is usually the most stable path (OpenWebUI runs tools internally).
TOOL_INVOCATION_STREAM="${TOOL_INVOCATION_STREAM:-false}"
TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS="${TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS:-120}"
TOOL_INVOCATION_ASYNC_POLL_INTERVAL_SECONDS="${TOOL_INVOCATION_ASYNC_POLL_INTERVAL_SECONDS:-0.3}"

TOOL_INVOCATION_REQUIRE_WEB_SCRAPER="${TOOL_INVOCATION_REQUIRE_WEB_SCRAPER:-true}"
TOOL_INVOCATION_REQUIRE_LLM_COUNCIL="${TOOL_INVOCATION_REQUIRE_LLM_COUNCIL:-true}"
TOOL_INVOCATION_REQUIRE_RUN_CODE="${TOOL_INVOCATION_REQUIRE_RUN_CODE:-true}"

TOOL_INVOCATION_ASSERT_METHOD="${TOOL_INVOCATION_ASSERT_METHOD:-db}"

# Stub model mode:
# - Optional, used to make tool invocation deterministic without calling a real upstream LLM.
# - Disabled by default because OpenWebUI may persist provider config in its data volume,
#   so env overrides are not always respected as a routing mechanism.
TOOL_INVOCATION_USE_STUB_MODEL="${TOOL_INVOCATION_USE_STUB_MODEL:-false}"
TOOL_INVOCATION_STUB_MODEL_ID="${TOOL_INVOCATION_STUB_MODEL_ID:-openwebui-tool-caller-stub}"
TOOL_INVOCATION_STUB_REWIRE_OPENWEBUI="${TOOL_INVOCATION_STUB_REWIRE_OPENWEBUI:-true}"
TOOL_INVOCATION_STUB_REWIRE_WAIT_SECONDS="${TOOL_INVOCATION_STUB_REWIRE_WAIT_SECONDS:-120}"

OUTPUT_DIR="${OUTPUT_DIR:-/tmp/openwebui_tool_invocation_$(date +%Y%m%d_%H%M%S)}"
VERBOSE="${VERBOSE:-false}"

# Runtime state
TMP_DIR=""
OPENWEBUI_AUTH_TOKEN=""
SELECTED_MODEL=""
declare -a CREATED_CHAT_IDS=()

OPENWEBUI_REWIRED=false

SERVER_PORT=""
SERVER_PID=""
SERVER_NONCE=""
LAST_TOOL_ERROR=""
LLM_COUNCIL_OLD_VALVES_FILE=""
JUPYTER_NONCE_FILE=""

# Counters
TESTS_TOTAL=0
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_WARNINGS=0
TESTS_SKIPPED=0

declare -a FAILURES=()
declare -a WARNINGS=()

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     OpenWebUI Tool Invocation Integration Test             ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_section() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}" >&2
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${CYAN}ℹ $1${NC}"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

is_true() {
    local value
    value="$(printf "%s" "${1:-}" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

is_native_function_calling() {
    local value
    value="$(printf "%s" "${TOOL_INVOCATION_FUNCTION_CALLING:-}" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "native" ]
}

normalize_base_url() {
    local value="$1"
    printf "%s" "${value%/}"
}

normalize_path() {
    local value="$1"
    if [ -z "$value" ]; then
        printf "/"
        return
    fi
    if [[ $value != /* ]]; then
        value="/$value"
    fi
    printf "%s" "$value"
}

uri_encode() {
    local value="$1"
    jq -nr --arg v "$value" '$v|@uri'
}

start_test() {
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    print_step "$1"
}

pass_test() {
    TESTS_PASSED=$((TESTS_PASSED + 1))
    print_success "OK"
}

fail_test() {
    local reason="$1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILURES+=("$reason")
    print_error "$reason"
}

warn_test() {
    local reason="$1"
    TESTS_WARNINGS=$((TESTS_WARNINGS + 1))
    WARNINGS+=("$reason")
    print_warning "$reason"
}

skip_test() {
    local reason="$1"
    TESTS_SKIPPED=$((TESTS_SKIPPED + 1))
    print_info "Skipped: $reason"
}

read_dot007_key() {
    local wanted="$1"
    local dotfile="$HOME/.007"

    [ -f "$dotfile" ] || return 1

    awk -v k="$wanted" '
        /^[[:space:]]*#/ { next }
        /^[[:space:]]*$/ { next }
        {
            line=$0
            sub(/^[[:space:]]*export[[:space:]]+/, "", line)
            sep = index(line, "=") ? "=" : (index(line, ":") ? ":" : "")
            if (sep == "") next
            split(line, parts, sep)
            key=parts[1]
            sub(/^[[:space:]]+/, "", key); sub(/[[:space:]]+$/, "", key)
            if (key != k) next
            val = substr(line, length(parts[1]) + 2)
            sub(/^[[:space:]]+/, "", val); sub(/[[:space:]]+$/, "", val)
            if ((val ~ /^".*"$/) || (val ~ /^\x27.*\x27$/)) {
                val = substr(val, 2, length(val)-2)
            }
            print val
            exit 0
        }
    ' "$dotfile"
}

require_dependencies() {
    local missing=0

    for cmd in curl jq; do
        if ! command_exists "$cmd"; then
            print_error "Missing dependency: $cmd"
            missing=1
        fi
    done

    if ! command_exists python3; then
        print_error "Missing dependency: python3 (required for local test server)"
        missing=1
    fi

    if [ "$TOOL_INVOCATION_ASSERT_METHOD" = "db" ]; then
        if ! command_exists docker; then
            print_error "Missing dependency: docker (required for DB-backed invocation assertions)"
            missing=1
        fi
    fi

    if is_true "$TOOL_INVOCATION_USE_STUB_MODEL" && is_true "$TOOL_INVOCATION_STUB_REWIRE_OPENWEBUI"; then
        if ! command_exists docker; then
            print_error "Missing dependency: docker (required to rewire OpenWebUI to the stub model)"
            missing=1
        fi
        if ! init_compose_cmd >/dev/null 2>&1; then
            print_error "Missing dependency: docker compose (required to rewire OpenWebUI to the stub model)"
            missing=1
        fi
    fi

    if [ "$missing" -ne 0 ]; then
        exit 1
    fi
}

init_compose_cmd() {
    if docker compose version >/dev/null 2>&1; then
        echo "docker compose"
        return 0
    fi

    if command_exists docker-compose; then
        echo "docker-compose"
        return 0
    fi

    return 1
}

docker_compose() {
    local compose_cmd
    compose_cmd="$(init_compose_cmd)" || return 1
    # shellcheck disable=SC2086
    $compose_cmd -f "${COMPOSE_FILE:-docker-compose.yml}" "$@"
}

docker_compose_ps_q() {
    local service="$1"
    docker_compose ps -q "$service" 2>/dev/null || true
}

get_openwebui_container_id() {
    local cid=""
    cid="$(docker_compose_ps_q "$OPENWEBUI_DOCKER_SERVICE" | tr -d '\r' | head -n 1)"
    if [ -z "$cid" ]; then
        cid="$(docker ps -q --filter "name=^/${OPENWEBUI_DOCKER_SERVICE}$" 2>/dev/null | head -n 1 || true)"
    fi
    printf "%s" "$cid"
}

wait_for_chat_message_done_db() {
    local chat_id="$1"
    local message_id="$2"
    local label="$3"
    local status_file="$OUTPUT_DIR/${label}.status.json"

    local cid=""
    cid="$(get_openwebui_container_id)"
    if [ -z "$cid" ]; then
        LAST_TOOL_ERROR="Could not resolve OpenWebUI container id for DB polling"
        return 1
    fi

    local started elapsed
    started="$(date +%s)"

    while true; do
        # Query the assistant message row (OpenWebUI stores chat_message.id as "${chat_id}-${message_id}").
        CHAT_ID="$chat_id" MESSAGE_ID="$message_id" docker exec -i "$cid" env CHAT_ID="$chat_id" MESSAGE_ID="$message_id" python3 - <<'PY' >"$status_file" 2>/dev/null || true
import json
import os
import sqlite3

chat_id = (os.environ.get("CHAT_ID") or "").strip()
message_id = (os.environ.get("MESSAGE_ID") or "").strip()
row_id = f"{chat_id}-{message_id}" if chat_id and message_id else ""

def is_empty_error(err) -> bool:
    if err is None:
        return True
    if isinstance(err, (dict, list)):
        return len(err) == 0
    if isinstance(err, str):
        s = err.strip()
        if s in ("", "null", "None"):
            return True
        try:
            parsed = json.loads(s)
            return parsed is None or parsed == {} or parsed == []
        except Exception:
            return False
    return False

con = sqlite3.connect("/app/backend/data/webui.db")
cur = con.cursor()
cur.execute("SELECT done, error, role, model_id, created_at, updated_at, content, output FROM chat_message WHERE id=?", (row_id,))
row = cur.fetchone()
if not row:
    print(json.dumps({"exists": False, "id": row_id}))
    raise SystemExit(0)

done, error, role, model_id, created_at, updated_at, content_raw, output_raw = row
err_present = not is_empty_error(error)
err_summary = None
if err_present:
    err_summary = (error if isinstance(error, str) else str(error))[:500]

def _as_str(val) -> str:
    if val is None:
        return ""
    if isinstance(val, (dict, list)):
        try:
            return json.dumps(val, ensure_ascii=False)
        except Exception:
            return str(val)
    if isinstance(val, str):
        return val
    return str(val)

output_present = False
for raw in (content_raw, output_raw):
    s = _as_str(raw).strip()
    if s and s not in ("null", "None", "[]", "{}"):
        output_present = True
        break

print(
    json.dumps(
        {
            "exists": True,
            "id": row_id,
            "done": bool(done),
            "output_present": bool(output_present),
            "error_present": bool(err_present),
            "error_summary": err_summary,
            "role": role,
            "model_id": model_id,
            "created_at": created_at,
            "updated_at": updated_at,
        },
        ensure_ascii=False,
    )
)
PY

        if jq -e '.exists == true and (.error_present != true) and (.done == true or .output_present == true)' "$status_file" >/dev/null 2>&1; then
            return 0
        fi

        if jq -e '.exists == true and (.error_present == true)' "$status_file" >/dev/null 2>&1; then
            LAST_TOOL_ERROR="Chat message error: $(jq -r '.error_summary // \"unknown\"' "$status_file" 2>/dev/null || true)"
            return 1
        fi

        elapsed=$(($(date +%s) - started))
        if [ "$elapsed" -ge "$TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS" ]; then
            LAST_TOOL_ERROR="Timed out waiting for chat message to complete (db, ${TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS}s)"
            return 1
        fi

        sleep "$TOOL_INVOCATION_ASYNC_POLL_INTERVAL_SECONDS"
    done
}

wait_for_chat_message_contains_nonce_db() {
    local chat_id="$1"
    local message_id="$2"
    local label="$3"
    local nonce="$4"

    local cid=""
    cid="$(get_openwebui_container_id)"
    if [ -z "$cid" ]; then
        LAST_TOOL_ERROR="Could not resolve OpenWebUI container id for DB polling"
        return 1
    fi

    local started elapsed
    started="$(date +%s)"

    while true; do
        CHAT_ID="$chat_id" MESSAGE_ID="$message_id" NONCE="$nonce" docker exec -i "$cid" env CHAT_ID="$chat_id" MESSAGE_ID="$message_id" NONCE="$nonce" python3 - <<'PY' >"$OUTPUT_DIR/${label}.nonce.status.json" 2>/dev/null || true
import json
import os
import sqlite3

chat_id = (os.environ.get("CHAT_ID") or "").strip()
message_id = (os.environ.get("MESSAGE_ID") or "").strip()
nonce = (os.environ.get("NONCE") or "").strip()
row_id = f"{chat_id}-{message_id}" if chat_id and message_id else ""

def is_empty_error(err) -> bool:
    if err is None:
        return True
    if isinstance(err, (dict, list)):
        return len(err) == 0
    if isinstance(err, str):
        s = err.strip()
        if s in ("", "null", "None"):
            return True
        try:
            parsed = json.loads(s)
            return parsed is None or parsed == {} or parsed == []
        except Exception:
            return False
    return False

def parse_json(val):
    if val is None:
        return None
    if isinstance(val, (dict, list)):
        return val
    if isinstance(val, str):
        s = val.strip()
        if not s:
            return ""
        try:
            return json.loads(s)
        except Exception:
            return s
    return val

def extract_text(node) -> str:
    texts = []
    def walk(n):
        if isinstance(n, dict):
            t = n.get("type")
            if t in ("output_text", "input_text") and isinstance(n.get("text"), str):
                texts.append(n["text"])
            for v in n.values():
                walk(v)
        elif isinstance(n, list):
            for it in n:
                walk(it)
    walk(node)
    return "".join(texts)

con = sqlite3.connect("/app/backend/data/webui.db")
cur = con.cursor()
cur.execute("SELECT done, error, content, output FROM chat_message WHERE id=?", (row_id,))
row = cur.fetchone()
if not row:
    print(json.dumps({"exists": False, "id": row_id}))
    raise SystemExit(0)

done, error, content_raw, output_raw = row
err_present = not is_empty_error(error)
err_summary = None
if err_present:
    err_summary = (error if isinstance(error, str) else str(error))[:500]

content = parse_json(content_raw)
output = parse_json(output_raw)

text = extract_text(output)
if not text:
    if isinstance(content, str):
        text = content
    else:
        try:
            text = json.dumps(content, ensure_ascii=False)
        except Exception:
            text = str(content)

found = bool(nonce and (nonce in text))

print(
    json.dumps(
        {
            "exists": True,
            "id": row_id,
            "done": bool(done),
            "error_present": bool(err_present),
            "error_summary": err_summary,
            "found": bool(found),
        },
        ensure_ascii=False,
    )
)
PY

        if jq -e '.exists == true and (.error_present == true)' "$OUTPUT_DIR/${label}.nonce.status.json" >/dev/null 2>&1; then
            LAST_TOOL_ERROR="Chat message error: $(jq -r '.error_summary // \"unknown\"' "$OUTPUT_DIR/${label}.nonce.status.json" 2>/dev/null || true)"
            return 1
        fi

        if jq -e '.exists == true and .found == true' "$OUTPUT_DIR/${label}.nonce.status.json" >/dev/null 2>&1; then
            return 0
        fi

        elapsed=$(($(date +%s) - started))
        if [ "$elapsed" -ge "$TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS" ]; then
            LAST_TOOL_ERROR="Timed out waiting for nonce to appear in chat message (db, ${TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS}s)"
            return 1
        fi

        sleep "$TOOL_INVOCATION_ASYNC_POLL_INTERVAL_SECONDS"
    done
}

dump_chat_message_output_text_db() {
    local chat_id="$1"
    local message_id="$2"
    local dest_file="$3"

    local cid=""
    cid="$(get_openwebui_container_id)"
    if [ -z "$cid" ]; then
        return 1
    fi

    docker exec -i "$cid" env CHAT_ID="$chat_id" MESSAGE_ID="$message_id" python3 - <<'PY' >"$dest_file" 2>/dev/null || true
import json
import os
import sqlite3

chat_id = (os.environ.get("CHAT_ID") or "").strip()
message_id = (os.environ.get("MESSAGE_ID") or "").strip()
row_id = f"{chat_id}-{message_id}" if chat_id and message_id else ""

def parse_json(val):
    if val is None:
        return None
    if isinstance(val, (dict, list)):
        return val
    if isinstance(val, str):
        s = val.strip()
        if not s:
            return ""
        try:
            return json.loads(s)
        except Exception:
            return s
    return val

def extract_output_text(node) -> str:
    texts = []
    def walk(n):
        if isinstance(n, dict):
            # OpenWebUI often stores assistant output under output_text, but tool outputs
            # can land under function_call_output as input_text depending on tool calling mode.
            if n.get("type") in ("output_text", "input_text") and isinstance(n.get("text"), str):
                texts.append(n["text"])
            for v in n.values():
                walk(v)
        elif isinstance(n, list):
            for it in n:
                walk(it)
    walk(node)
    return "".join(texts)

con = sqlite3.connect("/app/backend/data/webui.db")
cur = con.cursor()
cur.execute("SELECT content, output FROM chat_message WHERE id=?", (row_id,))
row = cur.fetchone()
if not row:
    raise SystemExit(0)

content_raw, output_raw = row
content = parse_json(content_raw)
output = parse_json(output_raw)

text = extract_output_text(output)
if not text:
    if isinstance(content, str):
        text = content
    else:
        text = json.dumps(content, ensure_ascii=False)

print(text, end="")
PY
}

request_get() {
    local url="$1"
    local output_file="$2"
    local auth_key="$3"

    if [ -n "$auth_key" ]; then
        curl -sS -m "$TOOL_INVOCATION_REQUEST_TIMEOUT" \
            -H "Authorization: Bearer $auth_key" \
            -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null || true
        return
    fi

    curl -sS -m "$TOOL_INVOCATION_REQUEST_TIMEOUT" \
        -o "$output_file" -w "%{http_code}" "$url" 2>/dev/null || true
}

try_mint_openwebui_admin_jwt() {
    if ! command_exists python3 || ! command_exists docker; then
        return 1
    fi

    if [ -z "${WEBUI_SECRET_KEY:-}" ]; then
        return 1
    fi

    local cid=""
    cid="$(docker_compose_ps_q "$OPENWEBUI_DOCKER_SERVICE" | tr -d '\r' | head -n 1)"
    if [ -z "$cid" ]; then
        cid="$(docker ps -q --filter "name=^/${OPENWEBUI_DOCKER_SERVICE}$" 2>/dev/null | head -n 1 || true)"
    fi
    if [ -z "$cid" ]; then
        return 1
    fi

    local admin_id=""
    admin_id="$(docker exec "$cid" python3 -c "import sqlite3; con=sqlite3.connect('/app/backend/data/webui.db'); cur=con.cursor(); cur.execute(\"SELECT id FROM user WHERE role='admin' ORDER BY created_at ASC LIMIT 1\"); row=cur.fetchone(); print(row[0] if row else '')" 2>/dev/null | tr -d '\r' | head -n 1 || true)"
    if [ -z "$admin_id" ]; then
        return 1
    fi

    local token=""
    token="$(
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
    "exp": int(time.time()) + 60 * 60 * 24 * 30,
    "jti": str(uuid.uuid4()),
}

header_b64 = b64url(json.dumps(header, separators=(",", ":"), sort_keys=True).encode("utf-8"))
payload_b64 = b64url(json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8"))
signing_input = f"{header_b64}.{payload_b64}".encode("utf-8")
sig = hmac.new(secret.encode("utf-8"), signing_input, hashlib.sha256).digest()
token = f"{header_b64}.{payload_b64}.{b64url(sig)}"
print(token)
PY
    )"

    if [ -z "$token" ]; then
        return 1
    fi

    OPENWEBUI_AUTH_TOKEN="$token"
    print_info "Authenticated via locally minted admin JWT (WEBUI_SECRET_KEY)"
    return 0
}

wait_for_openwebui_health() {
    local url="${OPENWEBUI_URL%/}${OPENWEBUI_HEALTH_PATH}"
    local started elapsed code
    started="$(date +%s)"

    while true; do
        code="$(curl -sS -m 5 -o /dev/null -w "%{http_code}" "$url" 2>/dev/null || true)"
        if [ "$code" = "200" ]; then
            return 0
        fi

        elapsed=$(($(date +%s) - started))
        if [ "$elapsed" -ge "$TOOL_INVOCATION_STUB_REWIRE_WAIT_SECONDS" ]; then
            return 1
        fi
        sleep 1
    done
}

rewire_openwebui_to_stub_model() {
    if ! is_true "$TOOL_INVOCATION_USE_STUB_MODEL" || ! is_true "$TOOL_INVOCATION_STUB_REWIRE_OPENWEBUI"; then
        return 0
    fi

    if ! command_exists docker; then
        return 1
    fi
    if ! init_compose_cmd >/dev/null 2>&1; then
        return 1
    fi

    local stub_base="http://host.docker.internal:${SERVER_PORT}/v1"
    print_step "Recreating OpenWebUI with stub model provider (${stub_base})"

    OPENAI_API_BASE_URL="$stub_base" OPENAI_API_BASE_URLS="$stub_base" \
        docker_compose up -d --force-recreate "$OPENWEBUI_DOCKER_SERVICE" >/dev/null 2>&1 || return 1

    OPENWEBUI_REWIRED=true

    if ! wait_for_openwebui_health; then
        return 1
    fi

    return 0
}

restore_openwebui_after_stub() {
    if [ "$OPENWEBUI_REWIRED" != "true" ]; then
        return 0
    fi

    if ! command_exists docker; then
        return 0
    fi
    if ! init_compose_cmd >/dev/null 2>&1; then
        return 0
    fi

    print_step "Restoring OpenWebUI container (default provider config)"
    docker_compose up -d --force-recreate "$OPENWEBUI_DOCKER_SERVICE" >/dev/null 2>&1 || true
    wait_for_openwebui_health || true
    OPENWEBUI_REWIRED=false
}

request_api() {
    local method="$1"
    local path="$2"
    local payload_file="$3"
    local output_file="$4"
    local url="${OPENWEBUI_URL%/}${path}"

    local -a args=(
        -sS
        -m "$TOOL_INVOCATION_REQUEST_TIMEOUT"
        -X "$method"
        "$url"
        -H "Accept: application/json"
        -o "$output_file"
        -w "%{http_code}"
    )

    if [ -n "$OPENWEBUI_AUTH_TOKEN" ]; then
        args+=(-H "Authorization: Bearer $OPENWEBUI_AUTH_TOKEN")
    fi
    if [ -n "$payload_file" ]; then
        args+=(-H "Content-Type: application/json" --data @"$payload_file")
    fi

    curl "${args[@]}" 2>/dev/null || true
}

openwebui_signin() {
    local output_file="$1"
    local signin_url="${OPENWEBUI_URL%/}${OPENWEBUI_SIGNIN_PATH}"
    local payload_file="$TMP_DIR/openwebui_signin_payload.json"

    jq -n \
        --arg email "$OPENWEBUI_SIGNIN_EMAIL" \
        --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
        '{email: $email, password: $password}' >"$payload_file"

    curl -sS -m "$TOOL_INVOCATION_REQUEST_TIMEOUT" \
        -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" \
        "$signin_url" -d @"$payload_file" 2>/dev/null || true
}

ensure_openwebui_auth_token() {
    if [ -n "$OPENWEBUI_TOKEN" ]; then
        OPENWEBUI_AUTH_TOKEN="$OPENWEBUI_TOKEN"
        return 0
    fi

    if [ -n "$OPENWEBUI_API_KEY" ]; then
        OPENWEBUI_AUTH_TOKEN="$OPENWEBUI_API_KEY"
        return 0
    fi

    # Optional: allow storing local OpenWebUI token under ~/.007.
    local dot_token
    dot_token="$(read_dot007_key "OPENWEBUI_TOKEN" 2>/dev/null || true)"
    if [ -z "$dot_token" ]; then
        dot_token="$(read_dot007_key "OPENWEBUI_API_KEY" 2>/dev/null || true)"
    fi
    if [ -z "$dot_token" ]; then
        dot_token="$(read_dot007_key "OPENWEBUI_LOCAL_API_KEY" 2>/dev/null || true)"
    fi
    if [ -n "$dot_token" ]; then
        OPENWEBUI_AUTH_TOKEN="$dot_token"
        return 0
    fi

    if ! is_true "$OPENWEBUI_AUTO_AUTH"; then
        try_mint_openwebui_admin_jwt >/dev/null 2>&1 || true
        return 0
    fi

    local signin_file="$TMP_DIR/openwebui_signin.json"
    local code
    code="$(openwebui_signin "$signin_file")"
    if [ "$code" = "200" ]; then
        OPENWEBUI_AUTH_TOKEN="$(jq -r '.token // empty' "$signin_file" 2>/dev/null || true)"
    fi

    if [ -z "$OPENWEBUI_AUTH_TOKEN" ]; then
        try_mint_openwebui_admin_jwt >/dev/null 2>&1 || true
    fi
    return 0
}

extract_model_ids() {
    local catalog_file="$1"
    jq -r '
        if type == "array" then
            .
        elif (type == "object" and (.data | type == "array")) then
            .data
        else
            []
        end
        | .[]?
        | .id // empty
    ' "$catalog_file" 2>/dev/null || true
}

pick_model_from_catalog() {
    local models_file="$1"
    local preferred
    preferred="$(
        cat <<'EOF'
openwebui-tool-caller-stub
gpt-5.3-codex
openai-codex
glm-5
minimax/chat-elite
minmax
qwen3-coder-flash
EOF
    )"

    local id
    for id in $preferred; do
        if extract_model_ids "$models_file" | grep -Fxq "$id"; then
            printf "%s" "$id"
            return 0
        fi
    done

    extract_model_ids "$models_file" | sed '/^$/d' | head -n 1
}

ensure_model_selected() {
    if [ -n "$TOOL_INVOCATION_MODEL" ]; then
        SELECTED_MODEL="$TOOL_INVOCATION_MODEL"
        return 0
    fi

    local models_file="$TMP_DIR/openwebui_models.json"
    local code
    code="$(request_get "${OPENWEBUI_URL%/}${OPENWEBUI_MODELS_PATH}" "$models_file" "$OPENWEBUI_AUTH_TOKEN")"
    if [ "$code" != "200" ]; then
        return 1
    fi

    SELECTED_MODEL="$(pick_model_from_catalog "$models_file" | tr -d '\r' || true)"
    [ -n "$SELECTED_MODEL" ]
}

generate_nonce() {
    python3 - <<'PY'
import secrets
print(secrets.token_hex(16))
PY
}

pick_free_port() {
    python3 - <<'PY'
import socket
s = socket.socket()
s.bind(("0.0.0.0", 0))
print(s.getsockname()[1])
s.close()
PY
}

write_local_server() {
    local server_file="$TMP_DIR/openwebui_tool_test_server.py"
    cat >"$server_file" <<'PY'
from __future__ import annotations

import json
import os
import re
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

NONCE = os.environ.get("OPENWEBUI_TOOL_NONCE", "").strip()
MODEL_ID = os.environ.get("OPENWEBUI_TOOL_MODEL_ID", "council-test-model").strip()
PORT = int(os.environ.get("PORT", "0"))

def _as_text(content) -> str:
    if isinstance(content, str):
        return content
    try:
        return json.dumps(content, ensure_ascii=False)
    except Exception:
        return str(content)


def _extract_user_text(messages) -> str:
    parts = []
    for m in messages or []:
        if not isinstance(m, dict):
            continue
        if m.get("role") not in ("system", "user", "assistant"):
            continue
        c = m.get("content")
        if isinstance(c, str) and c:
            parts.append(c)
    return "\n".join(parts)


def _extract_tool_outputs(messages) -> list[str]:
    outs: list[str] = []
    for m in messages or []:
        if not isinstance(m, dict):
            continue
        if m.get("role") != "tool":
            continue
        outs.append(_as_text(m.get("content", "")))
    return outs


def _available_tool_names(req: dict) -> list[str]:
    tools = req.get("tools") or []
    names: list[str] = []
    for t in tools:
        if not isinstance(t, dict):
            continue
        if t.get("type") != "function":
            continue
        fn = t.get("function")
        if not isinstance(fn, dict):
            continue
        name = fn.get("name")
        if isinstance(name, str) and name:
            names.append(name)
    return names


def _tool_choice_name(req: dict) -> str | None:
    tc = req.get("tool_choice")
    if not isinstance(tc, dict):
        return None
    fn = tc.get("function")
    if not isinstance(fn, dict):
        return None
    name = fn.get("name")
    if isinstance(name, str) and name:
        return name
    return None


def _pick_function_name(req: dict) -> str:
    names = _available_tool_names(req)
    chosen = _tool_choice_name(req)
    if chosen and (not names or chosen in names):
        return chosen

    text = _extract_user_text(req.get("messages"))
    if text and names:
        for name in names:
            if re.search(rf"\\b{re.escape(name)}\\b", text):
                return name

    if names:
        return names[0]

    if text:
        for known in ("scrape_url", "search_web", "council", "run_code"):
            if re.search(rf"\\b{known}\\b", text):
                return known

    return "council"


def _build_arguments(fn_name: str, req: dict) -> dict:
    text = _extract_user_text(req.get("messages"))

    if fn_name == "scrape_url":
        m = re.search(r"https?://[^\\s\"'<>]+", text)
        url = m.group(0) if m else f"http://host.docker.internal:{PORT}/scrape.html"
        return {"url": url}

    if fn_name == "search_web":
        query = ""
        m = re.search(r"query\\s*=\\s*[\"']([^\"']+)[\"']", text)
        if m:
            query = m.group(1).strip()
        if not query:
            m = re.search(r"query\\s*[:=]\\s*([^\\n]+)", text)
            if m:
                query = m.group(1).strip().strip("\"' ")

        num_results = 3
        m = re.search(r"num_results\\s*=\\s*(\\d+)", text)
        if m:
            try:
                num_results = int(m.group(1))
            except Exception:
                num_results = 3

        return {"query": query or "openwebui tool invocation test", "num_results": num_results}

    if fn_name == "run_code":
        path = ""
        # Prefer explicit absolute paths (the test includes one path on its own line).
        m = re.search(r"^(/[^\\s]+)\\s*$", text, flags=re.M)
        if m:
            path = m.group(1).strip()
        if not path:
            m = re.search(r"(/home/[^\\s]+)", text)
            if m:
                path = m.group(1).strip()
        if not path:
            path = "/home/jovyan/work/nonce.txt"

        code = f"with open(r'{path}', 'r') as f:\\n    print(f.read(), end='')\\n"
        return {"language": "python", "code": code}

    # Default: llm_council.council
    question = ""
    m = re.search(r"question\\s*[:=]\\s*[\"']([^\"']+)[\"']", text)
    if m:
        question = m.group(1).strip()
    return {"question": question or "Say council-ready."}


class Handler(BaseHTTPRequestHandler):
    server_version = "OpenWebUI-ToolTestServer/1.0"

    def log_message(self, fmt: str, *args) -> None:
        # Keep server logs quiet unless explicitly enabled.
        if os.environ.get("TOOL_SERVER_LOG", "0") in ("1", "true", "yes"):
            super().log_message(fmt, *args)

    def _send(self, status: int, content_type: str, body: bytes) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        path = self.path.split("?", 1)[0]

        if path in ("/v1/models", "/models"):
            payload = {
                "object": "list",
                "data": [{"id": MODEL_ID, "object": "model"}],
            }
            body = json.dumps(payload).encode("utf-8")
            self._send(200, "application/json; charset=utf-8", body)
            return

        if path == "/scrape.html":
            html = f"""<!doctype html>
<html lang="en">
  <head><meta charset="utf-8"><title>Scrape Test</title></head>
  <body>
    <h1>OpenWebUI Tool Invocation Test</h1>
    <p>nonce: {NONCE}</p>
  </body>
</html>
"""
            body = html.encode("utf-8")
            self._send(200, "text/html; charset=utf-8", body)
            return

        if path == "/nonce.txt":
            body = (NONCE + "\n").encode("utf-8")
            self._send(200, "text/plain; charset=utf-8", body)
            return

        self._send(404, "text/plain; charset=utf-8", b"not found\n")

    def do_POST(self) -> None:
        path = self.path.split("?", 1)[0]

        if path in ("/v1/chat/completions", "/chat/completions"):
            length = int(self.headers.get("Content-Length", "0") or "0")
            raw = self.rfile.read(length) if length else b"{}"
            try:
                req = json.loads(raw.decode("utf-8") or "{}")
            except Exception:
                req = {}

            tool_outputs = _extract_tool_outputs(req.get("messages"))
            if tool_outputs:
                content = "\n\n".join(tool_outputs)
                payload = {
                    "id": "cmpl-tool-test",
                    "object": "chat.completion",
                    "created": int(time.time()),
                    "model": MODEL_ID,
                    "choices": [
                        {
                            "index": 0,
                            "message": {"role": "assistant", "content": content},
                            "finish_reason": "stop",
                        }
                    ],
                }
                body = json.dumps(payload).encode("utf-8")
                self._send(200, "application/json; charset=utf-8", body)
                return

            tools = req.get("tools")
            if isinstance(tools, list) and tools:
                fn_name = _pick_function_name(req)
                args = _build_arguments(fn_name, req)
                call_id = "call_" + uuid.uuid4().hex[:24]
                payload = {
                    "id": "cmpl-tool-test",
                    "object": "chat.completion",
                    "created": int(time.time()),
                    "model": MODEL_ID,
                    "choices": [
                        {
                            "index": 0,
                            "message": {
                                "role": "assistant",
                                "content": None,
                                "tool_calls": [
                                    {
                                        "id": call_id,
                                        "type": "function",
                                        "function": {
                                            "name": fn_name,
                                            "arguments": json.dumps(args, separators=(",", ":")),
                                        },
                                    }
                                ],
                            },
                            "finish_reason": "tool_calls",
                        }
                    ],
                }
                body = json.dumps(payload).encode("utf-8")
                self._send(200, "application/json; charset=utf-8", body)
                return

            content = f"council-server-ok nonce={NONCE} model={MODEL_ID}"
            payload = {
                "id": "cmpl-tool-test",
                "object": "chat.completion",
                "created": int(time.time()),
                "model": MODEL_ID,
                "choices": [
                    {
                        "index": 0,
                        "message": {"role": "assistant", "content": content},
                        "finish_reason": "stop",
                    }
                ],
            }
            body = json.dumps(payload).encode("utf-8")
            self._send(200, "application/json; charset=utf-8", body)
            return

        self._send(404, "text/plain; charset=utf-8", b"not found\n")


def main() -> None:
    if not NONCE:
        raise SystemExit("OPENWEBUI_TOOL_NONCE missing")

    server = ThreadingHTTPServer(("0.0.0.0", PORT), Handler)
    server.serve_forever()


if __name__ == "__main__":
    main()
PY

    echo "$server_file"
}

start_local_server() {
    SERVER_NONCE="$(generate_nonce)"
    SERVER_PORT="$(pick_free_port)"

    local server_file
    server_file="$(write_local_server)"

    OPENWEBUI_TOOL_NONCE="$SERVER_NONCE" OPENWEBUI_TOOL_MODEL_ID="$TOOL_INVOCATION_STUB_MODEL_ID" PORT="$SERVER_PORT" \
        python3 "$server_file" >"$OUTPUT_DIR/tool_test_server.log" 2>&1 &
    SERVER_PID="$!"

    # Wait until the server is reachable from the host.
    local attempt=1
    while [ "$attempt" -le 30 ]; do
        if curl -fsS "http://127.0.0.1:${SERVER_PORT}/nonce.txt" >/dev/null 2>&1; then
            return 0
        fi
        sleep 0.1
        attempt=$((attempt + 1))
    done

    return 1
}

stop_local_server() {
    if [ -n "$SERVER_PID" ]; then
        kill "$SERVER_PID" >/dev/null 2>&1 || true
        wait "$SERVER_PID" >/dev/null 2>&1 || true
        SERVER_PID=""
    fi
}

restore_llm_council_valves() {
    if [ -z "$LLM_COUNCIL_OLD_VALVES_FILE" ] || [ ! -f "$LLM_COUNCIL_OLD_VALVES_FILE" ]; then
        return 0
    fi

    local tool_id_encoded
    tool_id_encoded="$(uri_encode "llm_council")"
    request_api "POST" "/api/v1/tools/id/${tool_id_encoded}/valves/update" "$LLM_COUNCIL_OLD_VALVES_FILE" \
        "$OUTPUT_DIR/llm_council_valves_restore.json" >/dev/null || true
}

cleanup() {
    stop_local_server

    if [ -n "$JUPYTER_NONCE_FILE" ]; then
        docker_compose exec -T jupyter sh -lc "rm -f '$JUPYTER_NONCE_FILE'" >/dev/null 2>&1 || true
        JUPYTER_NONCE_FILE=""
    fi

    if [ "${#CREATED_CHAT_IDS[@]}" -gt 0 ]; then
        local chat_id chat_id_encoded
        for chat_id in "${CREATED_CHAT_IDS[@]}"; do
            chat_id_encoded="$(uri_encode "$chat_id")"
            request_api "DELETE" "${OPENWEBUI_CHAT_BY_ID_PREFIX}/${chat_id_encoded}" "" \
                "$OUTPUT_DIR/openwebui_chat_cleanup_${chat_id}.json" >/dev/null || true
        done
    fi

    restore_llm_council_valves

    restore_openwebui_after_stub

    if [ -n "$TMP_DIR" ] && [ -d "$TMP_DIR" ]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT

openwebui_chat_completion() {
    local chat_id="$1"
    local message_id="$2"
    local session_id="$3"
    local tool_id="$4"
    local prompt="$5"
    local output_file="$6"
    local payload_file="$7"
    local function_name="${8:-}"
    local url="${OPENWEBUI_URL%/}${OPENWEBUI_CHAT_PATH}"

    local stream_value="false"
    if is_true "$TOOL_INVOCATION_STREAM"; then
        stream_value="true"
    fi

    jq -n \
        --arg model "$SELECTED_MODEL" \
        --arg prompt "$prompt" \
        --arg tool_id "$tool_id" \
        --arg chat_id "$chat_id" \
        --arg message_id "$message_id" \
        --arg session_id "$session_id" \
        --arg function_calling "$TOOL_INVOCATION_FUNCTION_CALLING" \
        --arg function_name "$function_name" \
        --arg force_tool_call "$TOOL_INVOCATION_FORCE_TOOL_CALL" \
        --argjson stream "$stream_value" \
        --argjson max_tokens "$TOOL_INVOCATION_MAX_TOKENS" \
        '{
            model: $model,
            messages: [{role: "user", content: $prompt}],
            max_tokens: $max_tokens,
            temperature: 0,
            stream: $stream,
            tool_ids: [$tool_id],
            chat_id: $chat_id,
            id: $message_id,
            session_id: $session_id,
            params: { function_calling: $function_calling }
        }
        + (
            # Only send OpenAI-style tool_choice when OpenWebUI is configured to use upstream
            # native tool calling (params.function_calling=native). In "default" mode OpenWebUI
            # executes tools internally and forwarding tool_choice upstream can break providers.
            if (
                ($function_name != "")
                and (($force_tool_call|ascii_downcase) | IN("true","1","yes"))
                and (($function_calling|ascii_downcase) == "native")
            )
            then {tool_choice: {type: "function", function: {name: $function_name}}}
            else {}
            end
        )' >"$payload_file"

    if is_true "$VERBOSE"; then
        jq '.tool_ids' "$payload_file" >/dev/null 2>&1 || true
    fi

    if [ -n "$OPENWEBUI_AUTH_TOKEN" ]; then
        curl -sS -m "$TOOL_INVOCATION_REQUEST_TIMEOUT" \
            -H "Authorization: Bearer $OPENWEBUI_AUTH_TOKEN" \
            -H "Content-Type: application/json" \
            -o "$output_file" -w "%{http_code}" \
            "$url" -d @"$payload_file" 2>/dev/null || true
        return
    fi

    curl -sS -m "$TOOL_INVOCATION_REQUEST_TIMEOUT" \
        -H "Content-Type: application/json" \
        -o "$output_file" -w "%{http_code}" \
        "$url" -d @"$payload_file" 2>/dev/null || true
}

create_ephemeral_chat() {
    local title="$1"
    local output_file="$2"

    local payload_file="$TMP_DIR/openwebui_create_chat.payload.json"

    jq -n --arg title "$title" '{chat:{title:$title,history:{messages:{}}}}' >"$payload_file"

    local code
    code="$(request_api "POST" "$OPENWEBUI_CHATS_NEW_PATH" "$payload_file" "$output_file")"
    if [ "$code" != "200" ]; then
        return 1
    fi

    local chat_id
    chat_id="$(jq -r '.id // empty' "$output_file" 2>/dev/null || true)"
    if [ -z "$chat_id" ]; then
        return 1
    fi

    CREATED_CHAT_IDS+=("$chat_id")
    printf "%s" "$chat_id"
}

delete_chat_best_effort() {
    local chat_id="$1"
    [ -n "$chat_id" ] || return 0

    local chat_id_encoded
    chat_id_encoded="$(uri_encode "$chat_id")"
    request_api "DELETE" "${OPENWEBUI_CHAT_BY_ID_PREFIX}/${chat_id_encoded}" "" \
        "$OUTPUT_DIR/openwebui_delete_chat_${chat_id}.json" >/dev/null || true
}

wait_for_chat_task_completion() {
    local chat_id="$1"
    local task_id="$2"

    local chat_id_encoded
    chat_id_encoded="$(uri_encode "$chat_id")"

    local tasks_file="$OUTPUT_DIR/openwebui_tasks_${chat_id}.json"
    local code elapsed started
    started="$(date +%s)"

    while true; do
        code="$(request_api "GET" "${OPENWEBUI_TASKS_BY_CHAT_PREFIX}/${chat_id_encoded}" "" "$tasks_file")"
        if [ "$code" != "200" ]; then
            LAST_TOOL_ERROR="Task poll failed (HTTP $code)"
            return 1
        fi

        if [ -z "$task_id" ]; then
            # Best effort: consider completed when there are no active tasks for the chat.
            if jq -e '(.task_ids // []) | length == 0' "$tasks_file" >/dev/null 2>&1; then
                return 0
            fi
        else
            if jq -e --arg tid "$task_id" '((.task_ids // []) | index($tid)) == null' "$tasks_file" >/dev/null 2>&1; then
                return 0
            fi
        fi

        elapsed=$(($(date +%s) - started))
        if [ "$elapsed" -ge "$TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS" ]; then
            LAST_TOOL_ERROR="Timed out waiting for background chat task to complete (${TOOL_INVOCATION_ASYNC_WAIT_TIMEOUT_SECONDS}s)"
            return 1
        fi

        sleep "$TOOL_INVOCATION_ASYNC_POLL_INTERVAL_SECONDS"
    done
}

run_tool_chat_and_assert_nonce() {
    local tool_id="$1"
    local prompt="$2"
    local label="$3"
    local nonce="$4"
    local function_name="${5:-}"

    local out_file="$OUTPUT_DIR/${label}.response.json"
    local payload_file="$OUTPUT_DIR/${label}.payload.json"
    local assistant_text_file="$OUTPUT_DIR/${label}.assistant.txt"
    local attempt=1
    local code
    local err=""
    local chat_id=""

    LAST_TOOL_ERROR=""

    while [ "$attempt" -le "$TOOL_INVOCATION_CHAT_RETRIES" ]; do
        chat_id="$(create_ephemeral_chat "Tool Invocation Test :: ${label}" "$OUTPUT_DIR/${label}.chat.new.json" 2>/dev/null || true)"
        if [ -z "$chat_id" ]; then
            err="failed to create chat (missing auth token or /api/v1/chats/new error)"
        else
            local message_id session_id
            message_id="tool-msg-$(generate_nonce)"
            session_id="tool-session-$(generate_nonce)"

            code="$(openwebui_chat_completion "$chat_id" "$message_id" "$session_id" "$tool_id" "$prompt" "$out_file" "$payload_file" "$function_name")"
            if [ "$code" = "200" ]; then
                if [ "$TOOL_INVOCATION_ASSERT_METHOD" = "db" ]; then
                    local task_id=""
                    task_id="$(jq -r '.task_id // empty' "$out_file" 2>/dev/null || true)"
                    if ! wait_for_chat_task_completion "$chat_id" "$task_id"; then
                        err="${LAST_TOOL_ERROR:-failed waiting for background chat task completion}"
                    elif ! wait_for_chat_message_contains_nonce_db "$chat_id" "$message_id" "$label" "$nonce"; then
                        err="${LAST_TOOL_ERROR:-nonce did not appear in chat message}"
                    else
                        dump_chat_message_output_text_db "$chat_id" "$message_id" "$assistant_text_file" || true
                        delete_chat_best_effort "$chat_id"
                        LAST_TOOL_ERROR=""
                        return 0
                    fi
                else
                    err="unsupported TOOL_INVOCATION_ASSERT_METHOD=${TOOL_INVOCATION_ASSERT_METHOD}"
                fi
            elif [ "$code" = "429" ] && is_true "$TOOL_INVOCATION_TREAT_429_AS_WARN"; then
                delete_chat_best_effort "$chat_id"
                LAST_TOOL_ERROR="HTTP 429 rate limit"
                return 2
            else
                err="$(jq -r '.error.message // .detail // .message // "unknown error"' "$out_file" 2>/dev/null || true)"
            fi
        fi

        delete_chat_best_effort "$chat_id"

        if [ "$attempt" -lt "$TOOL_INVOCATION_CHAT_RETRIES" ]; then
            sleep "$TOOL_INVOCATION_CHAT_RETRY_SLEEP_SECONDS"
        fi
        attempt=$((attempt + 1))
    done

    LAST_TOOL_ERROR="HTTP ${code:-unknown}: ${err:-unknown error}"
    return 1
}

test_web_scraper_scrape_url() {
    local url="http://host.docker.internal:${SERVER_PORT}/scrape.html"
    local prompt
    prompt="$(
        cat <<EOF
Call the tool function scrape_url on this URL: ${url}
Return ONLY the scraped text. Do not add commentary, do not summarize.
EOF
    )"

    start_test "web_scraper.scrape_url returns server nonce (proves tool executed)"
    if run_tool_chat_and_assert_nonce "web_scraper" "$prompt" "web_scraper_scrape_url" "$SERVER_NONCE" "scrape_url"; then
        pass_test
        return 0
    fi

    local status=$?
    if [ "$status" -eq 2 ]; then
        warn_test "web_scraper.scrape_url hit HTTP 429 rate limit"
        return 0
    fi

    fail_test "web_scraper.scrape_url failed (${LAST_TOOL_ERROR:-unknown error})"
    return 0
}

test_web_scraper_search_web_best_effort() {
    local prompt
    prompt="$(
        cat <<'EOF'
Call the tool function search_web with query="openwebui tool invocation test" and num_results=3.
Return ONLY the raw JSON returned by the tool.
EOF
    )"

    start_test "web_scraper.search_web best-effort (expects JSON with source/results)"
    local out_file="$OUTPUT_DIR/web_scraper_search_web.response.json"
    local payload_file="$OUTPUT_DIR/web_scraper_search_web.payload.json"
    local code
    local chat_id message_id session_id task_id chat_id_encoded
    # shellcheck disable=SC2034
    local chat_file
    chat_id="$(create_ephemeral_chat "Tool Invocation Test :: web_scraper_search_web" "$OUTPUT_DIR/web_scraper_search_web.chat.new.json" 2>/dev/null || true)"
    if [ -z "$chat_id" ]; then
        warn_test "web_scraper.search_web could not create a chat (missing auth or /api/v1/chats/new failed)"
        return 0
    fi

    message_id="tool-msg-$(generate_nonce)"
    session_id="tool-session-$(generate_nonce)"

    code="$(openwebui_chat_completion "$chat_id" "$message_id" "$session_id" "web_scraper" "$prompt" "$out_file" "$payload_file" "search_web")"
    if [ "$code" = "429" ] && is_true "$TOOL_INVOCATION_TREAT_429_AS_WARN"; then
        delete_chat_best_effort "$chat_id"
        warn_test "web_scraper.search_web hit HTTP 429 rate limit"
        return 0
    fi
    if [ "$code" != "200" ]; then
        delete_chat_best_effort "$chat_id"
        warn_test "web_scraper.search_web returned HTTP $code (non-fatal)"
        return 0
    fi

    task_id="$(jq -r '.task_id // empty' "$out_file" 2>/dev/null || true)"
    if [ "$TOOL_INVOCATION_ASSERT_METHOD" != "db" ]; then
        delete_chat_best_effort "$chat_id"
        warn_test "web_scraper.search_web skipped (requires TOOL_INVOCATION_ASSERT_METHOD=db)"
        return 0
    fi

    if ! wait_for_chat_task_completion "$chat_id" "$task_id"; then
        delete_chat_best_effort "$chat_id"
        warn_test "web_scraper.search_web task did not complete in time (${LAST_TOOL_ERROR:-timeout})"
        return 0
    fi

    if ! wait_for_chat_message_done_db "$chat_id" "$message_id" "web_scraper_search_web"; then
        delete_chat_best_effort "$chat_id"
        warn_test "web_scraper.search_web did not complete in time (${LAST_TOOL_ERROR:-timeout})"
        return 0
    fi

    local assistant_text_file="$OUTPUT_DIR/web_scraper_search_web.assistant.txt"
    dump_chat_message_output_text_db "$chat_id" "$message_id" "$assistant_text_file" || true
    if [ ! -s "$assistant_text_file" ]; then
        delete_chat_best_effort "$chat_id"
        warn_test "web_scraper.search_web produced empty assistant content"
        return 0
    fi

    local content
    content="$(cat "$assistant_text_file" 2>/dev/null || true)"

    json_ok="$(
        python3 - "$content" <<'PY'
import json
import re
import sys

text = sys.argv[1]
text = text.strip()
if text.startswith("```"):
    text = re.sub(r"^```[a-zA-Z0-9_-]*\\s*", "", text)
    text = re.sub(r"\\s*```\\s*$", "", text)

# Try strict parse first.
try:
    obj = json.loads(text)
except Exception:
    # Heuristic: find the first {...} block.
    m = re.search(r"\\{.*\\}", text, flags=re.S)
    if not m:
        print("no-json")
        raise SystemExit(0)
    try:
        obj = json.loads(m.group(0))
    except Exception:
        print("bad-json")
        raise SystemExit(0)

source = obj.get("source")
results = obj.get("results")
if isinstance(source, str) and isinstance(results, list):
    print("ok")
else:
    print("shape")
PY
    )"

    case "$json_ok" in
        ok)
            pass_test
            ;;
        *)
            warn_test "web_scraper.search_web returned content but JSON parse/shape failed (${json_ok})"
            ;;
    esac

    delete_chat_best_effort "$chat_id"
}

test_llm_council_deterministic() {
    local tool_id_encoded
    tool_id_encoded="$(uri_encode "llm_council")"

    LLM_COUNCIL_OLD_VALVES_FILE="$OUTPUT_DIR/llm_council_valves.old.json"
    request_api "GET" "/api/v1/tools/id/${tool_id_encoded}/valves" "" "$LLM_COUNCIL_OLD_VALVES_FILE" >/dev/null || true

    local new_valves_file="$OUTPUT_DIR/llm_council_valves.override.json"
    jq -n \
        --arg models "$TOOL_INVOCATION_STUB_MODEL_ID" \
        --arg base_url "http://host.docker.internal:${SERVER_PORT}/v1" \
        --arg api_keys "" \
        --argjson timeout 30 \
        '{MODELS: $models, API_BASE_URLS: $base_url, API_KEYS: $api_keys, TIMEOUT: $timeout}' >"$new_valves_file"

    start_test "llm_council.council returns server nonce (valves override proves tool executed)"
    local update_file="$OUTPUT_DIR/llm_council_valves.update.response.json"
    local update_code
    update_code="$(request_api "POST" "/api/v1/tools/id/${tool_id_encoded}/valves/update" "$new_valves_file" "$update_file")"
    if [ "$update_code" != "200" ]; then
        fail_test "llm_council valves update failed (HTTP $update_code)"
        return 0
    fi

    local prompt
    prompt="$(
        cat <<'EOF'
Call the tool function council with question: "Say council-ready."
Return ONLY the tool result, verbatim. Do not add any text.
EOF
    )"

    # shellcheck disable=SC2034
    local result
    if run_tool_chat_and_assert_nonce "llm_council" "$prompt" "llm_council_council" "$SERVER_NONCE" "council"; then
        pass_test
        return 0
    fi

    local status=$?
    if [ "$status" -eq 2 ]; then
        warn_test "llm_council.council hit HTTP 429 rate limit"
        return 0
    fi

    fail_test "llm_council.council failed (${LAST_TOOL_ERROR:-unknown error})"
    return 0
}

test_run_code_nonce_file() {
    if ! command_exists docker; then
        if is_true "$TOOL_INVOCATION_REQUIRE_RUN_CODE"; then
            fail_test "run_code requires docker to create a nonce file in Jupyter, but docker is missing"
            return 0
        fi
        skip_test "run_code test skipped (docker missing)"
        return 0
    fi
    if ! init_compose_cmd >/dev/null 2>&1; then
        if is_true "$TOOL_INVOCATION_REQUIRE_RUN_CODE"; then
            fail_test "run_code requires docker compose to reach Jupyter, but compose is missing"
            return 0
        fi
        skip_test "run_code test skipped (docker compose missing)"
        return 0
    fi

    local run_id
    run_id="$(generate_nonce)"
    JUPYTER_NONCE_FILE="/home/jovyan/work/openwebui_tool_nonce_${run_id}.txt"

    if ! docker_compose ps --status running >/dev/null 2>&1; then
        if is_true "$TOOL_INVOCATION_REQUIRE_RUN_CODE"; then
            fail_test "docker compose stack is not running"
            return 0
        fi
        skip_test "run_code test skipped (compose stack not running)"
        return 0
    fi

    docker_compose exec -T jupyter sh -lc "printf '%s\n' '${SERVER_NONCE}' > '${JUPYTER_NONCE_FILE}'" >/dev/null 2>&1 || true

    local prompt
    prompt="$(
        cat <<EOF
Call the tool function run_code with Python code that prints the contents of this file exactly:
${JUPYTER_NONCE_FILE}
Return ONLY what the tool prints. No markdown, no extra text.
EOF
    )"

    start_test "run_code.run_code returns Jupyter nonce (proves code execution ran)"
    if run_tool_chat_and_assert_nonce "run_code" "$prompt" "run_code_nonce_file" "$SERVER_NONCE" "run_code"; then
        pass_test
        return 0
    fi

    local status=$?
    if [ "$status" -eq 2 ]; then
        warn_test "run_code.run_code hit HTTP 429 rate limit"
        return 0
    fi

    fail_test "run_code.run_code failed (${LAST_TOOL_ERROR:-unknown error})"
    return 0
}

write_summary() {
    local summary_file="$OUTPUT_DIR/summary.txt"

    {
        echo "OpenWebUI Tool Invocation Integration Test"
        echo "Timestamp: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        echo "Target: $OPENWEBUI_URL"
        echo "Model: ${SELECTED_MODEL:-<unset>}"
        echo "Output Directory: $OUTPUT_DIR"
        echo
        echo "Total:    $TESTS_TOTAL"
        echo "Passed:   $TESTS_PASSED"
        echo "Failed:   $TESTS_FAILED"
        echo "Warnings: $TESTS_WARNINGS"
        echo "Skipped:  $TESTS_SKIPPED"
        echo
        if [ "${#FAILURES[@]}" -gt 0 ]; then
            echo "Failures:"
            printf -- "- %s\n" "${FAILURES[@]}"
            echo
        fi
        if [ "${#WARNINGS[@]}" -gt 0 ]; then
            echo "Warnings:"
            printf -- "- %s\n" "${WARNINGS[@]}"
        fi
    } >"$summary_file"

    print_section "Summary"
    print_info "Output directory: $OUTPUT_DIR"
    print_info "Total=$TESTS_TOTAL Passed=$TESTS_PASSED Failed=$TESTS_FAILED Warnings=$TESTS_WARNINGS Skipped=$TESTS_SKIPPED"

    if [ "${#FAILURES[@]}" -gt 0 ]; then
        local item
        for item in "${FAILURES[@]}"; do
            print_error "$item"
        done
    fi
    if [ "${#WARNINGS[@]}" -gt 0 ]; then
        local item
        for item in "${WARNINGS[@]}"; do
            print_warning "$item"
        done
    fi
}

main() {
    OPENWEBUI_URL="$(normalize_base_url "$OPENWEBUI_URL")"
    OPENWEBUI_CHAT_PATH="$(normalize_path "$OPENWEBUI_CHAT_PATH")"
    OPENWEBUI_MODELS_PATH="$(normalize_path "$OPENWEBUI_MODELS_PATH")"
    OPENWEBUI_HEALTH_PATH="$(normalize_path "$OPENWEBUI_HEALTH_PATH")"
    OPENWEBUI_SIGNIN_PATH="$(normalize_path "$OPENWEBUI_SIGNIN_PATH")"
    OPENWEBUI_CHATS_NEW_PATH="$(normalize_path "$OPENWEBUI_CHATS_NEW_PATH")"
    OPENWEBUI_CHAT_BY_ID_PREFIX="$(normalize_path "$OPENWEBUI_CHAT_BY_ID_PREFIX")"
    OPENWEBUI_TASKS_BY_CHAT_PREFIX="$(normalize_path "$OPENWEBUI_TASKS_BY_CHAT_PREFIX")"

    mkdir -p "$OUTPUT_DIR"
    TMP_DIR="$(mktemp -d)"

    print_header
    print_info "Target OpenWebUI URL: $OPENWEBUI_URL"

    require_dependencies
    ensure_openwebui_auth_token
    # OpenWebUI's native tool_calls execution path requires streaming responses.
    if is_native_function_calling && ! is_true "$TOOL_INVOCATION_STREAM"; then
        TOOL_INVOCATION_STREAM="true"
        print_info "function_calling=native requires stream=true; forcing TOOL_INVOCATION_STREAM=true"
    fi

    print_section "Preflight"
    start_test "OpenWebUI health endpoint returns HTTP 200"
    local health_file="$OUTPUT_DIR/openwebui_health.json"
    local health_code
    health_code="$(request_get "${OPENWEBUI_URL}${OPENWEBUI_HEALTH_PATH}" "$health_file" "$OPENWEBUI_AUTH_TOKEN")"
    if [ "$health_code" = "200" ]; then
        pass_test
    else
        fail_test "OpenWebUI health check failed (HTTP $health_code)"
        write_summary
        exit 1
    fi

    start_test "Start local tool test server (nonce-backed)"
    if start_local_server; then
        pass_test
        print_info "Local test server: http://127.0.0.1:${SERVER_PORT}"
    else
        fail_test "Failed to start local test server"
        write_summary
        exit 1
    fi

    if is_true "$TOOL_INVOCATION_USE_STUB_MODEL" && is_true "$TOOL_INVOCATION_STUB_REWIRE_OPENWEBUI"; then
        start_test "Rewire OpenWebUI OpenAI provider to local stub model"
        if rewire_openwebui_to_stub_model; then
            pass_test
        else
            fail_test "Failed to rewire OpenWebUI to stub model provider (docker/compose or health check failed)"
            write_summary
            exit 1
        fi
    fi

    start_test "Select a chat model from /api/models"
    if ensure_model_selected; then
        pass_test
        print_info "Selected model: $SELECTED_MODEL"
    else
        fail_test "Failed to select model (set TOOL_INVOCATION_MODEL or fix /api/models auth)"
        write_summary
        exit 1
    fi

    print_section "Tool Invocation"

    if is_true "$TOOL_INVOCATION_REQUIRE_WEB_SCRAPER"; then
        test_web_scraper_scrape_url
        test_web_scraper_search_web_best_effort
    else
        skip_test "web_scraper invocation tests disabled (TOOL_INVOCATION_REQUIRE_WEB_SCRAPER=false)"
    fi

    if is_true "$TOOL_INVOCATION_REQUIRE_LLM_COUNCIL"; then
        test_llm_council_deterministic
        restore_llm_council_valves
        LLM_COUNCIL_OLD_VALVES_FILE=""
    else
        skip_test "llm_council invocation tests disabled (TOOL_INVOCATION_REQUIRE_LLM_COUNCIL=false)"
    fi

    if is_true "$TOOL_INVOCATION_REQUIRE_RUN_CODE"; then
        test_run_code_nonce_file
    else
        skip_test "run_code invocation tests disabled (TOOL_INVOCATION_REQUIRE_RUN_CODE=false)"
    fi

    write_summary

    if [ "$TESTS_FAILED" -gt 0 ]; then
        exit 1
    fi
}

main "$@"
