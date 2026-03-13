#!/bin/bash

###############################################################################
# Enable Indigo Service tool in OpenWebUI
# 1) Start Indigo sidecar
# 2) Upsert Indigo tool row into OpenWebUI webui.db
# 3) Restart OpenWebUI to reload tool cache
# 4) Verify DB registration + container reachability
###############################################################################

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_SERVICE="${OPENWEBUI_SERVICE:-openwebui}"
OPENWEBUI_CONTAINER_NAME="${OPENWEBUI_CONTAINER_NAME:-$OPENWEBUI_SERVICE}"
INDIGO_TOOL_ID="${INDIGO_TOOL_ID:-indigo_chemistry}"
INDIGO_TOOL_NAME="${INDIGO_TOOL_NAME:-Indigo Chemistry}"
INDIGO_TOOL_SOURCE="${INDIGO_TOOL_SOURCE:-${SCRIPT_DIR}/local/plugins/indigo_chemistry_tool.py}"
INDIGO_TOOL_CATEGORY="${INDIGO_TOOL_CATEGORY:-chemistry}"
INDIGO_TOOL_API_BASE_URL="${INDIGO_TOOL_API_BASE_URL:-http://indigo-service/v2/indigo}"
INDIGO_TOOL_TIMEOUT_SECONDS="${INDIGO_TOOL_TIMEOUT_SECONDS:-20}"
NO_RESTART="${NO_RESTART:-false}"

usage() {
    cat <<'USAGE'
Usage: ./scripts/indigo/enable-indigo-live.sh [--no-restart]

Options:
  --no-restart  Do not restart openwebui after DB upsert
USAGE
}

while [ $# -gt 0 ]; do
    case "$1" in
        --no-restart)
            NO_RESTART=true
            shift
            ;;
        -h | --help)
            usage
            exit 0
            ;;
        *)
            print_error "Unknown argument: $1"
            usage
            exit 2
            ;;
    esac
done

main() {
    print_header "Enable Indigo Tool"

    if [ ! -f "$INDIGO_TOOL_SOURCE" ]; then
        print_error "Tool source file not found: $INDIGO_TOOL_SOURCE"
        exit 1
    fi

    if ! command_exists docker; then
        print_error "Docker is required"
        exit 1
    fi

    if ! docker info >/dev/null 2>&1; then
        print_error "Docker daemon is not running"
        exit 1
    fi

    if ! init_compose_cmd; then
        print_error "Docker Compose command unavailable"
        exit 1
    fi

    print_step "Starting Indigo Service sidecar"
    "${SELF_DIR}/start-indigo-service.sh"
    "${SELF_DIR}/check-indigo-service.sh" --quiet
    print_success "Indigo Service is healthy"

    print_step "Ensuring OpenWebUI container is running"
    if ! compose_service_running "$OPENWEBUI_SERVICE"; then
        print_error "OpenWebUI service '$OPENWEBUI_SERVICE' is not running"
        exit 1
    fi

    local openwebui_container_id
    openwebui_container_id="$(resolve_container_id "$OPENWEBUI_SERVICE" "$OPENWEBUI_CONTAINER_NAME" || true)"
    if [ -z "$openwebui_container_id" ]; then
        print_error "Unable to resolve running container for service '$OPENWEBUI_SERVICE'"
        exit 1
    fi

    print_step "Upserting Indigo tool row in OpenWebUI DB"
    local content_b64
    content_b64="$(base64 -w 0 "$INDIGO_TOOL_SOURCE")"

    docker_compose exec -T "$OPENWEBUI_SERVICE" sh -lc '
if command -v python3 >/dev/null 2>&1; then
    exec python3 - "$@"
fi
exec python - "$@"
' sh \
        "$INDIGO_TOOL_ID" \
        "$INDIGO_TOOL_NAME" \
        "$INDIGO_TOOL_API_BASE_URL" \
        "$INDIGO_TOOL_CATEGORY" \
        "$INDIGO_TOOL_TIMEOUT_SECONDS" \
        "$content_b64" <<'PY'
import base64
import json
import re
import sqlite3
import sys
import time

from open_webui.utils.plugin import load_tool_module_by_id, replace_imports
from open_webui.utils.tools import get_tool_specs

if len(sys.argv) != 7:
    raise RuntimeError(
        "Expected 6 args: tool_id tool_name api_base_url tool_category timeout_seconds content_b64"
    )

tool_id = sys.argv[1]
tool_name = sys.argv[2]
api_base_url = sys.argv[3]
tool_category = sys.argv[4]
timeout_seconds = int(sys.argv[5])
content = base64.b64decode(sys.argv[6]).decode("utf-8")

# Keep default valve URL aligned with deployment profile.
content = re.sub(
    r'default="http://indigo-service/v2/indigo"',
    f'default="{api_base_url}"',
    content,
    count=1,
)
content = replace_imports(content)

module, _ = load_tool_module_by_id(tool_id, content=content)
specs = get_tool_specs(module)
if not isinstance(specs, list) or not specs:
    raise RuntimeError("Generated Indigo tool specs are empty")

spec_names = [item.get("name") for item in specs if isinstance(item, dict)]

now = int(time.time())
meta = json.dumps({"version": "1.0.0", "category": tool_category}, ensure_ascii=True)
valves = json.dumps(
    {
        "INDIGO_API_BASE_URL": api_base_url,
        "TIMEOUT_SECONDS": timeout_seconds,
    },
    ensure_ascii=True,
)
specs_json = json.dumps(specs, ensure_ascii=True)

db_path = "/app/backend/data/webui.db"
conn = sqlite3.connect(db_path)
conn.row_factory = sqlite3.Row
cur = conn.cursor()

cur.execute("SELECT id FROM user WHERE role='admin' ORDER BY created_at LIMIT 1")
admin = cur.fetchone()
if admin is None:
    cur.execute("SELECT id FROM user ORDER BY created_at LIMIT 1")
    admin = cur.fetchone()
if admin is None:
    raise RuntimeError("Cannot resolve owner user_id for tool row")

user_id = admin["id"]
cur.execute("SELECT id, created_at FROM tool WHERE id = ?", (tool_id,))
existing = cur.fetchone()

if existing is None:
    cur.execute(
        """
        INSERT INTO tool (id, user_id, name, content, specs, meta, created_at, updated_at, valves, is_active)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 1)
        """,
        (tool_id, user_id, tool_name, content, specs_json, meta, now, now, valves),
    )
    action = "inserted"
else:
    created_at = int(existing["created_at"]) if existing["created_at"] else now
    cur.execute(
        """
        UPDATE tool
        SET user_id = ?, name = ?, content = ?, specs = ?, meta = ?, valves = ?, updated_at = ?, is_active = 1
        WHERE id = ?
        """,
        (user_id, tool_name, content, specs_json, meta, valves, now, tool_id),
    )
    cur.execute("UPDATE tool SET created_at = ? WHERE id = ?", (created_at, tool_id))
    action = "updated"

conn.commit()
conn.close()

print(f"{action}:{tool_id}")
print("specs:" + ",".join([name for name in spec_names if isinstance(name, str)]))
print("api_base_url:" + api_base_url)
PY

    print_success "Indigo tool row upserted"

    if [ "$NO_RESTART" != "true" ]; then
        print_step "Restarting openwebui to reload tool cache"
        docker_compose restart "$OPENWEBUI_SERVICE" >/dev/null
        print_success "${OPENWEBUI_SERVICE} restarted"
    else
        print_warning "Skipping openwebui restart (--no-restart)"
    fi

    print_step "Verifying tool registration"
    docker_compose exec -T "$OPENWEBUI_SERVICE" sh -lc '
if command -v python3 >/dev/null 2>&1; then
    exec python3 - "$@"
fi
exec python - "$@"
' sh \
        "$INDIGO_TOOL_ID" <<'PY'
import sqlite3
import sys

tool_id = sys.argv[1]
conn = sqlite3.connect("/app/backend/data/webui.db")
cur = conn.cursor()
cur.execute("SELECT id, name, is_active, valves FROM tool WHERE id = ?", (tool_id,))
row = cur.fetchone()
conn.close()

if not row:
    raise SystemExit(f"Tool {tool_id} missing from DB")
print(f"tool:{row[0]} name:{row[1]} active:{row[2]} valves:{row[3]}")
PY

    print_step "Verifying Indigo reachability from OpenWebUI container"
    local info_url http_code
    info_url="${INDIGO_TOOL_API_BASE_URL%/}/info"
    http_code="$(docker exec "$openwebui_container_id" sh -lc "curl -sS --connect-timeout 3 --max-time 8 -o /tmp/indigo_tool_info.json -w '%{http_code}' '${info_url}'" || true)"
    if [ "$http_code" != "200" ]; then
        print_error "OpenWebUI container cannot reach Indigo API at ${info_url} (HTTP ${http_code})"
        exit 1
    fi

    print_success "Indigo tool integration is ready"
    print_info "Tool ID: ${INDIGO_TOOL_ID}"
    print_info "API URL: ${INDIGO_TOOL_API_BASE_URL}"
}

main "$@"
