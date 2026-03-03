#!/bin/bash

###############################################################################
# Enable live SMILES retrieval in OpenWebUI
# 1) Start structure-search API
# 2) Upsert OpenWebUI tool in webui.db
# 3) Restart OpenWebUI to reload tools
# 4) Verify live connectivity and registration
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

SMILES_STRUCTURE_API_PORT="${SMILES_STRUCTURE_API_PORT:-8011}"
SMILES_TOOL_API_URL="${SMILES_TOOL_API_URL:-http://smiles-structure-api:${SMILES_STRUCTURE_API_PORT}/v1/structure}"
SMILES_TOOL_ID="${SMILES_TOOL_ID:-structure_search}"
SMILES_TOOL_NAME="${SMILES_TOOL_NAME:-Structure Search}"
SMILES_TOOL_SOURCE="${SMILES_TOOL_SOURCE:-${SCRIPT_DIR}/smiles-pipeline/plugins/structure_search_tool.py}"
SMILES_TOOL_CATEGORY="${SMILES_TOOL_CATEGORY:-chemistry}"
NO_RESTART="${NO_RESTART:-false}"

usage() {
    cat <<'EOF'
Usage: ./enable-smiles-retrieval-live.sh [--no-restart]

Options:
  --no-restart  Do not restart openwebui after DB upsert
EOF
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
    print_header "Enable Live SMILES Retrieval"

    if [ ! -f "$SMILES_TOOL_SOURCE" ]; then
        print_error "Tool source file not found: $SMILES_TOOL_SOURCE"
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

    print_step "Starting SMILES Structure API service"
    "${SCRIPT_DIR}/start-smiles-structure-api.sh"
    "${SCRIPT_DIR}/check-smiles-structure-api.sh" --quiet
    print_success "SMILES API is healthy"

    print_step "Ensuring OpenWebUI is running"
    if ! docker ps --format '{{.Names}}' | grep -qx 'openwebui'; then
        print_error "openwebui container is not running"
        exit 1
    fi

    print_step "Upserting '${SMILES_TOOL_NAME}' tool in OpenWebUI DB"
    local content_b64
    content_b64="$(base64 -w 0 "$SMILES_TOOL_SOURCE")"

    docker_compose exec -T openwebui python - "$SMILES_TOOL_ID" "$SMILES_TOOL_NAME" "$SMILES_TOOL_API_URL" "$SMILES_TOOL_CATEGORY" "$content_b64" <<'PY'
import base64
import json
import os
import re
import sqlite3
import time
import sys

from open_webui.utils.plugin import load_tool_module_by_id, replace_imports
from open_webui.utils.tools import get_tool_specs

db_path = "/app/backend/data/webui.db"
if len(sys.argv) != 6:
    raise RuntimeError("Expected 5 arguments: tool_id tool_name tool_api_url tool_category content_b64")

tool_id = sys.argv[1]
tool_name = sys.argv[2]
tool_api_url = sys.argv[3]
tool_category = sys.argv[4]
content = base64.b64decode(sys.argv[5]).decode("utf-8")

content = re.sub(
    r'default="http://[^"]+/v1/structure"',
    f'default="{tool_api_url}"',
    content,
    count=1,
)
content = replace_imports(content)

tool_module, _ = load_tool_module_by_id(tool_id, content=content)
specs = get_tool_specs(tool_module)
if not isinstance(specs, list) or not specs:
    raise RuntimeError("Generated tool specs are empty")

now = int(time.time())
meta = json.dumps({"version": "1.0.0", "category": tool_category}, ensure_ascii=True)
valves = json.dumps({"STRUCTURE_SEARCH_API_URL": tool_api_url}, ensure_ascii=True)
specs_json = json.dumps(specs, ensure_ascii=True)

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

spec_names = [item.get("name") for item in specs if isinstance(item, dict)]
print(f"{action}:{tool_id}")
print("specs:" + ",".join([name for name in spec_names if isinstance(name, str)]))
print("api_url:" + tool_api_url)
PY

    print_success "OpenWebUI tool row upserted"

    if [ "$NO_RESTART" != "true" ]; then
        print_step "Restarting openwebui to reload tool cache"
        docker_compose restart openwebui >/dev/null
        print_success "openwebui restarted"
    else
        print_warning "Skipping openwebui restart (--no-restart)"
    fi

    print_step "Verifying tool registration and runtime connectivity"
    docker_compose exec -T openwebui python - <<'PY'
import sqlite3

conn = sqlite3.connect("/app/backend/data/webui.db")
cur = conn.cursor()
cur.execute("SELECT id,name,is_active,valves FROM tool WHERE id = 'structure_search'")
row = cur.fetchone()
if not row:
    raise SystemExit("Tool structure_search missing from DB")
print(f"tool:{row[0]} name:{row[1]} active:{row[2]} valves:{row[3]}")
conn.close()
PY

    local health_url http_code
    health_url="${SMILES_TOOL_API_URL%/v1/structure}/health"
    http_code="$(docker exec openwebui sh -lc "curl -sS --connect-timeout 3 --max-time 8 -o /tmp/smiles_live_health.json -w '%{http_code}' '${health_url}'" || true)"
    if [ "$http_code" != "200" ]; then
        print_error "OpenWebUI container cannot reach SMILES API at ${health_url} (HTTP ${http_code})"
        exit 1
    fi

    print_success "Live SMILES retrieval is wired"
    print_info "Tool API URL: ${SMILES_TOOL_API_URL}"
    print_info "Try in chat with tool: search_similar_structures"
}

main "$@"
