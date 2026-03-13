#!/bin/bash
###############################################################################
# Repair OpenWebUI local tools (SQLite DB)
#
# Fixes two common issues:
# 1) /api/v1/tools/id/{id}/valves/spec returns null
#    Root cause: OpenWebUI loads Tools() *instances* (not modules). If a tool
#    defines `class Valves(BaseModel)` at module-level, the instance must expose
#    `self.Valves = Valves` for the router to find it.
#
# 2) tool.specs drift (e.g. specs contain ["pubmed"] but code defines search_*).
#    Root cause: old/incorrect specs stored in webui.db. OpenWebUI executes tools
#    based on stored specs, so drift can break the tool at runtime.
#
# This script patches tool content in-place inside the running OpenWebUI
# container, recomputes specs using OpenWebUI's own spec generator, and (by
# default) restarts the container to clear in-memory caches.
#
# Safety:
# - Creates a timestamped DB backup via ./scripts/admin/backup-openwebui-db.sh
# - Never prints secrets (no tokens/keys required)
###############################################################################

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "${SELF_DIR}/../.." && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

OPENWEBUI_SERVICE="${OPENWEBUI_SERVICE:-openwebui}"
OPENWEBUI_CONTAINER_NAME="${OPENWEBUI_CONTAINER_NAME:-$OPENWEBUI_SERVICE}"
NO_RESTART="false"
TOOL_ID_FILTERS=()

usage() {
    cat <<'EOF'
Usage:
  ./scripts/admin/repair-openwebui-tools.sh [--no-restart] [--tool-id TOOL_ID...]

Options:
  --no-restart     Do not restart the openwebui container after patching.
  --tool-id ID     Repair only the specified tool ID (repeatable).
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --no-restart)
            NO_RESTART="true"
            shift 1
            ;;
        --tool-id)
            if [ $# -lt 2 ] || [ -z "${2:-}" ]; then
                print_error "--tool-id requires a value"
                usage
                exit 2
            fi
            TOOL_ID_FILTERS+=("$2")
            shift 2
            ;;
        --tool-id=*)
            TOOL_ID_FILTERS+=("${1#*=}")
            shift 1
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

print_header "Repair OpenWebUI Tools (DB + Specs + Valves)"

if ! docker info >/dev/null 2>&1; then
    print_error "Docker daemon is not running."
    exit 1
fi

if ! init_compose_cmd; then
    print_error "Docker Compose command unavailable"
    exit 1
fi

if ! compose_service_running "$OPENWEBUI_SERVICE"; then
    print_error "OpenWebUI service '$OPENWEBUI_SERVICE' is not running."
    print_error "Start it first: docker compose up -d ${OPENWEBUI_SERVICE}"
    exit 1
fi

print_section "Backup"
print_step "Creating DB backup"
./scripts/admin/backup-openwebui-db.sh >/dev/null
print_success "DB backup created under ./backups"

print_section "Patch Tools + Specs"
print_step "Patching tools in /app/backend/data/webui.db (inside container)"

docker_compose exec -T "$OPENWEBUI_SERVICE" sh -lc '
if command -v python3 >/dev/null 2>&1; then
    exec python3 - "$@"
fi
exec python - "$@"
' sh "${TOOL_ID_FILTERS[@]}" <<'PY'
from __future__ import annotations

import json
import re
import sqlite3
import time
import sys
from typing import List, Tuple

from open_webui.utils.plugin import load_tool_module_by_id, replace_imports
from open_webui.utils.tools import get_tool_specs


DB_PATH = "/app/backend/data/webui.db"
TOOL_ID_FILTERS = [arg for arg in sys.argv[1:] if arg]


def _spec_names(specs_json: str) -> List[str]:
    try:
        specs = json.loads(specs_json or "[]")
        if not isinstance(specs, list):
            return []
        names = []
        for s in specs:
            if isinstance(s, dict) and isinstance(s.get("name"), str):
                names.append(s["name"])
        return names
    except Exception:
        return []


def ensure_instance_exposes_valves(content: str) -> Tuple[str, bool]:
    if not content:
        return content, False

    if not re.search(r"^\s*class\s+Valves\b", content, flags=re.M):
        return content, False

    # If already present, do nothing.
    if re.search(r"^\s*self\.Valves\s*=\s*Valves\s*$", content, flags=re.M):
        return content, False

    # Prefer inserting right after `self.valves = Valves()`.
    m = re.search(r"^(?P<indent>\s*)self\.valves\s*=\s*Valves\s*\([^\n]*\)\s*$", content, flags=re.M)
    if m:
        indent = m.group("indent")
        insert = f"{indent}self.Valves = Valves"
        # Insert after the matched line.
        end = m.end()
        content = content[:end] + "\n" + insert + content[end:]
        return content, True

    # Otherwise, try inserting inside Tools.__init__ right after the def line.
    lines = content.splitlines()
    out = []
    in_init = False
    init_indent = ""
    inserted = False
    for line in lines:
        out.append(line)
        if not in_init and re.match(r"^\s*def\s+__init__\s*\(\s*self\s*\)\s*:\s*$", line):
            in_init = True
            init_indent = re.match(r"^\s*", line).group(0)  # type: ignore
            continue
        if in_init and not inserted:
            # Insert at the first body line (or right away if next line is blank/comment).
            body_indent = init_indent + " " * 4
            if line.strip() == "" or line.lstrip().startswith("#"):
                continue
            out.insert(len(out) - 1, f"{body_indent}self.Valves = Valves")
            inserted = True
        if in_init:
            # Exit __init__ when indentation returns to <= init indent and line is not blank.
            if line.strip() and not line.startswith(init_indent + " " * 4):
                in_init = False
    if inserted:
        return "\n".join(out) + ("\n" if content.endswith("\n") else ""), True
    return content, False


conn = sqlite3.connect(DB_PATH)
conn.row_factory = sqlite3.Row
cur = conn.execute(
    """
    SELECT id, name, content, specs
    FROM tool
    WHERE content IS NOT NULL
      AND TRIM(content) != ''
      AND id NOT LIKE 'server:%'
    ORDER BY id
    """
)
rows = {r["id"]: r for r in cur.fetchall()}

errors: List[str] = []
now = int(time.time())
selected_tool_ids: List[str] = []

if TOOL_ID_FILTERS:
    seen = set()
    for tool_id in TOOL_ID_FILTERS:
        if tool_id in seen:
            continue
        seen.add(tool_id)
        selected_tool_ids.append(tool_id)
else:
    selected_tool_ids = list(rows.keys())

if not selected_tool_ids:
    raise SystemExit("No local OpenWebUI tools found in webui.db")

repaired_count = 0

for tool_id in selected_tool_ids:
    row = rows.get(tool_id)
    if row is None:
        errors.append(f"missing tool row: {tool_id}")
        continue

    old_specs = _spec_names(row["specs"] or "")
    content = row["content"] or ""

    patched, added_valves_attr = ensure_instance_exposes_valves(content)
    patched = replace_imports(patched)

    try:
        tool_module, _frontmatter = load_tool_module_by_id(tool_id, content=patched)
        specs = get_tool_specs(tool_module)
    except Exception as exc:
        errors.append(f"{tool_id}: failed to load and/or compute specs: {exc}")
        continue

    new_specs = [s.get("name") for s in specs if isinstance(s, dict)]
    repaired_count += 1

    conn.execute(
        "UPDATE tool SET content=?, specs=?, updated_at=? WHERE id=?",
        (patched, json.dumps(specs), now, tool_id),
    )

    valves_note = " (+Valves)" if added_valves_attr else ""
    print(f"{tool_id}: specs {old_specs} -> {new_specs}{valves_note}")

conn.commit()
conn.close()

print(f"repaired={repaired_count} selected={len(selected_tool_ids)}")

if errors:
    raise SystemExit("Tool repair completed with errors:\n- " + "\n- ".join(errors))

print("OK")
PY

print_success "Tools updated and specs regenerated"

if [ "$NO_RESTART" = "true" ]; then
    print_warning "Skipping container restart (--no-restart)"
    exit 0
fi

print_section "Restart"
print_step "Restarting openwebui container to clear cached tool modules/specs"
docker_compose restart "$OPENWEBUI_SERVICE" >/dev/null
print_success "${OPENWEBUI_SERVICE} restarted"
