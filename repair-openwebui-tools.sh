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
# - Creates a timestamped DB backup via ./backup-openwebui-db.sh
# - Never prints secrets (no tokens/keys required)
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

NO_RESTART="false"

usage() {
    cat <<'EOF'
Usage:
  ./repair-openwebui-tools.sh [--no-restart]

Options:
  --no-restart   Do not restart the openwebui container after patching.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --no-restart)
            NO_RESTART="true"
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

if ! docker ps --format '{{.Names}}' | grep -qx 'openwebui'; then
    print_error "Container 'openwebui' is not running."
    print_error "Start it first: docker compose up -d openwebui"
    exit 1
fi

print_section "Backup"
print_step "Creating DB backup"
./backup-openwebui-db.sh >/dev/null
print_success "DB backup created under ./backups"

print_section "Patch Tools + Specs"
print_step "Patching tools in /app/backend/data/webui.db (inside container)"

docker_compose exec -T openwebui python - <<'PY'
from __future__ import annotations

import json
import re
import sqlite3
import time
from typing import List, Tuple

from open_webui.utils.plugin import load_tool_module_by_id, replace_imports
from open_webui.utils.tools import get_tool_specs


DB_PATH = "/app/backend/data/webui.db"
TOOL_IDS = [
    "web_scraper",
    "llm_council",
    "run_code",
    "tavily",
    "pubmed",
    "arxiv",
    "github",
]


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

    # If already present, do nothing.
    if re.search(r"^\s*self\.Valves\s*=\s*Valves\s*$", content, flags=re.M):
        return content, False

    # Prefer inserting right after `self.valves = Valves()`.
    m = re.search(r"^(?P<indent>\s*)self\.valves\s*=\s*Valves\(\)\s*$", content, flags=re.M)
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
cur = conn.execute("SELECT id, name, content, specs FROM tool ORDER BY id")
rows = {r["id"]: r for r in cur.fetchall()}

errors: List[str] = []
now = int(time.time())

for tool_id in TOOL_IDS:
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

    conn.execute(
        "UPDATE tool SET content=?, specs=?, updated_at=? WHERE id=?",
        (patched, json.dumps(specs), now, tool_id),
    )

    valves_note = " (+Valves)" if added_valves_attr else ""
    print(f"{tool_id}: specs {old_specs} -> {new_specs}{valves_note}")

conn.commit()
conn.close()

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
docker_compose restart openwebui >/dev/null
print_success "openwebui restarted"
