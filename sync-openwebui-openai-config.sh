#!/bin/bash
#
# Sync OpenWebUI "openai" persistent config (webui.db) with local .env values.
#
# Why:
# - OpenWebUI stores connection settings in its internal SQLite config (webui.db).
# - If OPENAI_API_KEY or OPENAI_API_BASE_URL are changed in .env but webui.db
#   still has old values, model listing and routing silently break.
#
# What it does:
# - Updates `.env`: ensures OPENAI_API_KEYS matches OPENAI_API_KEY (single key).
# - Patches OpenWebUI SQLite config (inside container): sets openai.api_base_urls,
#   openai.api_keys and enables the provider.
# - Restarts the OpenWebUI container to reload config.
#
# Safety:
# - Never prints secret values.
# - Creates a DB backup via ./backup-openwebui-db.sh
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

main() {
    load_env_defaults
    print_header

    if ! command_exists docker; then
        print_error "docker is required"
        exit 1
    fi
    if ! docker info >/dev/null 2>&1; then
        print_error "Docker daemon is not running"
        exit 1
    fi
    if ! command_exists python3; then
        print_error "python3 is required"
        exit 1
    fi

    local env_file="${ENV_FILE:-$SCRIPT_DIR/.env}"
    if [ ! -f "$env_file" ]; then
        print_error "Missing env file: $env_file"
        exit 1
    fi

    local openai_key="${OPENAI_API_KEY:-}"
    if [ -z "$openai_key" ]; then
        print_error "OPENAI_API_KEY is required in .env"
        exit 1
    fi

    if ! docker ps --format '{{.Names}}' | grep -Fxq "openwebui"; then
        print_error "Container 'openwebui' is not running"
        print_error "Start it first: docker compose up -d openwebui"
        exit 1
    fi

    print_step "Backing up OpenWebUI DB"
    if [ -x "./backup-openwebui-db.sh" ]; then
        ./backup-openwebui-db.sh >/dev/null
        print_success "DB backup created under ./backups"
    else
        print_warning "backup-openwebui-db.sh not found or not executable; skipping DB backup"
    fi

    print_step "Normalizing .env (OPENAI_API_KEYS)"
    OPENAI_KEY="$openai_key" ENV_PATH="$env_file" python3 - <<'PY'
import os
import re
from pathlib import Path

env_path = Path(os.environ["ENV_PATH"])
openai_key = os.environ.get("OPENAI_KEY", "")
if not openai_key:
    raise SystemExit("OPENAI_KEY is empty")

def upsert_env(path: Path, key: str, value: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines(True)
    out = []
    found = False
    key_re = re.compile(rf"^{re.escape(key)}=")
    for line in lines:
        if key_re.match(line):
            out.append(f"{key}={value}\n")
            found = True
        else:
            out.append(line)
    if not found:
        if out and not out[-1].endswith("\n"):
            out[-1] += "\n"
        out.append(f"{key}={value}\n")
    path.write_text("".join(out), encoding="utf-8")

# Ensure OpenWebUI does not read a stale OPENAI_API_KEYS list.
upsert_env(env_path, "OPENAI_API_KEYS", openai_key)
PY
    print_success ".env normalized"

    local openai_base_url="${OPENAI_API_BASE_URL:-}"

    print_step "Patching OpenWebUI persistent config (openai.api_keys, openai.api_base_urls)"
    docker exec -i \
        -e "OPENAI_KEY=$openai_key" \
        -e "OPENAI_BASE_URL=$openai_base_url" \
        openwebui python3 - <<'PY'
from __future__ import annotations

import json
import os
import sqlite3

DB = "/app/backend/data/webui.db"
openai_key = os.environ.get("OPENAI_KEY", "").strip()
openai_base_url = os.environ.get("OPENAI_BASE_URL", "").strip()
if not openai_key:
    raise SystemExit("OPENAI_KEY is required")

con = sqlite3.connect(DB)
cur = con.cursor()
row = cur.execute("SELECT data, version FROM config WHERE id = 1").fetchone()
if not row:
    raise SystemExit("config row not found (id=1)")

data, version = row
obj = json.loads(data)
openai = obj.get("openai")
if not isinstance(openai, dict):
    openai = {}

openai["api_keys"] = [openai_key]
openai["enable"] = True
if openai_base_url:
    openai["api_base_urls"] = [openai_base_url]

obj["openai"] = openai
new_data = json.dumps(obj, ensure_ascii=False, separators=(",", ":"))
cur.execute("UPDATE config SET data = ?, version = ? WHERE id = 1", (new_data, version))
con.commit()
print("patched: api_keys=1, api_base_urls=%s" % (1 if openai_base_url else "unchanged"))
PY
    print_success "OpenWebUI config updated (keys + base_urls synced)"

    print_step "Restarting OpenWebUI container"
    COMPOSE_FILE="${COMPOSE_FILE:-docker-compose.yml}" docker_compose up -d --force-recreate openwebui >/dev/null
    print_success "OpenWebUI restarted"

    print_step "Done"
    print_success "OpenWebUI OpenAI config is synced with .env"
}

main "$@"
