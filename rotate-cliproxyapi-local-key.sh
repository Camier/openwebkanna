#!/bin/bash
#
# Rotate the *local* API key used between OpenWebUI <-> CLIProxyAPI.
#
# Why:
# - The default placeholder key (e.g. "layra-cliproxyapi-key") is intentionally flagged as weak.
# - A strong local key avoids accidental exposure if the port binding is widened later.
#
# What this does:
# - Generates a new random key (never printed).
# - Updates `.env`: OPENAI_API_KEY + CLIPROXYAPI_API_KEY
# - Updates the mounted CLIProxyAPI config (`CLIPROXYAPI_CONFIG`, default: cliproxyapi/config.local.yaml)
# - Restarts `cliproxyapi` and `openwebui` containers so changes take effect.
#
# Safety:
# - Does NOT print the key.
# - Only edits local, gitignored files (`.env`, `cliproxyapi/config.local.yaml` by default).
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

COMPOSE_FILE="${COMPOSE_FILE:-${SCRIPT_DIR}/docker-compose.yml}"
ENV_FILE="${ENV_FILE:-${SCRIPT_DIR}/.env}"

main() {
    print_header "Rotate CLIProxyAPI Local Key"

    if [ ! -f "$ENV_FILE" ]; then
        print_error "Missing .env file: $ENV_FILE"
        print_error "Create it first (copy from .env.example) and set required variables."
        exit 1
    fi

    local config_path="${CLIPROXYAPI_CONFIG:-./cliproxyapi/config.local.yaml}"
    config_path="$(resolve_path "$config_path")"

    if [ ! -f "$config_path" ]; then
        local template="${SCRIPT_DIR}/cliproxyapi/config.yaml"
        if [ ! -f "$template" ]; then
            print_error "Missing CLIProxyAPI config template: $template"
            exit 1
        fi

        print_step "Creating local CLIProxyAPI config"
        cp "$template" "$config_path"
        chmod 600 "$config_path" 2>/dev/null || true
        print_success "Created: ${config_path}"
    fi

    if ! command_exists openssl; then
        print_error "openssl is required to generate a strong key"
        exit 1
    fi

    local new_key=""
    new_key="$(openssl rand -hex 32)"

    print_step "Updating local secret values (key not displayed)"
    NEW_KEY="$new_key" ENV_FILE="$ENV_FILE" CONFIG_FILE="$config_path" python3 - <<'PY'
import os
import re
import sys
from pathlib import Path

new_key = os.environ.get("NEW_KEY", "")
env_file = Path(os.environ["ENV_FILE"])
config_file = Path(os.environ["CONFIG_FILE"])

if not new_key:
    sys.exit("NEW_KEY is empty")

def upsert_env(path: Path, key: str, value: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines(True)
    found = False
    out = []
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

upsert_env(env_file, "OPENAI_API_KEY", new_key)
upsert_env(env_file, "CLIPROXYAPI_API_KEY", new_key)
upsert_env(env_file, "OPENAI_API_KEYS", new_key)

text = config_file.read_text(encoding="utf-8")
pattern = r"^api-keys:\n(?:[ \t]*-[^\n]*\n)+"
replacement = f"api-keys:\n  - {new_key}\n"
new_text, n = re.subn(pattern, replacement, text, flags=re.M)
if n == 0:
    sys.exit("Could not find api-keys block in CLIProxyAPI config; aborting without changes.")
config_file.write_text(new_text, encoding="utf-8")
PY
    print_success "Updated .env and CLIProxyAPI config"

    # OpenWebUI also stores some OpenAI connection values in its persistent config DB.
    # If we rotate the local key, keep that DB config in sync too so /api/models and
    # /api/chat/completions continue to work without manual intervention.
    if command_exists docker && docker ps --format '{{.Names}}' | grep -Fxq "openwebui"; then
        print_step "Syncing OpenWebUI persistent OpenAI config (key not displayed)"
        docker exec -i -e "OPENAI_KEY=$new_key" openwebui python3 - <<'PY'
from __future__ import annotations

import json
import os
import sqlite3

DB = "/app/backend/data/webui.db"
openai_key = os.environ.get("OPENAI_KEY", "").strip()
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
obj["openai"] = openai

new_data = json.dumps(obj, ensure_ascii=False, separators=(",", ":"))
cur.execute("UPDATE config SET data = ?, version = ? WHERE id = 1", (new_data, version))
con.commit()
print("patched_openai_keys_count=1")
PY
        print_success "OpenWebUI persistent OpenAI config synced"
    else
        print_warning "OpenWebUI container not running; skipping persistent OpenAI config sync"
    fi

    print_step "Restarting containers"
    # Restart cliproxyapi first so OpenWebUI reconnects with the new key.
    docker_compose up -d --force-recreate cliproxyapi >/dev/null
    docker_compose up -d --force-recreate openwebui >/dev/null
    print_success "Containers restarted"

    if [ -x "${SCRIPT_DIR}/check-cliproxyapi.sh" ]; then
        print_step "Validating CLIProxyAPI health"
        if CLIPROXYAPI_CHECK_CHAT_COMPLETION=false "${SCRIPT_DIR}/check-cliproxyapi.sh" --quiet >/dev/null 2>&1; then
            print_success "CLIProxyAPI health check passed"
        else
            print_warning "CLIProxyAPI health check failed; inspect logs with: ./logs.sh"
        fi
    fi

    print_step "Done"
    print_success "Local key rotated successfully"
}

main "$@"
