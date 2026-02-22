#!/bin/bash
###############################################################################
# OpenWebUI user recovery helper (role + activation + optional password reset)
#
# Why this exists:
# - If WEBUI_AUTH is enabled, /api/models requires auth; if you're stuck in a
#   broken session you'll see "no models" in UI because the browser gets 401.
# - Sometimes an account can end up in a "pending" role. This script promotes
#   a user to "admin" (or "user") and ensures auth.active=1.
#
# Safety:
# - This script edits OpenWebUI's SQLite DB in-place. It always creates a backup
#   into ./backups first.
#
# Usage:
#   ./openwebui-user-admin.sh --email you@example.com
#   ./openwebui-user-admin.sh --email you@example.com --role user
#   ./openwebui-user-admin.sh --email you@example.com --reset-password
#   ./openwebui-user-admin.sh --email you@example.com --reset-password --no-restart
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults
cd "$SCRIPT_DIR"

usage() {
    cat <<'EOF'
Usage:
  ./openwebui-user-admin.sh --email you@example.com [--role admin|user] [--reset-password] [--no-restart]

Notes:
  - Always clears "pending" by setting role explicitly.
  - Ensures auth.active=1.
  - If --reset-password is set, you will be prompted for a new password (hidden input).
EOF
}

EMAIL=""
ROLE="admin"
RESET_PASSWORD="false"
NO_RESTART="false"

while [ $# -gt 0 ]; do
    case "$1" in
        --email)
            EMAIL="${2:-}"
            shift 2
            ;;
        --role)
            ROLE="${2:-}"
            shift 2
            ;;
        --reset-password)
            RESET_PASSWORD="true"
            shift 1
            ;;
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

if [ -z "$EMAIL" ]; then
    print_error "--email is required"
    usage
    exit 2
fi

if [ "$ROLE" != "admin" ] && [ "$ROLE" != "user" ]; then
    print_error "--role must be 'admin' or 'user' (got: $ROLE)"
    exit 2
fi

print_header "OpenWebUI User Recovery"

if ! docker info >/dev/null 2>&1; then
    print_error "Docker daemon is not running."
    exit 1
fi

if ! docker ps --format '{{.Names}}' | grep -qx 'openwebui'; then
    print_error "Container 'openwebui' is not running."
    print_error "Start it first: docker compose up -d openwebui"
    exit 1
fi

print_step "Backing up webui.db"
./backup-openwebui-db.sh >/dev/null
print_success "DB backup written under ./backups"

NEW_PASSWORD=""
if [ "$RESET_PASSWORD" = "true" ]; then
    print_step "Reading new password (hidden input)"
    read -r -s -p "New password: " NEW_PASSWORD
    echo
    read -r -s -p "Confirm password: " NEW_PASSWORD_CONFIRM
    echo
    if [ -z "$NEW_PASSWORD" ]; then
        print_error "Password cannot be empty."
        exit 2
    fi
    if [ "$NEW_PASSWORD" != "$NEW_PASSWORD_CONFIRM" ]; then
        print_error "Password confirmation does not match."
        exit 2
    fi
fi

print_step "Updating user role/activation in SQLite (inside container)"
if [ "$RESET_PASSWORD" = "true" ]; then
    # Pass password via stdin and Python code via -c to avoid stdin collisions.
    RESET_PASSWORD_PYTHON="$(
        cat <<'PY'
from __future__ import annotations

import os
import sqlite3
import sys
import time

import bcrypt

DB = "/app/backend/data/webui.db"
email = os.environ.get("OPENWEBUI_TARGET_EMAIL", "").strip()
role = os.environ.get("OPENWEBUI_TARGET_ROLE", "admin").strip()
pw = sys.stdin.read()

if not email:
    raise SystemExit("missing OPENWEBUI_TARGET_EMAIL")
if role not in {"admin", "user"}:
    raise SystemExit(f"invalid role: {role}")

pw = pw.rstrip("\n")
if not pw:
    raise SystemExit("empty password provided")

now = int(time.time())
hashed = bcrypt.hashpw(pw.encode("utf-8"), bcrypt.gensalt()).decode("utf-8")

con = sqlite3.connect(DB)
cur = con.cursor()

cur.execute("SELECT id FROM auth WHERE email = ?", (email,))
row = cur.fetchone()
if not row:
    raise SystemExit("no auth row for this email")
user_id = row[0]

cur.execute("UPDATE auth SET active = 1, password = ? WHERE id = ?", (hashed, user_id))
cur.execute("UPDATE user SET role = ?, updated_at = ? WHERE id = ?", (role, now, user_id))
con.commit()

cur.execute("SELECT id, name, email, role FROM user WHERE id = ?", (user_id,))
print("updated_user", cur.fetchone())
PY
    )"
    printf "%s" "$NEW_PASSWORD" | docker exec -i \
        -e "OPENWEBUI_TARGET_EMAIL=$EMAIL" \
        -e "OPENWEBUI_TARGET_ROLE=$ROLE" \
        openwebui python3 -c "$RESET_PASSWORD_PYTHON"
else
    docker exec -i \
        -e "OPENWEBUI_TARGET_EMAIL=$EMAIL" \
        -e "OPENWEBUI_TARGET_ROLE=$ROLE" \
        openwebui python3 - <<'PY'
from __future__ import annotations

import os
import sqlite3
import time

DB = "/app/backend/data/webui.db"
email = os.environ.get("OPENWEBUI_TARGET_EMAIL", "").strip()
role = os.environ.get("OPENWEBUI_TARGET_ROLE", "admin").strip()

if not email:
    raise SystemExit("missing OPENWEBUI_TARGET_EMAIL")
if role not in {"admin", "user"}:
    raise SystemExit(f"invalid role: {role}")

now = int(time.time())

con = sqlite3.connect(DB)
cur = con.cursor()

cur.execute("SELECT id FROM auth WHERE email = ?", (email,))
row = cur.fetchone()
if not row:
    raise SystemExit("no auth row for this email")
user_id = row[0]

cur.execute("UPDATE auth SET active = 1 WHERE id = ?", (user_id,))
cur.execute("UPDATE user SET role = ?, updated_at = ? WHERE id = ?", (role, now, user_id))
con.commit()

cur.execute("SELECT id, name, email, role FROM user WHERE id = ?", (user_id,))
print("updated_user", cur.fetchone())
PY
fi

print_success "User updated"
unset NEW_PASSWORD NEW_PASSWORD_CONFIRM

if [ "$NO_RESTART" = "true" ]; then
    print_warning "Skipping restart (--no-restart)."
    exit 0
fi

print_step "Restarting OpenWebUI container"
docker_compose restart openwebui >/dev/null
print_success "OpenWebUI restarted"

print_step "Next steps (browser)"
OPENWEBUI_PORT="${WEBUI_PORT:-3000}"
echo "  1. If the UI still shows no models: clear site data for http://127.0.0.1:${OPENWEBUI_PORT}"
echo "  2. Sign in again and reload."
