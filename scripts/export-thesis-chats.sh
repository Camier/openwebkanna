#!/bin/bash
###############################################################################
# Thesis Chat Export Script - Export Conversations for Academic Documentation
#
# Usage:
#   ./export-thesis-chats.sh [chat_id]
#
#   Without chat_id: Lists all chats with IDs
#   With chat_id: Exports specific chat to markdown
#
# Examples:
#   ./export-thesis-chats.sh                    # List all chats
#   ./export-thesis-chats.sh abc123             # Export chat abc123
#   ./export-thesis-chats.sh --all              # Export all chats
#
# Output: Markdown files in thesis-exports/ or logs/thesis-exports/
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXPORT_DIR="${PROJECT_ROOT}/thesis-exports"
FALLBACK_EXPORT_DIR="${PROJECT_ROOT}/logs/thesis-exports"
EFFECTIVE_EXPORT_DIR=""

# Load environment
source "${PROJECT_ROOT}/lib/init.sh" 2>/dev/null || {
    echo "Warning: Could not load lib/init.sh"
}

DB_PATH="${DB_PATH:-/app/backend/data/webui.db}"

resolve_export_dir() {
    if [[ -n ${EFFECTIVE_EXPORT_DIR} ]]; then
        return 0
    fi

    if mkdir -p "${EXPORT_DIR}" 2>/dev/null && [[ -w ${EXPORT_DIR} ]]; then
        EFFECTIVE_EXPORT_DIR="${EXPORT_DIR}"
        return 0
    fi

    mkdir -p "${FALLBACK_EXPORT_DIR}"
    EFFECTIVE_EXPORT_DIR="${FALLBACK_EXPORT_DIR}"
    echo "Warning: ${EXPORT_DIR} is not writable; exporting to ${EFFECTIVE_EXPORT_DIR}"
}

check_container() {
    if ! docker compose ps --services --status running | grep -qx "openwebui"; then
        echo "Error: OpenWebUI container is not running"
        echo "Start with: docker compose up -d"
        exit 1
    fi
}

list_chats() {
    echo "=== Available Chats ==="
    echo ""
    docker compose exec -T -e DB_PATH="${DB_PATH}" openwebui python3 - <<'PY'
from datetime import datetime, timezone
import os
import sqlite3

def normalize_timestamp(value):
    if value in (None, ""):
        return ""
    ts = int(value)
    if ts > 10**18:
        ts = ts / 1_000_000_000
    elif ts > 10**15:
        ts = ts / 1_000_000
    elif ts > 10**12:
        ts = ts / 1_000
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

db_path = os.environ["DB_PATH"]
conn = sqlite3.connect(db_path)
rows = conn.execute(
    """
    SELECT
        c.id,
        c.created_at,
        COALESCE(c.title, '(untitled)') AS title,
        COUNT(m.id) AS messages
    FROM chat c
    LEFT JOIN chat_message m ON c.id = m.chat_id
    WHERE c.user_id IS NOT NULL
    GROUP BY c.id
    ORDER BY c.created_at DESC
    LIMIT 50
    """
).fetchall()

if not rows:
    print("No chats found.")
else:
    print("CHAT ID\tCREATED\tMESSAGES\tTITLE")
    for chat_id, created, title, messages in rows:
        print(f"{chat_id}\t{normalize_timestamp(created)}\t{messages}\t{title}")
PY
    echo ""
    echo "To export a specific chat, run:"
    echo "  ./export-thesis-chats.sh <chat_id>"
}

export_chat() {
    local chat_id="$1"
    local output_file
    resolve_export_dir
    output_file="${EFFECTIVE_EXPORT_DIR}/chat-${chat_id}-$(date +%Y%m%d).md"

    echo "Exporting chat ${chat_id}..."

    docker compose exec -T \
        -e DB_PATH="${DB_PATH}" \
        -e CHAT_ID="${chat_id}" \
        openwebui python3 - >"${output_file}" <<'PY'
import os
import sqlite3
from datetime import datetime, timezone

def normalize_timestamp(value):
    if value in (None, ""):
        return ""
    ts = int(value)
    if ts > 10**18:
        ts = ts / 1_000_000_000
    elif ts > 10**15:
        ts = ts / 1_000_000
    elif ts > 10**12:
        ts = ts / 1_000
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

db_path = os.environ["DB_PATH"]
chat_id = os.environ["CHAT_ID"]
conn = sqlite3.connect(db_path)

chat = conn.execute(
    """
    SELECT
        COALESCE(title, '(untitled)') AS title,
        created_at
    FROM chat
    WHERE id = ?
    """,
    (chat_id,),
).fetchone()

if chat is None:
    raise SystemExit(f"Error: Chat {chat_id} not found")

messages = conn.execute(
    """
    SELECT
        created_at,
        COALESCE(role, 'unknown') AS role,
        COALESCE(content, '') AS content
    FROM chat_message
    WHERE chat_id = ?
    ORDER BY created_at ASC
    """,
    (chat_id,),
).fetchall()

title, created = chat
print("# Chat Export")
print()
print(f"**Chat ID:** {chat_id}")
print(f"**Title:** {title}")
print(f"**Created:** {normalize_timestamp(created)}")
print(f"**Exported:** {datetime.now(timezone.utc).isoformat()}")
print()
print("---")
print()

role_labels = {
    "user": "User",
    "assistant": "Assistant",
}

for timestamp, role, content in messages:
    label = role_labels.get(role, role)
    print(f"## {normalize_timestamp(timestamp)} | {label}")
    print()
    print(content.rstrip())
    print()
    print("---")
    print()
PY

    cat "${output_file}"

    echo ""
    echo "✓ Exported to: ${output_file}"
    echo ""
    echo "Citation format for thesis:"
    echo "  OpenWebUI conversation (Chat ID: ${chat_id}), exported $(date +%Y-%m-%d)."
}

export_all_chats() {
    local chat_ids=()
    local chat_id=""

    resolve_export_dir

    echo "Exporting all chats..."

    mapfile -t chat_ids < <(
        docker compose exec -T -e DB_PATH="${DB_PATH}" openwebui python3 - <<'PY'
import os
import sqlite3

db_path = os.environ["DB_PATH"]
conn = sqlite3.connect(db_path)
for (chat_id,) in conn.execute(
    "SELECT id FROM chat WHERE user_id IS NOT NULL ORDER BY created_at DESC"
):
    print(chat_id)
PY
    )

    for chat_id in "${chat_ids[@]}"; do
        export_chat "${chat_id}"
        echo "---"
    done

    echo ""
    echo "✓ All chats exported to: ${EFFECTIVE_EXPORT_DIR}/"
}

main() {
    case "${1:-}" in
        "")
            check_container
            list_chats
            ;;
        "--all" | "-a")
            check_container
            export_all_chats
            ;;
        "--help" | "-h")
            cat <<'EOF'
Thesis Chat Export Script

Usage:
  ./export-thesis-chats.sh [OPTIONS] [CHAT_ID]

Options:
  --all, -a       Export all chats
  --help, -h      Show this help

Examples:
  ./export-thesis-chats.sh              # List all chats
  ./export-thesis-chats.sh abc123       # Export chat abc123
  ./export-thesis-chats.sh --all        # Export all chats

Output:
  Markdown files in thesis-exports/ directory
  Falls back to logs/thesis-exports/ if thesis-exports/ is not writable
EOF
            ;;
        *)
            check_container
            export_chat "$1"
            ;;
    esac
}

main "$@"
