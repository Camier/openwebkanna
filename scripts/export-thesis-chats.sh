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
# Output: Markdown files in thesis-exports/ directory
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXPORT_DIR="${PROJECT_ROOT}/thesis-exports"

# Load environment
source "${PROJECT_ROOT}/lib/init.sh" 2>/dev/null || {
    echo "Warning: Could not load lib/init.sh"
}

DB_PATH="${DB_PATH:-/app/backend/data/webui.db}"

check_container() {
    if ! docker compose ps | grep -q "openwebui.*running"; then
        echo "Error: OpenWebUI container is not running"
        echo "Start with: docker compose up -d"
        exit 1
    fi
}

list_chats() {
    echo "=== Available Chats ==="
    echo ""
    docker compose exec -T openwebui sqlite3 "${DB_PATH}" <<'EOF'
SELECT
    c.id,
    datetime(c.created_at/1000000000, 'unixepoch') as created,
    c.title,
    COUNT(m.id) as messages
FROM chat c
LEFT JOIN chat_message m ON c.id = m.chat_id
WHERE c.user_id IS NOT NULL
GROUP BY c.id
ORDER BY c.created_at DESC
LIMIT 50;
EOF
    echo ""
    echo "To export a specific chat, run:"
    echo "  ./export-thesis-chats.sh <chat_id>"
}

export_chat() {
    local chat_id="$1"
    local output_file="${EXPORT_DIR}/chat-${chat_id}-$(date +%Y%m%d).md"

    mkdir -p "${EXPORT_DIR}"

    echo "Exporting chat ${chat_id}..."

    # Export chat metadata and messages
    docker compose exec -T openwebui sqlite3 "${DB_PATH}" <<EOF | tee "${output_file}"
.mode markdown
.header on
.print # Chat Export
.print
.print **Chat ID:** ${chat_id}
.print **Exported:** $(date -Iseconds)
.print
---

SELECT
    datetime(m.created_at/1000000000, 'unixepoch') as timestamp,
    CASE m.role
        WHEN 'user' THEN '👤 User'
        WHEN 'assistant' THEN '🤖 Assistant'
        ELSE m.role
    END as role,
    m.content
FROM chat_message m
WHERE m.chat_id = '${chat_id}'
ORDER BY m.created_at ASC;
EOF

    echo ""
    echo "✓ Exported to: ${output_file}"
    echo ""
    echo "Citation format for thesis:"
    echo "  OpenWebUI conversation (Chat ID: ${chat_id}), exported $(date +%Y-%m-%d)."
}

export_all_chats() {
    mkdir -p "${EXPORT_DIR}"

    echo "Exporting all chats..."

    # Get all chat IDs
    chat_ids=$(docker compose exec -T openwebui sqlite3 "${DB_PATH}" \
        "SELECT id FROM chat WHERE user_id IS NOT NULL ORDER BY created_at DESC;")

    for chat_id in ${chat_ids}; do
        export_chat "${chat_id}"
        echo "---"
    done

    echo ""
    echo "✓ All chats exported to: ${EXPORT_DIR}/"
}

main() {
    check_container

    case "${1:-}" in
        "")
            list_chats
            ;;
        "--all" | "-a")
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
EOF
            ;;
        *)
            export_chat "$1"
            ;;
    esac
}

main "$@"
