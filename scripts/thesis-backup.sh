#!/bin/bash
###############################################################################
# Thesis Backup Script - Daily Automated Backup for Academic Research
#
# Usage:
#   1. Save as ~/thesis-backup.sh
#   2. Make executable: chmod +x ~/thesis-backup.sh
#   3. Add to crontab: crontab -e
#      0 2 * * * /home/YOUR_USERNAME/thesis-backup.sh
#
# This script backs up:
#   - OpenWebUI database (chats, settings, knowledge base)
#   - Extracted paper data
#   - Configuration files
#
# Backups are stored in multiple locations for redundancy.
###############################################################################

set -euo pipefail

# Configuration
DATE=$(date +%Y%m%d_%H%M%S)
THESIS_BACKUP_DIR="${HOME}/thesis-backups"
OPENWEBUI_DIR="/LAB/@thesis/openwebui"
RETENTION_DAYS=30

# Ensure backup directory exists
mkdir -p "${THESIS_BACKUP_DIR}"

echo "=== Thesis Backup Started: ${DATE} ==="

# 1. OpenWebUI Database Backup
echo "[1/4] Backing up OpenWebUI database..."
cd "${OPENWEBUI_DIR}"
if [[ -x ./backup-openwebui-db.sh ]]; then
    ./backup-openwebui-db.sh
    # Copy latest backup to thesis directory
    LATEST_BACKUP=$(ls -t backups/webui.db.backup_* 2>/dev/null | head -1)
    if [[ -n ${LATEST_BACKUP} ]]; then
        cp "${LATEST_BACKUP}" "${THESIS_BACKUP_DIR}/"
        echo "  ✓ Database backed up: $(basename ${LATEST_BACKUP})"
    fi
else
    echo "  ✗ backup-openwebui-db.sh not found or not executable"
    exit 1
fi

# 2. Knowledge Base Export
echo "[2/4] Exporting knowledge base metadata..."
if [[ -d "${OPENWEBUI_DIR}/data/corpus" ]]; then
    tar czf "${THESIS_BACKUP_DIR}/corpus-${DATE}.tar.gz" \
        -C "${OPENWEBUI_DIR}/data" \
        corpus/ 2>/dev/null || echo "  ! Corpus backup skipped (may not exist)"
fi

# 3. Configuration Backup
echo "[3/4] Backing up configuration..."
tar czf "${THESIS_BACKUP_DIR}/config-${DATE}.tar.gz" \
    -C "${OPENWEBUP_DIR}" \
    .env cliproxyapi/config.yaml docker-compose.yml 2>/dev/null ||
    echo "  ! Some config files may be missing"

# 4. Cleanup old backups
echo "[4/4] Cleaning up backups older than ${RETENTION_DAYS} days..."
find "${THESIS_BACKUP_DIR}" -name "*.backup_*" -mtime +${RETENTION_DAYS} -delete 2>/dev/null || true
find "${THESIS_BACKUP_DIR}" -name "corpus-*.tar.gz" -mtime +${RETENTION_DAYS} -delete 2>/dev/null || true
find "${THESIS_BACKUP_DIR}" -name "config-*.tar.gz" -mtime +${RETENTION_DAYS} -delete 2>/dev/null || true

# Summary
echo ""
echo "=== Backup Complete ==="
echo "Backup location: ${THESIS_BACKUP_DIR}"
echo "Recent backups:"
ls -lh "${THESIS_BACKUP_DIR}" | tail -5
echo ""
echo "Next steps:"
echo "  1. Copy backups to external storage (USB, cloud)"
echo "  2. Verify backup integrity: tar tzf ${THESIS_BACKUP_DIR}/corpus-${DATE}.tar.gz"
echo "  3. For restore: ./restore-openwebui-db.sh ${THESIS_BACKUP_DIR}/webui.db.backup_YYYYMMDD_HHMMSS"
