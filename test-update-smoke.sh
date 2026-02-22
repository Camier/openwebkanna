#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Lightweight smoke tests for update.sh argument and restore flows.
# Uses command stubs so no real Docker/Compose calls are made.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"

TMP_DIR="$(mktemp -d)"
PROJECT_DIR="$TMP_DIR/project"
FAKEBIN_DIR="$TMP_DIR/fakebin"

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

mkdir -p "$PROJECT_DIR" "$FAKEBIN_DIR"

log() {
    echo "[$1] $2"
}

run_case() {
    local name=$1
    local should_succeed=$2
    shift 2
    local logfile="$TMP_DIR/${name}.log"

    set +e
    "$@" >"$logfile" 2>&1
    local status=$?
    set -e

    if ((should_succeed == 0)) && ((status != 0)); then
        log "FAIL" "$name (exit=$status)"
        cat "$logfile"
        exit 1
    fi
    if ((should_succeed == 1)) && ((status == 0)); then
        log "FAIL" "$name unexpectedly succeeded"
        cat "$logfile"
        exit 1
    fi

}

assert_contains() {
    local needle=$1
    local file=$2

    if ! grep -Fq "$needle" "$file"; then
        log "FAIL" "Expected output missing: $needle"
        cat "$file"
        exit 1
    fi
}

# Stub docker so restore mode can run without a real Docker daemon.
cat >"$FAKEBIN_DIR/docker" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

if [[ "${1-}" == "run" ]]; then
    # Simulate successful container execution.
    exit 0
fi

exit 0
EOF
chmod +x "$FAKEBIN_DIR/docker"

# Stub docker-compose to prevent touching real services.
cat >"$FAKEBIN_DIR/docker-compose" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

exit 0
EOF
chmod +x "$FAKEBIN_DIR/docker-compose"

# Use a fixture copy of update.sh and lib/ to keep test artifacts isolated from repo data.
cp "$SCRIPT_DIR/update.sh" "$PROJECT_DIR/update.sh"
chmod +x "$PROJECT_DIR/update.sh"
cp -r "$SCRIPT_DIR/lib" "$PROJECT_DIR/lib"

cd "$PROJECT_DIR"
export PATH="$FAKEBIN_DIR:$PATH"

run_case "restore_without_backup_arg_fails" 1 ./update.sh --restore
assert_contains "Please specify a backup file to restore" "${TMP_DIR}/restore_without_backup_arg_fails.log"

run_case "list_backups_without_files" 0 ./update.sh --list-backups
assert_contains "No backups found" "${TMP_DIR}/list_backups_without_files.log"

mkdir -p "$PROJECT_DIR/backups"
touch "$PROJECT_DIR/backups/openwebui_backup_20260101_120000.tar.gz"

run_case "restore_valid_backup_succeeds" 0 ./update.sh --force --restore "$PROJECT_DIR/backups/openwebui_backup_20260101_120000.tar.gz"
assert_contains "Backup restored successfully" "${TMP_DIR}/restore_valid_backup_succeeds.log"

log "PASS" "All smoke tests passed"
