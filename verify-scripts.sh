#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Script quality verification for shell entrypoints.
# 1) Syntax validation with bash -n
# 2) Lint with shellcheck (if installed)
# 3) Targeted smoke test for update.sh
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"

cd "$SCRIPT_DIR"

FAILED=0
STRICT_ALL="${VERIFY_ALL_STRICT:-0}"
FAILED_TARGET=0

mapfile -t SCRIPTS < <(rg --files -g '*.sh' "$SCRIPT_DIR" | sort)
mapfile -t TARGET_SCRIPTS < <(printf "%s\n" \
    "$SCRIPT_DIR/update.sh" \
    "$SCRIPT_DIR/test-update-smoke.sh" \
    "$SCRIPT_DIR/audit-openwebui-plugins.sh" \
    "$SCRIPT_DIR/verify-scripts.sh")

is_target_script() {
    local candidate=$1
    local target
    for target in "${TARGET_SCRIPTS[@]}"; do
        if [ "$candidate" = "$target" ]; then
            return 0
        fi
    done
    return 1
}

if [ ${#SCRIPTS[@]} -eq 0 ]; then
    print_info "No shell scripts found"
    exit 0
fi

print_info "Running bash syntax checks..."
for script in "${SCRIPTS[@]}"; do
    if ! bash -n "$script"; then
        print_error "bash -n failed: $script"
        FAILED=1
    else
        print_success "bash -n: $script"
    fi
done

if command -v shellcheck >/dev/null 2>&1; then
    print_info "Running shellcheck..."
    for script in "${SCRIPTS[@]}"; do
        if ! is_target_script "$script"; then
            SCOPED_FLAG="(informational)"
        else
            SCOPED_FLAG="(critical)"
        fi
        if ! shellcheck -x -s bash "$script"; then
            print_error "shellcheck failed: $script"
            FAILED=1
            if is_target_script "$script"; then
                FAILED_TARGET=1
            fi
            if ((STRICT_ALL == 1)); then
                FAILED_TARGET=1
            fi
        else
            print_success "shellcheck $SCOPED_FLAG: $script"
        fi
    done
else
    print_warning "shellcheck not installed; skipping"
    # Keep as a warning so verify stays informative without failing in environments
    # where shellcheck is unavailable.
fi

print_info "Running update.sh smoke tests..."
if ! "$SCRIPT_DIR/test-update-smoke.sh"; then
    print_error "update.sh smoke test failed"
    FAILED=1
fi

if ((STRICT_ALL == 1)); then
    print_info "STRICT mode enabled; shellcheck failures block verification"
fi

if ((FAILED == 0)); then
    print_success "verify-scripts completed successfully"
    exit 0
fi

if ((FAILED_TARGET != 0)); then
    print_error "verify-scripts found issues in critical scripts"
    exit 1
fi

if ((FAILED != 0)); then
    print_warning "verify-scripts found shellcheck issues outside critical paths"
    if ((STRICT_ALL == 1)); then
        print_error "strict mode requested; treat as failure"
        exit 1
    fi
fi

print_success "verify-scripts completed successfully"
exit 0
