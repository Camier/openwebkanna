#!/bin/bash

###############################################################################
# Real Integration Guard
# Fails if runtime/config files include mock-like implementation markers.
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"

cd "$SCRIPT_DIR"

main() {
    local patterns
    local targets
    local results

    if ! command -v rg >/dev/null 2>&1; then
        print_error "ripgrep (rg) is required for audit-no-mock.sh"
        exit 1
    fi

    print_info "Running real-integration audit (runtime + config files)"

    patterns='\\b(mock|dummy|fake|stub|test[[:space:]_-]*double|monkey[[:space:]_-]*patch|simulated?[[:space:]_-]*(provider|backend|response))\\b'

    targets="$(git ls-files '*.sh' '*.yml' '*.yaml' '*.toml' '*.json' 2>/dev/null || true)"
    if [ -z "$targets" ]; then
        print_warning "No tracked runtime/config files found to audit"
        exit 0
    fi

    # Exclude this script from keyword matching to avoid self-reference noise.
    results="$(printf '%s\n' "$targets" | rg -v '^audit-no-mock\.sh$' | xargs -r rg -n -i -e "$patterns" || true)"

    if [ -n "$results" ]; then
        print_error "Mock-like markers detected in tracked runtime/config files:"
        echo "$results"
        exit 1
    fi

    print_success "No mock-like implementation markers detected"
}

main "$@"
