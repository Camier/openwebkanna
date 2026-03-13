#!/bin/bash

###############################################################################
# Sync root compatibility copies from canonical config/ sources.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# shellcheck disable=SC1091
source "${PROJECT_ROOT}/lib/init.sh"

cd "${PROJECT_ROOT}"

MANIFEST="config/compatibility-copies.txt"
MODE="apply"

show_help() {
    cat <<'EOF'
Usage: ./scripts/sync-compatibility-copies.sh [OPTIONS]

Options:
  --dry-run            Show planned copy actions without writing files.
  -h, --help           Show this help.

Behavior:
  Copies each root compatibility file from its canonical config/ source as
  declared in config/compatibility-copies.txt.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            MODE="dry-run"
            shift
            ;;
        -h | --help)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown argument: $1"
            show_help
            exit 1
            ;;
    esac
done

if [[ ! -f ${MANIFEST} ]]; then
    print_error "Missing manifest: ${MANIFEST}"
    exit 1
fi

print_header "Sync Compatibility Copies"
print_info "Mode: ${MODE}"

synced=0
unchanged=0

while IFS= read -r entry || [[ -n ${entry} ]]; do
    entry="${entry%$'\r'}"
    entry="${entry#"${entry%%[![:space:]]*}"}"
    entry="${entry%"${entry##*[![:space:]]}"}"

    [[ -z ${entry} ]] && continue
    [[ ${entry} == \#* ]] && continue

    if [[ ${entry} != *:* ]]; then
        print_error "Invalid manifest entry: ${entry}"
        exit 1
    fi

    root_path="${entry%%:*}"
    canonical_path="${entry#*:}"

    if [[ ! -f ${canonical_path} ]]; then
        print_error "Missing canonical file: ${canonical_path}"
        exit 1
    fi

    mkdir -p "$(dirname "${root_path}")"

    if [[ -f ${root_path} ]] && cmp -s "${root_path}" "${canonical_path}"; then
        print_success "Aligned: ${root_path}"
        unchanged=$((unchanged + 1))
        continue
    fi

    if [[ ${MODE} == "dry-run" ]]; then
        print_info "Would sync ${root_path} <- ${canonical_path}"
    else
        cp "${canonical_path}" "${root_path}"
        print_success "Synced ${root_path} <- ${canonical_path}"
    fi

    synced=$((synced + 1))
done <"${MANIFEST}"

print_section "Summary"
print_info "Synced: ${synced}"
print_info "Already aligned: ${unchanged}"
