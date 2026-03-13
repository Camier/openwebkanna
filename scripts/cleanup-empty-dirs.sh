#!/bin/bash

###############################################################################
# Empty Directory Cleanup
# Identifies empty directories while excluding runtime/system directories.
# These are typically runtime artifacts that should not be tracked in git.
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

# shellcheck disable=SC1091
source "${SCRIPT_DIR}/../lib/init.sh"
load_env_defaults

DRY_RUN=false
SHOW_HELP=false

show_help() {
    cat <<'EOF'
Usage: ./scripts/cleanup-empty-dirs.sh [OPTIONS]

Options:
  --dry-run    Show what would be removed without deleting
  -h, --help   Show this help message

Description:
  Finds and optionally removes truly empty directories.

  Excluded from cleanup:
  - .git/ (git metadata)
  - .venvs/, .conda/ (local environment directories)
  - archive/ (deprecated but kept for reference)
  - data/milvus/ (runtime state)
  - data/processing/prod_max/ and data/processing/prod_max_multimodal/ (processing outputs)
  - logs/, backups/, certs/ (runtime directories)
  - Any directory containing .gitkeep files

  Targeted cleanup inside processing outputs:
  - Empty per-paper `data/processing/prod_max/*/images/` directories
  - Empty per-paper `data/processing/prod_max_multimodal/*/images/` directories

Examples:
  # Show empty directories without removing
  ./scripts/cleanup-empty-dirs.sh --dry-run

  # Remove empty directories
  ./scripts/cleanup-empty-dirs.sh
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h | --help)
            SHOW_HELP=true
            shift
            ;;
        *)
            print_error "Unknown argument: $1"
            SHOW_HELP=true
            shift 1
            ;;
    esac
done

if [ "$SHOW_HELP" = true ]; then
    show_help
    exit 0
fi

print_header "Empty Directory Cleanup"

# Directories to always exclude
export EXCLUDE_PATTERN="\.git|\.venvs|\.conda|archive|data/milvus|logs|backups|certs|data/pdfs|data/extractions|data/corpus|data/processing|jupyter|cliproxyapi|searxng|thesis-exports"

print_step "Finding empty directories..."
print_info "Excluded: .git, .venvs/.conda, archive, data/milvus, data/processing, logs, backups, certs, runtime dirs"
print_info ""

# Find empty directories, excluding special dirs. Processing outputs stay excluded
# except for per-paper empty images/ placeholders, which are safe generated residue.
general_empty_dirs=$(find . -type d -empty \
    ! -path "./.git*" \
    ! -path "./.venvs*" \
    ! -path "./.conda*" \
    ! -path "./archive*" \
    ! -path "./data/milvus*" \
    ! -path "./logs*" \
    ! -path "./backups*" \
    ! -path "./certs*" \
    ! -path "./data/pdfs*" \
    ! -path "./data/extractions*" \
    ! -path "./data/corpus*" \
    ! -path "./data/processing*" \
    ! -path "./jupyter*" \
    ! -path "./cliproxyapi*" \
    ! -path "./searxng*" \
    ! -path "./thesis-exports*" \
    2>/dev/null | sort)

targeted_processing_empty_dirs=$(find ./data/processing \
    \( -path "./data/processing/prod_max/*/images" -o -path "./data/processing/prod_max_multimodal/*/images" \) \
    -type d -empty 2>/dev/null | sort)

empty_dirs=$(printf '%s\n%s\n' "${general_empty_dirs}" "${targeted_processing_empty_dirs}" | sed '/^$/d' | sort -u)

if [ -z "$empty_dirs" ]; then
    print_success "No empty directories found"
    exit 0
fi

count=$(echo "$empty_dirs" | wc -l | tr -d '[:space:]')
print_info "Found $count empty directorie(s):"
print_info ""

echo "$empty_dirs" | while read -r dir; do
    echo "  - $dir"
done

print_info ""

if [ "$DRY_RUN" = true ]; then
    print_warning "Dry run - no directories removed"
    print_info ""
    print_info "To remove these directories, run without --dry-run"
    exit 0
fi

# Confirm removal
print_warning "This will remove $count empty directorie(s)"
read -r -p "Continue? [y/N] " confirm

if [[ ! $confirm =~ ^[Yy]$ ]]; then
    print_info "Cancelled"
    exit 0
fi

# Remove empty directories
removed=0
failed=0

while IFS= read -r dir; do
    if rmdir "$dir" 2>/dev/null; then
        print_success "Removed: $dir"
        removed=$((removed + 1))
    else
        print_error "Failed to remove: $dir"
        failed=$((failed + 1))
    fi
done <<<"$empty_dirs"

print_info ""
print_info "Removed: $removed, Failed: $failed"

if [ "$failed" -gt 0 ]; then
    print_warning "Some directories could not be removed (may have files or permission issues)"
fi

print_info ""
print_success "Cleanup complete"
