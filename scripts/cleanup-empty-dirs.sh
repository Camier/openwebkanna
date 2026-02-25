#!/bin/bash

###############################################################################
# Empty Directory Cleanup
# Identifies empty directories (excluding .git, archive, data/milvus, and prod_max)
# These are typically runtime artifacts that should not be tracked in git.
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

source "${SCRIPT_DIR}/lib/init.sh"
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
  - archive/ (deprecated but kept for reference)
  - data/milvus/ (runtime state)
  - prod_max/ and prod_max_multimodal/ (processing outputs)
  - logs/, backups/, certs/ (runtime directories)
  - Any directory containing .gitkeep files

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
export EXCLUDE_PATTERN="\.git|archive|data/milvus|prod_max|prod_max_multimodal|logs|backups|certs|data/pdfs|data/extractions|data/corpus|jupyter|cliproxyapi|searxng|thesis-exports"

print_step "Finding empty directories..."
print_info "Excluded: .git, archive, data/milvus, prod_max*, logs, backups, certs, runtime dirs"
print_info ""

# Find empty directories, excluding special dirs
empty_dirs=$(find . -type d -empty \
    ! -path "./.git*" \
    ! -path "./archive*" \
    ! -path "./data/milvus*" \
    ! -path "./prod_max*" \
    ! -path "./prod_max_multimodal*" \
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
read -p "Continue? [y/N] " confirm

if [[ ! $confirm =~ ^[Yy]$ ]]; then
    print_info "Cancelled"
    exit 0
fi

# Remove empty directories
removed=0
failed=0

echo "$empty_dirs" | while read -r dir; do
    if rmdir "$dir" 2>/dev/null; then
        print_success "Removed: $dir"
        removed=$((removed + 1))
    else
        print_error "Failed to remove: $dir"
        failed=$((failed + 1))
    fi
done

print_info ""
print_info "Removed: $removed, Failed: $failed"

if [ $failed -gt 0 ]; then
    print_warning "Some directories could not be removed (may have files or permission issues)"
fi

print_info ""
print_success "Cleanup complete"
