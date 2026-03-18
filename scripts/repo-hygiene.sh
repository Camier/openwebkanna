#!/bin/bash

###############################################################################
# Repository Hygiene
# Removes cache directories and quarantines runtime drift artifacts.
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_ROOT"

# shellcheck disable=SC1091
source "${SCRIPT_DIR}/../lib/init.sh"
load_env_defaults

MODE="dry-run"
STAMP="$(date -u +%F)"
QUARANTINE_BASE="artifacts/repo-hygiene"

show_help() {
    cat <<'EOF'
Usage: ./scripts/repo-hygiene.sh [OPTIONS]

Options:
  --apply                 Apply cleanup and quarantine actions.
  --dry-run               Show planned actions only (default).
  --stamp YYYY-MM-DD      Override report/quarantine date stamp.
  -h, --help              Show this help message.

Behavior:
  1) Removes cache directories:
     - __pycache__
     - .pytest_cache
     - .mypy_cache
     - .ruff_cache
  2) Quarantines known runtime drift artifacts:
     - config/searxng/settings.yml.new
     - logs/*.pid
  3) Quarantines stale timestamped rollup artifacts, keeping only the newest file per family:
     - artifacts/data-quality/curation_stats.*.json
     - data/metadata/quarantine/backups/chunks_corpus.jsonl.pre_curate.*.bak
     - data/metadata/quarantine/chunks/chunks_removed.*.jsonl
  4) Quarantines broken symlinks found in the repo.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)
            MODE="apply"
            shift
            ;;
        --dry-run)
            MODE="dry-run"
            shift
            ;;
        --stamp)
            STAMP="${2:-}"
            if [[ -z ${STAMP} ]]; then
                print_error "--stamp requires a value"
                exit 1
            fi
            shift 2
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

QUARANTINE_DIR="${QUARANTINE_BASE}/${STAMP}/quarantine"
REPORT_DIR="${QUARANTINE_BASE}/${STAMP}/reports"
REPORT_FILE="${REPORT_DIR}/cleanup-manifest.txt"

print_header "Repository Hygiene"
print_info "Mode: ${MODE}"
print_info "Stamp: ${STAMP}"

mapfile -t cache_dirs < <(
    find . -type d \
        \( -name "__pycache__" -o -name ".pytest_cache" -o -name ".mypy_cache" -o -name ".ruff_cache" \) \
        ! -path "./.git/*" \
        -print 2>/dev/null | sed 's#^\./##' | sort -u
)

declare -a quarantine_candidates=()
declare -a known_globs=(
    "config/searxng/settings.yml.new"
    "searxng/settings.yml.new"
    "logs/*.pid"
)
declare -a rolling_keep_specs=(
    "artifacts/data-quality/curation_stats.*.json:1"
    "data/metadata/quarantine/backups/chunks_corpus.jsonl.pre_curate.*.bak:1"
    "data/metadata/quarantine/chunks/chunks_removed.*.jsonl:1"
)

for glob in "${known_globs[@]}"; do
    while IFS= read -r match; do
        [[ -n ${match} ]] && quarantine_candidates+=("${match}")
    done < <(compgen -G "${glob}" || true)
done

for spec in "${rolling_keep_specs[@]}"; do
    glob="${spec%%:*}"
    keep_count="${spec##*:}"

    mapfile -t rollup_matches < <(compgen -G "${glob}" | sort || true)
    if [[ ${#rollup_matches[@]} -le ${keep_count} ]]; then
        continue
    fi

    prune_count=$((${#rollup_matches[@]} - keep_count))
    for ((i = 0; i < prune_count; i++)); do
        [[ -n ${rollup_matches[$i]} ]] && quarantine_candidates+=("${rollup_matches[$i]}")
    done
done

while IFS= read -r broken_link; do
    [[ -n ${broken_link} ]] && quarantine_candidates+=("${broken_link#./}")
done < <(find . -xtype l -print 2>/dev/null)

if [[ ${#quarantine_candidates[@]} -gt 0 ]]; then
    mapfile -t quarantine_candidates < <(printf '%s\n' "${quarantine_candidates[@]}" | sed '/^$/d' | sort -u)
fi

print_section "Detected Artifacts"
print_info "Cache directories: ${#cache_dirs[@]}"
print_info "Quarantine candidates: ${#quarantine_candidates[@]}"

if [[ ${#cache_dirs[@]} -gt 0 ]]; then
    print_info "Cache directories:"
    printf '  - %s\n' "${cache_dirs[@]}"
fi

if [[ ${#quarantine_candidates[@]} -gt 0 ]]; then
    print_info "Quarantine candidates:"
    printf '  - %s\n' "${quarantine_candidates[@]}"
fi

if [[ ${MODE} == "dry-run" ]]; then
    print_warning "Dry-run only. Re-run with --apply to execute."
    exit 0
fi

mkdir -p "${QUARANTINE_DIR}" "${REPORT_DIR}"
{
    echo "# Repository hygiene manifest"
    echo "timestamp_utc=$(date -u +%FT%TZ)"
    echo "mode=${MODE}"
    echo
} >"${REPORT_FILE}"

removed_cache_count=0
for dir in "${cache_dirs[@]}"; do
    if [[ -d ${dir} ]]; then
        rm -rf "${dir}"
        removed_cache_count=$((removed_cache_count + 1))
        printf 'removed_cache\t%s\n' "${dir}" >>"${REPORT_FILE}"
    fi
done

quarantine_count=0
for path in "${quarantine_candidates[@]}"; do
    if [[ -e ${path} || -L ${path} ]]; then
        dest="${QUARANTINE_DIR}/${path}"
        mkdir -p "$(dirname "${dest}")"
        mv "${path}" "${dest}"
        quarantine_count=$((quarantine_count + 1))
        printf 'quarantined\t%s\t%s\n' "${path}" "${dest}" >>"${REPORT_FILE}"
    fi
done

print_section "Summary"
print_success "Removed cache directories: ${removed_cache_count}"
print_success "Quarantined artifacts: ${quarantine_count}"
print_info "Manifest: ${REPORT_FILE}"
