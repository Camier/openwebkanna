#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${PROJECT_ROOT}/lib/init.sh"

cd "${PROJECT_ROOT}"

PRIMARY_DOCS=(
    "README.md"
    "docs/ssot/stack.md"
    "docs/ssot/stack.yaml"
    "config/README.md"
    "scripts/README.md"
    "artifacts/README.md"
    "cliproxyapi/README.md"
    "research/README.md"
    "docs/runbooks/SETUP.md"
    "docs/runbooks/OPERATIONS.md"
    "docs/runbooks/PREREQUISITES.md"
    "docs/runbooks/TROUBLESHOOTING.md"
    "MCP_INTEGRATION_GUIDE.md"
)

VERSION_ALIGNMENT_DOCS=(
    "README.md"
    "docs/runbooks/SETUP.md"
    "docs/runbooks/OPERATIONS.md"
    "docs/runbooks/PREREQUISITES.md"
    "docs/runbooks/TROUBLESHOOTING.md"
    "MCP_INTEGRATION_GUIDE.md"
    "SECURITY.md"
    "AGENTS.md"
    "docs/status/WEBSEARCH_DEPLOYMENT_STATUS.md"
)

LINK_DOCS=(
    "README.md"
    "CLAUDE.md"
    "GEMINI.md"
    "MCP_INTEGRATION_GUIDE.md"
    "docs/README.md"
    "docs/ssot/stack.md"
    "docs/legacy/README.md"
    "docs/runbooks/SETUP.md"
    "docs/runbooks/OPERATIONS.md"
    "docs/runbooks/PREREQUISITES.md"
    "docs/runbooks/TROUBLESHOOTING.md"
    "docs/status/WEBSEARCH_DEPLOYMENT_STATUS.md"
    "config/README.md"
    "scripts/README.md"
    "artifacts/README.md"
    "cliproxyapi/README.md"
    "research/README.md"
    "docs/reference/openwebui/OPENWEBUI_PIPELINES_GUIDE.md"
    "docs/reference/openwebui/OPENWEBUI_WEBSEARCH_DEEP_DIVE.md"
    "docs/reference/openwebui/OPENWEBUI_MASTER_REFERENCE.md"
    "docs/reference/openwebui/openwebui_documentation.md"
    "docs/reference/openwebui/openwebui_env_reference.md"
    "docs/reference/openwebui/openwebui_rag_technical_reference.md"
    "docs/reference/openwebui/openwebui-tools-functions-guide.md"
)

INVALID_TUNE_FLAGS=(
    "--top-k"
    "--chunk-size"
    "--chunk-overlap"
    "--snapshot-only"
)

FAILURES=0

COMPAT_COPY_MANIFEST="config/compatibility-copies.txt"
COMPAT_COPY_PAIRS=()

current_openwebui_version() {
    local image
    image=$(sed -n 's/^OPENWEBUI_IMAGE=.*:\(v[0-9][0-9.]*\)$/\1/p' config/env/.env.example | head -n 1)
    if [ -z "$image" ]; then
        record_failure "Unable to derive current OpenWebUI version from config/env/.env.example"
        return 1
    fi
    printf '%s\n' "$image"
}

record_failure() {
    local message=$1
    print_error "$message"
    FAILURES=$((FAILURES + 1))
}

require_file() {
    local file=$1
    if [ ! -f "$file" ]; then
        record_failure "Missing required documentation file: $file"
        return 1
    fi
    return 0
}

check_required_docs_exist() {
    local file
    print_section "Required Doc Presence Check"

    for file in "${PRIMARY_DOCS[@]}"; do
        if require_file "$file"; then
            print_success "Required doc present: $file"
        fi
    done
}

load_compatibility_copy_pairs() {
    if ! require_file "$COMPAT_COPY_MANIFEST"; then
        return 1
    fi

    COMPAT_COPY_PAIRS=()

    local entry=""
    while IFS= read -r entry || [ -n "$entry" ]; do
        entry="${entry%$'\r'}"
        entry="${entry#"${entry%%[![:space:]]*}"}"
        entry="${entry%"${entry##*[![:space:]]}"}"

        [ -z "$entry" ] && continue
        [[ $entry == \#* ]] && continue

        if [[ $entry != *:* ]]; then
            record_failure "Invalid compatibility copy manifest entry: $entry"
            continue
        fi

        COMPAT_COPY_PAIRS+=("$entry")
    done <"$COMPAT_COPY_MANIFEST"

    if ((${#COMPAT_COPY_PAIRS[@]} == 0)); then
        record_failure "No compatibility copy pairs declared in $COMPAT_COPY_MANIFEST"
        return 1
    fi

    print_success "Loaded ${#COMPAT_COPY_PAIRS[@]} compatibility copy pair(s) from $COMPAT_COPY_MANIFEST"
    return 0
}

check_compatibility_copies() {
    local entry
    print_section "Compatibility Copy Check"

    if ! load_compatibility_copy_pairs; then
        return
    fi

    for entry in "${COMPAT_COPY_PAIRS[@]}"; do
        local root_path="${entry%%:*}"
        local canonical_path="${entry#*:}"

        if ! require_file "$root_path"; then
            continue
        fi
        if ! require_file "$canonical_path"; then
            continue
        fi

        if cmp -s "$root_path" "$canonical_path"; then
            print_success "Compatibility copy aligned: $root_path"
        else
            record_failure "Compatibility copy drift: $root_path != $canonical_path"
        fi
    done
}

normalize_link_target() {
    local target=$1

    target="${target#"${target%%[![:space:]]*}"}"
    target="${target%"${target##*[![:space:]]}"}"

    if [[ $target == \<*\> ]]; then
        target="${target#<}"
        target="${target%>}"
    else
        target="${target%%[[:space:]]*}"
    fi

    target="${target#\"}"
    target="${target%\"}"
    target="${target#\'}"
    target="${target%\'}"

    printf '%s\n' "$target"
}

is_relative_md_link() {
    local target=$1
    case "$target" in
        "" | http://* | https://* | mailto:* | \#* | /*)
            return 1
            ;;
    esac

    local path_without_ref="${target%%\#*}"
    path_without_ref="${path_without_ref%%\?*}"

    [[ $path_without_ref == *.md ]]
}

extract_links() {
    local file=$1
    grep -Eo '\[[^]]+\]\([^)]*\)' "$file" | sed -E 's/^[^[]*\[[^]]+\]\(([^)]*)\)$/\1/' || true
    grep -Eo '^[[:space:]]*\[[^]]+\]:[[:space:]]+[^[:space:]]+' "$file" | sed -E 's/^[[:space:]]*\[[^]]+\]:[[:space:]]+//' || true
}

check_stale_version_markers() {
    local current_version
    current_version=$(current_openwebui_version) || return
    local marker="v0.8.3"
    local file
    print_section "Stale Version Marker Check"
    for file in "${VERSION_ALIGNMENT_DOCS[@]}"; do
        if ! require_file "$file"; then
            continue
        fi

        local matches
        matches=$(rg -n --fixed-strings -- "$marker" "$file" | rg -v "scan baseline" || true)
        if [ -n "$matches" ]; then
            record_failure "Stale version marker '$marker' found in $file"
            while IFS= read -r line; do
                print_info "$file:$line"
            done <<<"$matches"
        else
            print_success "No stale marker in $file (current OpenWebUI version: $current_version)"
        fi
    done
}

check_stale_port_markers() {
    local file
    local marker="http://localhost:3010"
    print_section "Stale Port Marker Check"

    for file in "${VERSION_ALIGNMENT_DOCS[@]}"; do
        if ! require_file "$file"; then
            continue
        fi

        if [ "$file" = "docs/status/WEBSEARCH_DEPLOYMENT_STATUS.md" ]; then
            print_success "Live-status doc may intentionally reference localhost:3010: $file"
            continue
        fi

        local matches
        matches=$(rg -n --fixed-strings -- "$marker" "$file" || true)
        if [ -n "$matches" ]; then
            record_failure "Stale OpenWebUI port marker '$marker' found in $file"
            while IFS= read -r line; do
                print_info "$file:$line"
            done <<<"$matches"
        else
            print_success "No stale localhost:3010 marker in $file"
        fi
    done
}

check_stale_websearch_topology_markers() {
    local file="docs/status/WEBSEARCH_DEPLOYMENT_STATUS.md"
    local marker
    print_section "Web Search Topology Check"
    if ! require_file "$file"; then
        return
    fi

    for marker in "/LAB/@ai_hub/searxng" "systemd service" "host.docker.internal:8888/search"; do
        local matches
        matches=$(rg -n --fixed-strings -- "$marker" "$file" || true)
        if [ -n "$matches" ]; then
            record_failure "Stale web search topology marker '$marker' found in $file"
            while IFS= read -r line; do
                print_info "$file:$line"
            done <<<"$matches"
        fi
    done
}

check_invalid_operations_flags() {
    local file="docs/runbooks/OPERATIONS.md"
    local flag
    local found_invalid=0
    print_section "Operations Flag Check"
    if ! require_file "$file"; then
        return
    fi

    for flag in "${INVALID_TUNE_FLAGS[@]}"; do
        local matches
        matches=$(rg -n --fixed-strings -- "$flag" "$file" || true)
        if [ -n "$matches" ]; then
            found_invalid=1
            record_failure "Invalid tune flag '$flag' found in $file"
            while IFS= read -r line; do
                print_info "$file:$line"
            done <<<"$matches"
        fi
    done

    if ((found_invalid == 0)); then
        print_success "No invalid tune CLI flags found in $file"
    fi
}

check_mcp_phrase_consistency() {
    local file="MCP_INTEGRATION_GUIDE.md"
    local phrase="not OpenAPI"
    print_section "MCP Phrase Consistency Check"
    if ! require_file "$file"; then
        return
    fi

    local matches
    matches=$(rg -n -i --fixed-strings -- "$phrase" "$file" || true)
    if [ -n "$matches" ]; then
        record_failure "Contradictory phrase '$phrase' found in $file"
        while IFS= read -r line; do
            print_info "$file:$line"
        done <<<"$matches"
    else
        print_success "No contradictory phrase found in $file"
    fi
}

check_relative_markdown_links() {
    local file
    print_section "Relative Markdown Link Check"

    for file in "${LINK_DOCS[@]}"; do
        if ! require_file "$file"; then
            continue
        fi

        local file_failures=0
        while IFS= read -r raw_target; do
            [ -n "$raw_target" ] || continue

            local target
            target=$(normalize_link_target "$raw_target")
            if ! is_relative_md_link "$target"; then
                continue
            fi

            local path_without_ref
            path_without_ref="${target%%\#*}"
            path_without_ref="${path_without_ref%%\?*}"

            local resolved_path
            resolved_path="$(dirname "$file")/$path_without_ref"
            if [ ! -f "$resolved_path" ]; then
                record_failure "Broken markdown link in $file -> $target"
                file_failures=$((file_failures + 1))
            fi
        done < <(extract_links "$file")

        if ((file_failures == 0)); then
            print_success "Relative markdown links valid in $file"
        fi
    done
}

main() {
    print_header "Doc Consistency Check"

    check_stale_version_markers
    check_stale_port_markers
    check_stale_websearch_topology_markers
    check_invalid_operations_flags
    check_mcp_phrase_consistency
    check_required_docs_exist
    check_compatibility_copies
    check_relative_markdown_links

    if ((FAILURES == 0)); then
        print_success "Doc consistency checks passed"
        exit 0
    fi

    print_error "Doc consistency checks failed with $FAILURES issue(s)"
    exit 1
}

main "$@"
