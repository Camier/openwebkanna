#!/usr/bin/env bash

###############################################################################
# OpenWebUI plugin code audit
# Verifies that tool/function Python code stored in webui.db compiles and imports cleanly.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

TARGET_DIR="${TARGET_DIR:-$SCRIPT_DIR}"

OPENWEBUI_SERVICE="${OPENWEBUI_SERVICE:-openwebui}"
COMPOSE_FILE="${COMPOSE_FILE:-docker-compose.yml}"
FOCUS="${PLUGIN_AUDIT_FOCUS:-all}"
PLUGIN_AUDIT_TIMEOUT_SECONDS="${PLUGIN_AUDIT_TIMEOUT_SECONDS:-60}"
PLUGIN_AUDIT_CHECK_PIPELINES_DIR="${PLUGIN_AUDIT_CHECK_PIPELINES_DIR:-false}"
PLUGIN_AUDIT_IMPORT_CHECK="${PLUGIN_AUDIT_IMPORT_CHECK:-true}"
PIPELINES_DIR="${PIPELINES_DIR:-}"
# shellcheck disable=SC2034
QUIET=false
# shellcheck disable=SC2034
COMPOSE_CMD=()

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --focus all|tool|function   Scope audit to a specific plugin type (default: all)
  --target PATH               Target directory containing compose file (default: script dir)
  --compose-file PATH         Compose file path relative to target dir (default: docker-compose.yml)
  --check-pipelines-dir       Compile Python files from PIPELINES_DIR after DB audit
  --quiet                     Only print failures and final status
  -h, --help                  Show this help

Environment:
  PLUGIN_AUDIT_IMPORT_CHECK=true|false   Validate that imported modules are installed (default: true)
EOF
}

validate_focus() {
    case "$FOCUS" in
        all | tool | function) return 0 ;;
    esac
    print_error "Invalid focus: $FOCUS (expected: all, tool, or function)"
    return 1
}

is_openwebui_running() {
    local services=""
    init_compose_cmd || return 1
    services="$("${COMPOSE_CMD[@]}" -f "$COMPOSE_FILE" ps --status running --services 2>/dev/null || true)"
    printf "%s\n" "$services" | grep -qx "$OPENWEBUI_SERVICE"
}

# Detect timeout command at runtime (Linux: timeout, macOS: gtimeout)
get_timeout_cmd() {
    if command_exists timeout; then
        printf "timeout"
    elif command_exists gtimeout; then
        printf "gtimeout"
    else
        printf ""
    fi
}

run_python_audit() {
    init_compose_cmd || return 1

    local audit_script="${SCRIPT_DIR}/lib/audit_plugins.py"
    if [ ! -f "$audit_script" ]; then
        print_error "Audit script not found: $audit_script"
        return 1
    fi

    local timeout_cmd
    timeout_cmd="$(get_timeout_cmd)"

    local -a compose_exec_cmd=("${COMPOSE_CMD[@]}" -f "$COMPOSE_FILE" exec -T "$OPENWEBUI_SERVICE" python "$audit_script" "$FOCUS" "$PLUGIN_AUDIT_IMPORT_CHECK")

    if [ -n "$timeout_cmd" ]; then
        "$timeout_cmd" "${PLUGIN_AUDIT_TIMEOUT_SECONDS}s" "${compose_exec_cmd[@]}"
        return $?
    fi

    # No timeout available, run directly
    "${compose_exec_cmd[@]}"
}

run_filesystem_pipeline_audit() {
    local pipelines_root="$PIPELINES_DIR"
    local failed=0
    local file=""

    if ! is_true "$PLUGIN_AUDIT_CHECK_PIPELINES_DIR"; then
        return 0
    fi

    if [ -z "$pipelines_root" ]; then
        print_warning "PIPELINES_DIR is not set; skipping filesystem pipeline audit"
        return 0
    fi

    if [[ $pipelines_root != /* ]]; then
        pipelines_root="${TARGET_DIR}/${pipelines_root}"
    fi

    if [ ! -d "$pipelines_root" ]; then
        print_error "PIPELINES_DIR does not exist: $pipelines_root"
        return 1
    fi

    if ! command_exists python; then
        print_error "python is required for filesystem pipeline audit"
        return 1
    fi

    print_step "Compiling pipeline .py files from PIPELINES_DIR"
    while IFS= read -r file; do
        if ! python -m py_compile "$file" >/dev/null 2>&1; then
            print_error "Syntax error in pipeline file: $file"
            failed=1
        fi
    done < <(find "$pipelines_root" -type f -name '*.py' | sort)

    if [ "$failed" -ne 0 ]; then
        return 1
    fi

    print_success "Filesystem pipeline audit passed (${pipelines_root})"
    return 0
}

main() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --focus)
                if [[ $# -lt 2 ]]; then
                    print_error "--focus requires a value"
                    exit 1
                fi
                FOCUS="$2"
                shift 2
                ;;
            --focus=*)
                FOCUS="${1#*=}"
                shift
                ;;
            focus=*)
                FOCUS="${1#*=}"
                shift
                ;;
            --target)
                if [[ $# -lt 2 ]]; then
                    print_error "--target requires a path"
                    exit 1
                fi
                TARGET_DIR="$2"
                shift 2
                ;;
            --compose-file)
                if [[ $# -lt 2 ]]; then
                    print_error "--compose-file requires a path"
                    exit 1
                fi
                COMPOSE_FILE="$2"
                shift 2
                ;;
            --check-pipelines-dir)
                PLUGIN_AUDIT_CHECK_PIPELINES_DIR=true
                shift
                ;;
            --quiet)
                # shellcheck disable=SC2034
                QUIET=true
                shift
                ;;
            -h | --help)
                usage
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done

    validate_focus || exit 1

    if [ ! -d "$TARGET_DIR" ]; then
        print_error "Target directory not found: $TARGET_DIR"
        exit 1
    fi
    cd "$TARGET_DIR"

    if ! command -v docker >/dev/null 2>&1; then
        print_error "docker is required"
        exit 1
    fi
    if ! init_compose_cmd; then
        exit 1
    fi

    print_header "OpenWebUI Plugin Code Audit"
    print_step "Checking OpenWebUI runtime availability"

    if ! is_openwebui_running; then
        print_error "OpenWebUI service '$OPENWEBUI_SERVICE' is not running"
        print_warning "Start services first (e.g., ./deploy.sh) then rerun this audit"
        exit 1
    fi

    print_step "Auditing plugin code from webui.db (focus: $FOCUS)"
    set +e
    audit_output="$(run_python_audit 2>&1)"
    audit_rc=$?
    set -e

    checked_count="$(printf "%s\n" "$audit_output" | sed -n 's/^AUDIT_CHECKED=//p' | tail -n1)"
    issue_count="$(printf "%s\n" "$audit_output" | sed -n 's/^AUDIT_ISSUES=//p' | tail -n1)"

    if [ -z "$checked_count" ]; then
        print_error "Audit failed before producing a result"
        printf "%s\n" "$audit_output" >&2
        exit 1
    fi

    if [ "$audit_rc" -ne 0 ]; then
        if [ "$audit_rc" -eq 124 ] || [ "$audit_rc" -eq 137 ]; then
            print_error "Plugin code audit timed out after ${PLUGIN_AUDIT_TIMEOUT_SECONDS}s"
            exit 1
        fi
        print_error "Plugin code audit failed (${issue_count:-unknown} issues)"
        printf "%s\n" "$audit_output" | grep '^ISSUE|' | while IFS='|' read -r _ table id name line _offset msg code; do
            echo "  - ${table}:${id} (${name}) line ${line}: ${msg}"
            echo "    code: ${code}"
        done
        exit 1
    fi

    print_success "Plugin code audit passed (checked ${checked_count} entries)"

    if ! run_filesystem_pipeline_audit; then
        exit 1
    fi
}

main "$@"
