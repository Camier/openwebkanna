#!/bin/bash

###############################################################################
# CLIProxyAPI Official Release Installer
# Downloads and installs upstream CLIProxyAPI release binaries from GitHub
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

load_env_defaults() {
    local env_file=""
    local line=""
    local key=""
    local value=""

    if [ -f "$SCRIPT_DIR/.env" ]; then
        env_file="$SCRIPT_DIR/.env"
    elif [ -f ".env" ]; then
        env_file=".env"
    fi

    if [ -z "$env_file" ]; then
        return 0
    fi

    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%$'\r'}"
        line="${line#"${line%%[![:space:]]*}"}"
        [ -z "$line" ] && continue
        [[ "$line" = \#* ]] && continue
        [[ "$line" != *=* ]] && continue

        key="${line%%=*}"
        value="${line#*=}"
        key="$(printf "%s" "$key" | tr -d '[:space:]')"

        [[ "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        if [ -n "${!key+x}" ]; then
            continue
        fi

        if [[ "$value" == \"*\" ]] && [[ "$value" == *\" ]]; then
            value="${value:1:${#value}-2}"
        elif [[ "$value" == \'*\' ]] && [[ "$value" == *\' ]]; then
            value="${value:1:${#value}-2}"
        fi

        printf -v "$key" "%s" "$value"
        export "$key"
    done < "$env_file"
}

load_env_defaults

CLIPROXYAPI_REPO="${CLIPROXYAPI_REPO:-router-for-me/CLIProxyAPI}"
CLIPROXYAPI_VERSION="${CLIPROXYAPI_VERSION:-latest}"
CLIPROXYAPI_RELEASE_API="https://api.github.com/repos/${CLIPROXYAPI_REPO}/releases/latest"
CLIPROXYAPI_RELEASE_BASE_URL="https://github.com/${CLIPROXYAPI_REPO}/releases/download"
CLIPROXYAPI_INSTALL_PATH="${CLIPROXYAPI_INSTALL_PATH:-./bin/cliproxyapi}"
CLIPROXYAPI_UPSTREAM_CONFIG_SNAPSHOT="${CLIPROXYAPI_UPSTREAM_CONFIG_SNAPSHOT:-./cliproxyapi/config.upstream.example.yaml}"
CLIPROXYAPI_RELEASE_METADATA_FILE="${CLIPROXYAPI_RELEASE_METADATA_FILE:-./cliproxyapi/upstream-release.txt}"
CLIPROXYAPI_FORCE_REINSTALL="${CLIPROXYAPI_FORCE_REINSTALL:-false}"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

print_header() {
    echo -e "${CYAN}${BOLD}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║     CLIProxyAPI Official Installer                        ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_step() {
    echo -e "\n${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}" >&2
}

print_info() {
    echo -e "${CYAN}ℹ $1${NC}"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ]
}

resolve_path() {
    local candidate="$1"
    if [[ "$candidate" = /* ]]; then
        printf "%s" "$candidate"
    else
        printf "%s/%s" "$SCRIPT_DIR" "$candidate"
    fi
}

detect_platform_asset_suffix() {
    local os_name=""
    local arch_name=""
    os_name="$(uname -s | tr '[:upper:]' '[:lower:]')"
    arch_name="$(uname -m)"

    case "$os_name" in
        linux|darwin)
            ;;
        *)
            print_error "Unsupported OS for upstream binary install: $os_name"
            return 1
            ;;
    esac

    case "$arch_name" in
        x86_64|amd64)
            arch_name="amd64"
            ;;
        aarch64|arm64)
            arch_name="arm64"
            ;;
        *)
            print_error "Unsupported architecture for upstream binary install: $arch_name"
            return 1
            ;;
    esac

    printf "%s_%s" "$os_name" "$arch_name"
}

fetch_latest_release_metadata() {
    local output_file="$1"
    curl -sSfL "$CLIPROXYAPI_RELEASE_API" -o "$output_file"
}

json_value() {
    local file="$1"
    local expr="$2"
    local value=""

    if command_exists jq; then
        value="$(jq -r "$expr" "$file" 2>/dev/null || true)"
    else
        case "$expr" in
            ".tag_name")
                value="$(sed -n 's/.*"tag_name"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$file" | head -n 1)"
                ;;
            ".published_at")
                value="$(sed -n 's/.*"published_at"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$file" | head -n 1)"
                ;;
        esac
    fi

    if [ "$value" = "null" ]; then
        value=""
    fi

    printf "%s" "$value"
}

verify_checksums() {
    local checksums_path="$1"
    local archive_name="$2"
    local archive_path="$3"
    local expected_hash=""
    local actual_hash=""

    expected_hash="$(grep -E "[[:space:]]${archive_name}$" "$checksums_path" | awk '{print $1}' | head -n 1)"
    if [ -z "$expected_hash" ]; then
        print_warning "Could not find checksum entry for ${archive_name}; skipping checksum verification"
        return 0
    fi

    if command_exists sha256sum; then
        actual_hash="$(sha256sum "$archive_path" | awk '{print $1}')"
    elif command_exists shasum; then
        actual_hash="$(shasum -a 256 "$archive_path" | awk '{print $1}')"
    else
        print_warning "No sha256 tool found; skipping checksum verification"
        return 0
    fi

    if [ "$actual_hash" != "$expected_hash" ]; then
        print_error "Checksum verification failed for ${archive_name}"
        print_info "Expected: ${expected_hash}"
        print_info "Actual:   ${actual_hash}"
        return 1
    fi

    print_success "Checksum verified for ${archive_name}"
    return 0
}

main() {
    local install_path=""
    local install_dir=""
    local install_basename=""
    local canonical_link_path=""
    local tmp_dir=""
    local release_tag=""
    local release_published_at=""
    local version_number=""
    local platform_suffix=""
    local archive_name=""
    local release_url=""
    local archive_path=""
    local checksums_path=""
    local metadata_path=""
    local extracted_binary=""
    local extracted_config=""
    local snapshot_path=""
    local metadata_file_path=""
    local help_output=""
    local version_line=""

    print_header

    if ! command_exists curl; then
        print_error "curl is required"
        return 1
    fi
    if ! command_exists tar; then
        print_error "tar is required"
        return 1
    fi

    install_path="$(resolve_path "$CLIPROXYAPI_INSTALL_PATH")"
    install_dir="$(dirname "$install_path")"
    install_basename="$(basename "$install_path")"
    canonical_link_path="${install_dir}/cli-proxy-api"
    snapshot_path="$(resolve_path "$CLIPROXYAPI_UPSTREAM_CONFIG_SNAPSHOT")"
    metadata_file_path="$(resolve_path "$CLIPROXYAPI_RELEASE_METADATA_FILE")"

    if [ -x "$install_path" ] && ! is_true "$CLIPROXYAPI_FORCE_REINSTALL"; then
        print_step "CLIProxyAPI binary already present"
        print_info "Path: $install_path"
        print_info "Set CLIPROXYAPI_FORCE_REINSTALL=true to refresh from upstream"
        return 0
    fi

    print_step "Resolving release metadata"
    tmp_dir="$(mktemp -d)"
    trap 'rm -rf "$tmp_dir"' EXIT
    metadata_path="${tmp_dir}/release.json"

    if [ "$CLIPROXYAPI_VERSION" = "latest" ]; then
        fetch_latest_release_metadata "$metadata_path"
        release_tag="$(json_value "$metadata_path" ".tag_name")"
        release_published_at="$(json_value "$metadata_path" ".published_at")"
    else
        if [[ "$CLIPROXYAPI_VERSION" == v* ]]; then
            release_tag="$CLIPROXYAPI_VERSION"
        else
            release_tag="v${CLIPROXYAPI_VERSION}"
        fi
        release_published_at="unknown"
    fi

    if [ -z "$release_tag" ]; then
        print_error "Unable to resolve CLIProxyAPI release tag"
        return 1
    fi

    version_number="${release_tag#v}"
    platform_suffix="$(detect_platform_asset_suffix)"
    archive_name="CLIProxyAPI_${version_number}_${platform_suffix}.tar.gz"
    release_url="${CLIPROXYAPI_RELEASE_BASE_URL}/${release_tag}/${archive_name}"
    archive_path="${tmp_dir}/${archive_name}"
    checksums_path="${tmp_dir}/checksums.txt"

    print_info "Release: ${release_tag}"
    print_info "Published: ${release_published_at}"
    print_info "Asset: ${archive_name}"

    print_step "Downloading official release asset"
    curl -sSfL "$release_url" -o "$archive_path"
    print_success "Downloaded ${archive_name}"

    if curl -sSfL "${CLIPROXYAPI_RELEASE_BASE_URL}/${release_tag}/checksums.txt" -o "$checksums_path"; then
        verify_checksums "$checksums_path" "$archive_name" "$archive_path"
    else
        print_warning "Unable to download checksums.txt; proceeding without checksum verification"
    fi

    print_step "Installing binary"
    tar -xzf "$archive_path" -C "$tmp_dir"
    extracted_binary="${tmp_dir}/cli-proxy-api"
    extracted_config="${tmp_dir}/config.example.yaml"

    if [ ! -f "$extracted_binary" ]; then
        print_error "Release archive does not contain cli-proxy-api binary"
        return 1
    fi

    mkdir -p "$install_dir"
    install -m 755 "$extracted_binary" "$install_path"
    ln -sf "$install_basename" "$canonical_link_path"
    print_success "Installed binary to $install_path"

    if [ -f "$extracted_config" ]; then
        mkdir -p "$(dirname "$snapshot_path")"
        cp "$extracted_config" "$snapshot_path"
        print_success "Saved upstream config snapshot to $snapshot_path"
    fi

    mkdir -p "$(dirname "$metadata_file_path")"
    {
        printf "repo=%s\n" "$CLIPROXYAPI_REPO"
        printf "tag=%s\n" "$release_tag"
        printf "published_at=%s\n" "$release_published_at"
        printf "asset=%s\n" "$archive_name"
        printf "installed_at_utc=%s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
    } > "$metadata_file_path"
    print_success "Wrote release metadata to $metadata_file_path"

    help_output="$("$install_path" --help 2>&1 || true)"
    version_line="$(printf "%s\n" "$help_output" | grep -m1 "CLIProxyAPI Version:" || true)"
    if [ -n "$version_line" ]; then
        print_info "$version_line"
    fi

    print_success "CLIProxyAPI official setup complete"
    return 0
}

main "$@"
