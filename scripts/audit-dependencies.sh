#!/bin/bash

###############################################################################
# audit-dependencies.sh - Security audit for Docker image dependencies
#
# Checks Docker image versions and runs Trivy vulnerability scans for
# OpenWebUI, CLIProxyAPI, pgvector, and Jupyter containers.
#
# Usage:
#   ./scripts/audit-dependencies.sh [--skip-scan]
#
# Options:
#   --skip-scan    Skip Trivy vulnerability scanning (version check only)
#
# Exit codes:
#   0 - No CRITICAL CVEs found
#   1 - CRITICAL CVEs found or error occurred
###############################################################################

set -euo pipefail

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LIB_DIR="${PROJECT_ROOT}/lib"

# Source shared libraries if they exist
if [[ -f "${LIB_DIR}/colors.sh" ]]; then
    source "${LIB_DIR}/colors.sh"
fi

if [[ -f "${LIB_DIR}/print-utils.sh" ]]; then
    source "${LIB_DIR}/print-utils.sh"
fi

if [[ -f "${LIB_DIR}/docker-helpers.sh" ]]; then
    source "${LIB_DIR}/docker-helpers.sh"
fi

# Fallback color definitions if colors.sh not available
if [[ -z ${RED:-} ]]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    CYAN='\033[0;36m'
    NC='\033[0m'
    BOLD='\033[1m'
fi

# Fallback print functions if print-utils.sh not available
if ! declare -f print_header >/dev/null 2>&1; then
    print_header() {
        echo -e "${CYAN}${BOLD}"
        echo "+------------------------------------------------------------+"
        echo "|                    $1                    |"
        echo "+------------------------------------------------------------+"
        echo -e "${NC}"
    }
fi

if ! declare -f print_step >/dev/null 2>&1; then
    print_step() { echo -e "\n${BLUE}${BOLD}>> $1${NC}"; }
fi

if ! declare -f print_success >/dev/null 2>&1; then
    print_success() { echo -e "${GREEN}[OK] $1${NC}"; }
fi

if ! declare -f print_error >/dev/null 2>&1; then
    print_error() { echo -e "${RED}[ERROR] $1${NC}" >&2; }
fi

if ! declare -f print_warning >/dev/null 2>&1; then
    print_warning() { echo -e "${YELLOW}[WARN] $1${NC}"; }
fi

if ! declare -f print_info >/dev/null 2>&1; then
    print_info() { echo -e "${CYAN}[INFO] $1${NC}"; }
fi

if ! declare -f print_section >/dev/null 2>&1; then
    print_section() { echo -e "\n${BOLD}--- $1 ---${NC}\n"; }
fi

###############################################################################
# Configuration
###############################################################################

# Images to audit (name:service_name format)
# service_name corresponds to docker-compose service names
declare -A IMAGES=(
    ["openwebui"]="openwebui"
    ["cliproxyapi"]="cliproxyapi"
    ["pgvector"]="pgvector"
    ["jupyter"]="jupyter"
)

# Track CVE counts
declare -A CVE_COUNTS=(
    ["CRITICAL"]=0
    ["HIGH"]=0
    ["MEDIUM"]=0
    ["LOW"]=0
    ["UNKNOWN"]=0
)

SKIP_SCAN=false
CRITICAL_FOUND=false

###############################################################################
# Helper Functions
###############################################################################

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Security audit for Docker image dependencies.

Options:
    --skip-scan    Skip Trivy vulnerability scanning (version check only)
    -h, --help     Show this help message

Examples:
    $(basename "$0")              # Full audit with Trivy scans
    $(basename "$0") --skip-scan  # Version check only

Exit codes:
    0 - No CRITICAL CVEs found
    1 - CRITICAL CVEs found or error occurred
EOF
}

check_trivy_installed() {
    if ! command -v trivy >/dev/null 2>&1; then
        print_warning "Trivy is not installed. Install with:"
        echo "  curl -sf https://raw.githubusercontent.com/aquasecurity/trivy/main/contrib/install.sh | sh -s -- -b /usr/local/bin"
        echo "  Or: brew install trivy"
        return 1
    fi
    return 0
}

get_image_info() {
    local service_name="$1"
    local compose_file="${COMPOSE_FILE:-${PROJECT_ROOT}/docker-compose.yml}"

    # Try to get image from docker-compose config
    if [[ -f $compose_file ]]; then
        local image
        image=$(docker compose -f "$compose_file" config 2>/dev/null |
            grep -A2 "^[[:space:]]*${service_name}:" |
            grep "image:" |
            head -1 |
            awk '{print $2}' || true)
        if [[ -n $image ]]; then
            echo "$image"
            return 0
        fi
    fi

    # Fallback: try running container
    local running_image
    running_image=$(docker ps --filter "name=${service_name}" --format "{{.Image}}" 2>/dev/null | head -1 || true)
    if [[ -n $running_image ]]; then
        echo "$running_image"
        return 0
    fi

    return 1
}

get_image_version() {
    local image="$1"

    # Get image ID and inspect for version info
    local image_id
    image_id=$(docker images --filter "reference=${image}" --format "{{.ID}}" 2>/dev/null | head -1 || true)

    if [[ -z $image_id ]]; then
        echo "not-pulled"
        return 1
    fi

    # Try to get version from labels or repo tags
    local version_info
    version_info=$(docker inspect "$image_id" --format '{{.Config.Labels.version}}' 2>/dev/null || true)

    if [[ -n $version_info && $version_info != "<nil>" ]]; then
        echo "$version_info"
        return 0
    fi

    # Fallback to image tag
    local tag
    tag=$(echo "$image" | cut -d: -f2)
    if [[ -n $tag && $tag != "$image" ]]; then
        echo "$tag"
        return 0
    fi

    echo "unknown"
    return 0
}

check_image_freshness() {
    local image="$1"
    local image_id
    image_id=$(docker images --filter "reference=${image}" --format "{{.ID}}" 2>/dev/null | head -1 || true)

    if [[ -z $image_id ]]; then
        print_warning "Image not pulled locally: ${image}"
        return 1
    fi

    # Check when image was created
    local created
    created=$(docker inspect "$image_id" --format '{{.Created}}' 2>/dev/null || true)

    if [[ -n $created ]]; then
        local created_ts
        created_ts=$(date -d "$created" +%s 2>/dev/null || date -j -f "%Y-%m-%dT%H:%M:%S" "$created" +%s 2>/dev/null || echo "0")
        local now_ts
        now_ts=$(date +%s)
        local age_days=$(((now_ts - created_ts) / 86400))

        if [[ $age_days -gt 90 ]]; then
            print_warning "Image is ${age_days} days old - consider updating"
        elif [[ $age_days -gt 30 ]]; then
            print_info "Image is ${age_days} days old"
        else
            print_info "Image age: ${age_days} days"
        fi
    fi

    return 0
}

run_trivy_scan() {
    local service_name="$1"
    local image="$2"
    local output_file="${PROJECT_ROOT}/.trivy-scan-${service_name}.json"

    print_step "Running Trivy scan for ${service_name}..."

    # Check if image exists locally
    if ! docker image inspect "$image" >/dev/null 2>&1; then
        print_warning "Image not found locally, pulling ${image}..."
        if ! docker pull "$image" 2>/dev/null; then
            print_error "Failed to pull image: ${image}"
            return 1
        fi
    fi

    # Run Trivy scan
    if ! trivy image \
        --quiet \
        --format json \
        --output "$output_file" \
        --severity "CRITICAL,HIGH,MEDIUM,LOW,UNKNOWN" \
        "$image" 2>/dev/null; then
        print_warning "Trivy scan completed with warnings for ${service_name}"
    fi

    # Parse results
    if [[ -f $output_file ]]; then
        parse_trivy_results "$service_name" "$output_file"
    else
        print_error "Trivy output file not created"
        return 1
    fi

    return 0
}

parse_trivy_results() {
    local service_name="$1"
    local output_file="$2"

    print_section "Vulnerabilities for ${service_name}"

    if [[ ! -s $output_file ]]; then
        print_info "No vulnerabilities found"
        return 0
    fi

    # Count CVEs by severity using jq if available
    local has_jq=true
    if ! command -v jq >/dev/null 2>&1; then
        has_jq=false
        print_warning "jq not installed - using basic parsing"
    fi

    if $has_jq; then
        # Extract Results array and count vulnerabilities
        local results
        results=$(jq -r '.Results // []' "$output_file" 2>/dev/null)

        # Count by severity
        local critical high medium low unknown
        critical=$(echo "$results" | jq '[.[]?.Vulnerabilities // [] | .[]? | select(.Severity=="CRITICAL")] | length' 2>/dev/null || echo "0")
        high=$(echo "$results" | jq '[.[]?.Vulnerabilities // [] | .[]? | select(.Severity=="HIGH")] | length' 2>/dev/null || echo "0")
        medium=$(echo "$results" | jq '[.[]?.Vulnerabilities // [] | .[]? | select(.Severity=="MEDIUM")] | length' 2>/dev/null || echo "0")
        low=$(echo "$results" | jq '[.[]?.Vulnerabilities // [] | .[]? | select(.Severity=="LOW")] | length' 2>/dev/null || echo "0")
        unknown=$(echo "$results" | jq '[.[]?.Vulnerabilities // [] | .[]? | select(.Severity=="UNKNOWN")] | length' 2>/dev/null || echo "0")

        # Accumulate totals
        CVE_COUNTS["CRITICAL"]=$((CVE_COUNTS["CRITICAL"] + critical))
        CVE_COUNTS["HIGH"]=$((CVE_COUNTS["HIGH"] + high))
        CVE_COUNTS["MEDIUM"]=$((CVE_COUNTS["MEDIUM"] + medium))
        CVE_COUNTS["LOW"]=$((CVE_COUNTS["LOW"] + low))
        CVE_COUNTS["UNKNOWN"]=$((CVE_COUNTS["UNKNOWN"] + unknown))

        # Print summary for this image
        echo ""
        echo "  CRITICAL: ${RED}${critical}${NC}"
        echo "  HIGH:     ${YELLOW}${high}${NC}"
        echo "  MEDIUM:   ${BLUE}${medium}${NC}"
        echo "  LOW:      ${GREEN}${low}${NC}"
        echo "  UNKNOWN:  ${CYAN}${unknown}${NC}"

        # List CRITICAL CVEs
        if [[ $critical -gt 0 ]]; then
            CRITICAL_FOUND=true
            echo ""
            echo -e "${RED}CRITICAL CVEs:${NC}"
            echo "$results" | jq -r '.[]?.Vulnerabilities // [] | .[]? | select(.Severity=="CRITICAL") | "  - \(.VulnerabilityID): \(.Title // .PkgName)"' 2>/dev/null | head -10
            if [[ $(echo "$results" | jq '[.[]?.Vulnerabilities // [] | .[]? | select(.Severity=="CRITICAL")] | length' 2>/dev/null) -gt 10 ]]; then
                echo "  ... (showing first 10)"
            fi
        fi
    else
        # Fallback: simple grep-based counting
        local critical high medium low unknown
        critical=$(grep -c '"Severity":"CRITICAL"' "$output_file" 2>/dev/null || echo "0")
        high=$(grep -c '"Severity":"HIGH"' "$output_file" 2>/dev/null || echo "0")
        medium=$(grep -c '"Severity":"MEDIUM"' "$output_file" 2>/dev/null || echo "0")
        low=$(grep -c '"Severity":"LOW"' "$output_file" 2>/dev/null || echo "0")
        unknown=$(grep -c '"Severity":"UNKNOWN"' "$output_file" 2>/dev/null || echo "0")

        CVE_COUNTS["CRITICAL"]=$((CVE_COUNTS["CRITICAL"] + critical))
        CVE_COUNTS["HIGH"]=$((CVE_COUNTS["HIGH"] + high))
        CVE_COUNTS["MEDIUM"]=$((CVE_COUNTS["MEDIUM"] + medium))
        CVE_COUNTS["LOW"]=$((CVE_COUNTS["LOW"] + low))
        CVE_COUNTS["UNKNOWN"]=$((CVE_COUNTS["UNKNOWN"] + unknown))

        echo ""
        echo "  CRITICAL: ${RED}${critical}${NC}"
        echo "  HIGH:     ${YELLOW}${high}${NC}"
        echo "  MEDIUM:   ${BLUE}${medium}${NC}"
        echo "  LOW:      ${GREEN}${low}${NC}"
        echo "  UNKNOWN:  ${CYAN}${unknown}${NC}"

        if [[ $critical -gt 0 ]]; then
            CRITICAL_FOUND=true
        fi
    fi

    return 0
}

print_summary_report() {
    print_header "Audit Summary"

    echo ""
    echo "CVE Counts by Severity:"
    echo "  CRITICAL: ${RED}${CVE_COUNTS["CRITICAL"]}${NC}"
    echo "  HIGH:     ${YELLOW}${CVE_COUNTS["HIGH"]}${NC}"
    echo "  MEDIUM:   ${BLUE}${CVE_COUNTS["MEDIUM"]}${NC}"
    echo "  LOW:      ${GREEN}${CVE_COUNTS["LOW"]}${NC}"
    echo "  UNKNOWN:  ${CYAN}${CVE_COUNTS["UNKNOWN"]}${NC}"
    echo ""

    local total=0
    for severity in CRITICAL HIGH MEDIUM LOW UNKNOWN; do
        total=$((total + CVE_COUNTS[$severity]))
    done

    echo "Total CVEs: ${BOLD}${total}${NC}"
    echo ""

    if $CRITICAL_FOUND; then
        print_error "CRITICAL vulnerabilities detected!"
        echo ""
        echo "Remediation steps:"
        echo "  1. Update base images to latest versions"
        echo "  2. Review and apply security patches"
        echo "  3. Re-run this audit after updates"
        echo ""
        echo "Detailed scan results saved in:"
        echo "  ${PROJECT_ROOT}/.trivy-scan-*.json"
        return 1
    elif [[ ${CVE_COUNTS["HIGH"]} -gt 0 ]]; then
        print_warning "HIGH severity vulnerabilities found - review recommended"
        return 0
    else
        print_success "No CRITICAL or HIGH vulnerabilities found"
        return 0
    fi
}

###############################################################################
# Main
###############################################################################

main() {
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --skip-scan)
                SKIP_SCAN=true
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

    print_header "Dependency Security Audit"

    # Check Docker availability
    if ! command -v docker >/dev/null 2>&1; then
        print_error "Docker is not installed or not in PATH"
        exit 1
    fi

    # Check Trivy unless skipped
    if ! $SKIP_SCAN; then
        if ! check_trivy_installed; then
            print_warning "Skipping vulnerability scans (Trivy not available)"
            SKIP_SCAN=true
        fi
    fi

    # Check each image
    print_section "Image Version Check"

    declare -A IMAGE_VERSIONS=()

    for service_name in "${IMAGES[@]}"; do
        print_step "Checking ${service_name}..."

        local image
        image=$(get_image_info "$service_name")

        if [[ -z $image ]]; then
            print_warning "Could not determine image for ${service_name}"
            IMAGE_VERSIONS["${service_name}"]="unknown"
            continue
        fi

        print_info "Image: ${image}"

        local version
        version=$(get_image_version "$image")
        IMAGE_VERSIONS["${service_name}"]="${version}"
        print_info "Version: ${version}"

        check_image_freshness "$image"

        if ! $SKIP_SCAN; then
            run_trivy_scan "$service_name" "$image"
        fi
    done

    # Print version summary
    print_section "Version Summary"
    for service_name in "${IMAGES[@]}"; do
        printf "  %-20s %s\n" "${service_name}:" "${IMAGE_VERSIONS[${service_name}]:-unknown}"
    done

    # Print final summary
    if ! $SKIP_SCAN; then
        print_summary_report
        exit_code=$?
    else
        print_info "Vulnerability scanning was skipped"
        exit_code=0
    fi

    exit $exit_code
}

# Run main
main "$@"
