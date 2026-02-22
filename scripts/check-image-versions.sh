#!/bin/bash

###############################################################################
# check-image-versions.sh - Check Docker image versions against latest
#
# Usage:
#   ./scripts/check-image-versions.sh [--json] [--quiet]
#
# Exit codes:
#   0 - All images up to date
#   1 - One or more images outdated (>= 90 days old or newer version available)
#   2 - Script error (missing dependencies, network issues, etc.)
###############################################################################

set -euo pipefail

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Source print utilities
source "${PROJECT_ROOT}/lib/print-utils.sh"

# Configuration
AGE_WARNING_DAYS=90
DOCKER_COMPOSE_FILE="${PROJECT_ROOT}/docker-compose.yml"

# Output options
OUTPUT_JSON=false
QUIET_MODE=false

# Tracking for exit code
ANY_OUTDATED=false
ANY_ERROR=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --json)
            OUTPUT_JSON=true
            shift
            ;;
        --quiet | -q)
            QUIET_MODE=true
            shift
            ;;
        --help | -h)
            echo "Usage: $0 [--json] [--quiet]"
            echo ""
            echo "Options:"
            echo "  --json    Output results as JSON"
            echo "  --quiet   Only show outdated images"
            echo "  --help    Show this help message"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 2
            ;;
    esac
done

# Check dependencies
check_dependencies() {
    local missing=()
    for cmd in curl jq docker; do
        if ! command -v "$cmd" &>/dev/null; then
            missing+=("$cmd")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        print_error "Missing required dependencies: ${missing[*]}"
        exit 2
    fi

    if [[ ! -f $DOCKER_COMPOSE_FILE ]]; then
        print_error "docker-compose.yml not found at $DOCKER_COMPOSE_FILE"
        exit 2
    fi
}

# Extract image info from docker-compose.yml
# Returns: image_name|tag_or_digest
extract_image_from_compose() {
    local service="$1"
    local image_line

    # Get the image line for the service (handles env var defaults)
    image_line=$(grep -A 50 "^  ${service}:" "$DOCKER_COMPOSE_FILE" | grep -m 1 "image:" | sed 's/.*image:[[:space:]]*//')

    if [[ -z $image_line ]]; then
        return 1
    fi

    echo "$image_line"
}

# Parse image string into registry, repo, tag
# Handles: registry/repo:tag, registry/repo@sha256:digest, repo:tag
parse_image() {
    local image="$1"
    local registry=""
    local repo=""
    local tag=""
    local digest=""

    # Extract registry if present
    if [[ $image =~ ^([a-zA-Z0-9.-]+)/(.+)$ ]]; then
        # Check if it's a registry (contains . or :) or just org/repo
        local first_part="${BASH_REMATCH[1]}"
        local rest="${BASH_REMATCH[2]}"

        if [[ $first_part =~ \. || $first_part =~ : ]]; then
            registry="$first_part"
            image="$rest"
        elif [[ $first_part == "ghcr.io" || $first_part == "docker.io" || $first_part == "quay.io" ]]; then
            registry="$first_part"
            image="$rest"
        else
            # It's org/repo format, reset
            image="${first_part}/${rest}"
        fi
    fi

    # Handle digest format (@sha256:...)
    if [[ $image =~ ^(.+)@(sha256:[a-fA-F0-9]+)$ ]]; then
        repo="${BASH_REMATCH[1]}"
        digest="${BASH_REMATCH[2]}"
        tag=""
    # Handle tag format
    elif [[ $image =~ ^(.+):([^:]+)$ ]]; then
        repo="${BASH_REMATCH[1]}"
        tag="${BASH_REMATCH[2]}"
    else
        repo="$image"
        tag="latest"
    fi

    # Default registry
    if [[ -z $registry ]]; then
        if [[ $repo =~ ^[a-z0-9-]+/[a-z0-9-]+ ]]; then
            registry="docker.io"
        fi
    fi

    echo "${registry}|${repo}|${tag}|${digest}"
}

# Fetch latest tag from Docker Hub API
fetch_docker_hub_latest() {
    local repo="$1"
    local current_tag="$2"

    # Docker Hub API URL
    local api_url="https://hub.docker.com/v2/repositories/${repo}/tags?page_size=100"

    local response
    response=$(curl -sf "$api_url" 2>/dev/null) || return 1

    # Parse tags based on current tag pattern
    local latest_tag
    if [[ $current_tag =~ ^v[0-9] ]]; then
        # v-prefixed version tag - find highest version with same prefix pattern
        latest_tag=$(echo "$response" | jq -r '
            .results[]
            | select(.name | test("^v[0-9]"))
            | .name
        ' 2>/dev/null | sort -V | tail -1)
    elif [[ $current_tag =~ ^[0-9] ]]; then
        # Numeric version tag - find highest version
        latest_tag=$(echo "$response" | jq -r '
            .results[]
            | select(.name | test("^[0-9]"))
            | .name
        ' 2>/dev/null | sort -V | tail -1)
    elif [[ $current_tag =~ ^pg[0-9]+$ ]]; then
        # PostgreSQL version tag like pg16 - find matching pattern
        latest_tag=$(echo "$response" | jq -r '
            .results[]
            | select(.name | test("^pg[0-9]+$"))
            | .name
        ' 2>/dev/null | sort -V | tail -1)
    elif [[ $current_tag == "latest" ]]; then
        # latest tag - find most recent non-latest
        latest_tag=$(echo "$response" | jq -r '
            .results[]
            | select(.name != "latest")
            | .name
        ' 2>/dev/null | head -1)
    else
        # Unknown pattern - find exact match or most recent
        latest_tag=$(echo "$response" | jq -r '
            .results
            | sort_by(.last_updated)
            | reverse
            | .[0]
            | .name
        ' 2>/dev/null)
    fi

    echo "${latest_tag:-unknown}"
}

# Fetch latest tag from GitHub Container Registry
# For GHCR images, we use GitHub Releases API since versioned tags are published there
fetch_ghcr_latest() {
    local repo="$1"
    local current_tag="$2"

    # Try GitHub Releases API first (for projects that publish versioned releases)
    local github_repo="${repo}"
    local release_url="https://api.github.com/repos/${github_repo}/releases/latest"

    local response
    response=$(curl -sf "$release_url" 2>/dev/null) || {
        # Fall back to GHCR API if GitHub releases don't exist
        local token
        token=$(curl -sf "https://ghcr.io/token?scope=repository:${repo}:pull" 2>/dev/null | jq -r '.token') || return 1
        response=$(curl -sfH "Authorization: Bearer $token" "https://ghcr.io/v2/${repo}/tags/list" 2>/dev/null) || return 1

        # Return last tag from GHCR
        echo "$response" | jq -r '.tags[-1]' 2>/dev/null || echo "unknown"
        return
    }

    # Extract version tag from GitHub release
    local latest_tag
    latest_tag=$(echo "$response" | jq -r '.tag_name' 2>/dev/null)

    if [[ -z $latest_tag || $latest_tag == "null" ]]; then
        echo "unknown"
    else
        echo "$latest_tag"
    fi
}

# Get image last updated date from Docker Hub
fetch_docker_hub_image_date() {
    local repo="$1"
    local tag="$2"

    local api_url="https://hub.docker.com/v2/repositories/${repo}/tags/${tag}"
    local response

    response=$(curl -sf "$api_url" 2>/dev/null) || return 1

    echo "$response" | jq -r '.last_updated' 2>/dev/null
}

# Get image last updated date from GHCR (via GitHub Releases API)
fetch_ghcr_image_date() {
    local repo="$1"
    local tag="$2"

    # Try GitHub Releases API for date (for versioned releases)
    local release_url="https://api.github.com/repos/${repo}/releases/tags/${tag}"
    local response

    response=$(curl -sf "$release_url" 2>/dev/null) || {
        # Fall back to GHCR manifest inspection
        local token
        token=$(curl -sf "https://ghcr.io/token?scope=repository:${repo}:pull" 2>/dev/null | jq -r '.token') || return 1
        local manifest
        manifest=$(curl -sfH "Authorization: Bearer $token" \
            -H "Accept: application/vnd.docker.distribution.manifest.v1+json" \
            "https://ghcr.io/v2/${repo}/manifests/${tag}" 2>/dev/null) || return 1

        # Extract history v1 compatibility info for date
        echo "$manifest" | jq -r '.history[0].v1Compatibility' 2>/dev/null | jq -r '.created' 2>/dev/null
        return
    }

    # Get published_at from release
    echo "$response" | jq -r '.published_at' 2>/dev/null
}

# Calculate age in days from ISO date string
calculate_age_days() {
    local date_str="$1"

    if [[ -z $date_str || $date_str == "null" ]]; then
        echo "unknown"
        return
    fi

    local image_epoch current_epoch age_seconds

    # Convert ISO date to epoch
    image_epoch=$(date -d "$date_str" +%s 2>/dev/null) || {
        echo "unknown"
        return
    }
    current_epoch=$(date +%s)

    age_seconds=$((current_epoch - image_epoch))
    echo $((age_seconds / 86400))
}

# Compare version strings (handles v prefix, pg prefix, and semver)
compare_versions() {
    local v1="$1"
    local v2="$2"

    # Exact match first
    if [[ $v1 == "$v2" ]]; then
        echo "equal"
        return
    fi

    # Strip v prefix for comparison
    v1="${v1#v}"
    v2="${v2#v}"

    # Handle pg-prefixed versions (pg16, pg17, etc.)
    if [[ $v1 == pg* && $v2 == pg* ]]; then
        local n1="${v1#pg}"
        local n2="${v2#pg}"
        if [[ $n1 == "$n2" ]]; then
            echo "equal"
        elif [[ $n1 -lt $n2 ]] 2>/dev/null; then
            echo "older"
        else
            echo "newer"
        fi
        return
    fi

    # Handle latest tag
    if [[ $v1 == "latest" ]]; then
        echo "equal" # Assume latest is current
        return
    fi

    # Use sort -V to compare version strings
    if [[ "$(echo -e "$v1\n$v2" | sort -V | head -1)" == "$v1" ]]; then
        echo "older"
    else
        echo "newer"
    fi
}

# Check a single image
check_image() {
    local service="$1"
    local image_spec="$2"

    local parsed registry repo tag digest
    parsed=$(parse_image "$image_spec")
    IFS='|' read -r registry repo tag digest <<<"$parsed"

    local current_version="${tag:-${digest:0:19}}"
    local latest_version="unknown"
    local age_days="unknown"
    local status="unknown"
    local status_code=0

    # Determine if we have a digest (pinned) or tag
    local is_pinned=false
    if [[ -n $digest ]]; then
        is_pinned=true
        current_version="${digest:0:19}..."
    fi

    # Fetch latest version based on registry
    if [[ $registry == "ghcr.io" ]]; then
        latest_version=$(fetch_ghcr_latest "$repo" "$tag") || latest_version="error"

        if [[ $is_pinned == "false" && $latest_version != "error" && $latest_version != "unknown" ]]; then
            local image_date
            image_date=$(fetch_ghcr_image_date "$repo" "$tag") || image_date=""
            age_days=$(calculate_age_days "$image_date")
        fi
    else
        # Docker Hub or other
        latest_version=$(fetch_docker_hub_latest "$repo" "$tag") || latest_version="error"

        if [[ $is_pinned == "false" && $latest_version != "error" && $latest_version != "unknown" ]]; then
            local image_date
            image_date=$(fetch_docker_hub_image_date "$repo" "$tag") || image_date=""
            age_days=$(calculate_age_days "$image_date")
        fi
    fi

    # Determine status
    if [[ $latest_version == "error" ]]; then
        status="error"
        status_code=2
        ANY_ERROR=true
    elif [[ $is_pinned == "true" ]]; then
        status="pinned"
        status_code=0
    elif [[ $latest_version == "unknown" ]]; then
        status="unknown"
        status_code=0
    else
        local comparison
        comparison=$(compare_versions "$tag" "$latest_version")

        if [[ $comparison == "older" ]]; then
            status="outdated"
            status_code=1
            ANY_OUTDATED=true
        else
            status="current"
            status_code=0
        fi

        # Check age
        if [[ $age_days != "unknown" && $age_days -ge $AGE_WARNING_DAYS ]]; then
            if [[ $status != "outdated" ]]; then
                status="old"
            fi
            status_code=1
            ANY_OUTDATED=true
        fi
    fi

    # Format output
    local full_image="${registry:+${registry}/}${repo}"
    local age_display="${age_days}"
    if [[ $age_days != "unknown" ]]; then
        age_display="${age_days}d"
    fi

    # Output based on mode
    if [[ $OUTPUT_JSON == "true" ]]; then
        echo "{\"service\":\"$service\",\"image\":\"$full_image\",\"current\":\"$current_version\",\"latest\":\"$latest_version\",\"age_days\":${age_days:-null},\"status\":\"$status\"}"
    else
        # Table format - skip if quiet and not outdated
        if [[ $QUIET_MODE == "true" && $status_code -eq 0 ]]; then
            return $status_code
        fi

        # Colorize status (use echo -e for proper escape sequence interpretation)
        local status_colored="$status"
        case "$status" in
            current) status_colored="${GREEN}current${NC}" ;;
            outdated) status_colored="${RED}outdated${NC}" ;;
            old) status_colored="${YELLOW}old (${age_days}d)${NC}" ;;
            pinned) status_colored="${CYAN}pinned${NC}" ;;
            error) status_colored="${RED}error${NC}" ;;
        esac

        printf "%-15s %-40s %-20s %-20s %-10s " \
            "$service" \
            "$full_image" \
            "$current_version" \
            "$latest_version" \
            "$age_display"
        echo -e "$status_colored"
    fi

    return $status_code
}

# Main function
main() {
    check_dependencies

    # Define images to check
    declare -A IMAGES=(
        ["openwebui"]='${OPENWEBUI_IMAGE:-ghcr.io/open-webui/open-webui:v0.8.3}'
        ["cliproxyapi"]='${CLIPROXYAPI_IMAGE:-eceasy/cli-proxy-api:v6.8.18}'
        ["postgres"]="pgvector/pgvector:pg16"
        ["jupyter"]='${JUPYTER_IMAGE:-jupyter/scipy-notebook:latest}'
    )

    # Resolve env vars in image specs
    declare -A RESOLVED_IMAGES=()

    # Extract only image-related vars from .env (avoid sourcing entire file which may have issues)
    if [[ -f "${PROJECT_ROOT}/.env" ]]; then
        OPENWEBUI_IMAGE=$(grep -E "^OPENWEBUI_IMAGE=" "${PROJECT_ROOT}/.env" 2>/dev/null | head -1 | cut -d'=' -f2- || true)
        CLIPROXYAPI_IMAGE=$(grep -E "^CLIPROXYAPI_IMAGE=" "${PROJECT_ROOT}/.env" 2>/dev/null | head -1 | cut -d'=' -f2- || true)
        JUPYTER_IMAGE=$(grep -E "^JUPYTER_IMAGE=" "${PROJECT_ROOT}/.env" 2>/dev/null | head -1 | cut -d'=' -f2- || true)
    fi

    RESOLVED_IMAGES["openwebui"]="${OPENWEBUI_IMAGE:-ghcr.io/open-webui/open-webui:v0.8.3}"
    RESOLVED_IMAGES["cliproxyapi"]="${CLIPROXYAPI_IMAGE:-eceasy/cli-proxy-api:v6.8.18}"
    RESOLVED_IMAGES["postgres"]="pgvector/pgvector:pg16"
    RESOLVED_IMAGES["jupyter"]="${JUPYTER_IMAGE:-jupyter/scipy-notebook:latest}"

    if [[ $OUTPUT_JSON == "true" ]]; then
        echo "["
        local first=true
        for service in "${!RESOLVED_IMAGES[@]}"; do
            if [[ $first == "true" ]]; then
                first=false
            else
                echo ","
            fi
            check_image "$service" "${RESOLVED_IMAGES[$service]}" || true
        done
        echo "]"
    else
        print_header "Docker Image Version Check"

        if [[ $QUIET_MODE != "true" ]]; then
            echo ""
            printf "%-15s %-40s %-20s %-20s %-10s %s\n" \
                "SERVICE" "IMAGE" "CURRENT" "LATEST" "AGE" "STATUS"
            printf "%-15s %-40s %-20s %-20s %-10s %s\n" \
                "-------" "-----" "-------" "------" "---" "------"
        fi

        for service in openwebui cliproxyapi postgres jupyter; do
            check_image "$service" "${RESOLVED_IMAGES[$service]}" || true
        done

        echo ""

        if [[ $ANY_OUTDATED == "true" ]]; then
            print_warning "One or more images are outdated or >= ${AGE_WARNING_DAYS} days old"
            print_info "Run 'docker compose pull' to update images"
        elif [[ $ANY_ERROR == "true" ]]; then
            print_warning "Some version checks failed (network issues may be the cause)"
        else
            print_success "All images are up to date"
        fi
    fi

    # Return appropriate exit code
    if [[ $ANY_ERROR == "true" ]]; then
        exit 2
    elif [[ $ANY_OUTDATED == "true" ]]; then
        exit 1
    else
        exit 0
    fi
}

main "$@"
