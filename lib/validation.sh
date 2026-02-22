#!/bin/bash

###############################################################################
# validation.sh - Shared validation and path utility functions
#
# Provides common validation and normalization utilities for shell scripts.
#
# Functions:
#   is_true()           - Check if a value represents boolean true
#   command_exists()    - Check if a command is available in PATH
#   normalize_base_url() - Remove trailing slash from URL
#   normalize_path()    - Ensure path starts with / and remove trailing slash
#   resolve_path()      - Resolve relative paths to absolute using SCRIPT_DIR
#
# Usage:
#   source /path/to/lib/validation.sh
#   if is_true "$SOME_FLAG"; then ...
#   URL="$(normalize_base_url "$URL")"
###############################################################################

###############################################################################
# is_true
# Check if a value represents boolean true.
#
# Arguments:
#   $1 - Value to check
#
# Returns:
#   0 (true) if value is "true", "1", "yes", or "on" (case-insensitive)
#   1 (false) otherwise
#
# Example:
#   if is_true "$ENABLED"; then
#       echo "Feature is enabled"
#   fi
###############################################################################
is_true() {
    local value
    value="$(printf "%s" "$1" | tr '[:upper:]' '[:lower:]')"
    [ "$value" = "true" ] || [ "$value" = "1" ] || [ "$value" = "yes" ] || [ "$value" = "on" ]
}

###############################################################################
# command_exists
# Check if a command is available in PATH.
#
# Arguments:
#   $1 - Command name to check
#
# Returns:
#   0 if command exists, 1 otherwise
#
# Example:
#   if ! command_exists curl; then
#       echo "curl is required" >&2
#       exit 1
#   fi
###############################################################################
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

###############################################################################
# normalize_base_url
# Remove trailing slash from a base URL.
#
# Arguments:
#   $1 - URL to normalize
#
# Output:
#   URL without trailing slash (stdout)
#
# Example:
#   BASE_URL="$(normalize_base_url "$BASE_URL")"
#   # "http://example.com/" -> "http://example.com"
###############################################################################
normalize_base_url() {
    local value="$1"
    printf "%s" "${value%/}"
}

###############################################################################
# normalize_path
# Ensure a path starts with / and has no trailing slash.
#
# Arguments:
#   $1 - Path to normalize
#
# Output:
#   Normalized path (stdout), defaults to "/" if empty
#
# Example:
#   API_PATH="$(normalize_path "$API_PATH")"
#   # "v1/models" -> "/v1/models"
#   # "/v1/models/" -> "/v1/models"
###############################################################################
normalize_path() {
    local value="$1"
    if [ -z "$value" ]; then
        printf "/"
        return
    fi
    if [[ $value != /* ]]; then
        value="/$value"
    fi
    printf "%s" "${value%/}"
}

###############################################################################
# resolve_path
# Resolve a path to absolute, using SCRIPT_DIR for relative paths.
#
# Requires:
#   SCRIPT_DIR must be set before calling (typically at script start)
#
# Arguments:
#   $1 - Path to resolve (absolute or relative)
#
# Output:
#   Absolute path (stdout)
#
# Example:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   CONFIG_FILE="$(resolve_path "$CONFIG_FILE")"
#   # "./config.yaml" -> "/path/to/script/config.yaml"
###############################################################################
resolve_path() {
    local candidate="$1"
    if [[ $candidate == /* ]]; then
        printf "%s" "$candidate"
    else
        printf "%s/%s" "$SCRIPT_DIR" "$candidate"
    fi
}
