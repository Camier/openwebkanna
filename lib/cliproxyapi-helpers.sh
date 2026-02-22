#!/bin/bash

###############################################################################
# cliproxyapi-helpers.sh - CLIProxyAPI configuration utilities
#
# Provides functions for parsing and manipulating CLIProxyAPI config.yaml files.
#
# Functions:
#   extract_first_config_api_key() - Get first API key from config
#   rewrite_first_config_api_key() - Replace first API key in config
#   extract_config_aliases()       - List all model aliases from config
#   ensure_cliproxyapi_config_api_key_alignment() - Sync env key with config
#
# Usage:
#   source /path/to/lib/cliproxyapi-helpers.sh
#   key=$(extract_first_config_api_key "$CONFIG_FILE")
###############################################################################

###############################################################################
# extract_first_config_api_key
# Extract the first API key from a CLIProxyAPI config.yaml file.
#
# Arguments:
#   $1 - Path to config.yaml file
#
# Output:
#   First API key value (stdout), empty if not found
#
# Returns:
#   0 on success (even if no key found), non-zero on file error
#
# Example:
#   key=$(extract_first_config_api_key "./cliproxyapi/config.yaml")
###############################################################################
extract_first_config_api_key() {
    local config_file="$1"

    [ -f "$config_file" ] || return 0

    awk '
        /^[[:space:]]*api-keys:[[:space:]]*$/ {
            in_keys = 1
            next
        }
        in_keys == 1 && /^[[:space:]]*-[[:space:]]*/ {
            line = $0
            sub(/^[[:space:]]*-[[:space:]]*/, "", line)
            sub(/[[:space:]]+#.*/, "", line)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
            if ((line ~ /^".*"$/) || (line ~ /^'"'"'.*'"'"'$/)) {
                line = substr(line, 2, length(line) - 2)
            }
            print line
            exit
        }
        in_keys == 1 && $0 !~ /^[[:space:]]/ {
            in_keys = 0
        }
    ' "$config_file"
}

###############################################################################
# rewrite_first_config_api_key
# Replace the first API key in a CLIProxyAPI config.yaml file.
#
# Arguments:
#   $1 - Path to config.yaml file
#   $2 - New API key value
#
# Returns:
#   0 on success, 1 on awk error, 2 if no key was found to replace
#
# Example:
#   rewrite_first_config_api_key "./cliproxyapi/config.yaml" "new-secret-key"
###############################################################################
rewrite_first_config_api_key() {
    local config_file="$1"
    local new_key="$2"
    local tmp_file=""

    tmp_file="$(mktemp)"
    if ! awk -v new_key="$new_key" '
        BEGIN {
            in_keys = 0
            replaced = 0
        }
        {
            if ($0 ~ /^[[:space:]]*api-keys:[[:space:]]*$/) {
                in_keys = 1
                print
                next
            }

            if (in_keys == 1 && replaced == 0 && $0 ~ /^[[:space:]]*-[[:space:]]*/) {
                match($0, /^[[:space:]]*/)
                indent = substr($0, RSTART, RLENGTH)
                key = new_key
                gsub(/"/, "\\\"", key)
                print indent "- \"" key "\""
                replaced = 1
                in_keys = 0
                next
            }

            if (in_keys == 1 && $0 !~ /^[[:space:]]/) {
                in_keys = 0
            }

            print
        }
        END {
            if (replaced == 0) {
                exit 2
            }
        }
    ' "$config_file" >"$tmp_file"; then
        rm -f "$tmp_file"
        return 1
    fi

    mv "$tmp_file" "$config_file"
}

###############################################################################
# extract_config_aliases
# Extract all model aliases from a CLIProxyAPI config.yaml file.
#
# Arguments:
#   $1 - Path to config.yaml file
#
# Output:
#   List of alias values, one per line, sorted and unique (stdout)
#
# Returns:
#   0 on success
#
# Example:
#   aliases=$(extract_config_aliases "./cliproxyapi/config.yaml")
###############################################################################
extract_config_aliases() {
    local config_file="$1"

    [ -f "$config_file" ] || return 0

    awk '
        /^[[:space:]]*alias:[[:space:]]*/ {
            line = $0
            sub(/^[[:space:]]*alias:[[:space:]]*/, "", line)
            sub(/[[:space:]]+#.*/, "", line)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
            if ((line ~ /^".*"$/) || (line ~ /^'"'"'.*'"'"'$/)) {
                line = substr(line, 2, length(line) - 2)
            }
            if (line != "") {
                print line
            }
        }
    ' "$config_file" | sort -u
}

###############################################################################
# ensure_cliproxyapi_config_api_key_alignment
# Ensure CLIProxyAPI config file API key matches the environment variable.
#
# This function checks if the first API key in the config matches
# CLIPROXYAPI_API_KEY. If not, it can either sync the config or fail
# depending on configuration.
#
# Required globals (set before calling):
#   CLIPROXYAPI_CONFIG           - Path to config.yaml
#   CLIPROXYAPI_API_KEY          - Expected API key (optional)
#   CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH - Fail on mismatch (default: true)
#   CLIPROXYAPI_SYNC_CONFIG_API_KEY - Sync config on mismatch (default: true)
#
# Returns:
#   0 on success or no enforcement needed
#   1 on failure (config not found, mismatch when enforced, sync failed)
#
# Example:
#   CLIPROXYAPI_CONFIG="./cliproxyapi/config.yaml"
#   CLIPROXYAPI_API_KEY="my-key"
#   ensure_cliproxyapi_config_api_key_alignment
###############################################################################
ensure_cliproxyapi_config_api_key_alignment() {
    local config_api_key=""

    if [ ! -f "$CLIPROXYAPI_CONFIG" ]; then
        print_error "CLIProxyAPI config file not found: $CLIPROXYAPI_CONFIG"
        return 1
    fi

    if [ -z "${CLIPROXYAPI_API_KEY:-}" ]; then
        return 0
    fi

    config_api_key="$(extract_first_config_api_key "$CLIPROXYAPI_CONFIG" || true)"
    if [ -z "$config_api_key" ]; then
        print_warning "No api-keys entry found in $CLIPROXYAPI_CONFIG"
        if is_true "${CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH:-true}"; then
            print_error "Cannot validate authenticated readiness without api-keys in config"
            return 1
        fi
        return 0
    fi

    if [ "$config_api_key" = "$CLIPROXYAPI_API_KEY" ]; then
        return 0
    fi

    if is_true "${CLIPROXYAPI_SYNC_CONFIG_API_KEY:-true}"; then
        if rewrite_first_config_api_key "$CLIPROXYAPI_CONFIG" "$CLIPROXYAPI_API_KEY"; then
            print_warning "Synced first api-keys entry in $CLIPROXYAPI_CONFIG to match CLIPROXYAPI_API_KEY"
            return 0
        fi
        print_error "Failed to sync api-keys entry in $CLIPROXYAPI_CONFIG"
        return 1
    fi

    print_error "CLIPROXYAPI_API_KEY does not match first api-keys entry in $CLIPROXYAPI_CONFIG"
    print_info "Set CLIPROXYAPI_SYNC_CONFIG_API_KEY=true or update config.yaml api-keys to match OPENAI_API_KEY"
    if is_true "${CLIPROXYAPI_ENFORCE_CONFIG_API_KEY_MATCH:-true}"; then
        return 1
    fi
    return 0
}

###############################################################################
# resolved_expected_aliases
# Get list of expected CLIProxyAPI model aliases.
#
# Uses CLIPROXYAPI_EXPECT_ALIASES if set, otherwise extracts from config.
#
# Required globals (set before calling):
#   CLIPROXYAPI_EXPECT_ALIASES - Space-separated list of expected aliases (optional)
#   CLIPROXYAPI_CONFIG         - Path to config.yaml (used if EXPECT_ALIASES empty)
#
# Output:
#   List of expected aliases, one per line (stdout)
#
# Example:
#   for alias in $(resolved_expected_aliases); do ...
###############################################################################
resolved_expected_aliases() {
    if [ -n "${CLIPROXYAPI_EXPECT_ALIASES:-}" ]; then
        printf "%s\n" $CLIPROXYAPI_EXPECT_ALIASES
        return 0
    fi

    extract_config_aliases "${CLIPROXYAPI_CONFIG:-./cliproxyapi/config.yaml}"
}
