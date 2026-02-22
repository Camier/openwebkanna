#!/bin/bash

###############################################################################
# env-loader.sh - Shared environment variable loader
#
# Provides load_env_defaults() function for sourcing .env files with:
# - Automatic .env file discovery (SCRIPT_DIR/.env or ./.env)
# - Comment and empty line handling
# - Quoted value stripping (both single and double quotes)
# - Only exports variables not already set (respects existing environment)
# - Safe key validation (alphanumeric + underscore, starting with letter/_)
#
# Usage:
#   source /path/to/lib/env-loader.sh
#   load_env_defaults
###############################################################################

load_env_defaults() {
    local env_file=""
    local line=""
    local key=""
    local value=""

    # Find .env file - prefer SCRIPT_DIR/.env, fallback to ./.env
    if [ -f "$SCRIPT_DIR/.env" ]; then
        env_file="$SCRIPT_DIR/.env"
    elif [ -f ".env" ]; then
        env_file=".env"
    fi

    # Return silently if no .env file found
    if [ -z "$env_file" ]; then
        return 0
    fi

    # Read and parse .env file line by line
    while IFS= read -r line || [ -n "$line" ]; do
        # Strip carriage return (Windows compatibility)
        line="${line%$'\r'}"

        # Trim leading whitespace
        line="${line#"${line%%[![:space:]]*}"}"

        # Skip empty lines
        [ -z "$line" ] && continue

        # Skip comment lines
        [[ $line == \#* ]] && continue

        # Skip lines without assignment
        [[ $line != *=* ]] && continue

        # Extract key and value
        key="${line%%=*}"
        value="${line#*=}"

        # Trim whitespace from key
        key="$(printf "%s" "$key" | tr -d '[:space:]')"

        # Validate key format (must be valid shell variable name)
        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue

        # Skip if variable is already set in environment
        if [ -n "${!key+x}" ]; then
            continue
        fi

        # Strip surrounding double quotes
        if [[ $value == \"*\" ]] && [[ $value == *\" ]]; then
            value="${value:1:${#value}-2}"
        # Strip surrounding single quotes
        elif [[ $value == \'*\' ]] && [[ $value == *\' ]]; then
            value="${value:1:${#value}-2}"
        fi

        # Export the variable
        printf -v "$key" "%s" "$value"
        export "${key?}"
    done <"$env_file"
}
