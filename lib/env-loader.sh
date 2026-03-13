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

        # Trim only outer whitespace from value so common "KEY = value" edits work.
        value="${value#"${value%%[![:space:]]*}"}"
        value="${value%"${value##*[![:space:]]}"}"

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

###############################################################################
# read_dot007_key
# Read a key from ~/.007, supporting KEY=value, export KEY=value, and KEY: value.
#
# Arguments:
#   $1 - Key name to read
#
# Output:
#   Value to stdout
#
# Returns:
#   0 if the key was found, 1 otherwise
###############################################################################
read_dot007_key() {
    local wanted="$1"
    local dotfile="$HOME/.007"

    [ -f "$dotfile" ] || return 1

    awk -v k="$wanted" '
        /^[[:space:]]*#/ { next }
        /^[[:space:]]*$/ { next }
        {
            line=$0
            sub(/^[[:space:]]*export[[:space:]]+/, "", line)
            sep = index(line, "=") ? "=" : (index(line, ":") ? ":" : "")
            if (sep == "") next
            split(line, parts, sep)
            key=parts[1]
            sub(/^[[:space:]]+/, "", key); sub(/[[:space:]]+$/, "", key)
            if (key != k) next
            val = substr(line, length(parts[1]) + 2)
            sub(/^[[:space:]]+/, "", val); sub(/[[:space:]]+$/, "", val)
            if ((val ~ /^".*"$/) || (val ~ /^\x27.*\x27$/)) {
                val = substr(val, 2, length(val)-2)
            }
            print val
            exit 0
        }
    ' "$dotfile"
}
