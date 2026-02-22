#!/bin/bash

###############################################################################
# init.sh - Single-source library loader for OpenWebUI deployment scripts
#
# Sources all library modules in correct dependency order:
#   1. colors.sh        - ANSI color definitions
#   2. validation.sh    - Boolean checks, path normalization
#   3. print-utils.sh   - Formatted output functions (depends on colors.sh)
#   4. env-loader.sh    - Environment variable loading
#   5. docker-helpers.sh - Docker Compose utilities
#   6. cliproxyapi-helpers.sh - CLIProxyAPI config parsing
#   7. http-helpers.sh  - HTTP requests and model extraction
#
# Usage:
#   #!/bin/bash
#   set -e
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SCRIPT_DIR}/lib/init.sh"
#   load_env_defaults
#   cd "$SCRIPT_DIR"
#
# This replaces inline copies of:
#   - load_env_defaults() function
#   - Color variable definitions (RED, GREEN, etc.)
#   - print_header, print_step, print_success, print_error, etc.
#   - init_compose_cmd, docker_compose functions
#   - is_true function
#   - CLIProxyAPI config parsing functions
#   - HTTP and model extraction helpers
###############################################################################

# Resolve lib directory
_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source libraries in dependency order
source "${_LIB_DIR}/colors.sh"
source "${_LIB_DIR}/validation.sh"
source "${_LIB_DIR}/print-utils.sh"
source "${_LIB_DIR}/env-loader.sh"
source "${_LIB_DIR}/docker-helpers.sh"
source "${_LIB_DIR}/cliproxyapi-helpers.sh"
source "${_LIB_DIR}/http-helpers.sh"

# Unset temporary variable
unset _LIB_DIR
