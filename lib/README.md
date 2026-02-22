# OpenWebUI Shell Library Reference Guide

This document provides comprehensive documentation for the shared shell libraries used across the OpenWebUI deployment project.

## Overview

The `lib/` directory contains reusable shell script components that enforce consistent coding patterns across all operational scripts. These libraries provide:

- **Standardized output formatting** with color support
- **Docker Compose abstraction** for version compatibility
- **Environment file parsing** with security validation
- **Common validation utilities** for paths and booleans
- **Python plugin auditing** for OpenWebUI

---

## Quick Start

```bash
#!/bin/bash
set -e

# Standard header pattern
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source required libraries
source "${SCRIPT_DIR}/lib/print-utils.sh"
source "${SCRIPT_DIR}/lib/docker-helpers.sh"
source "${SCRIPT_DIR}/lib/validation.sh"
source "${SCRIPT_DIR}/lib/env-loader.sh"

# Load environment defaults
load_env_defaults

# Your script logic here
print_header "My Script"
print_step "Doing something"
print_success "Done"
```

---

## Library API Reference

### 1. colors.sh - ANSI Color Definitions

**File:** `lib/colors.sh`
**Purpose:** Standard ANSI color codes for terminal output

#### Variables

| Variable | Value | Usage |
|----------|-------|-------|
| `RED` | `\033[0;31m` | Error messages |
| `GREEN` | `\033[0;32m` | Success messages |
| `YELLOW` | `\033[1;33m` | Warning messages |
| `BLUE` | `\033[0;34m` | Step indicators |
| `CYAN` | `\033[0;36m` | Headers and info |
| `NC` | `\033[0m` | Reset color |
| `BOLD` | `\033[1m` | Emphasis |

#### Usage

```bash
source "${SCRIPT_DIR}/lib/colors.sh"
echo -e "${GREEN}Success${NC}"
echo -e "${RED}Error${NC}" >&2
```

#### Notes
- Always use `echo -e` to enable escape sequence interpretation
- Always end colored output with `${NC}` to reset
- Error messages should be redirected to `>&2`

---

### 2. print-utils.sh - Standard Print Functions

**File:** `lib/print-utils.sh`
**Purpose:** Consistent formatted output functions
**Dependencies:** `colors.sh` (auto-sourced)

#### Functions

##### `print_header "Title Text"`
Prints a fancy bordered header, centered in a 60-character box.

```bash
print_header "Deployment Script"
# Output:
# +------------------------------------------------------------+
# |                   Deployment Script                        |
# +------------------------------------------------------------+
```

##### `print_step "Step description"`
Prints a step indicator with blue arrow prefix.

```bash
print_step "Starting services"
# Output:
# >> Starting services
```

##### `print_success "Message"`
Prints a green checkmark-style success message.

```bash
print_success "Services started"
# Output: [OK] Services started
```

##### `print_error "Message"`
Prints a red error message to stderr.

```bash
print_error "Failed to start"
# Output (to stderr): [ERROR] Failed to start
```

##### `print_warning "Message"`
Prints a yellow warning message.

```bash
print_warning "Using fallback"
# Output: [WARN] Using fallback
```

##### `print_info "Message"`
Prints a cyan informational message.

```bash
print_info "Processing..."
# Output: [INFO] Processing...
```

##### `print_section "Section Name"`
Prints a section divider.

```bash
print_section "Configuration"
# Output:
# --- Configuration ---
```

#### Usage

```bash
source "${SCRIPT_DIR}/lib/print-utils.sh"

print_header "My Script"
print_step "Initializing"
print_info "Loading config"
print_success "Ready"
```

---

### 3. docker-helpers.sh - Docker Compose Abstraction

**File:** `lib/docker-helpers.sh`
**Purpose:** Handle Docker Compose v1 vs v2 differences automatically

#### Global Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `COMPOSE_CMD` | `()` | Array holding detected compose command |
| `COMPOSE_FILE` | `docker-compose.yml` | Path to compose file |

#### Functions

##### `init_compose_cmd()`
Detects and initializes the `COMPOSE_CMD` array. Prefers `docker compose` (v2 plugin) over `docker-compose` (v1 standalone).

**Returns:**
- `0` on success
- `1` if Docker Compose is not available

```bash
if init_compose_cmd; then
    echo "Using: ${COMPOSE_CMD[*]}"
fi
```

##### `docker_compose [args...]`
Wrapper that invokes the detected compose command with the configured compose file.

**Arguments:** All arguments passed through to docker compose

**Returns:** Exit code from docker compose, or `1` if init fails

```bash
# Automatically detects and uses correct compose command
docker_compose up -d
docker_compose logs -f openwebui
docker_compose exec -T openwebui python --version
```

#### Usage Pattern

```bash
source "${SCRIPT_DIR}/lib/docker-helpers.sh"

# Set custom compose file if needed
COMPOSE_FILE="docker-compose.prod.yml"

# Use wrapper (init_compose_cmd called automatically)
docker_compose up -d
docker_compose ps
```

---

### 4. env-loader.sh - Environment File Parser

**File:** `lib/env-loader.sh`
**Purpose:** Safe .env file loading with validation

#### Functions

##### `load_env_defaults()`
Parses and loads environment variables from `.env` file with:
- Automatic file discovery (`SCRIPT_DIR/.env` or `./.env`)
- Comment and empty line filtering
- Quote stripping (single and double)
- Variable name validation (alphanumeric + underscore)
- Environment variable precedence (existing vars not overwritten)
- Windows line ending compatibility

**Requires:** `SCRIPT_DIR` to be set before calling

**Returns:** Always `0` (silently succeeds if no .env found)

```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/env-loader.sh"
load_env_defaults
```

#### Parsing Behavior

| Input | Result |
|-------|--------|
| `KEY=value` | `KEY=value` |
| `KEY="value"` | `KEY=value` (quotes stripped) |
| `KEY='value'` | `KEY=value` (quotes stripped) |
| `# comment` | Ignored |
| `KEY=value # comment` | Not supported (inline comments) |
| `  KEY  =  value  ` | `KEY=value` (whitespace trimmed) |
| `123_KEY=value` | Ignored (invalid name) |

#### Usage Pattern

```bash
#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/env-loader.sh"

# Load defaults from .env
load_env_defaults

# Override with script-specific defaults if needed
MY_VAR="${MY_VAR:-default_value}"
```

---

### 5. validation.sh - Validation Utilities

**File:** `lib/validation.sh`
**Purpose:** Common validation and normalization functions

#### Functions

##### `is_true "value"`
Checks if a value represents boolean true.

**Arguments:**
- `$1` - Value to check

**Returns:**
- `0` (true) if value is: `true`, `1`, `yes`, `on` (case-insensitive)
- `1` (false) otherwise

```bash
if is_true "$ENABLED"; then
    echo "Feature is enabled"
fi
```

##### `command_exists "command"`
Checks if a command is available in PATH.

**Arguments:**
- `$1` - Command name to check

**Returns:**
- `0` if command exists
- `1` otherwise

```bash
if ! command_exists curl; then
    echo "curl is required" >&2
    exit 1
fi
```

##### `normalize_base_url "url"`
Removes trailing slash from a URL.

**Arguments:**
- `$1` - URL to normalize

**Output:** URL without trailing slash (stdout)

```bash
BASE_URL="$(normalize_base_url "$BASE_URL")"
# "http://example.com/" -> "http://example.com"
```

##### `normalize_path "path"`
Ensures a path starts with `/` and has no trailing slash.

**Arguments:**
- `$1` - Path to normalize

**Output:** Normalized path (stdout), defaults to "/" if empty

```bash
API_PATH="$(normalize_path "$API_PATH")"
# "v1/models" -> "/v1/models"
# "/v1/models/" -> "/v1/models"
```

##### `resolve_path "path"`
Resolves a relative path to absolute using `SCRIPT_DIR`.

**Requires:** `SCRIPT_DIR` must be set before calling

**Arguments:**
- `$1` - Path to resolve (absolute or relative)

**Output:** Absolute path (stdout)

```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$(resolve_path "$CONFIG_FILE")"
# "./config.yaml" -> "/path/to/script/config.yaml"
```

#### Usage Pattern

```bash
source "${SCRIPT_DIR}/lib/validation.sh"

# Check dependencies
for cmd in curl jq docker; do
    if ! command_exists "$cmd"; then
        echo "Missing: $cmd" >&2
        exit 1
    fi
done

# Normalize URLs
API_BASE="$(normalize_base_url "$API_BASE")"
```

---

### 6. audit_plugins.py - Python Plugin Auditor

**File:** `lib/audit_plugins.py`
**Purpose:** Audit OpenWebUI tool/function Python code stored in webui.db

#### Overview
Verifies that Python code stored in the OpenWebUI database:
1. Compiles without syntax errors
2. Can import required dependencies (optional)

#### Functions

##### `find_missing_imports(source, import_check=True)`
Finds modules that are imported but not available in the Python environment.

**Arguments:**
- `source` (str): Python source code
- `import_check` (bool): Whether to perform import checking

**Returns:** List of tuples `(module_name, line_number)` for missing imports

##### `audit_plugins(focus, import_check, db_path)`
Audits plugin code from webui.db.

**Arguments:**
- `focus` (str): Scope - `"all"`, `"tool"`, or `"function"`
- `import_check` (bool): Whether to check import dependencies
- `db_path` (str): Path to webui.db file

**Returns:** Tuple `(checked_count, issues_list)`

#### CLI Usage

```bash
# Audit all plugins (tools and functions)
python3 lib/audit_plugins.py all

# Audit only tools
python3 lib/audit_plugins.py tool

# Audit without import checking
python3 lib/audit_plugins.py all false

# Use custom database path
python3 lib/audit_plugins.py all --db-path /path/to/webui.db
```

#### Output Format

```
AUDIT_CHECKED=5
AUDIT_ISSUES=1
ISSUE|tool|123|My Tool|15|0|invalid syntax|def foo(
```

**Issue fields (pipe-separated):**
1. `ISSUE` - Fixed prefix
2. Table name (`tool` or `function`)
3. Row ID
4. Plugin name
5. Line number
6. Column offset
7. Error message
8. Code snippet

**Exit codes:**
- `0` - No issues found
- `1` - One or more issues found

---

## Library Dependencies

```
colors.sh
    ↑
print-utils.sh

docker-helpers.sh (standalone)

env-loader.sh (standalone)

validation.sh (standalone)

audit_plugins.py (standalone, Python)
```

**Important:** `print-utils.sh` automatically sources `colors.sh`. When using `print-utils.sh`, you don't need to source `colors.sh` separately.

---

## Usage Patterns

### Pattern 1: Root-Level Scripts (Full Library Stack)

For scripts in the project root that need full functionality:

```bash
#!/bin/bash

###############################################################################
# Script Name - Brief description
###############################################################################

set -e  # Exit on error

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source libraries
source "${SCRIPT_DIR}/lib/print-utils.sh"
source "${SCRIPT_DIR}/lib/docker-helpers.sh"
source "${SCRIPT_DIR}/lib/validation.sh"
source "${SCRIPT_DIR}/lib/env-loader.sh"

# Load environment
load_env_defaults
cd "$SCRIPT_DIR"

# Configuration with defaults
MY_PORT="${MY_PORT:-8080}"
VERBOSE="${VERBOSE:-false}"

# Main logic
main() {
    print_header "Script Title"

    # Check dependencies
    for cmd in curl jq docker; do
        if ! command_exists "$cmd"; then
            print_error "Missing dependency: $cmd"
            exit 1
        fi
    done

    print_step "Starting operations"

    # Use docker compose wrapper
    docker_compose up -d

    if is_true "$VERBOSE"; then
        print_info "Verbose mode enabled"
    fi

    print_success "Complete"
}

main "$@"
```

### Pattern 2: Scripts in `scripts/` Directory

For scripts in the `scripts/` subdirectory:

```bash
#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LIB_DIR="${PROJECT_ROOT}/lib"

# Source with existence checks for portability
if [[ -f "${LIB_DIR}/print-utils.sh" ]]; then
    source "${LIB_DIR}/print-utils.sh"
fi

# Fallback definitions if libraries unavailable
if ! declare -f print_header >/dev/null 2>&1; then
    print_header() { echo "=== $1 ==="; }
fi
```

### Pattern 3: Minimal Scripts (Single Library)

For simple scripts needing only one library:

```bash
#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/validation.sh"

# Just use what you need
if ! command_exists docker; then
    echo "Docker required" >&2
    exit 1
fi
```

---

## Best Practices

### 1. Script Header Template

Every script should start with:

```bash
#!/bin/bash

###############################################################################
# Script Name - One-line description
#
# Usage:
#   ./script-name.sh [options]
#
# Options:
#   -h, --help    Show help message
###############################################################################

set -e  # Exit on error (add -u and -o pipefail for stricter mode)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
```

### 2. Error Handling

```bash
# Always redirect errors to stderr
print_error "Something went wrong"  # Already does >&2
echo "Error message" >&2             # Manual redirect

# Check exit codes
if ! docker_compose up -d; then
    print_error "Failed to start services"
    exit 1
fi

# Use trap for cleanup
cleanup() {
    print_info "Cleaning up..."
}
trap cleanup EXIT
```

### 3. Color Usage

```bash
# Use library functions (preferred)
print_success "Done"
print_error "Failed"

# Or use colors directly when needed
echo -e "${GREEN}Custom success format${NC}"
echo -e "${YELLOW}Warning: ${NC}details here"

# Always reset colors
printf "${CYAN}Loading...${NC}\n"  # printf doesn't need -e
```

### 4. Environment Variable Defaults

```bash
# Use :- for defaults, := for setting if unset
PORT="${PORT:-8080}"           # Default if unset or empty
DEBUG="${DEBUG:-false}"        # Boolean default
NAME="${NAME:=default}"        # Set and export default

# Check boolean flags
if is_true "$DEBUG"; then
    VERBOSE=true
fi
```

### 5. Path Handling

```bash
# Always use SCRIPT_DIR for relative paths
CONFIG="${SCRIPT_DIR}/config.yaml"

# Use resolve_path for user-provided paths
USER_PATH="$(resolve_path "$1")"

# Normalize URLs and API paths
API_URL="$(normalize_base_url "$API_URL")"
API_PATH="$(normalize_path "$API_PATH")"
```

### 6. Library Loading with Fallbacks

For scripts that may run outside the project:

```bash
# Try to load libraries, provide fallbacks
if [[ -f "${SCRIPT_DIR}/lib/print-utils.sh" ]]; then
    source "${SCRIPT_DIR}/lib/print-utils.sh"
else
    # Minimal fallback implementations
    print_error() { echo "ERROR: $1" >&2; }
    print_success() { echo "OK: $1"; }
fi
```

---

## Extending the Libraries

### Adding a New Print Function

Edit `lib/print-utils.sh`:

```bash
# Print a debug message (only when DEBUG=true)
print_debug() {
    if [[ "${DEBUG:-false}" == "true" ]]; then
        echo -e "${CYAN}[DEBUG] $1${NC}" >&2
    fi
}
```

### Adding a New Validation Function

Edit `lib/validation.sh`:

```bash
###############################################################################
# is_integer
# Check if a value is a valid integer.
#
# Arguments:
#   $1 - Value to check
#
# Returns:
#   0 if integer, 1 otherwise
###############################################################################
is_integer() {
    [[ "$1" =~ ^-?[0-9]+$ ]]
}
```

### Adding a New Docker Helper

Edit `lib/docker-helpers.sh`:

```bash
###############################################################################
# docker_health_check
# Check if a container is healthy.
#
# Arguments:
#   $1 - Container name
#
# Returns:
#   0 if healthy, 1 otherwise
###############################################################################
docker_health_check() {
    local container="$1"
    local status
    status=$(docker inspect --format='{{.State.Health.Status}}' "$container" 2>/dev/null)
    [[ "$status" == "healthy" ]]
}
```

### Creating a New Library

1. Create file: `lib/my-library.sh`
2. Add header with description and usage
3. Define functions with proper documentation
4. Follow existing patterns for consistency

Example template:

```bash
#!/bin/bash

###############################################################################
# my-library.sh - Brief description
#
# Usage:
#   source /path/to/lib/my-library.sh
#   my_function "argument"
###############################################################################

# Dependencies (if any)
_MY_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_MY_LIB_DIR}/colors.sh"

###############################################################################
# my_function
# Description of what this does.
#
# Arguments:
#   $1 - First argument description
#   $2 - Second argument description
#
# Returns:
#   0 on success, 1 on failure
#
# Example:
#   my_function "arg1" "arg2"
###############################################################################
my_function() {
    local arg1="$1"
    local arg2="$2"

    # Implementation
    return 0
}
```

---

## Migration Guide

### Migrating Scripts to Use Libraries

**Before (inline definitions):**

```bash
#!/bin/bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

echo -e "${GREEN}Starting...${NC}"
```

**After (using libraries):**

```bash
#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/print-utils.sh"

print_step "Starting..."
```

### Replacing Inline load_env_defaults

Many root scripts have inline copies of `load_env_defaults()`. To use the library version:

1. Remove the inline function definition
2. Add `source "${SCRIPT_DIR}/lib/env-loader.sh"`
3. Keep the `load_env_defaults` call

---

## Testing Library Changes

Always test library changes with the verification script:

```bash
./verify-scripts.sh
```

Or test individual components:

```bash
# Test colors
source lib/colors.sh
echo -e "${GREEN}Green${NC} ${RED}Red${NC}"

# Test print utils
source lib/print-utils.sh
print_header "Test"
print_success "OK"

# Test docker helpers
source lib/docker-helpers.sh
docker_compose ps

# Test validation
source lib/validation.sh
is_true "yes" && echo "Yes is true"
command_exists "bash" && echo "bash exists"
```

---

## Summary

| Library | Use When | Key Functions |
|---------|----------|---------------|
| `colors.sh` | Direct color control | Color variables |
| `print-utils.sh` | Formatted output | `print_*` functions |
| `docker-helpers.sh` | Docker Compose ops | `docker_compose()` |
| `env-loader.sh` | .env file loading | `load_env_defaults()` |
| `validation.sh` | Input validation | `is_true()`, `command_exists()` |
| `audit_plugins.py` | Plugin auditing | CLI tool for webui.db |

For questions or issues, refer to existing scripts that use these libraries for working examples.
