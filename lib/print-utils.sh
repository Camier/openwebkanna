#!/bin/bash

###############################################################################
# print-utils.sh - Shared print functions for shell scripts
#
# Usage:
#   source /path/to/lib/print-utils.sh
#   print_header "My Script Title"
#   print_step "Doing something"
#   print_success "It worked"
#   print_error "Something failed"
###############################################################################

# Get the directory where this script is located (resolves symlinks)
_PRINT_UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source color definitions
source "${_PRINT_UTILS_DIR}/colors.sh"

# Print a fancy bordered header
# Usage: print_header "Title Text"
print_header() {
    local title="$1"
    local width=60
    local padding
    local left_pad=""
    local right_pad=""
    padding=$(((width - ${#title}) / 2))

    # Calculate padding for centering
    for ((i = 0; i < padding; i++)); do
        left_pad+=" "
    done
    right_pad="$left_pad"
    # Add extra space if title length is odd
    if [ $((${#title} % 2)) -ne 0 ]; then
        right_pad+=" "
    fi

    echo -e "${CYAN}${BOLD}"
    echo "+------------------------------------------------------------+"
    echo "|${left_pad}${title}${right_pad}|"
    echo "+------------------------------------------------------------+"
    echo -e "${NC}"
}

# Print a step with arrow prefix
# Usage: print_step "Step description"
print_step() {
    echo -e "\n${BLUE}${BOLD}>> $1${NC}"
}

# Print a green checkmark success message
# Usage: print_success "Operation completed"
print_success() {
    echo -e "${GREEN}[OK] $1${NC}"
}

# Print a red error message to stderr
# Usage: print_error "Something went wrong"
print_error() {
    echo -e "${RED}[ERROR] $1${NC}" >&2
}

# Print a yellow warning message
# Usage: print_warning "This is a warning"
print_warning() {
    echo -e "${YELLOW}[WARN] $1${NC}"
}

# Print a magenta info message
# Usage: print_info "Informational note"
print_info() {
    echo -e "${MAGENTA}[INFO] $1${NC}"
}

# Print a section divider
# Usage: print_section "Section Name"
print_section() {
    echo -e "\n${BOLD}--- $1 ---${NC}\n"
}
