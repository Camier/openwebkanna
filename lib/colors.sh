#!/bin/bash

###############################################################################
# colors.sh - Shared ANSI color definitions for shell scripts
#
# Usage:
#   source /path/to/lib/colors.sh
#   echo -e "${GREEN}Success${NC}"
#
# shellcheck disable=SC2034
# ^ Variables are intentionally exported for sourcing scripts to use
###############################################################################

# Standard ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No color (reset)
BOLD='\033[1m'
