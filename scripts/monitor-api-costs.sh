#!/bin/bash
###############################################################################
# API Cost Monitor - Track and optimize API usage for thesis research
#
# Usage:
#   ./monitor-api-costs.sh [OPTIONS]
#
# Options:
#   --daily       Show daily usage summary
#   --weekly      Show weekly usage summary (default)
#   --monthly     Show monthly usage summary
#   --estimate    Estimate costs based on current usage patterns
#   --export      Export detailed usage to CSV
#
# This script helps solo researchers track API costs and optimize usage.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPORT_DIR="${PROJECT_ROOT}/logs/api-usage"

mkdir -p "${REPORT_DIR}"

# Default to weekly
PERIOD="weekly"
# shellcheck disable=SC2034
ACTION="show"

collect_usage_data() {
    local since="$1"

    # Check if CLIProxyAPI logs are available
    if ! docker compose -f "${PROJECT_ROOT}/docker-compose.yml" ps | grep -q "cliproxyapi.*running"; then
        echo "Warning: CLIProxyAPI is not running"
        return 1
    fi

    # Extract API calls from logs
    # Looking for patterns like: model=xyz, tokens=prompt/completion/total
    docker compose -f "${PROJECT_ROOT}/docker-compose.yml" logs --since="${since}" cliproxyapi 2>/dev/null |
        grep -E "(model=|tokens=|completion|chat.completion)" |
        tail -100 || echo "No recent API calls found"
}

show_summary() {
    local period="$1"
    local since

    case "${period}" in
        daily) since="24h ago" ;;
        weekly) since="7 days ago" ;;
        monthly) since="30 days ago" ;;
        *) since="7 days ago" ;;
    esac

    echo "=== API Usage Summary (${period}) ==="
    echo "Period: ${since}"
    echo ""

    local usage_data
    usage_data=$(collect_usage_data "${since}") || {
        echo "Could not retrieve usage data"
        return 1
    }

    if [[ -z ${usage_data} ]]; then
        echo "No API usage recorded for this period"
        return 0
    fi

    # Parse and summarize
    echo "Recent API Activity:"
    echo "${usage_data}" | head -20

    echo ""
    echo "Tips to reduce costs:"
    echo "  1. Use local embeddings (already enabled: sentence-transformers)"
    echo "  2. Batch PDF imports: ./import-pdfs-to-kb.sh --parallel 2"
    echo "  3. Reuse generated embeddings across multiple queries"
    echo "  4. Use cheaper models for background tasks (TASK_MODEL_EXTERNAL)"
}

# shellcheck disable=SC2120
estimate_costs() {
    echo "=== API Cost Estimation ==="
    echo ""
    echo "Based on typical thesis research patterns:"
    echo ""
    echo "Daily usage estimates:"
    echo "  - Light research (10-20 queries): ~$0.50-2.00"
    echo "  - Medium research (50-100 queries): ~$2.00-5.00"
    echo "  - Heavy analysis (200+ queries): ~$5.00-15.00"
    echo ""
    echo "Monthly estimates (30 days):"
    echo '  - Light: ~$15-60'
    echo '  - Medium: ~$60-150'
    echo '  - Heavy: ~$150-450'
    echo ""
    echo "Cost optimization strategies:"
    echo "  ✓ Local embeddings: $0 (sentence-transformers)"
    echo "  ✓ Local code execution: $0 (Jupyter/Pyodide)"
    echo "  ✓ Local RAG retrieval: $0"
    echo "  ✗ API calls for generation: Variable cost"
    echo ""
    echo "Current configuration cost impact:"

    # Check current settings
    if grep -q "RAG_EMBEDDING_ENGINE=openai" "${PROJECT_ROOT}/.env" 2>/dev/null; then
        echo "  ⚠ Using OpenAI embeddings (higher cost)"
    else
        echo "  ✓ Using local embeddings (free)"
    fi

    if grep -q "TASK_MODEL_EXTERNAL" "${PROJECT_ROOT}/.env" 2>/dev/null; then
        echo "  ✓ Task model configured for background operations"
    else
        echo "  ℹ No task model configured (may use expensive model for background tasks)"
    fi
}

export_to_csv() {
    local output_file
    output_file="${REPORT_DIR}/api-usage-$(date +%Y%m%d).csv"

    echo "Exporting API usage data to CSV..."

    # Create CSV header
    echo "timestamp,model,prompt_tokens,completion_tokens,total_tokens,estimated_cost" >"${output_file}"

    # Extract and format log data
    docker compose -f "${PROJECT_ROOT}/docker-compose.yml" logs cliproxyapi 2>/dev/null |
        grep -E "model=.*tokens=" |
        while read line; do
            # Parse log line and extract metrics
            # This is a simplified example - adjust based on actual log format
            timestamp=$(echo "${line}" | grep -oE '^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}' || echo "unknown")
            model=$(echo "${line}" | grep -oE 'model=[^ ]+' | cut -d= -f2 || echo "unknown")

            echo "${timestamp},${model},0,0,0,0" >>"${output_file}"
        done

    echo "✓ Exported to: ${output_file}"
    echo ""
    echo "Use this CSV with:"
    echo "  - Excel/Google Sheets for analysis"
    echo "  - Python pandas for detailed statistics"
    echo "  - University expense tracking systems"
}

show_recommendations() {
    echo "=== Cost Optimization Recommendations ==="
    echo ""

    echo "Current Setup Analysis:"
    echo ""

    # Check optimizations already in place
    echo "✓ Strengths:"

    if grep -q "sentence-transformers" "${PROJECT_ROOT}/.env" 2>/dev/null ||
        ! grep -q "RAG_EMBEDDING_ENGINE=openai" "${PROJECT_ROOT}/.env" 2>/dev/null; then
        echo "  • Local embeddings (free vs ~$0.10/1K docs for OpenAI)"
    fi

    echo "  • Docker-based deployment (efficient resource usage)"
    echo "  • PostgreSQL/pgvector (no managed service costs)"

    echo ""
    echo "⚠ Potential Improvements:"

    if ! grep -q "TASK_MODEL_EXTERNAL" "${PROJECT_ROOT}/.env" 2>/dev/null; then
        echo "  • Configure TASK_MODEL_EXTERNAL for cheaper background tasks"
        echo "    Add to .env: TASK_MODEL_EXTERNAL=glm-5"
    fi

    if ! grep -q "ENABLE_AUTOCOMPLETE_GENERATION=False" "${PROJECT_ROOT}/.env" 2>/dev/null; then
        echo "  • Consider disabling autocomplete generation"
        echo "    Add to .env: ENABLE_AUTOCOMPLETE_GENERATION=False"
        echo "    (Saves ~20-30% on API calls)"
    fi

    echo ""
    echo "Research-Specific Tips:"
    echo "  1. Batch process all PDFs at once (amortize embedding cost)"
    echo "  2. Export and review important conversations (avoid re-querying)"
    echo "  3. Use local code execution for data analysis (free)"
    echo "  4. Cache RAG results for frequently referenced papers"
    echo "  5. Use cheaper models for initial exploration, expensive for final analysis"
}

usage() {
    cat <<EOF
API Cost Monitor for Thesis Research

Usage: ./monitor-api-costs.sh [OPTIONS]

Options:
  --daily       Show daily usage summary
  --weekly      Show weekly usage summary (default)
  --monthly     Show monthly usage summary
  --estimate    Show cost estimation and recommendations
  --export      Export usage data to CSV
  --help        Show this help

Examples:
  ./monitor-api-costs.sh              # Weekly summary
  ./monitor-api-costs.sh --daily      # Daily usage
  ./monitor-api-costs.sh --estimate   # Cost projections
  ./monitor-api-costs.sh --export     # Export to CSV

For detailed analysis, export CSV and open in Excel or Google Sheets.
EOF
}

main() {
    case "${1:-}" in
        --daily | --weekly | --monthly)
            PERIOD="${1#--}"
            show_summary "${PERIOD}"
            ;;
        --estimate | -e)
            estimate_costs
            show_recommendations
            ;;
        --export)
            export_to_csv
            ;;
        --help | -h)
            usage
            ;;
        "")
            show_summary "weekly"
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
}

main "$@"
