#!/bin/bash
###############################################################################
# Thesis RAG Workflow Helper - Optimize research queries for academic writing
#
# Usage:
#   ./thesis-rag-helper.sh [COMMAND]
#
# Commands:
#   query-templates    Show effective query templates for thesis research
#   document-tips      Show document organization best practices
#   citation-guide     Show how to cite RAG sources in your thesis
#   export-methodology Export research methodology documentation
#
# This script helps optimize RAG usage for ethnopharmacological research.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

show_query_templates() {
    cat <<'EOF'
=== Effective RAG Query Templates for Thesis Research ===

1. MECHANISM/MODE OF ACTION QUERIES
   (Best for pharmacology chapters)

   ❌ Ineffective: "Tell me about Sceletium"

   ✅ Effective: "What are the reported serotonin reuptake inhibition
      mechanisms of mesembrine in Sceletium tortuosum? Include specific
      binding affinities if mentioned."

   Template: "What are the [mechanism] of [compound] in [plant]?
              Include [specific metrics] if mentioned."

2. COMPARATIVE ANALYSIS QUERIES
   (Best for literature review chapters)

   ❌ Ineffective: "Compare different studies"

   ✅ Effective: "Compare the alkaloid profiles reported by Smith 2019
      and Jones 2021 for Sceletium tortuosum. What are the key differences
      in extraction methodology and detected compounds?"

   Template: "Compare [finding] reported by [Author Year] and [Author Year]
              for [subject]. What are the key differences in [aspect]?"

3. TRADITIONAL USE/ETHNOBOTANY QUERIES
   (Best for traditional knowledge chapters)

   ✅ Effective: "Based on the papers in the knowledge base, what are the
      traditional preparation methods for kanna (Sceletium tortuosum)
      documented by indigenous communities? Include specific citations
      with page numbers where available."

   Template: "Based on [document set], what are the [traditional practices]
              documented by [group]? Include [citation format]."

4. STRUCTURED DATA EXTRACTION
   (Best for creating summary tables)

   ✅ Effective: "Create a table of all pharmacological studies on Sceletium
      tortuosum including: Author/Year, Compound Tested, Methodology,
      Sample Size, Key Findings, and Limitations noted by authors."

   Template: "Create a [format] of all [studies] on [topic] including:
              [field1], [field2], [field3], etc."

5. GAP ANALYSIS QUERIES
   (Best for research proposal/methodology)

   ✅ Effective: "Based on the reviewed literature, what aspects of
      Sceletium tortuosum pharmacokinetics remain understudied?
      What methodological limitations are consistently mentioned
      across multiple studies?"

   Template: "Based on [literature set], what [gaps] remain in [field]?
              What [limitations] are consistently mentioned?"

6. SAFETY/TOXICOLOGY QUERIES
   (Best for safety chapters)

   ✅ Effective: "What toxicological data is available for Sceletium
      tortuosum alkaloids? Include LD50 values, reported adverse effects,
      contraindications, and safety thresholds from clinical studies."

   Template: "What [safety data] is available for [subject]? Include
              [specific metrics] and [safety parameters]."

TIPS FOR EFFECTIVE QUERIES:
  • Be specific - name exact compounds, authors, or study types
  • Ask for citations - "Include specific citations" triggers source tracking
  • Request structured output - tables, lists, comparisons
  • Mention page numbers - helps verify claims against original papers
  • Use follow-up questions - drill deeper into interesting findings

EOF
}

show_document_tips() {
    cat <<'EOF'
=== Document Organization Best Practices ===

1. FILE NAMING CONVENTION
   Format: Year_Author_Plant-Part_Preparation.pdf

   Examples:
   ✓ 2019_Smith_Sceletium-whole-plant_Aqueous-extract.pdf
   ✓ 2021_Jones_Sceletium-root_Ethanol-extract.pdf
   ✗ smith paper.pdf (too vague)
   ✗ sceletium.pdf (no year or author)

2. DIRECTORY STRUCTURE
   Organize PDFs by topic for easier batch processing:

   data/pdfs/
   ├── sceletium-pharmacology/
   │   ├── 2019_Smith_Mesembrine_SRI-study.pdf
   │   └── 2021_Jones_Alkaloid-profile_HPLC.pdf
   ├── sceletium-traditional-use/
   │   ├── 2018_VanWyk_Ethnobotany_Khoisan.pdf
   │   └── 2020_Louw_Preparation-methods.pdf
   ├── sceletium-toxicology/
   │   ├── 2017_Harvey_Safety-assessment.pdf
   │   └── 2022_Patel_LD50-study.pdf
   └── other-medicinal-plants/
       └── (other species for comparative analysis)

3. BATCH PROCESSING STRATEGY

   Step 1: Organize PDFs into topic folders
   Step 2: Process one topic at a time:

   ./import-pdfs-to-kb.sh --parallel 2 data/pdfs/sceletium-pharmacology/
   ./import-pdfs-to-kb.sh --parallel 2 data/pdfs/sceletium-traditional-use/

   Step 3: Verify embeddings in OpenWebUI admin panel

4. METADATA ENRICHMENT

   Before importing, ensure PDFs have:
   ✓ Clear titles (visible in OpenWebUI document list)
   ✓ Author/year in filename (helps citation tracking)
   ✓ OCR text layer (required for text extraction)

   Check OCR: pdfinfo document.pdf | grep "PDF version"
   If needed, OCR with: ocrmypdf input.pdf output.pdf

5. KNOWLEDGE BASE MAINTENANCE

   Monthly tasks:
   • Review and remove duplicate papers
   • Update with newly published research
   • Re-embed updated versions of papers
   • Export citations to reference manager

EOF
}

show_citation_guide() {
    cat <<'EOF'
=== Citing RAG Sources in Your Thesis ===

1. ACADEMIC INTEGRITY GUIDELINES

   RAG as Research Assistant (Not Primary Source):
   • Use RAG to FIND relevant papers
   • ALWAYS verify claims against original PDFs
   • Cite the original papers, not the RAG system
   • Document RAG use in methodology chapter

2. CITATION FORMAT FOR THESIS

   In-text citation (correct):
   "Mesembrine demonstrates serotonin reuptake inhibition (Smith et al.,
   2019, p. 45; Jones, 2021, p. 112)."

   In-text citation (incorrect):
   "OpenWebUI says mesembrine inhibits serotonin reuptake."

   Methodology section (example):
   "Literature was analyzed using OpenWebUI v0.8.3 with a locally-hosted
   RAG system. The knowledge base contained 47 papers on Sceletium tortuosum
   (see Appendix A for complete list). Queries were structured to extract
   specific pharmacological data with citation verification against original
   sources."

3. SOURCE VERIFICATION CHECKLIST

   Before including any claim in your thesis:

   □ Locate the cited paper in your PDF collection
   □ Open the paper to the cited page
   □ Verify the exact claim appears in the original
   □ Check if context changes the interpretation
   □ Note any limitations the original authors mentioned
   □ Add to reference manager (Zotero/Mendeley)

4. DOCUMENTING YOUR RESEARCH PROCESS

   Export methodology documentation:

   ./export-thesis-chats.sh --all

   This creates:
   • Complete query history
   • RAG responses with citations
   • Your follow-up questions
   • Final verified conclusions

   Store in: thesis/methodology/rag-research-log/

5. REFERENCE MANAGER INTEGRATION

   Recommended workflow:
   1. Find paper via RAG query
   2. Export citation from OpenWebUI
   3. Import to Zotero/Mendeley
   4. Add notes about key findings
   5. Tag with thesis chapter/section

6. ETHICAL CONSIDERATIONS

   For Traditional Knowledge:
   • Follow Nagoya Protocol if applicable
   • Acknowledge indigenous communities
   • Respect data sovereignty principles
   • Document consent and collaboration

   Example acknowledgment:
   "Traditional use information was shared with permission from [Community]
   elders during [Date/Project]. This research respects indigenous data
   sovereignty principles as outlined in [Agreement]."

EOF
}

export_methodology() {
    local export_dir="${PROJECT_ROOT}/thesis-exports/methodology"
    mkdir -p "${export_dir}"

    local timestamp
    timestamp=$(date +%Y%m%d_%H%M%S)
    local output_file="${export_dir}/rag-methodology-${timestamp}.md"

    cat >"${output_file}" <<EOF
# RAG Research Methodology Documentation

**Generated:** $(date -Iseconds)
**Thesis Topic:** [Your thesis title here]
**Researcher:** [Your name]

## System Configuration

### Software Stack
- OpenWebUI v0.8.3
- PostgreSQL with pgvector extension
- CLIProxyAPI for model access
- Sentence-transformers for embeddings

### Knowledge Base Composition
EOF

    # Count documents in knowledge base if available
    if [[ -f "${PROJECT_ROOT}/data/corpus/biblio_corpus.jsonl" ]]; then
        local doc_count
        doc_count=$(wc -l <"${PROJECT_ROOT}/data/corpus/biblio_corpus.jsonl" 2>/dev/null || echo "unknown")
        echo "- Total Documents: ${doc_count}" >>"${output_file}"
    fi

    # List PDF directories
    echo "" >>"${output_file}"
    echo "### Document Collections" >>"${output_file}"
    for dir in "${PROJECT_ROOT}/data/pdfs"/*/; do
        if [[ -d ${dir} ]]; then
            local count
            count=$(find "${dir}" -name "*.pdf" | wc -l)
            echo "- $(basename "${dir}"): ${count} PDFs" >>"${output_file}"
        fi
    done

    cat >>"${output_file}" <<'EOF'

## RAG Configuration

### Embedding Settings
- Model: sentence-transformers/all-MiniLM-L6-v2
- Chunk Size: 1500 characters
- Chunk Overlap: 150 characters
- Top-K Retrieval: 5 chunks

### Query Methodology
- Structured queries with specific parameters
- Citation verification against original sources
- Iterative refinement with follow-up questions
- Export and review of significant findings

## Research Log

[Document your specific research queries and findings here]

### Chapter: [Chapter Name]

**Research Questions:**
- [Question 1]
- [Question 2]

**Key Findings:**
- [Finding with verified citation]
- [Finding with verified citation]

**Papers Consulted:**
- Author (Year). Title. Journal.

## Quality Assurance

### Verification Process
- [ ] All claims verified against original papers
- [ ] Citations checked for accuracy
- [ ] Limitations noted for each study
- [ ] Conflicting findings documented
- [ ] Reference manager updated

## Appendices

### Appendix A: Complete Bibliography
[List all papers in knowledge base]

### Appendix B: Query History
[Export from ./export-thesis-chats.sh --all]

---

**Methodology Template Version:** 1.0
**Last Updated:** $(date +%Y-%m-%d)
EOF

    echo "✓ Methodology template exported to:"
    echo "  ${output_file}"
    echo ""
    echo "Next steps:"
    echo "  1. Fill in your thesis title and name"
    echo "  2. Document your research questions"
    echo "  3. Record key findings with citations"
    echo "  4. Run ./export-thesis-chats.sh --all"
    echo "  5. Attach complete bibliography as Appendix A"
}

usage() {
    cat <<'EOF'
Thesis RAG Workflow Helper

Usage: ./thesis-rag-helper.sh [COMMAND]

Commands:
  query-templates      Show effective query templates for thesis research
  document-tips        Show document organization best practices
  citation-guide       Show how to cite RAG sources in your thesis
  export-methodology   Export methodology documentation template
  help                 Show this help

Examples:
  ./thesis-rag-helper.sh query-templates    # Get query templates
  ./thesis-rag-helper.sh citation-guide     # Learn citation best practices
  ./thesis-rag-helper.sh export-methodology # Generate methodology doc

For solo academic research on ethnopharmacology (Sceletium tortuosum).
EOF
}

main() {
    case "${1:-}" in
        query-templates)
            show_query_templates
            ;;
        document-tips)
            show_document_tips
            ;;
        citation-guide)
            show_citation_guide
            ;;
        export-methodology)
            export_methodology
            ;;
        help | --help | -h)
            usage
            ;;
        "")
            echo "Thesis RAG Workflow Helper"
            echo ""
            echo "Available commands:"
            echo "  query-templates      Query templates for effective research"
            echo "  document-tips        Document organization best practices"
            echo "  citation-guide       Academic citation guidelines"
            echo "  export-methodology   Generate methodology documentation"
            echo ""
            echo "Run './thesis-rag-helper.sh help' for detailed usage."
            ;;
        *)
            echo "Unknown command: $1"
            usage
            exit 1
            ;;
    esac
}

main "$@"
