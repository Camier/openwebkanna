# OpenWebUI Plugins, Tools & Addons for Thesis Research

## Executive Summary

For your ethnopharmacology thesis on Sceletium tortuosum, here are the most valuable additions to your OpenWebUI setup, categorized by research phase.

---

## Must-Have: Enable Native Mode (Agentic Mode)

**Priority: CRITICAL**

This unlocks all built-in tools. Without it, you're missing most functionality.

```bash
# In your .env file or docker-compose.yml environment:
ENABLE_NATIVE_TOOL_CALLS=True
ENABLE_NATIVE_API_TOOL_CALLS=True
ENABLE_NATIVE_API_PIPE=True
```

**Verify it's working:**
- Go to OpenWebUI Settings → Models
- Select a model → Enable "Tools"
- Look for tools icon (🔧) in chat interface

---

## 1. Built-in Tools (Already Available)

Once Native Mode is enabled, these work immediately:

### Knowledge Base Tools (Essential for Your RAG Setup)
- **`query_knowledge_files`** - Search inside your uploaded PDFs
- **`view_knowledge_file`** - Read full paper content
- **`list_knowledge_bases`** - See all your document collections

**How to use:**
```
"Search my knowledge base for traditional Khoisan preparation
 methods of Sceletium tortuosum"
```

### Web Research Tools
- **`search_web`** - Search the web (uses your SearXNG)
- **`fetch_url`** - Extract content from specific URLs
- **`query_knowledge_bases`** - Find relevant documents

### Note-Taking Tools (Great for Thesis Organization)
- **`write_note`** - Save research findings
- **`search_notes`** - Find previous notes
- **`add_memory`** - Store facts about your research

### Code Execution
- **`execute_code`** - Run Python/R for data analysis

---

## 2. Recommended Community Plugins

### A. Zotero MCP Integration ⭐ TOP PRIORITY

**What:** Connect your Zotero library directly to OpenWebUI

**Why you need it:**
- Your thesis citations are probably in Zotero
- Semantic search across your entire library
- AI can find papers you forgot you had
- Direct BibTeX export for LaTeX

**Installation:**

```bash
# 1. Install Zotero MCP server
pip install zotero-mcp

# 2. Configure in OpenWebUI:
# Settings → Admin Settings → Tools → Add Tool
# Or use Functions section if Tools not available

# 3. Set environment variables in .env:
ZOTERO_USER_ID=your_zotero_user_id
ZOTERO_API_KEY=your_zotero_api_key
ZOTERO_LIBRARY_TYPE=user  # or 'group' for shared libraries
```

**Get your credentials:**
- Zotero User ID: https://www.zotero.org/settings/keys (look for "Your userID")
- API Key: https://www.zotero.org/settings/keys/new

**Use case examples:**
```
"Find all papers in my Zotero library about Sceletium alkaloids"
"Export BibTeX citations for my pharmacology chapter"
"What papers do I have about serotonin reuptake inhibition?"
```

---

### B. arXiv Search Tool

**What:** Search latest research on arXiv

**Why:** Stay current with new Sceletium/pharmacology research

**Installation:**
```bash
# Via OpenWebUI interface:
# 1. Go to Workspace → Tools
# 2. Click "Get Tools" or "Community"
# 3. Search for "arxiv"
# 4. Install and enable
```

**Or manual install:**
```bash
# Download from: https://github.com/Haervwe/open-webui-tools
cd /LAB/@thesis/openwebui
# Place in appropriate tools directory
```

**Use case:**
```
"Search arXiv for recent papers on Sceletium tortuosum pharmacology"
"Find new research on mesembrine mechanisms published this year"
```

---

### C. Doc Builder (MD/PDF Export) ⭐ ESSENTIAL

**What:** Export chats as formatted PDFs

**Why:** Turn research conversations into thesis draft sections

**Installation:**
```bash
# Via OpenWebUI Hub:
# 1. Go to Workspace → Tools
# 2. Click "Get Tools"
# 3. Search: "Doc Builder"
# 4. Install "Doc Builder MD/PDF v1.8"
```

**Use case:**
```
# After a productive research chat:
"Export this conversation as PDF for my literature review chapter"
```

---

### D. Data Analyst Tool

**What:** Natural language data analysis with pandas/visualization

**Why:** Analyze phytochemical data, create figures for thesis

**Installation (Docker):**
```yaml
# Add to docker-compose.yml
services:
  openwebui-data-analyst:
    image: yoontheopen/openwebui-data-analyst:latest
    ports:
      - "8081:8080"
    volumes:
      - ./data/analyst:/app/data
    environment:
      - API_KEY=your-openwebui-api-key
```

**Or use your existing Jupyter:**
Your Jupyter setup already handles this - just use `execute_code` tool.

---

## 3. Configuration Optimizations

### Enable All Tools by Default

```bash
# In .env
DEFAULT_MODELS_TOOL_CHOICE=auto  # Auto-select tools
ENABLE_TOOLS=True
ENABLE_FUNCTION_CALLING=True
```

### Optimize RAG for Academic Papers

```bash
# In .env - already good, but verify:
RAG_EMBEDDING_MODEL=sentence-transformers/all-MiniLM-L6-v2
RAG_TOP_K=5
CHUNK_SIZE=1500
CHUNK_OVERLAP=150
ENABLE_RAG_HYBRID_SEARCH=True  # BM25 + vector search
ENABLE_RAG_WEB_SEARCH=False     # Don't mix web results with your papers
```

### Configure Memory for Long Research Sessions

```bash
# Enable memory to remember context across chats
ENABLE_MEMORY=True
MEMORY_EXPIRY_DAYS=365  # Keep thesis context for entire project
```

---

## 4. External Tool Integrations

### Zotero (Already covered above)
**Integration:** MCP protocol
**Use:** Citation management

### Obsidian (Optional but powerful)
**What:** Link OpenWebUI with Obsidian note-taking

**Why:** Many academics use Obsidian for research notes

**Setup:**
```bash
# Use OpenWebUI's note export + Obsidian folder sync
# Or use API to sync directly (advanced)
```

### LaTeX Integration
**For thesis writing:**

1. Export chats as Markdown (Doc Builder)
2. Convert to LaTeX with pandoc:
   ```bash
   pandoc research-chat.md -o chapter.tex --biblatex
   ```

3. Or use custom pipe for direct LaTeX export (see below)

---

## 5. Custom Functions You Should Create

### A. Thesis Citation Formatter

Create a custom function that formats citations for your thesis:

```python
# Save to: functions/thesis_citation.py
"""
Format citations for ethnopharmacology thesis
"""

def format_citation(author, year, title, journal, doi=None):
    """Format citation in Vancouver style (common for pharmacology)"""
    citation = f"{author}. {title}. {journal}. {year};"
    if doi:
        citation += f" doi:{doi}"
    return citation

def export_to_bibtex(papers):
    """Export paper list to BibTeX"""
    # Implementation here
    pass
```

### B. Phytochemical Data Extractor

Custom tool to extract structured data from your papers:

```python
# Save to: functions/phytochemical_extractor.py
"""
Extract phytochemical data from papers
"""

def extract_alkaloid_data(paper_content):
    """
    Extract: compound name, concentration, extraction method,
    biological activity, etc.
    """
    # Use LLM to parse paper and extract structured data
    pass
```

---

## 6. SearXNG Enhancements

You already have SearXNG. Optimize it for academic search:

```yaml
# In your SearXNG config (if accessible)
engines:
  - name: google scholar
    enabled: true
  - name: pubmed
    enabled: true
  - name: semantic scholar
    enabled: true
```

**Use with OpenWebUI:**
```
"Search the web for: mesembrine pharmacokinetics site:pubmed.ncbi.nlm.nih.gov"
```

---

## 7. Quick-Start Installation Script

Save this as `install-academic-tools.sh`:

```bash
#!/bin/bash
# Install academic research tools for OpenWebUI

echo "=== Installing Academic Tools ==="

# 1. Enable Native Mode
echo "[1/4] Enabling Native Mode..."
cat >> .env << 'EOF'

# Academic Research Tools
ENABLE_NATIVE_TOOL_CALLS=True
ENABLE_NATIVE_API_TOOL_CALLS=True
ENABLE_MEMORY=True
MEMORY_EXPIRY_DAYS=365
ENABLE_RAG_HYBRID_SEARCH=True
EOF

# 2. Install Zotero MCP
echo "[2/4] Installing Zotero MCP..."
pip3 install zotero-mcp --user || echo "Note: Install manually if pip fails"

echo ""
echo "Next steps:"
echo "  1. Add ZOTERO_USER_ID and ZOTERO_API_KEY to .env"
echo "  2. Restart OpenWebUI: docker compose restart"
echo "  3. Go to Workspace → Tools → Install Doc Builder"
echo "  4. Go to Workspace → Tools → Install arXiv Search"
echo "  5. Enable Native Mode in model settings"
echo ""
echo "Zotero credentials: https://www.zotero.org/settings/keys"
```

---

## 8. Recommended Workflow with Tools

### Phase 1: Literature Discovery
1. Use **arXiv Search** to find new papers
2. Use **Zotero MCP** to check if you already have them
3. Download to `data/pdfs/` and import to KB

### Phase 2: Deep Analysis
1. Upload PDFs to Knowledge Base
2. Use **`query_knowledge_files`** to ask specific questions
3. Use **`write_note`** to capture key findings
4. Use **`execute_code`** to analyze extracted data

### Phase 3: Writing
1. Research conversation → **Doc Builder** → PDF
2. Export notes with **export-thesis-chats.sh**
3. Zotero → BibTeX → LaTeX
4. Pandoc → Word (if advisor requires)

### Phase 4: Review
1. **`search_notes`** to find previous analysis
2. **`add_memory`** to store new insights
3. Cross-reference with knowledge base

---

## Summary: Priority Installation List

| Priority | Tool | Why |
|----------|------|-----|
| 🔴 Critical | Native Mode | Unlocks all built-in tools |
| 🔴 Critical | Knowledge Base Tools | Core RAG functionality |
| 🟡 High | Zotero MCP | Citation management |
| 🟡 High | Doc Builder | Export for thesis |
| 🟢 Medium | arXiv Search | Stay current |
| 🟢 Medium | Data Analyst | Statistical analysis |
| 🔵 Low | Custom Functions | Specialized workflows |

---

## Next Steps

1. **Immediate:** Enable Native Mode and restart
2. **Today:** Install Zotero MCP, configure credentials
3. **This week:** Install Doc Builder, arXiv Search
4. **Ongoing:** Create custom functions as needed

Your current setup with these additions will be a complete academic research platform!
