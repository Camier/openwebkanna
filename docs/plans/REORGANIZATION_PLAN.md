# Repository Reorganization Plan

Historical note:
- This document captures an earlier cleanup plan.
- Use `README.md`, `docs/REPO_MAP.md`, `docs/ssot/stack.md`, and `config/README.md` for the current canonical structure and ownership rules.

**Goal:** Organize for clarity, keep what's used, de-emphasize/remove what isn't, consistent naming, incremental behavior-preserving changes.

**Principles:**
1. **Incremental:** Each PR is independently mergeable, no breaking changes
2. **Behavior-preserving:** Scripts continue working after moves (update paths)
3. **Clear hierarchy:** Active vs. deprecated vs. reference material
4. **Consistent naming:** Noun-verb or verb-noun patterns, no mixed styles

---

## Current State Analysis

### Active Scripts (Keep in root)
- `deploy.sh` ✅
- `status.sh` ✅
- `cleanup.sh` ✅
- `logs.sh` ✅
- `import-pdfs-to-kb.sh` ✅
- `test-rag.sh` ✅
- `test-api.sh` ✅
- `tune-openwebui-documents.sh` ✅
- `llm-council.sh` ✅
- `backup-openwebui-db.sh` ✅
- `openwebui-user-admin.sh` ✅
- `audit-no-mock.sh` ✅
- `verify-scripts.sh` ✅
- `update.sh` ✅

### Documentation Overload
- **Core docs (keep in root):** README.md, SETUP.md, TROUBLESHOOTING.md, SECURITY.md, API_EXAMPLES.md
- **Move to docs/:** 15+ reference guides, technical deep dives

### Processing Artifacts
- `prod_max/` - 129 paper renderings → move to `data/processing/prod_max/`
- `prod_max_multimodal/` - multimodal chunks → move to `data/processing/prod_max_multimodal/`

---

## PR Units (Prioritized)

### PR 1: Documentation Hierarchy (High Priority)

**Rationale:** 100+ Markdown files clutter root. Create clear hierarchy.

**Changes:**
```bash
# Create structure
mkdir -p docs/{reference,guides,architecture,internals}

# Move reference docs
mv openwebui_rag_technical_reference.md docs/reference/rag-system.md
mv rag_technical_reference_review.md docs/reference/rag-review.md
mv openwebui_env_reference.md docs/reference/environment-vars.md
mv openwebui_documentation.md docs/reference/openwebui-overview.md
mv OPENWEBUI_PIPELINES_GUIDE.md docs/reference/pipelines.md
mv OPENWEBUI_WEBSEARCH_DEEP_DIVE.md docs/reference/websearch.md
mv OPENWEBUI_MASTER_REFERENCE.md docs/reference/master-reference.md
mv openwebui-tools-functions-guide.md docs/reference/tools-functions.md
mv MCP_INTEGRATION_GUIDE.md docs/reference/mcp-integration.md

# Move guides
mv ACADEMIC_TOOLS_GUIDE.md docs/guides/academic-tools.md
mv API_EXAMPLES.md docs/guides/api-examples.md  # Keep copy in root

# Keep in root (core docs)
# - README.md (main entry point)
# - SETUP.md (getting started)
# - TROUBLESHOOTING.md (problem resolution)
# - SECURITY.md (security policies)
# - PREREQUISITES.md (requirements)
```

**Update internal links:**
```bash
# Script to update relative links in moved files
find docs -name "*.md" -exec sed -i 's|](../|](../../|g' {} \;
find docs -name "*.md" -exec sed -i 's|](./|](.|g' {} \;
```

**Behavior preserved:** All docs accessible, links updated

---

### PR 3: Processing Artifacts Restructure (Medium Priority)

**Rationale:** `prod_max/` and `prod_max_multimodal/` are processing outputs, not configuration. Should live under `data/`.

**Changes:**
```bash
# Move processing outputs
mkdir -p data/processing
mv prod_max/ data/processing/prod_max/
mv prod_max_multimodal/ data/processing/prod_max_multimodal/

# Update scripts that reference these paths
# import-pdfs-to-kb.sh - no change (uses data/pdfs/)
# embed_images_in_chunks.py - update input/output paths
```

**Git consideration:** These directories are likely gitignored. If tracked, use `git mv` for history preservation.

**Behavior preserved:** Scripts updated to use new paths

---

### PR 4: Configuration Consolidation (Medium Priority)

**Rationale:** Keep `config/env/.env.example` and docs consistently aligned with the canonical LiteLLM-first baseline.

**Changes:**
```bash
# config/env/.env.example - Keep section order aligned with runtime defaults
# Keep LiteLLM section as primary OpenAI-compatible upstream

# Current order (target):
# 1. Authentication
# 2. Jupyter
# 3. Ollama
# 4. LiteLLM Proxy (PRIMARY)
# 5. MCPO
# 6. RAG
```

**Behavior preserved:** All configs still work, just reorganized

---

### PR 5: Create OPERATIONS.md Quick Reference (Medium Priority)

**Rationale:** README.md is 500+ lines. Operators need quick command reference without navigation.

**Changes:**
```bash
# Create new file: OPERATIONS.md
cat > OPERATIONS.md << 'EOF'
# OpenWebUI Operations Quick Reference

**For detailed setup:** See [README.md](README.md)
**For troubleshooting:** See [TROUBLESHOOTING.md](TROUBLESHOOTING.md)

---

## Daily Operations

### Start Stack
```bash
./deploy.sh --no-logs
```

### Check Status
```bash
./status.sh
```

### View Logs
```bash
./logs.sh                    # All services
./logs.sh openwebui          # Single service
docker compose logs -f       # Follow mode
```

### Restart Services
```bash
docker compose restart openwebui
```

### Stop Stack
```bash
docker compose down
./cleanup.sh
```

---

## RAG Operations

### Import PDFs (Serial)
```bash
./import-pdfs-to-kb.sh --dir data/pdfs --kb-name "My Papers"
```

### Import PDFs (Parallel, 3 workers)
```bash
IMPORT_PDFS_PARALLELISM=3 ./import-pdfs-to-kb.sh
```

### Tune RAG Settings
```bash
./scripts/admin/sync-openwebui-retrieval-config.sh
```

### Test RAG
```bash
./test-rag.sh --baseline     # Fast checks (2 min)
./test-rag.sh --full         # Full suite (15 min)
```

---

## Backup & Recovery

### Backup Database
```bash
./backup-openwebui-db.sh
```

### Restore from Backup
```bash
# Stop stack
docker compose down

# Restore volume
tar -xzf backups/openwebui-db-*.tar.gz -C /var/lib/docker/volumes/openwebui_data/

# Restart
docker compose up -d
```

---

## Quick Diagnostics

### Service Health
```bash
curl http://localhost:3000/health
curl http://localhost:8000/docs  # MCPO
```

### Model Availability
```bash
OPENWEBUI_API_KEY="<admin-bearer-token>" curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" http://localhost:3000/api/models | jq '.data[].id'
```

### RAG Settings
```bash
curl -s http://localhost:3000/api/v1/retrieval/ | jq
```

---

## Common Issues

| Issue | Command |
|-------|---------|
| No models in dropdown | `./openwebui-user-admin.sh --email you@example.com --role admin` |
| Pending activation | `./openwebui-user-admin.sh --email you@example.com --role user` |
| 502 errors | `./backup-openwebui-db.sh && docker compose restart openwebui` |
| SSL errors | Place CA at `certs/ca-bundle.pem`, set `REQUESTS_CA_BUNDLE=/certs/ca-bundle.pem` |
EOF
```

**Behavior preserved:** New file adds value, doesn't change existing behavior

---

### PR 6: Cleanup Empty Directories (Low Priority)

**Rationale:** Empty directories should have `.gitkeep` for explicit tracking, or be removed if truly unused.

**Changes:**
```bash
# Find empty directories
find /LAB/@thesis/openwebui -type d -empty

# Add .gitkeep to intentionally-empty dirs
touch data/metadata/.gitkeep
touch data/notes/.gitkeep
touch certs/.gitkeep
touch backups/.gitkeep
touch logs/.gitkeep

# Create cleanup script: scripts/cleanup-empty-dirs.sh
cat > scripts/cleanup-empty-dirs.sh << 'EOF'
#!/bin/bash
# Remove truly empty directories (not tracked via .gitkeep)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

echo "Finding empty directories..."
find . -type d -empty ! -path "./.git/*" -print

echo ""
echo "To remove, run:"
echo "  find . -type d -empty ! -path './.git/*' -delete"
EOF

chmod +x scripts/cleanup-empty-dirs.sh
```

---

### PR 7: Script Naming Consistency (Low Priority)

**Rationale:** Mixed naming patterns:
- `thesis-backup.sh` (topic-noun)
- `export-thesis-chats.sh` (verb-topic-noun)
- `thesis-rag-helper.sh` (topic-noun-noun)

**Proposed standard:** **verb-noun** for action scripts, **noun-helpers** for libraries

**Changes:**
```bash
# Rename inconsistent scripts (keep backward compat with symlinks)

# Action scripts (verb-noun)
mv scripts/thesis-backup.sh scripts/backup-thesis.sh
mv scripts/thesis-rag-helper.sh scripts/rag-helper-thesis.sh
mv scripts/export-thesis-chats.sh scripts/export-chats-thesis.sh

# Create symlinks for backward compatibility
ln -s backup-thesis.sh thesis-backup.sh
ln -s export-chats-thesis.sh export-thesis-chats.sh

# Library scripts (noun-helpers)
# Already consistent: docker-helpers.sh, http-helpers.sh
```

**Better approach:** Don't rename. Script names are muscle memory. Instead, document naming convention for NEW scripts only.

**Recommendation:** Skip this PR. Focus on functionality over naming perfection.

---

## Implementation Timeline

### Week 1 (High Priority)
- ✅ PR 1: remove legacy operator-path references
- ✅ PR 2: Documentation hierarchy
- ⚠️ Verify all links work after moves

### Week 2 (Medium Priority)
- ✅ PR 3: Processing artifacts
- ✅ PR 4: Configuration consolidation
- ✅ PR 5: OPERATIONS.md

### Week 3 (Low Priority)
- ✅ PR 6: Empty directory cleanup
- ❌ PR 7: Skip (not worth breaking muscle memory)

---

## Risk Mitigation

### Before Each PR
1. Run `./test-rag.sh --baseline` to establish baseline
2. Ensure no open work in progress
3. Create git branch: `git checkout -b reorg/pr-N`

### After Each PR
1. Run `./test-rag.sh --baseline` to verify no breakage
2. Test key workflows:
   ```bash
   ./deploy.sh --no-logs
   ./status.sh
   ./import-pdfs-to-kb.sh --limit 1  # Test single import
   ```
3. Verify all moved documentation links work

### Rollback Plan
Each PR is independently revertible:
```bash
git revert HEAD  # Undo last PR
git checkout master  # Return to clean state
```

---

## Post-Reorganization Structure

```
openwebui/
├── README.md                    # Main entry point (trimmed to 200 lines)
├── OPERATIONS.md                # NEW: Quick reference
├── SETUP.md                     # Getting started
├── TROUBLESHOOTING.md           # Problem resolution
├── SECURITY.md                  # Security policies
├── PREREQUISITES.md             # Requirements
├── API_EXAMPLES.md              # API usage examples
│
├── docs/                        # NEW: Documentation hierarchy
│   ├── reference/
│   │   ├── rag-system.md
│   │   ├── environment-vars.md
│   │   └── ...
│   ├── guides/
│   │   ├── academic-tools.md
│   │   └── api-examples.md
│   └── architecture/
│       └── ...
│
├── *.sh                         # Active scripts (40 files)
├── lib/                         # Shared libraries
├── scripts/                     # Utility scripts
│   ├── check-image-versions.sh
│   ├── audit-dependencies.sh
│   └── ...
│
├── data/
│   ├── pdfs/                    # Source PDFs
│   ├── extractions/             # Per-paper data
│   ├── corpus/                  # JSONL corpora
│   ├── processing/              # NEW: prod_max moved here
│   │   ├── prod_max/
│   │   └── prod_max_multimodal/
│   └── ...
│
├── config/
│   ├── compose/docker-compose.yml
│   └── env/.env.example
└── ...
```

---

## Success Metrics

- [ ] All scripts execute without path errors
- [ ] `./test-rag.sh --baseline` passes
- [ ] Documentation links work (no 404s)
- [ ] New users can find key commands in OPERATIONS.md
- [ ] Root directory has <50 files (currently ~100)
- [ ] Clear separation: active vs. archived vs. reference
