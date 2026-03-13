# OpenWebUI Operations Quick Reference

Last updated: 2026-03-13 (UTC)

Use this runbook when:

- the baseline stack already exists
- you need day-2 operator commands and validation loops

Use other docs for:

- first deploy: [SETUP.md](SETUP.md)
- machine and env readiness: [PREREQUISITES.md](PREREQUISITES.md)
- symptom-led recovery: [TROUBLESHOOTING.md](TROUBLESHOOTING.md)

Use this runbook only after the baseline stack has already been configured.
If you are still preparing `.env` or the first deploy, go back to `SETUP.md`.

## Entrypoint Policy

- Root `*.sh` scripts are the canonical operator surface.

Common daily commands:

```bash
./deploy.sh --no-logs
./status.sh
./logs.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

## Daily Operations

First-time environment preparation lives in `SETUP.md`.
Use this runbook for routine operator actions after the baseline stack has already been configured.

Normal post-change validation loop:

```bash
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

#### Start Stack
```bash
./deploy.sh --no-logs
```

#### Check Status
```bash
./status.sh
# Optional legacy sidecar check (only if CLIPROXYAPI_ENABLED=true):
./check-cliproxyapi.sh
```

#### View Logs
```bash
./logs.sh                    # All services
./logs.sh openwebui          # Single service
docker compose logs -f       # Follow mode
docker compose logs --tail=100 openwebui
```

#### Restart Services
```bash
docker compose restart openwebui
# Optional legacy sidecar restart (only if CLIPROXYAPI_ENABLED=true):
./restart-cliproxyapi.sh
```

#### Stop Stack
```bash
docker compose down
./cleanup.sh                 # Includes archived fallback cleanup if you explicitly enabled it
```

---

## RAG Operations

### Import PDFs

**Serial (default):**
```bash
./import-pdfs-to-kb.sh --dir data/pdfs --kb-name "My Papers"
```

**Parallel (faster for large batches):**
```bash
# 3 workers
IMPORT_PDFS_PARALLELISM=3 ./import-pdfs-to-kb.sh

# Or via command line
./import-pdfs-to-kb.sh --parallel 3
```

**Limit to first N PDFs:**
```bash
./import-pdfs-to-kb.sh --limit 5
```

**Custom directory:**
```bash
./import-pdfs-to-kb.sh --dir /path/to/pdfs --kb-name "Research Papers"
```

### Tune RAG Settings

Prefer `./test-rag.sh --baseline` for the normal health signal.
Use the raw retrieval endpoint only for advanced debugging when you already have an admin bearer token:

```bash
OPENWEBUI_API_KEY="<admin-bearer-token>" \
curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "http://127.0.0.1:${WEBUI_PORT:-3000}/api/v1/retrieval/" | jq
```

**Apply optimized settings (env overrides + `--apply`):**
```bash
# For focused factual retrieval
TUNE_TOP_K=8 \
TUNE_CHUNK_SIZE=1800 \
TUNE_CHUNK_OVERLAP=200 \
./tune-openwebui-documents.sh --apply

# For broader conceptual retrieval
TUNE_TOP_K=15 \
TUNE_CHUNK_SIZE=3000 \
TUNE_CHUNK_OVERLAP=600 \
./tune-openwebui-documents.sh --apply
```

These are intentional tuning examples, not hidden baseline defaults. The committed baseline remains `RAG_TOP_K=15`, `CHUNK_SIZE=3000`, `CHUNK_OVERLAP=600`, and `RAG_EMBEDDING_MODEL=pritamdeka/S-PubMedBert-MS-MARCO`.

**Snapshot current config to a custom directory:**
```bash
./tune-openwebui-documents.sh --snapshot-dir backups/rag-snapshots
```

**Restore retrieval/embedding settings from a snapshot:**
```bash
./tune-openwebui-documents.sh \
  --restore backups/rag-snapshots/openwebui-documents-snapshot-YYYYMMDD-HHMMSS.json
```

### Embedding Profiles (Reusable KB Use Cases)

Use the embedding profiles architecture for clean model lifecycle, KB bindings, and drift checks.
Source of truth:
- `docs/runbooks/EMBEDDING_PROFILES.md`

Daily operator flow:

```bash
# Inspect runtime and available lanes
./manage-openwebui-embedding-profiles.sh list
./manage-openwebui-embedding-profiles.sh lanes

# Activate biomedical text lane
./manage-openwebui-embedding-profiles.sh use-kb --lane sceletium --prewarm

# Diagnose compatibility/drift
./manage-openwebui-embedding-profiles.sh diagnose
```

### Test RAG

**Quick baseline (2 min):**
```bash
./test-rag.sh --baseline
```

**Full suite (15 min):**
```bash
./test-rag.sh --full
```

**Verbose output:**
```bash
./test-rag.sh --baseline --verbose
```

### Vector Database

**Check pgvector status (advanced debug):**
```bash
docker compose exec postgres psql -U openwebui -d openwebui -c \
  "SELECT COUNT(*) FROM collection_embeddings"
```

**Backup vector data:**
```bash
./backup-openwebui-db.sh
```

---

## LLM Council

### Quick Evaluation
```bash
./llm-council.sh --prompt "Compare RAG vs fine-tuning for academic papers"
```

### Multi-Model Benchmark
```bash
./llm-council.sh --models "glm-5 minimax/chat-elite" \
  --prompt "Explain the mechanism of action of Sceletium tortuosum"
```

### Batch Evaluation
```bash
./llm-council.sh --prompts-file data/notes/evaluation_prompts.txt \
  --models "glm-5 minimax/chat-elite" \
  --output-dir logs/llm-council/batch-1
```

### Anti-Position-Bias Testing
```bash
./llm-council.sh --models "glm-5 minimax/chat-elite" \
  --prompt "Which response is more accurate: A or B?" \
  --position-swap
```

---

## Legacy CLIProxyAPI Operations

### Start/Stop
```bash
./start-cliproxyapi.sh
./stop-cliproxyapi.sh
./restart-cliproxyapi.sh
```

### Health Check
```bash
./check-cliproxyapi.sh
./check-cliproxyapi.sh --models
```

### Configure Providers
```bash
# Generate config from ~/.007 secrets
./configure-cliproxyapi-providers.sh

# Setup OAuth (interactive)
./configure-cliproxyapi-oauth.sh
```

### Rotate API Key
```bash
./rotate-cliproxyapi-local-key.sh
```

### Test OAuth Aliases
```bash
./test-cliproxyapi-oauth.sh
./test-openwebui-cliproxy-routing.sh
```

---

## MCP Operations

### Configure MCP Servers
```bash
# Verify MCPO endpoints
./configure-mcpo-openapi-servers.sh --verify

# Apply configuration
./configure-mcpo-openapi-servers.sh
```

### Test MCP Integration
```bash
./test-mcp.sh
```

### Available MCP Servers
- `filesystem` - Access to `/data`, `/thesis-exports`
- `memory` - Session knowledge graph
- `fetch` - Web page fetching
- `zotero` - Zotero bibliographic database
- `time` - Time/date utilities

---

## Backup & Recovery

### Create Backup
```bash
./backup-openwebui-db.sh
# Creates: backups/webui.db.YYYYMMDD-HHMMSS.sqlite
```

### Restore from Backup
```bash
# 1. Select backup file to restore
BACKUP_FILE="backups/webui.db.YYYYMMDD-HHMMSS.sqlite"

# 2. Optional but recommended: take a pre-restore safety backup
./backup-openwebui-db.sh --output "backups/webui.db.pre-restore.$(date +%Y%m%d-%H%M%S).sqlite"

# 3. Stop OpenWebUI to avoid writes during restore
docker compose stop openwebui

# 4. Restore SQLite DB file inside the OpenWebUI container
docker cp "${BACKUP_FILE}" "$(docker compose ps -q openwebui)":/app/backend/data/webui.db

# 5. Start OpenWebUI and verify health
docker compose start openwebui
curl -fsS "http://127.0.0.1:${WEBUI_PORT:-3000}/health" >/dev/null && echo "openwebui healthy"
```

### Export Thesis Chats
```bash
./scripts/export-thesis-chats.sh
# Creates: thesis-exports/chat-<id>-YYYYMMDD.md
# Falls back to: logs/thesis-exports/chat-<id>-YYYYMMDD.md if thesis-exports/ is not writable
```

---

## Quick Diagnostics

### Service Health
```bash
# Prefer the script surface first
./status.sh
./test-api.sh --baseline

# Raw probes for advanced debugging
curl "http://127.0.0.1:${WEBUI_PORT:-3000}/health"

# CLIProxyAPI (only if enabled)
curl http://localhost:8317/health

# MCPO
curl http://localhost:8000/docs

# PostgreSQL
docker compose exec postgres pg_isready -U openwebui
```

### Model Availability
```bash
# From CLIProxyAPI (only if enabled)
./check-cliproxyapi.sh --models

# From OpenWebUI (advanced debug; requires an admin bearer token)
OPENWEBUI_API_KEY="<admin-bearer-token>" \
curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "http://127.0.0.1:${WEBUI_PORT:-3000}/api/models" | jq '.data[].id'
```

### Check RAG Settings
```bash
OPENWEBUI_API_KEY="<admin-bearer-token>" \
curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "http://127.0.0.1:${WEBUI_PORT:-3000}/api/v1/retrieval/" | jq '
  {
    model: .RAG_EMBEDDING_MODEL,
    top_k: .RAG_TOP_K,
    chunk_size: .CHUNK_SIZE,
    chunk_overlap: .CHUNK_OVERLAP
  }'
```

### Docker Resources
```bash
# Container status
docker compose ps

# Resource usage
docker stats --no-stream openwebui mcpo

# Disk usage
docker system df
```

---

## Common Issues

Symptom-led recovery now lives in `TROUBLESHOOTING.md`.
Use that runbook for auth drift, empty model lists, `502` recovery, SSL issues, vector-backend checks, and sidecar failures.

---

## User Management

### Promote to Admin
```bash
./openwebui-user-admin.sh --email user@example.com --role admin
```

### Activate User
```bash
./openwebui-user-admin.sh --email user@example.com --role user
```

### Reset User Password
```bash
./openwebui-user-admin.sh --email user@example.com --reset-password
```

### List Users (SQLite)
Advanced debug only. Prefer `./openwebui-user-admin.sh --email ...` for normal account management.

```bash
docker compose exec -T openwebui python3 - <<'PY'
import sqlite3

con = sqlite3.connect("/app/backend/data/webui.db")
for email, role, created_at in con.execute(
    'SELECT email, role, created_at FROM "user" ORDER BY created_at'
):
    print(f"{email}\t{role}\t{created_at}")
PY
```

---

## Maintenance

### Update Docker Images
```bash
# Check for updates
./scripts/check-image-versions.sh

# Pull and restart
docker compose pull && docker compose up -d

# Validate
./test-rag.sh --baseline
```

### Security Scanning
```bash
# Scan all images
./scripts/audit-dependencies.sh

# Scan single image
trivy image ghcr.io/open-webui/open-webui:v0.8.10
```

### Cleanup Resources
```bash
# Stop everything, clean Docker
./cleanup.sh

# Remove unused volumes
docker volume prune

# Remove unused images
docker image prune -a
```

### Verify Scripts
```bash
./verify-scripts.sh
./test-update-smoke.sh
```

---

## Config And Path References

Use the canonical sources instead of this runbook for static reference data:

- Env defaults and baseline values: `config/env/.env.example`
- Canonical config edit surface: `config/README.md`
- Runtime topology and service roles: `docs/ssot/stack.md`
- Repository layout and key paths: `docs/REPO_MAP.md`

---

## Support

- **Architecture:** [../../README.md](../../README.md)
- **API Examples:** [../../API_EXAMPLES.md](../../API_EXAMPLES.md)
- **Environment Reference:** [../reference/openwebui/openwebui_env_reference.md](../reference/openwebui/openwebui_env_reference.md)
- **Repository Map:** [../REPO_MAP.md](../REPO_MAP.md)
