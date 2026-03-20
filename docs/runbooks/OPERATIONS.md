# OpenWebUI Operations Quick Reference

Last updated: 2026-03-20 (UTC)

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

Canonical future RAG path during migration:

```bash
curl -fsS "${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}/health"
curl -fsS "${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}/ready"
curl -fsS -X POST "${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}/api/v1/retrieve" \
  -H 'content-type: application/json' \
  -d '{"query":"mesembrine structure","top_k":3}'
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
```

#### View Logs
```bash
./logs.sh                    # All services
./logs.sh openwebui          # Single service
docker compose --project-directory . -f config/compose/docker-compose.yml logs -f       # Follow mode
docker compose --project-directory . -f config/compose/docker-compose.yml logs --tail=100 openwebui
```

#### Restart Services
```bash
docker compose --project-directory . -f config/compose/docker-compose.yml restart openwebui
```

#### Stop Stack
```bash
docker compose --project-directory . -f config/compose/docker-compose.yml down
./cleanup.sh
```

---

## RAG Operations

Current baseline vs target:
- Current operator baseline still includes the compose-managed OpenWebUI retrieval path over `pgvector` as migration coverage.
- Canonical future RAG path is the host-native `multimodal_retrieval_api`.
- Operator scripts now validate the canonical retrieval service first, with legacy OpenWebUI retrieval kept as explicit migration coverage.

### Import PDFs

**Serial (default):**
```bash
./scripts/rag/import-pdfs-to-kb.sh --dir data/pdfs --kb-name "My Papers"
```

**Parallel (faster for large batches):**
```bash
# 3 workers
IMPORT_PDFS_PARALLELISM=3 ./scripts/rag/import-pdfs-to-kb.sh

# Or via command line
./scripts/rag/import-pdfs-to-kb.sh --parallel 3
```

**Limit to first N PDFs:**
```bash
./scripts/rag/import-pdfs-to-kb.sh --limit 5
```

**Custom directory:**
```bash
./scripts/rag/import-pdfs-to-kb.sh --dir /path/to/pdfs --kb-name "Research Papers"
```

### Tune RAG Settings

Prefer `./test-rag.sh --baseline` for the normal health signal.
For advanced debugging, inspect and update the live retrieval config directly through the OpenWebUI API when you already have an admin bearer token:

```bash
OPENWEBUI_API_KEY="<admin-bearer-token>"
curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "http://127.0.0.1:${WEBUI_PORT:-3000}/api/v1/retrieval/config" | jq
```

Prefer this deployment-sync path to keep `.env` and live retrieval DB in sync (including multimodal extraction defaults):

```bash
./scripts/admin/sync-openwebui-retrieval-config.sh
```

**Apply an explicit retrieval update:**
```bash
OPENWEBUI_API_KEY="<admin-bearer-token>"
curl -s -X POST \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  -H "Content-Type: application/json" \
  -d '{"TOP_K":8,"CHUNK_SIZE":1800,"CHUNK_OVERLAP":200,"CONTENT_EXTRACTION_ENGINE":"docling","PDF_EXTRACT_IMAGES":true}' \
  "http://127.0.0.1:${WEBUI_PORT:-3000}/api/v1/retrieval/config/update" | jq
```

The committed baseline remains `RAG_TOP_K=15`, `CHUNK_SIZE=3000`, `CHUNK_OVERLAP=600`, and `RAG_EMBEDDING_MODEL=pritamdeka/S-PubMedBert-MS-MARCO`.

### Operate The Canonical Future Retrieval API

Start or recover the host-native service with explicit settings from `.env` when you are validating or repairing the canonical RAG path:

```bash
PYTHONPATH=/LAB/@thesis/openwebui \
python3 -m uvicorn \
  services.multimodal_retrieval_api.app:app \
  --host "${MULTIMODAL_RETRIEVAL_API_HOST:-127.0.0.1}" \
  --port "${MULTIMODAL_RETRIEVAL_API_PORT:-8510}"
```

Required host-native settings:
- `MULTIMODAL_RETRIEVAL_API_QDRANT_URL`
- `MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH`

Useful direct checks:

```bash
curl -fsS "${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}/health" | jq
curl -fsS "${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}/ready" | jq
curl -fsS -X POST "${MULTIMODAL_RETRIEVAL_API_URL:-http://127.0.0.1:8510}/api/v1/retrieve" \
  -H 'content-type: application/json' \
  -d '{"query":"mesembrine structure","top_k":3}' | jq
```

Interpretation:
- Treat this as the canonical retrieval contract for the repo.
- Do not treat green legacy OpenWebUI retrieval checks as proof that the canonical service is healthy.

**Snapshot current retrieval + embedding config before changing it:**
```bash
OPENWEBUI_API_KEY="<admin-bearer-token>"
jq -n \
  --argjson retrieval "$(curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" "http://127.0.0.1:${WEBUI_PORT:-3000}/api/v1/retrieval/config")" \
  --argjson embedding "$(curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" "http://127.0.0.1:${WEBUI_PORT:-3000}/api/v1/retrieval/embedding")" \
  '{retrieval:$retrieval, embedding:$embedding}' \
  > "backups/openwebui-retrieval-$(date -u +%Y%m%dT%H%M%SZ).json"
```

**Check whether indexed chunks match the active embedding model:**
```bash
docker compose --project-directory . -f config/compose/docker-compose.yml exec postgres psql -U openwebui -d openwebui -c \
  "SELECT DISTINCT collection_name, vmetadata->>'embedding_config' AS embedding_config FROM document_chunk ORDER BY collection_name;"
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
docker compose --project-directory . -f config/compose/docker-compose.yml exec postgres psql -U openwebui -d openwebui -c \
  "SELECT COUNT(*) FROM collection_embeddings"
```

**Backup vector data:**
```bash
./scripts/admin/backup-openwebui-db.sh
```

## MCP Operations

### Configure MCP Servers
```bash
# Verify MCPO endpoints
./scripts/mcp/configure-mcpo-openapi-servers.sh --verify

# Apply configuration
./scripts/mcp/configure-mcpo-openapi-servers.sh
```

### Test MCP Integration
```bash
./scripts/mcp/test-mcp.sh
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
./scripts/admin/backup-openwebui-db.sh
# Creates: backups/webui.db.YYYYMMDD-HHMMSS.sqlite
```

### Restore from Backup
```bash
# 1. Select backup file to restore
BACKUP_FILE="backups/webui.db.YYYYMMDD-HHMMSS.sqlite"

# 2. Optional but recommended: take a pre-restore safety backup
./scripts/admin/backup-openwebui-db.sh --output "backups/webui.db.pre-restore.$(date +%Y%m%d-%H%M%S).sqlite"

# 3. Stop OpenWebUI to avoid writes during restore
docker compose --project-directory . -f config/compose/docker-compose.yml stop openwebui

# 4. Restore SQLite DB file inside the OpenWebUI container
docker cp "${BACKUP_FILE}" "$(docker compose --project-directory . -f config/compose/docker-compose.yml ps -q openwebui)":/app/backend/data/webui.db

# 5. Start OpenWebUI and verify health
docker compose --project-directory . -f config/compose/docker-compose.yml start openwebui
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

# MCPO
curl http://localhost:8000/docs

# PostgreSQL
docker compose --project-directory . -f config/compose/docker-compose.yml exec postgres pg_isready -U openwebui
```

### Model Availability
```bash
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
docker compose --project-directory . -f config/compose/docker-compose.yml ps

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
./scripts/admin/openwebui-user-admin.sh --email user@example.com --role admin
```

### Activate User
```bash
./scripts/admin/openwebui-user-admin.sh --email user@example.com --role user
```

### Reset User Password
```bash
./scripts/admin/openwebui-user-admin.sh --email user@example.com --reset-password
```

### List Users (SQLite)
Advanced debug only. Prefer `./scripts/admin/openwebui-user-admin.sh --email ...` for normal account management.

```bash
docker compose --project-directory . -f config/compose/docker-compose.yml exec -T openwebui python3 - <<'PY'
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
docker compose --project-directory . -f config/compose/docker-compose.yml pull && \
docker compose --project-directory . -f config/compose/docker-compose.yml up -d

# Validate
./test-rag.sh --baseline
```

### Security Scanning
```bash
# Scan all images with trivy
trivy image --severity CRITICAL,HIGH,MEDIUM,LOW,UNKNOWN \
  pgvector/pgvector:pg16 \
  ghcr.io/open-webui/open-webui:v0.8.10 \
  quay.io/jupyter/scipy-notebook:2026-02-19

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
./scripts/testing/verify-scripts.sh
```

---

## Config And Path References

Use the canonical sources instead of this runbook for static reference data:

- Env defaults and baseline values: `config/env/.env.example`
- Canonical config edit surface: `config/README.md`
- Runtime topology and service roles: `docs/ssot/stack.md`
- Repository layout and key paths: `docs/REPO_MAP.md`

---

## Deploy preflight recovery

`./deploy.sh` includes a Docker runtime preflight that checks cgroup/container startup before full stack launch.

If preflight fails with a cgroup error:
- confirm Docker daemon starts cleanly and `/etc/docker/daemon.json` is explicit about the host cgroup contract
- on this deployment path, keep this runtime stanza in daemon config for stability:
  - `"exec-opts": ["native.cgroupdriver=cgroupfs"]`
- restart Docker and rerun `./deploy.sh`
- then run:
  - `./status.sh`
  - `./test-api.sh --baseline`
  - `./test-rag.sh --baseline`

---

## Support

- **Architecture:** [../../README.md](../../README.md)
- **API Examples:** [../../API_EXAMPLES.md](../../API_EXAMPLES.md)
- **Environment Reference:** [../reference/openwebui/openwebui_env_reference.md](../reference/openwebui/openwebui_env_reference.md)
- **Repository Map:** [../REPO_MAP.md](../REPO_MAP.md)
