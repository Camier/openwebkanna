# OpenWebUI Operations Quick Reference

**For detailed setup:** See [README.md](README.md)
**For troubleshooting:** See [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
**For prerequisites:** See [PREREQUISITES.md](PREREQUISITES.md)

---

## Quick Start

### First-Time Setup
```bash
# 1. Copy environment template
cp .env.example .env

# 2. Edit .env - set required values:
#    - WEBUI_SECRET_KEY (>=32 bytes)
#    - OPENAI_API_KEY (or CLIProxyAPI key)
#    - JUPYTER_TOKEN

# 3. Deploy stack
./deploy.sh --no-logs

# 4. Check status
./status.sh

# 5. Open browser
#    http://localhost:3000
```

### Daily Operations

#### Start Stack
```bash
./deploy.sh --no-logs
```

#### Check Status
```bash
./status.sh
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
./restart-cliproxyapi.sh
```

#### Stop Stack
```bash
docker compose down
./cleanup.sh                 # Include vLLM (if enabled)
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

**View current settings:**
```bash
curl -s http://localhost:3000/api/v1/retrieval/ | jq
```

**Apply optimized settings:**
```bash
# For precise factual retrieval
./tune-openwebui-documents.sh --top-k 8 --chunk-size 1500 --chunk-overlap 300

# For conceptual/broad queries
./tune-openwebui-documents.sh --top-k 15 --chunk-size 3000 --chunk-overlap 600
```

**Snapshot current config:**
```bash
./tune-openwebui-documents.sh --snapshot-only
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

**Check pgvector status:**
```bash
docker exec openwebui_postgres psql -U openwebui -d openwebui -c \
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

## CLIProxyAPI Operations

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
# Creates: backups/openwebui-db-YYYYMMDD-HHMMSS.tar.gz
```

### Restore from Backup
```bash
# 1. Stop stack
docker compose down

# 2. Extract backup
tar -xzf backups/openwebui-db-*.tar.gz \
  -C /var/lib/docker/volumes/openwebui_data/

# 3. Restart
docker compose up -d

# 4. Verify
curl http://localhost:3000/health
```

### Export Thesis Chats
```bash
./scripts/export-thesis-chats.sh
# Creates: thesis-exports/chats-YYYYMMDD-HHMMSS.json
```

---

## Quick Diagnostics

### Service Health
```bash
# OpenWebUI
curl http://localhost:3000/health

# CLIProxyAPI
curl http://localhost:8317/health

# MCPO
curl http://localhost:8000/docs

# PostgreSQL
docker exec openwebui_postgres pg_isready -U openwebui
```

### Model Availability
```bash
# From CLIProxyAPI
./check-cliproxyapi.sh --models

# From OpenWebUI
curl -s -H "Authorization: Bearer $API_KEY" \
  http://localhost:3000/api/models | jq '.data[].id'
```

### Check RAG Settings
```bash
curl -s http://localhost:3000/api/v1/retrieval/ | jq '
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
docker stats --no-stream openwebui cliproxyapi mcpo

# Disk usage
docker system df
```

---

## Common Issues

| Problem | Solution |
|---------|----------|
| **No models in dropdown** | `./openwebui-user-admin.sh --email you@example.com --role admin` |
| **Pending activation** | `./openwebui-user-admin.sh --email you@example.com --role user` |
| **502 errors** | `./backup-openwebui-db.sh && docker compose restart openwebui` |
| **SSL certificate errors** | Place CA at `certs/ca-bundle.pem`, set `REQUESTS_CA_BUNDLE=/certs/ca-bundle.pem` in `.env` |
| **PDF stuck in processing** | Check Docling logs: `docker compose logs --tail=50 docling` |
| **Vector DB connection failed** | Verify pgurl: `docker exec openwebui_postgres pg_isready -U openwebui` |
| **CLIProxyAPI not reachable** | Check if running: `./check-cliproxyapi.sh` |
| **MCP tools not appearing** | Apply patches: `./apply-openwebui-tool-patches.sh` |

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

### List Users
```bash
docker exec openwebui psql -U openwebui -d openwebui -c \
  "SELECT email, role, created_at FROM \"user\" ORDER BY created_at"
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
trivy image ghcr.io/open-webui/open-webui:v0.8.3
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

## Environment Variables

### Required
```bash
WEBUI_SECRET_KEY=           # >=32 bytes for session signing
OPENAI_API_KEY=             # Or CLIProxyAPI key
JUPYTER_TOKEN=              # Jupyter authentication
POSTGRES_PASSWORD=          # Database password
```

### Recommended
```bash
CLIPROXYAPI_API_KEY=        # CLIProxyAPI local key
OPENWEBUI_SIGNIN_EMAIL=     # Admin email
OPENWEBUI_SIGNIN_PASSWORD=  # Admin password
```

### RAG Configuration
```bash
RAG_EMBEDDING_MODEL=sentence-transformers/all-MiniLM-L6-v2
RAG_TOP_K=15
CHUNK_SIZE=3000
CHUNK_OVERLAP=600
```

---

## File Locations

| Location | Contents |
|----------|----------|
| `data/pdfs/` | Source PDF files |
| `data/extractions/` | Per-paper extracted data |
| `data/corpus/` | JSONL corpora |
| `data/processing/` | Processed outputs (prod_max) |
| `cliproxyapi/config.yaml` | CLIProxyAPI configuration |
| `mcp/config.json` | MCP server definitions |
| `logs/` | Runtime logs |
| `backups/` | Database backups |
| `certs/` | CA bundles for SSL |

---

## Support

- **Architecture:** [docs/reference/rag-system.md](docs/reference/rag-system.md)
- **API Examples:** [API_EXAMPLES.md](API_EXAMPLES.md)
- **Environment Reference:** [docs/reference/environment-vars.md](docs/reference/environment-vars.md)
- **Repository Map:** [AGENTS.md](AGENTS.md)
