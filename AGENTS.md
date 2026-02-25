# Repository Map: OpenWebUI + LiteLLM RAG Deployment

## TL;DR
This repository is an **orchestration and deployment layer** for OpenWebUI + LiteLLM, optimized for ethnopharmacological research (specifically Sceletium tortuosum and related medicinal plants). OpenWebUI, LiteLLM, and supporting services run in Docker. Primary operations are shell-first, with supporting Python utilities for extraction, validation, and pipeline quality checks.

## Project Overview

### Purpose
Deploy and manage a complete RAG (Retrieval-Augmented Generation) system for academic research on ethnopharmacological literature. The system ingests PDFs, processes them into chunks, stores vectors, and enables semantic search and question-answering via LLMs through OpenWebUI.

### Architecture Diagram
```text
┌─────────────────┐      ┌──────────────────┐      ┌─────────────────────────┐
│   OpenWebUI     │──────▶  LiteLLM Proxy   │──────▶  Upstream Providers    │
│   (Docker)      │      │   (Docker)       │      │  - openrouter/*         │
│   Port 3000     │      │   Port 4000      │      │  - deepseek/*           │
└─────────────────┘      └──────────────────┘      │  - zai-glm-*            │
         │                                         │  - minimax-*            │
         │                                         │  - cliproxy/*           │
         ▼                                         │  - ollama-cloud/*       │
┌─────────────────┐                              └─────────────────────────┘
│  PostgreSQL     │
│  (pgvector)     │      ┌──────────────────┐
└─────────────────┘      │  CLIProxyAPI     │  (optional sidecar, port 8317)
         │               │  OAuth aliases   │
         ▼               └──────────────────┘
┌─────────────────┐
│     Jupyter     │
│ (code execution)│
└─────────────────┘
```

### Technology Stack

| Component | Technology | Version/Notes |
|-----------|-----------|---------------|
| Container Runtime | Docker + Docker Compose | v2 plugin preferred |
| Web UI | OpenWebUI | v0.8.3 (pinned) |
| API Proxy | LiteLLM | External Docker container, port 4000 |
| Vector Database | PostgreSQL + pgvector | pg16 |
| Code Execution | Jupyter (scipy-notebook) | Latest (configurable) |
| Scripting | Bash | All scripts use `set -e` |
| Optional sidecar | CLIProxyAPI | v6.8.18 (OAuth alias sidecar, not primary) |

## Directory Structure

```
openwebui/
├── *.sh                         # 40+ operational scripts (deploy, test, manage)
├── docker-compose.yml           # Container definitions (5 services)
├── docker-compose.rg.yml        # Alternative compose configuration
├── .env / .env.example          # Configuration (NEVER commit .env)
├── Dockerfile.openwebui-rg      # Custom OpenWebUI image build
│
├── cliproxyapi/                 # CLIProxyAPI configuration
│   ├── config.yaml              # Provider and alias mapping (tracked)
│   ├── config.local.yaml        # Local secrets (gitignored)
│   ├── config.upstream.example.yaml
│   ├── auth/                    # OAuth credentials (gitignored)
│   ├── logs/                    # Runtime logs (gitignored)
│   └── vendor/                  # Vendor files (gitignored)
│
├── lib/                         # Shared shell script libraries
│   ├── init.sh                 # Single-source loader (sources all below)
│   ├── colors.sh                # ANSI color definitions
│   ├── print-utils.sh           # Standard print functions
│   ├── docker-helpers.sh        # Docker Compose abstraction
│   ├── env-loader.sh            # .env file parser
│   ├── validation.sh            # Validation utilities
│   ├── cliproxyapi-helpers.sh   # CLIProxyAPI config parsing
│   ├── http-helpers.sh          # HTTP request helpers
│   └── audit_plugins.py         # Python plugin auditor
│
├── scripts/                     # Additional utilities
│   ├── audit-dependencies.sh    # Security scanning wrapper
│   └── check-image-versions.sh  # Version checking
│
├── data/                        # Research data (gitignored)
│   ├── corpus/                  # JSONL corpus files
│   │   ├── biblio_corpus.jsonl  # Bibliographic metadata
│   │   ├── chunks_corpus.jsonl  # RAG chunks (~40MB)
│   │   └── images_manifest.json
│   ├── extractions/{id}/        # Per-paper structured data
│   │   ├── raw/                 # Processing artifacts
│   │   ├── render/              # HTML/Markdown renditions
│   │   ├── rag/                 # Pre-computed RAG chunks
│   │   └── _assets/             # Extracted images
│   ├── pdfs/                    # Original source PDFs
│   ├── metadata/                # Supporting metadata
│   └── notes/                   # Research notes
│
├── bin/                         # CLIProxyAPI binary (symlink + executable)
├── mcp/                         # MCP server configuration
│   ├── config.json              # MCPO server definitions
│   └── config.research.optional.json
├── thesis-exports/              # Thesis export output directory
├── searxng/                     # SearXNG configuration (optional)
├── jupyter/                     # Jupyter configuration
│   └── jupyter_server_config.py
├── certs/                       # CA bundles for SSL (gitignored)
├── backups/                     # DB backups (gitignored)
├── logs/                        # Runtime logs (gitignored)
│   └── llm-council/             # LLM council evaluation outputs
└── *.md                         # Documentation files
    ├── README.md                # Main documentation
    ├── SETUP.md                 # Setup guide
    ├── PREREQUISITES.md         # Prerequisites checklist
    ├── TROUBLESHOOTING.md       # Problem resolution
    ├── SECURITY.md              # Security documentation
    └── API_EXAMPLES.md          # API usage examples
```

## Key Configuration Files

### `.env` (created from `.env.example`)
Primary configuration file. Key variables:

| Variable | Default | Purpose |
|----------|---------|---------|
| `WEBUI_PORT` | `3000` | Web interface port |
| `WEBUI_SECRET_KEY` | (required) | Session signing key (>=32 bytes) |
| `OPENAI_API_BASE_URL` | `http://cliproxyapi:8317/v1` | CLIProxyAPI connection (Docker DNS) |
| `OPENAI_API_KEY` | (required) | API key for CLIProxyAPI |
| `CLIPROXYAPI_ENABLED` | `true` | Enable CLIProxyAPI |
| `CLIPROXYAPI_DOCKER_MANAGED` | `true` | Docker-managed CLIProxyAPI |
| `CLIPROXYAPI_IMAGE` | `eceasy/cli-proxy-api:v6.8.18` | CLIProxyAPI image |
| `OPENWEBUI_IMAGE` | `ghcr.io/open-webui/open-webui:v0.8.3` | OpenWebUI image |
| `RAG_EMBEDDING_MODEL` | `sentence-transformers/all-MiniLM-L6-v2` | Embedding model |
| `RAG_TOP_K` | `5` | Chunks to retrieve |
| `CHUNK_SIZE` | `1500` | Text chunk size (chars) |
| `CHUNK_OVERLAP` | `150` | Overlap between chunks |
| `VECTOR_DB` | `chroma` | Vector store (chroma/pgvector) |
| `PGVECTOR_DB_URL` | (auto) | PostgreSQL connection string |
| `JUPYTER_TOKEN` | (required) | Jupyter authentication token |
| `ENABLE_CODE_EXECUTION` | `true` | Enable code execution |
| `ENABLE_WEB_SEARCH` | `false` | Enable web search via SearXNG |

### `docker-compose.yml`
Defines 5 services:
1. **postgres**: PostgreSQL 16 with pgvector extension
2. **jupyter**: Jupyter notebook for code execution
3. **cliproxyapi**: CLIProxyAPI OAuth proxy
4. **mcpo**: MCP-to-OpenAPI proxy (v0.0.19)

5. **openwebui**: OpenWebUI web interface
### `cliproxyapi/config.yaml`
CLIProxyAPI configuration. Contains:
- Server settings (host, port, TLS)
- Local API keys
- Provider configurations (Z.ai glm-5, MiniMax MiniMax-M2.5)
- OAuth model aliases:
  - codex: gpt-5-codex → openai-codex
  - qwen: qwen3-coder-plus → qwen-cli
  - kimi: kimi-k2.5 → kimi-cli
- Routing strategy: round-robin
**IMPORTANT**: This file contains placeholders. Real credentials go in `config.local.yaml` (gitignored).

## How To Run

### Quick Start (Complete Deployment)
```bash
# 1. Prepare environment
cp .env.example .env
# Edit .env to set required values (WEBUI_SECRET_KEY, OPENAI_API_KEY, JUPYTER_TOKEN)

# 2. Deploy everything
./deploy.sh --no-logs

# 3. Check status
./status.sh
./check-cliproxyapi.sh
```

### Service Management

#### CLIProxyAPI (Docker)
```bash
./start-cliproxyapi.sh       # Start CLIProxyAPI container
./check-cliproxyapi.sh       # Health check: process, port, API endpoints
./stop-cliproxyapi.sh        # Graceful shutdown
./restart-cliproxyapi.sh     # Stop + start
```

#### OpenWebUI (Docker)
```bash
docker compose up -d         # Start OpenWebUI stack
docker compose down          # Stop and remove containers
./status.sh                  # Combined status of all services
./logs.sh                    # View logs from all services
```

#### vLLM Server (Optional Fallback)
```bash
./start-vllm.sh              # Start vLLM server (background, logs to vllm.log)
./stop-vllm.sh               # Graceful shutdown
./restart-vllm.sh            # Stop + start
./check-vllm.sh              # Health check
```

### Maintenance
```bash
./update.sh                  # Pull latest images, backup, restart
./cleanup.sh                 # Stop everything, clean Docker resources
./backup-openwebui-db.sh     # Create DB snapshot
```

## How To Test

### API and Integration Tests
```bash
./test-api.sh                           # All API endpoint validation
./test-api.sh --baseline                # Fast baseline checks only

./test-rag.sh                           # End-to-end RAG flow test
./test-rag.sh --baseline                # Fast baseline RAG checks

./test-cliproxyapi-oauth.sh             # OAuth alias validation
./test-openwebui-cliproxy-routing.sh    # OpenWebUI-to-CLIProxyAPI E2E routing
./test-openwebui-tools-endpoints.sh     # Tool endpoint validation
./test-openwebui-tools-invocation.sh    # Tool execution verification
```

### Baseline Validation Workflow
```bash
# 1. Deploy baseline
./deploy.sh --no-logs

# 2. Validate runtime
./status.sh

# 3. Enforce real-integration guard
./audit-no-mock.sh

# 4. Run fast baseline checks
./test-rag.sh --baseline
./test-api.sh --baseline
```

## Entrypoints and URLs

| Service | URL | Notes |
|---------|-----|-------|
| Web UI | `http://localhost:3000` | Port configurable via `WEBUI_PORT` |
| CLIProxyAPI | `http://localhost:8317/v1` | OpenAI-compatible API |
| CLIProxyAPI Health | `http://localhost:8317/health` | Returns 200 when ready |
| OpenWebUI Health | `http://localhost:3000/health` | Container healthcheck |
| Jupyter | `http://localhost:8890` | Token authentication required |
| MCPO | `http://localhost:8000` | MCP-to-OpenAPI proxy |

## Code Style Guidelines

### Shell Scripts

**Structure** — every script follows this pattern:
```bash
#!/bin/bash
# Brief description

set -e  # ALWAYS — exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Configuration variables at top
PORT=8317
PID_FILE="vllm.pid"

# Source shared libraries
source "${SCRIPT_DIR}/lib/colors.sh"
source "${SCRIPT_DIR}/lib/print-utils.sh"

# Helper functions, then main logic
```

**Standard Colors** (from `lib/colors.sh`):
```bash
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'      # No color
BOLD='\033[1m'
```

**Helper Functions** (from `lib/print-utils.sh`):
- `print_header "Title"` — ASCII box for script title
- `print_step "Doing something"` — `>` prefix, bold blue
- `print_success "It worked"` — `[OK]` prefix, green
- `print_error "Failed"` — `[ERROR]` prefix, red to stderr
- `print_warning "Caution"` — `[WARN]` prefix, yellow
- `print_info "Note"` — `[INFO]` prefix, cyan
- `print_section "Name"` — `--- Section Name ---` divider

**Key Conventions**:
- Use `${VAR:-default}` for optional env vars with defaults
- Manage background processes via PID files; check for stale PIDs before starting
- Check port availability with `lsof -Pi :PORT` before binding
- Verify dependencies (`docker`, `curl`, `jq`) before main logic
- Use `nohup` with log redirection for background processes
- Section headers use `###...###` comment borders

### Docker / Compose
- Compose version `3.8`
- Name containers, volumes, and networks explicitly (no auto-generated names)
- Define healthchecks with `curl -f` against health endpoints
- Layer env vars: `env_file: .env` + explicit `environment:` with `${VAR:-default}`
- Use `extra_hosts` for `host.docker.internal` access

### Configuration Management
- `.env.example` is the canonical reference for all variables
- Group variables with `# ====` section headers and inline comments
- Never commit `.env` or any file containing real credentials
- OAuth credentials in `cliproxyapi/auth/` are gitignored
- Generate local configs from `~/.007` using `configure-cliproxyapi-providers.sh`

### Data Management
- Preserve `data/extractions/{id}/` structure (raw, render, rag, _assets)
- Corpus files are JSONL format with specific schema
- New papers: place PDFs in `data/pdfs/`, process through extraction pipeline
- The `data/` directory is gitignored; only commit structure documentation

## Testing Strategy

### Test Categories
1. **Unit/Script Tests**: `verify-scripts.sh`, `test-update-smoke.sh`
2. **Integration Tests**: `test-api.sh`, `test-rag.sh`
3. **OAuth Tests**: `test-cliproxyapi-oauth.sh`
4. **E2E Tests**: `test-openwebui-cliproxy-routing.sh`
5. **Tool Tests**: `test-openwebui-tools-endpoints.sh`, `test-openwebui-tools-invocation.sh`

### Test Patterns
- All test scripts use counter variables for pass/fail tracking
- Use `trap` for automatic cleanup on exit
- Set `VERBOSE=true` for detailed output
- Override URLs via environment variables (e.g., `OPENWEBUI_URL`)
- Baseline mode (`--baseline`) for fast, reproducible checks

### Continuous Validation
```bash
# Pre-commit validation
bash -n *.sh                    # Syntax check
./verify-scripts.sh             # Full script verification

# Post-deploy validation
./test-rag.sh --baseline
./test-api.sh --baseline
./audit-no-mock.sh
```

## Security Considerations

### Critical Security Rules
1. **Never commit `.env`** — Contains secrets and credentials
2. **Never commit `cliproxyapi/auth/`** — Contains OAuth tokens
3. **Use strong secrets**: `WEBUI_SECRET_KEY` must be >= 32 bytes
4. **Rotate keys periodically**: Use `rotate-cliproxyapi-local-key.sh`
5. **Scan images regularly**: Use `scripts/audit-dependencies.sh`

### Security Features
- API key rotation script: `rotate-cliproxyapi-local-key.sh`
- Container security scanning: `scripts/audit-dependencies.sh`
- SBOM generation and CVE tracking in `SECURITY.md`
- CA bundle mounting for corporate proxies: `./certs/`

### CVE Status (as of 2026-02-18)
| Image | Version | CRITICAL | HIGH | Status |
|-------|---------|----------|------|--------|
| OpenWebUI | v0.8.3 | 5 | 178 | Under Review |
| CLIProxyAPI | v6.8.18 | 0 | - | Updated |
| pgvector | pg16 | 2 | 6 | Under Review |

### Response Timeframes
| Severity | Maximum Response Time |
|----------|----------------------|
| CRITICAL | 7 days |
| HIGH | 30 days |
| MEDIUM | 90 days |
| LOW | Next maintenance window |

## Risky Areas and Common Issues

### Port Conflicts
- **Port 3000**: OpenWebUI (commonly conflicts with Grafana)
- **Port 8317**: CLIProxyAPI
- **Port 8888**: SearXNG (host-local)
- **Port 8890**: Jupyter
- **Port 8000**: MCPO (MCP-to-OpenAPI proxy)

Check with: `lsof -nP -iTCP:PORT -sTCP:LISTEN`

### Docker DNS
OpenWebUI must use `http://cliproxyapi:8317/v1` (Docker DNS), NOT `host.docker.internal`. Both containers must be on the same Docker network (`openwebui_network`).

### OAuth Credentials
CLIProxyAPI requires valid OAuth auth state in `cliproxyapi/auth/`. Manual OAuth setup is expected. If missing, run:
```bash
./configure-cliproxyapi-oauth.sh
```

### Data Volume
`openwebui_data` Docker volume holds all user data. Always backup before changes:
```bash
./backup-openwebui-db.sh
```

### vLLM Fallback
Only used when CLIProxyAPI is unhealthy AND `AUTO_START_VLLM_ON_CLIPROXYAPI_FAILURE=true`. Requires Conda env at `~/.conda/envs/vllm`.

## Development Workflow

### Standard Development Cycle
```bash
# 1. Verify services before changes
./status.sh

# 2. Make config changes
# Edit .env or cliproxyapi/config.yaml

# 3. Config changes → restart containers
docker compose down && docker compose up -d

# 4. CLIProxyAPI config changes
# Edit cliproxyapi/config.yaml
./restart-cliproxyapi.sh

# 5. vLLM changes (if using fallback)
# Edit start-vllm.sh
./restart-vllm.sh

# 6. Always run RAG tests after changes
./test-rag.sh --baseline
```

### Adding New Scripts
1. Follow the shell script structure template
2. Source shared libraries from `lib/`
3. Use standard print functions
4. Add to `verify-scripts.sh` validation
5. Document in this file

### Version Updates
```bash
# Check for updates
./scripts/check-image-versions.sh

# Update pinned version in .env
sed -i 's|^OPENWEBUI_IMAGE=.*|OPENWEBUI_IMAGE=ghcr.io/open-webui/open-webui:vX.Y.Z|' .env

# Pull and restart
docker compose pull && docker compose up -d

# Validate
./test-rag.sh --baseline
```

## Utility Scripts Reference

### User Management
- `openwebui-user-admin.sh` — Promote/activate/reset users in OpenWebUI DB

### Provider Configuration
- `configure-ollama-cloud.sh` — Configure OpenWebUI to use Ollama Cloud
- `configure-cliproxyapi-providers.sh` — Generate local config with real provider keys
- `configure-cliproxyapi-oauth.sh` — Interactive OAuth setup

### Key Management
- `rotate-cliproxyapi-local-key.sh` — Rotate API keys between OpenWebUI and CLIProxyAPI
- `sync-openwebui-openai-config.sh` — Sync DB config with .env values

### Tool Management
- `audit-openwebui-plugins.sh` — Audit tool/function code in webui.db
- `repair-openwebui-tools.sh` — Repair common tool issues
- `apply-openwebui-tool-patches.sh` — Apply patches via admin API

### Document Tuning
- `tune-openwebui-documents.sh` — Tune RAG settings with snapshot/restore
- `import-pdfs-to-kb.sh` — Import PDFs with optional parallelism

### LLM Council
- `llm-council.sh` — Multi-model evaluation through CLIProxyAPI

## Data Schema

### Corpus JSONL Format
```json
{
  "doc_id": "sha256_hash",
  "paper_id": "Year - Author - Title",
  "citekey": "AuthorYear_Title_Snippet",
  "title": "Full Paper Title",
  "authors": ["Author Name"],
  "year": 2024,
  "doi": "10.xxxx/...",
  "venue": "Journal Name",
  "publisher": "Publisher Name",
  "chunk_count": 100,
  "image_count": 5,
  "extraction_path": "data/extractions/..."
}
```

## Troubleshooting Quick Reference

### Models disappeared / Empty dropdown
```bash
# Clear browser cookies/localStorage for OpenWebUI origin
# OR server-side fix:
./openwebui-user-admin.sh --email you@example.com --role admin
```

### Pending activation
```bash
./openwebui-user-admin.sh --email you@example.com --role user
```

### 502 / Container crash-loop
```bash
./backup-openwebui-db.sh
docker compose logs --tail 200 openwebui
docker compose restart openwebui
```

### SSL certificate errors
1. Place CA bundle at `certs/ca-bundle.pem`
2. Set in `.env`: `REQUESTS_CA_BUNDLE=/certs/ca-bundle.pem`
3. Restart: `docker compose up -d --force-recreate openwebui`

### Web search fails from container
```bash
# Allow docker-to-host traffic (example for ufw)
sudo ufw allow from 172.16.0.0/12 to any port 8888
```

## External Dependencies

### Required
- Docker daemon
- Docker Compose plugin (v2)
- `curl`
- `jq`

### Optional
- vLLM Python environment (`~/.conda/envs/vllm`)
- SearXNG (host-local on port 8888)
- Trivy (for security scanning)

## License and Attribution

This repository contains deployment orchestration for:
- [OpenWebUI](https://github.com/open-webui/open-webui) — MIT License
- [CLIProxyAPI](https://hub.docker.com/r/eceasy/cli-proxy-api) — Commercial/OAuth

Research data in `data/` is for ethnopharmacological research on Sceletium tortuosum and related South African medicinal plants.

---

*Last updated: 2026-02-21 (UTC)*
*Repository: /LAB/@thesis/openwebui*
