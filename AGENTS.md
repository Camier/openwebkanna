# Repository Map: OpenWebUI + vLLM RAG Deployment

## TL;DR
This repository is an **orchestration and deployment layer** for OpenWebUI + vLLM, optimized for ethnopharmacological research. OpenWebUI runs in Docker; vLLM runs natively on the host for direct GPU access. No custom Python source code — all logic is in shell scripts and configuration files.

## How To Run

### Complete Deployment
```bash
./deploy.sh            # Starts both vLLM and OpenWebUI
```

### vLLM Server (Native — requires Conda env at ~/.conda/envs/vllm)
```bash
./start-vllm.sh        # Start vLLM server (background, logs to vllm.log)
./stop-vllm.sh         # Graceful shutdown (SIGTERM, then optional SIGKILL)
./restart-vllm.sh      # Stop + start
./check-vllm.sh        # Health check: process, port, API endpoints
```

### OpenWebUI (Docker)
```bash
docker compose up -d   # Start OpenWebUI container
docker compose down    # Stop and remove container
./status.sh            # Combined status of all services
./logs.sh              # View logs from all services
```

### Maintenance
```bash
./update.sh            # Pull latest OpenWebUI image, backup, restart
./cleanup.sh           # Stop everything, clean Docker resources
```

## How To Test
```bash
./test-api.sh          # All API endpoint validation (curl + jq based)
./test-rag.sh          # End-to-end RAG: embeddings -> vector DB -> retrieval -> generation
```
Both test scripts use counter variables and `trap` for automatic cleanup. Set `VERBOSE=true` for detailed output. Override URLs via `OPENWEBUI_URL` and `VLLM_URL` env vars.

## Repo Topology
```
openwebui/
├── *.sh                    # Operational scripts (deploy, test, manage)
├── docker-compose.yml      # OpenWebUI container definition
├── .env / .env.example     # Configuration (never commit .env)
├── vllm.log / vllm.pid     # Runtime artifacts (gitignored)
├── data/
│   ├── corpus/             # JSONL files: biblio_corpus.jsonl, chunks_corpus.jsonl
│   ├── extractions/{id}/   # Per-paper structured data:
│   │   ├── raw/            #   Processing artifacts
│   │   ├── render/         #   HTML/Markdown renditions
│   │   ├── rag/            #   Pre-computed RAG chunks
│   │   └── _assets/        #   Extracted images
│   ├── pdfs/               # Original source PDFs
│   ├── metadata/           # Supporting metadata
│   └── notes/              # Research notes
└── *.md                    # Documentation (README, SETUP, PREREQUISITES, etc.)
```

## Entrypoints
| Service | URL | Notes |
|---------|-----|-------|
| Web UI | `http://localhost:3000` | OpenWebUI (port configurable via `WEBUI_PORT`) |
| vLLM API | `http://localhost:8000/v1` | OpenAI-compatible API |
| vLLM Health | `http://localhost:8000/health` | Returns 200 when ready |
| OpenWebUI Health | `http://localhost:3000/health` | Container healthcheck |

## Architecture
- **vLLM**: Runs natively via Conda (`~/.conda/envs/vllm`). Model: `meta-llama/Llama-3.1-8B-Instruct`. GPU memory utilization: `0.25` (set in `start-vllm.sh`).
- **OpenWebUI**: Docker container (`ghcr.io/open-webui/open-webui:main`). Data persisted in `openwebui_data` volume.
- **Communication**: Container reaches host via `host.docker.internal:8000` (`extra_hosts` in compose).
- **Vector Store**: Chroma (default) or FAISS, managed within the OpenWebUI container.

## Config & Env Vars
Primary config is `.env` (created from `.env.example`). Key variables:

| Variable | Default | Purpose |
|----------|---------|---------|
| `WEBUI_PORT` | `3000` | Web interface port |
| `OPENAI_API_BASE_URL` | `http://host.docker.internal:8000/v1` | vLLM connection |
| `RAG_EMBEDDING_MODEL` | `sentence-transformers/all-MiniLM-L6-v2` | Embedding model |
| `RAG_TOP_K` | `5` | Chunks to retrieve |
| `CHUNK_SIZE` | `1500` | Text chunk size (chars) |
| `CHUNK_OVERLAP` | `150` | Overlap between chunks |
| `VECTOR_DB` | `chroma` | Vector store backend |
| `WEBUI_AUTH` | `false` | Enable authentication |

## Risky Areas
- **GPU Memory**: `GPU_MEMORY_UTILIZATION` in `start-vllm.sh`. Too high causes OOM. Currently 0.25.
- **Port Conflicts**: Ports 3000 and 8000 must be free. Grafana commonly conflicts on 3000.
- **host.docker.internal**: Requires `extra_hosts` mapping in `docker-compose.yml`. Without it, container can't reach vLLM.
- **Data Volume**: `openwebui_data` Docker volume holds all user data. Always backup via `./update.sh` before changes.

## Code Style Guidelines

### Shell Scripts (primary codebase)

**Structure** — every script follows this pattern:
```bash
#!/bin/bash
# Brief description

set -e  # ALWAYS — exit on error

# Configuration variables at top
PORT=8000
PID_FILE="vllm.pid"

# Color definitions (standard set — copy exactly)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'
BOLD='\033[1m'

# Helper functions, then main logic
```

**Helper Functions** (use these exact patterns):
- `print_header()` — ASCII box (`╔═══╗`) for script title
- `print_step()` / `print_section()` — `▶` prefix, bold blue
- `print_success()` — `✓` prefix, green
- `print_error()` — `✗` prefix, red, output to `>&2`
- `print_warning()` — yellow prefix
- `print_info()` — `ℹ` prefix, blue

**Conventions**:
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

### Configuration
- `.env.example` is canonical reference for all variables
- Group variables with `# ====` section headers and inline comments
- Never commit `.env`

### Data Management
- Preserve `data/extractions/{id}/` structure (raw, render, rag, _assets)
- Corpus files are JSONL format
- New papers: place PDFs in `data/pdfs/`, process through extraction pipeline

### Development Workflow
1. `./status.sh` — verify services before changes
2. Config changes → `docker compose down && docker compose up -d`
3. vLLM changes → edit `start-vllm.sh` → `./restart-vllm.sh`
4. Always run `./test-rag.sh` after changes
