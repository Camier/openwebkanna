# OpenWebUI + CLIProxyAPI RAG Deployment

Last updated: 2026-02-12 (UTC)

This repository is an orchestration layer for an academic-papers RAG workflow in OpenWebUI with a real CLIProxyAPI backend.

## Current architecture

- OpenWebUI runs in Docker (`openwebui` service).
- CLIProxyAPI runs in Docker (`cliproxyapi` service) and exposes an OpenAI-compatible API.
- OpenWebUI reaches CLIProxyAPI through Docker DNS: `http://cliproxyapi:8317/v1`.
- vLLM is optional host-side fallback and is not started when CLIProxyAPI is healthy.

```text
OpenWebUI (container) --> CLIProxyAPI (container) --> OAuth-backed providers
                                         \
                                          \--> optional host vLLM fallback
```

## What is production-relevant here

- No mock, dummy, or monkey fallback path for integration checks.
- OAuth aliases validated for:
1. `openai-codex`
2. `qwen-cli`
3. `kimi-cli`
- Manual OAuth login is expected and supported.

## Prerequisites

- Docker + Docker Compose plugin
- `curl`, `jq`
- Local OAuth credentials prepared manually for providers you use

Optional:
- vLLM Python environment (only if you intentionally enable fallback)

## Quick start

1. Copy env file:
```bash
cp .env.example .env
```

2. Ensure these values are set in `.env`:
```bash
OPENAI_API_BASE_URL=http://cliproxyapi:8317/v1
OPENAI_API_BASE_URLS=http://cliproxyapi:8317/v1
CLIPROXYAPI_DOCKER_MANAGED=true
CLIPROXYAPI_ENABLED=true
```

3. Start stack:
```bash
./deploy.sh --no-logs
```

4. Check health:
```bash
./status.sh
docker compose ps
./check-cliproxyapi.sh
```

## Troubleshooting

If the UI shows "no models", "pending activation", or you get 401/502 surprises, use:
- `TROUBLESHOOTING.md`
- `./backup-openwebui-db.sh` (DB snapshot)
- `./openwebui-user-admin.sh` (promote/activate/reset user)

## OpenWebUI baseline enhancements (as of 2026-02-17)

Before enabling advanced features, baseline defaults now include:

- Local SearXNG via host bridge (`host.docker.internal:${SEARXNG_PORT}`).
- Code execution + code-interpreter bootstrap enabled using Pyodide.
- OpenWebUI image aligned with current upstream release reference:
  `ghcr.io/open-webui/open-webui:v0.8.3`.

If you need a current release tag from upstream manually:

```bash
curl -s https://api.github.com/repos/open-webui/open-webui/releases/latest | jq -r '.tag_name'
```

Update `.env` / `.env.example` only after reviewing compatibility.

To pin your current env to that tag quickly (optional):

```bash
OPENWEBUI_TAG="$(curl -s https://api.github.com/repos/open-webui/open-webui/releases/latest | jq -r '.tag_name')"
sed -i "s|^OPENWEBUI_IMAGE=.*|OPENWEBUI_IMAGE=ghcr.io/open-webui/open-webui:${OPENWEBUI_TAG}|" .env
```

Use this baseline workflow before enabling advanced features:

1. Deploy baseline stack:
```bash
./deploy.sh --no-logs
```

2. Validate runtime baseline:
```bash
./status.sh
```

3. Enforce real-integration guard (no mock-like markers):
```bash
./audit-no-mock.sh
```

4. Run fast baseline RAG checks:
```bash
./test-rag.sh --baseline
```

5. Run fast baseline API checks:
```bash
./test-api.sh --baseline
```

`--baseline` keeps tests reproducible and quick by validating health, models, retrieval, and web-search wiring without running heavier full regression flows.

If local SearXNG is temporarily unavailable, baseline checks continue by default.
Host-side baseline probes automatically try `127.0.0.1` when `SEARXNG_QUERY_URL` uses `host.docker.internal`.
To enforce strict web-search readiness, run:

```bash
BASELINE_REQUIRE_WEB_SEARCH=true ./test-rag.sh --baseline
```

Note: OpenWebUI runs in Docker, so web-search must work from the container path too.
`./test-rag.sh --baseline` now probes both host-local and `docker exec openwebui` connectivity; strict mode enforces both.

## Dependency Management

### Checking for Updates

Check which Docker images have newer versions available:

```bash
./scripts/check-image-versions.sh
```

This compares locally pinned versions against upstream registry tags for:
- OpenWebUI (`ghcr.io/open-webui/open-webui`)
- CLIProxyAPI (if using a versioned image)
- Any other images defined in `docker-compose.yml`

To check a specific image manually:

```bash
docker run --rm quay.io/skopeo/stable list-tags docker://ghcr.io/open-webui/open-webui | jq '.Tags[-10:]'
```

### Security Scanning

Run vulnerability scans on container images:

```bash
./scripts/audit-dependencies.sh
```

This wrapper runs Trivy against all images in the compose stack. For a single image:

```bash
trivy image ghcr.io/open-webui/open-webui:v0.8.3
```

For a deeper SBOM-based scan:

```bash
trivy image --format spdx-json --output sbom.json ghcr.io/open-webui/open-webui:v0.8.3
trivy sbom sbom.json
```

Review CVE severity levels and prioritize:
- CRITICAL/HIGH: patch or upgrade promptly
- MEDIUM/LOW: assess exploitability in your context

### Updating Docker Images

1. Check for updates:
   ```bash
   ./scripts/check-image-versions.sh
   ```

2. Review changelog for breaking changes:
   - OpenWebUI: https://github.com/open-webui/open-webui/releases
   - CLIProxyAPI: check upstream release notes

3. Update image tags in `.env` or `docker-compose.yml`:
   ```bash
   # Example: pin to a specific OpenWebUI release
   sed -i 's|^OPENWEBUI_IMAGE=.*|OPENWEBUI_IMAGE=ghcr.io/open-webui/open-webui:v0.8.4|' .env
   ```

4. Pull and restart:
   ```bash
   docker compose pull && docker compose up -d
   ```

5. Validate baseline:
   ```bash
   ./status.sh
   ./test-api.sh --baseline
   ./test-rag.sh --baseline
   ```

6. If issues occur, rollback:
   ```bash
   # Revert to previous tag in .env, then:
   docker compose pull && docker compose up -d
   ```

### Version Pinning Policy

**When to pin versions:**
- Production deployments: always pin to specific tags (e.g., `v0.8.3`, not `latest`)
- After validated upgrades: update the pinned version in `.env.example`
- When reproducibility matters: CI/CD pipelines, shared environments

**When `latest` is acceptable:**
- Local development with frequent rebuilds
- Throwaway test environments
- When you explicitly want auto-updates on pull

**Recommended pattern in `.env`:**
```bash
# Pinned versions (update after validation)
OPENWEBUI_IMAGE=ghcr.io/open-webui/open-webui:v0.8.3
CLIPROXYAPI_IMAGE=your-registry/cliproxyapi:1.2.0
```

**Avoid:**
- Unpinned `latest` in production
- Digest pinning (`@sha256:...`) unless you have automated digest updates
- Mixing pinned and unpinned images in the same stack

## OAuth setup flow (manual)

If you already completed OAuth locally, keep using the existing auth state under `cliproxyapi/auth/`.

If you need to (re)configure:
```bash
./configure-cliproxyapi-oauth.sh
```

Then validate aliases and chat routing:
```bash
./test-cliproxyapi-oauth.sh
```

## OpenWebUI routing verification

Run end-to-end routing check from OpenWebUI to CLIProxyAPI:

```bash
./test-openwebui-cliproxy-routing.sh
```

This verifies:
- OpenWebUI can reach `http://cliproxyapi:8317/v1/models`
- OpenWebUI model list contains required aliases
- Chat completions route through each alias

## LLM Council

LLM Council is a multi-model evaluation system with candidate/judge pattern:
- **Candidates** generate responses to prompts
- **Judges** evaluate and vote on the best response
- Reports are generated in `logs/llm-council/`

### Configuration

Set defaults in `.env` (see `.env.example` for all options):

```bash
COUNCIL_MODELS=glm-5 minimax/chat-elite
COUNCIL_JUDGES=glm-5 minimax/chat-elite
COUNCIL_JUDGE_OUTPUT_FORMAT=json
COUNCIL_JUDGE_RETRIES=2
COUNCIL_POSITION_SWAP=false
```

### CLI Usage

Run multi-model council evaluation through CLIProxyAPI:

```bash
./llm-council.sh --models "glm-5 minimax/chat-elite" --judges "glm-5 minimax/chat-elite" \
  --prompt "Propose a robust rollout plan for OpenWebUI upgrades."
```

For reasoning-heavy judges (for example `minimax/chat-elite`), increase strict vote retries:

```bash
./llm-council.sh --models "glm-5 minimax/chat-elite" --judges "glm-5 minimax/chat-elite" \
  --prompt "Reply with exactly council-ok" \
  --judge-max-tokens 256 --judge-retries 3 --judge-force-max-tokens 96
```

Enable anti-position-bias voting (A/B then B/A consistency check):

```bash
./llm-council.sh --models "glm-5 minimax/chat-elite" --judges "glm-5 minimax/chat-elite" \
  --prompt "Compare two rollout plans and pick the better one." \
  --position-swap
```

Use strict JSON judge output with schema-like validation:

```bash
./llm-council.sh --models "glm-5 minimax/chat-elite" --judges "glm-5 minimax/chat-elite" \
  --prompt "Reply with exactly council-ok" \
  --judge-output-format json --judge-retries 3
```

Run a prompt set from file (one prompt per line):

```bash
./llm-council.sh --models "glm-5 minimax/chat-elite" --prompts-file ./data/notes/council_prompts.txt
```

### OpenWebUI Tool

To use LLM Council as an in-chat tool:

1. Apply the tool patch:
   ```bash
   ./apply-openwebui-tool-patches.sh
   ```

2. In chat, use the `council` function:
   ```
   Use the council tool to answer: What are the tradeoffs of RAG vs fine-tuning?
   ```

3. Configure via Valves in OpenWebUI admin:
   - `MODELS`: Comma-separated model IDs
   - `API_BASE_URLS`: Provider endpoints (semicolon-separated)
   - `TIMEOUT`: Request timeout in seconds

### Outputs

- Console summary with per-model vote totals
- Markdown report in `logs/llm-council/` with:
  - CLIProxyAPI version and configuration
  - Per-prompt votes and winners
  - Full candidate responses (optional)
  - Judge validation statistics

## Service commands

- Deploy everything:
```bash
./deploy.sh
```

- CLIProxyAPI lifecycle:
```bash
./start-cliproxyapi.sh
./check-cliproxyapi.sh
./stop-cliproxyapi.sh
./restart-cliproxyapi.sh
```

- vLLM lifecycle (optional fallback):
```bash
./start-vllm.sh
./check-vllm.sh
./stop-vllm.sh
./restart-vllm.sh
```

- Runtime visibility:
```bash
./status.sh
./logs.sh
```

## MCP integration (official)

This stack now includes `mcpo` (MCP-to-OpenAPI proxy) as a Docker service.

1. Ensure these env values are set in `.env` (or keep defaults from `.env.example`):
```bash
WEBUI_SECRET_KEY=<strong-random-secret-at-least-32-bytes>
MCPO_IMAGE=ghcr.io/open-webui/mcpo:main
MCPO_PORT=8000
MCPO_CONFIG=./mcp/config.json
ZOTERO_DATA_DIR=/home/miko/Zotero
```

2. Start/restart the stack:
```bash
./deploy.sh --no-logs
```

3. In OpenWebUI, add tool servers:
   - Admin Settings -> External Tools -> Add Server
   - For MCPO endpoints, use Type: `OpenAPI`
   - URL examples:
     - `http://mcpo:8000/filesystem`
     - `http://mcpo:8000/memory`
     - `http://mcpo:8000/fetch`
     - `http://mcpo:8000/time`
     - `http://mcpo:8000/zotero`
   - For native Streamable HTTP MCP servers (without MCPO), use Type: `MCP (Streamable HTTP)`

   **Automated configuration** (alternative to manual UI):
   ```bash
   # First, set your admin API key in .env:
   echo "OPENWEBUI_API_KEY=your-admin-key" >> .env

   # Verify MCPO endpoints are accessible:
   ./configure-mcpo-openapi-servers.sh --verify

   # Apply configuration (dry-run first):
   ./configure-mcpo-openapi-servers.sh --dry-run
   ./configure-mcpo-openapi-servers.sh
   ```

4. Verify MCPO endpoint:
```bash
curl -sS http://127.0.0.1:8000/docs >/dev/null && echo "mcpo up"
```

Important:
- Keep `WEBUI_SECRET_KEY` stable. Changing it invalidates encrypted stored credentials.
- Use `OpenAPI` type for MCPO proxied endpoints.
- Use `MCP (Streamable HTTP)` type for native MCP servers.
- For stdio/SSE MCP servers, configure them in `mcp/config.json` and expose through MCPO.

Preconfigured thesis MCP servers in `mcp/config.json`:
- `filesystem` (read/write under `/thesis-exports`, read-only access to `/data`)
- `memory` (session memory/knowledge graph helper)
- `fetch` (web page fetching for paper/URL extraction)
- `zotero` (local Zotero library via `/zotero/zotero.sqlite`, read-only)
- `time` (time/date utilities)

Optional external-research template:
- `mcp/config.research.optional.json` (example streamable-http endpoints for Zotero plugin MCP at `23120` and remote services)

- PDF import with optional bounded parallelism (default is serial):
```bash
./import-pdfs-to-kb.sh --parallel 3
# or via env:
IMPORT_PDFS_PARALLELISM=3 ./import-pdfs-to-kb.sh
```

## Integration and regression tests

- API coverage:
```bash
./test-api.sh
./test-api.sh --baseline
```

- RAG flow:
```bash
./test-rag.sh
./test-rag.sh --baseline
```

- OAuth alias regression:
```bash
./test-cliproxyapi-oauth.sh
```

- OpenWebUI routing regression:
```bash
./test-openwebui-cliproxy-routing.sh
```

## Important behavior in `deploy.sh`

- If CLIProxyAPI is healthy and `AUTO_SKIP_VLLM_IF_CLIPROXYAPI_HEALTHY=true`, vLLM is skipped.
- If CLIProxyAPI is unhealthy and `AUTO_START_VLLM_ON_CLIPROXYAPI_FAILURE=true`, fallback vLLM startup is allowed.
- Docker-managed CLIProxyAPI is controlled by:
1. `CLIPROXYAPI_DOCKER_MANAGED`
2. `CLIPROXYAPI_DOCKER_SERVICE`

## Utility scripts

### Provider configuration

- `configure-ollama-cloud.sh`: Configure OpenWebUI to use Ollama Cloud (https://ollama.com). Reads API key from `~/.007` (OLLAMA_CLOUD_API_KEY or OLLAMA_API_KEY), signs into OpenWebUI, and updates the Ollama connection via the admin API.

- `configure-cliproxyapi-providers.sh`: Generate a local `cliproxyapi/config.local.yaml` with real provider keys (Z.ai, MiniMax) sourced from `~/.007`. Keeps secrets out of the tracked config.

### Key rotation

- `rotate-cliproxyapi-local-key.sh`: Rotate the local API key used between OpenWebUI and CLIProxyAPI. Generates a new random key, updates `.env` and CLIProxyAPI config, syncs OpenWebUI's persistent DB config, and restarts containers.

### Tool management

- `audit-openwebui-plugins.sh`: Audit tool/function Python code stored in `webui.db`. Verifies syntax compilation and checks for missing import dependencies. Use `--focus tool|function|all` to scope the audit.

- `repair-openwebui-tools.sh`: Repair common tool issues (missing Valves exposure, drifted specs). Patches tool content in-place, recomputes specs, and restarts OpenWebUI to clear caches.

- `apply-openwebui-tool-patches.sh`: Apply patches to OpenWebUI tools via the admin API. Ensures only intended public functions are exposed (hides private `_` helpers from specs).

- `sync-openwebui-openai-config.sh`: Sync OpenWebUI's persistent OpenAI config (in `webui.db`) with local `.env` values. Use after rotating `OPENAI_API_KEY` to keep the DB in sync.

### Testing and verification

- `verify-scripts.sh`: Run script quality verification for all shell scripts. Validates bash syntax, runs shellcheck (if installed), and executes targeted smoke tests for `update.sh`.

- `test-update-smoke.sh`: Lightweight smoke tests for `update.sh` argument and restore flows. Uses command stubs so no real Docker/Compose calls are made.

- `test-openwebui-tools-endpoints.sh`: Tool-by-tool functional audit of OpenWebUI endpoints (tools, functions, code execution, web search). Validates valves/specs schemas and endpoint responses.

- `test-openwebui-tools-invocation.sh`: End-to-end integration test proving tools actually execute through `POST /api/chat/completions`. Uses local HTTP servers with nonces to verify real tool execution.

## Key files

- `docker-compose.yml`: OpenWebUI + CLIProxyAPI services (SearXNG is host-local)
- `.env.example`: canonical configuration template
- `OPENWEBUI_ADVANCED_FEATURES_PLAYBOOK.md`: phased advanced-feature rollout and hardening guide
- `cliproxyapi/config.yaml`: CLIProxyAPI provider and alias mapping
- `deploy.sh`: orchestrated startup logic with health gates
- `test-openwebui-cliproxy-routing.sh`: OpenWebUI-to-CLIProxyAPI E2E validation

## Known constraints

- `OPENWEBUI_URL` and `VLLM_URL` in `.env` are used by wrapper scripts, not Docker service discovery.
- `cli-proxy-api.sh` command namespace still uses `vllm` as the OpenAI-compatible service label.
  In this repository, that label usually points to CLIProxyAPI through `OPENAI_API_BASE_URL`.

## Troubleshooting

- If `cliproxyapi` is unhealthy in compose:
```bash
docker compose logs --tail=200 cliproxyapi
```

- If OpenWebUI does not list aliases:
```bash
curl -sS -H "Authorization: Bearer ${CLIPROXYAPI_API_KEY}" http://127.0.0.1:8317/v1/models | jq
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" http://localhost:${WEBUI_PORT:-3000}/api/models | jq
```

- If port 8317 is busy from an old host process:
```bash
lsof -nP -iTCP:8317 -sTCP:LISTEN
```

## Scope note

This README is the source of truth for the current deployment mode as of 2026-02-12.
Legacy vLLM-first documentation has been intentionally replaced with CLIProxyAPI-first operations.
