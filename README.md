# OpenWebUI + LiteLLM RAG Deployment

Last updated: 2026-03-13 (UTC)

This repository is an orchestration layer for an academic-papers RAG workflow in OpenWebUI with LiteLLM as the reference OpenAI-compatible backend.

Use this file when you need:
- the shortest path to a healthy local deploy
- the canonical daily command surface
- routing to the right source-of-truth document before editing or debugging

Do not use this file as the full procedure reference:
- use `docs/runbooks/*.md` for step-by-step operator flows
- use `docs/ssot/stack.md` for runtime truth
- use `config/README.md` before changing config files

Layout note:
- Root `docker-compose.yml`, `docker-compose.rg.yml`, `.env.example`, `mcp/`, `jupyter/`, and `searxng/` are compatibility copies of the canonical files under `config/`.
- Canonical config paths live under `config/`, canonical runbooks/references live under `docs/`, and canonical operator scripts live at the repo root.
- Local directory maps now exist for the main subtrees: `config/README.md`, `scripts/README.md`, `artifacts/README.md`, `cliproxyapi/README.md`, `research/README.md`, and `local/README.md`.
- Local-only helper assets now live under `local/` instead of cluttering the repo root. This includes sidecar binaries under `local/bin/`, optional tool payloads under `local/plugins/`, and the separate entity-maintenance workspace under `local/entities/`.
- `certs/` is reserved for local TLS material and placeholders such as `.gitkeep`; certificate files remain ignored by default.

Editing rule:
- Change `config/*` first for runtime behavior.
- Use the root copies as operator entrypoints and compatibility surfaces.
- Treat `docs/reference/openwebui/*` as upstream snapshot material, not live repo truth.

## Fast navigation

Use these first before exploring the rest of the tree:

- `docs/ssot/stack.md`: canonical runtime topology, run modes, service registry, and drift protocol.
- `docs/REPO_MAP.md`: repo-wide layout, operator entrypoints, config roots, and port inventory.
- `docs/README.md`: doc routing so runbooks, guides, plans, and reference snapshots do not compete.
- `config/README.md`: canonical config tree and root compatibility copies.
- `scripts/README.md`: secondary utilities and maintenance helpers.

## Canonical edit targets

When you need to change behavior, use these paths first:

- Compose topology: `config/compose/docker-compose.yml`
- Env defaults: `config/env/.env.example`
- MCP registry: `config/mcp/config.json`
- Jupyter config: `config/jupyter/jupyter_server_config.py`
- SearXNG config: `config/searxng/settings.yml`
- Embedding lanes/bindings: `config/embeddings/`

## Document ownership

- `README.md` is the operator front door: quick start, high-signal commands, and the main workflow.
- `docs/ssot/stack.md` owns runtime truth: run modes, critical services, critical paths, and drift gaps.
- `docs/REPO_MAP.md` owns repository shape: where config, scripts, subprojects, and entrypoints live.
- `config/README.md` owns the canonical edit surface under `config/` and the root compatibility-copy map.

## Runtime summary

- Default supported mode is Docker Compose with OpenWebUI in containers and LiteLLM as the primary host-assisted upstream at `http://host.docker.internal:4000/v1`.
- `searxng`, `cliproxyapi`, `open-terminal`, and `indigo-service` remain optional sidecars gated by profiles or env flags.
- The committed `.env.example` now pins the shipped `MCPO_IMAGE` and `INDIGO_SERVICE_IMAGE` by digest so fresh deploys do not drift when upstream tags move.
- For the full runtime model, service registry, and drift notes, use `docs/ssot/stack.md`.

```text
OpenWebUI (container) --> LiteLLM (host/container) --> providers
                     \
                      \--> optional sidecars:
                           - CLIProxyAPI (legacy OAuth workflows)
                           - Open Terminal (terminal/file integration)
                           - Indigo Service (chemistry REST APIs)
```

## Operational constraints

- No mock, dummy, or monkey fallback path for integration checks.
- OAuth aliases validated in this repo: `openai-codex`, `qwen-cli`, `kimi-cli`.
- Manual OAuth login is expected and supported.

## Prerequisites

- Docker + Docker Compose plugin
- `curl`, `jq`
- LiteLLM credentials and reachable LiteLLM endpoint

Not part of the baseline:
- CLIProxyAPI OAuth credentials
- Open Terminal API key
- Indigo Service enablement
- archived `vLLM` fallback tooling

## Quick start

Entrypoint note:
- Root `*.sh` scripts are the canonical operator surface.

1. Copy env file:
```bash
cp .env.example .env
```

2. Ensure these values are set in `.env` before the first deploy:
```bash
WEBUI_SECRET_KEY=<stable-random-secret-at-least-32-chars>
JUPYTER_TOKEN=<random-jupyter-token>
CODE_EXECUTION_JUPYTER_AUTH_TOKEN=<same-as-JUPYTER_TOKEN>
CODE_INTERPRETER_JUPYTER_AUTH_TOKEN=<same-as-JUPYTER_TOKEN>
POSTGRES_PASSWORD=<strong-local-password>
OPENAI_API_BASE_URL=http://host.docker.internal:4000/v1
OPENAI_API_BASE_URLS=http://host.docker.internal:4000/v1
OPENAI_API_KEY=<litellm-master-key>
VECTOR_DB=pgvector
CLIPROXYAPI_ENABLED=false
OPEN_TERMINAL_ENABLED=false
INDIGO_SERVICE_ENABLED=false
```

3. Sanity-check the config:
```bash
docker compose config >/dev/null
```

4. Start stack:
```bash
./deploy.sh --no-logs
```

5. Check health:
```bash
./status.sh
docker compose ps
./test-rag.sh --baseline
./test-api.sh --baseline
```

Optional legacy sidecar check (only if `CLIPROXYAPI_ENABLED=true`):
```bash
./check-cliproxyapi.sh
```

## Daily operations

Use the root scripts directly.

```bash
./deploy.sh --no-logs
./status.sh
./logs.sh
./audit-no-mock.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

High-signal local checks:

```bash
curl -fsS "http://127.0.0.1:${MCPO_PORT:-8000}/docs" >/dev/null && echo "mcpo up"
OPENWEBUI_TEST_MODEL="openrouter/openai/gpt-5-mini" ./test-api.sh --baseline
```

## Local workflow guides

Use these local documents instead of treating this file as a dump of every procedure:

- [docs/runbooks/PREREQUISITES.md](docs/runbooks/PREREQUISITES.md): first-time setup and local prerequisites.
- [docs/runbooks/OPERATIONS.md](docs/runbooks/OPERATIONS.md): day-to-day operations, updates, LLM Council, MCP setup, and maintenance flows.
- [docs/runbooks/TROUBLESHOOTING.md](docs/runbooks/TROUBLESHOOTING.md): account, model, RAG, proxy, and web-search triage.
- [docs/runbooks/EMBEDDING_PROFILES.md](docs/runbooks/EMBEDDING_PROFILES.md): embedding lanes, model switching, and KB bindings.
- [scripts/README.md](scripts/README.md): utility/audit helpers such as `check-image-versions.sh` and `audit-dependencies.sh`.

Advanced optional flows:

- CLIProxyAPI legacy OAuth sidecar: `setup-cliproxyapi.sh`, `start-cliproxyapi.sh`, `configure-cliproxyapi-oauth.sh`, `test-cliproxyapi-oauth.sh`, `import-qwen-auth.sh`, `cli-proxy-api.sh`, `local/bin/`, `docs/runbooks/OPERATIONS.md`
- Open Terminal smoke/integration flow: `test-openwebui-open-terminal.sh`
- Indigo sidecar and tool registration: `start-indigo-service.sh`, `check-indigo-service.sh`, `enable-indigo-live.sh`, local tool source under `local/plugins/`
- OpenWebUI tool repair/admin helpers: `audit-openwebui-plugins.sh`, `test-openwebui-tools-endpoints.sh`, `repair-openwebui-tools.sh`, `apply-openwebui-tool-patches.sh`
- Entity-maintenance workspace, when present locally: `local/entities/README.md`
- archived `vLLM` fallback scripts: `archive/`

## Focused local tasks

Embedding profile flow:

```bash
./manage-openwebui-embedding-profiles.sh list
./manage-openwebui-embedding-profiles.sh lanes
./manage-openwebui-embedding-profiles.sh use-kb --lane sceletium --prewarm
./manage-openwebui-embedding-profiles.sh diagnose
```

Validation loop:

```bash
./test-rag.sh --baseline
./test-api.sh --baseline
```

Legacy sidecar OAuth flow, only if you intentionally enable `cliproxyapi`:

```bash
./configure-cliproxyapi-oauth.sh
./test-cliproxyapi-oauth.sh
```

LLM Council quick run:

```bash
./llm-council.sh --prompt "Compare RAG vs fine-tuning for this local stack."
```

## Maintenance

Use the local maintenance surface, not this README, for detailed policy and edge cases:

```bash
./scripts/check-image-versions.sh
./scripts/audit-dependencies.sh
./backup-openwebui-db.sh
./update.sh
```

When updating runtime behavior:

- edit `config/` first, not the root compatibility copies
- run `./scripts/sync-compatibility-copies.sh` if the root compatibility copies need to be refreshed
- re-run `./scripts/check-doc-consistency.sh`
- re-run `./status.sh`, `./test-rag.sh --baseline`, and `./test-api.sh --baseline`

## Scope note

This README is intentionally the local front page, not the full operator manual.
Runtime truth lives in `docs/ssot/stack.md`; repo shape lives in `docs/REPO_MAP.md`; procedure detail lives in the linked runbooks and subproject docs.
