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
- Canonical config paths live under `config/`, canonical runbooks/references live under `docs/`, the small daily operator surface lives at the repo root, and secondary maintenance or optional workflows live under `scripts/`.
- Runtime service code that should not be reduced to shell wrappers lives under `services/`.
- Local directory maps now exist for the main subtrees: `config/README.md`, `scripts/README.md`, `artifacts/README.md`, `research/README.md`, and `local/README.md`.
- Local-only helper assets now live under `local/` instead of cluttering the repo root. This includes sidecar binaries under `local/bin/`, optional tool payloads under `local/plugins/`, and the separate entity-maintenance workspace under `local/entities/`.
- `certs/` is reserved for local TLS material and placeholders such as `.gitkeep`; certificate files remain ignored by default.

Editing rule:
- Change `config/*` first for runtime behavior.
- Use `scripts/*` for narrower admin, testing, ingest, and optional-sidecar workflows.
- Treat `docs/reference/openwebui/*` as upstream snapshot material, not live repo truth.

## Fast navigation

Use these first before exploring the rest of the tree:

- `docs/ssot/stack.md`: canonical runtime topology, run modes, service registry, and drift protocol.
- `docs/REPO_MAP.md`: repo-wide layout, operator entrypoints, config roots, and port inventory.
- `docs/README.md`: doc routing so runbooks, guides, plans, and reference snapshots do not compete.
- `config/README.md`: canonical config tree.
- `scripts/README.md`: secondary utilities and maintenance helpers.

## Canonical edit targets

When you need to change behavior, use these paths first:

- Compose topology: `config/compose/docker-compose.yml`
- Env defaults: `config/env/.env.example`
- MCP registry: `config/mcp/config.json`
- Jupyter config: `config/jupyter/jupyter_server_config.py`
- SearXNG config: `config/searxng/settings.yml`

## Document ownership

- `README.md` is the operator front door: quick start, high-signal commands, and the main workflow.
- `docs/ssot/stack.md` owns runtime truth: run modes, critical services, critical paths, and drift gaps.
- `docs/REPO_MAP.md` owns repository shape: where config, scripts, subprojects, and entrypoints live.
- `config/README.md` owns the canonical edit surface under `config/`.

## Runtime summary

- Default supported mode is Docker Compose with OpenWebUI in containers and LiteLLM as the primary host-assisted upstream at `http://host.docker.internal:4000/v1`.
- `searxng`, `open-terminal`, and `indigo-service` remain optional sidecars gated by profiles or env flags.
- The committed `config/env/.env.example` now pins the shipped `MCPO_IMAGE` and `INDIGO_SERVICE_IMAGE` by digest so fresh deploys do not drift when upstream tags move.
- For the full runtime model, service registry, and drift notes, use `docs/ssot/stack.md`.

```text
OpenWebUI (container) --> LiteLLM (host/container) --> providers
                     \
                      \--> optional sidecars:
                           - Open Terminal (terminal/file integration)
                           - Indigo Service (chemistry REST APIs)
```

## Operational constraints

- No mock, dummy, or monkey fallback path for integration checks.

## Prerequisites

- Docker + Docker Compose plugin
- `curl`, `jq`
- LiteLLM credentials and reachable LiteLLM endpoint

Not part of the baseline:
- Open Terminal API key
- Indigo Service enablement

## Quick start

Entrypoint note:
- Root `*.sh` scripts are the primary daily operator surface.
- Secondary workflows live under `scripts/`.

1. Copy env file:
```bash
cp config/env/.env.example .env
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

## Daily operations

Use the root scripts directly.

```bash
./deploy.sh --no-logs
./status.sh
./logs.sh
./scripts/testing/audit-no-mock.sh
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
- [docs/runbooks/OPERATIONS.md](docs/runbooks/OPERATIONS.md): day-to-day operations, retrieval, MCP setup, and maintenance flows.
- [docs/runbooks/TROUBLESHOOTING.md](docs/runbooks/TROUBLESHOOTING.md): account, model, RAG, proxy, and web-search triage.
- [scripts/README.md](scripts/README.md): utility/audit helpers such as `check-image-versions.sh`.

Advanced optional flows:

- Open Terminal smoke/integration flow: `scripts/open-terminal/test-openwebui-open-terminal.sh`
- Indigo sidecar and tool registration: `scripts/indigo/start-indigo-service.sh`, `scripts/indigo/check-indigo-service.sh`, `scripts/indigo/enable-indigo-live.sh`, local tool source under `local/plugins/`
- OpenWebUI tool repair/admin helpers: `scripts/testing/audit-openwebui-plugins.sh`, `scripts/testing/test-openwebui-tools-endpoints.sh`, `scripts/admin/repair-openwebui-tools.sh`, `scripts/admin/apply-openwebui-tool-patches.sh`
- Entity-maintenance workspace, when present locally: `local/entities/README.md`

## Focused local tasks

Validation loop:

```bash
./test-rag.sh --baseline
./test-api.sh --baseline
```

Multimodal retrieval artifact from an existing OpenWebUI knowledge base:

```bash
./scripts/rag/render-multimodal-answer.sh \
  --kb-name "Sceletium Research" \
  --query "What is the structure of mesembrine?"
```

Direct multimodal retrieval API against the native `rag_evidence` Qdrant backend:

```bash
PYTHONPATH=/LAB/@thesis/openwebui \
MULTIMODAL_RETRIEVAL_API_QDRANT_URL=http://localhost:6333 \
MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH=/path/to/nemotron-model \
MULTIMODAL_RETRIEVAL_API_TEXT_SPARSE_MODEL_NAME=Qdrant/bm25 \
/LAB/@thesis/wow/.venv/bin/python -m uvicorn \
  services.multimodal_retrieval_api.app:app \
  --host 127.0.0.1 \
  --port 8510

curl -fsS http://127.0.0.1:8510/health
curl -fsS http://127.0.0.1:8510/ready
curl -fsS -X POST http://127.0.0.1:8510/api/v1/retrieve \
  -H 'content-type: application/json' \
  -d '{"query":"mesembrine structure","top_k":3}'
```

One-shot migration/backfill into `rag_evidence`:

```bash
PYTHONPATH=/LAB/@thesis/openwebui \
python scripts/rag/migrate_rag_evidence.py \
  --manifest artifacts/rag/manifests/<run>.json \
  --output-report artifacts/rag/migrate_rag_evidence.report.json
```

Optional compatibility env file:

```bash
export MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE=/LAB/@thesis/wow/.env
```

Build and use the CLIP-backed multimodal figure index (V2):

```bash
python3 scripts/rag/build_multimodal_index.py
./scripts/rag/render-multimodal-answer-v2.sh \
  --kb-name "Sceletium Research" \
  --query "What is the structure of mesembrine?"
```

Direct retrieval config inspection:

```bash
./scripts/admin/sync-openwebui-retrieval-config.sh

OPENWEBUI_API_KEY="<admin-bearer-token>"
curl -s -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "http://127.0.0.1:${WEBUI_PORT:-3000}/api/v1/retrieval/config" | jq
```

## Maintenance

Use the local maintenance surface, not this README, for detailed policy and edge cases:

```bash
./scripts/check-image-versions.sh
./update.sh
```

When updating runtime behavior:

- edit `config/` first for runtime behavior
- re-run `./scripts/check-doc-consistency.sh`
- re-run `./status.sh`, `./test-rag.sh --baseline`, and `./test-api.sh --baseline`

## Scope note

This README is intentionally the local front page, not the full operator manual.
Runtime truth lives in `docs/ssot/stack.md`; repo shape lives in `docs/REPO_MAP.md`; procedure detail lives in the linked runbooks and subproject docs.
