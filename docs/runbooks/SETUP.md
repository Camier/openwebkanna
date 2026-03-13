# Setup Guide: OpenWebUI + LiteLLM

Last updated: 2026-03-13 (UTC)

This is the canonical setup path for this repository.
Default mode is LiteLLM-first with Docker-managed services.

Use this runbook when:
- you are setting up a fresh local clone
- you want to rebuild the baseline from scratch

Exit criteria for this runbook:
- `docker compose config >/dev/null` passes
- `./status.sh` shows the core services healthy
- `./test-rag.sh --baseline` passes
- `./test-api.sh --baseline` passes

## 1. Prepare environment

```bash
cd /LAB/@thesis/openwebui
cp .env.example .env
```

Edit `.env` before the first deploy, then ensure these values:

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
```

For the full committed baseline, use `config/env/.env.example`. Keep this runbook list as the minimum first-deploy contract, not the exhaustive env reference.

Optional sidecars stay off for the baseline unless you intentionally enable them later:
- `OPEN_TERMINAL_ENABLED=false`
- `INDIGO_SERVICE_ENABLED=false`

## 2. Configure LiteLLM access

- Start LiteLLM and confirm `/v1/models` is reachable from host.
- Ensure `OPENAI_API_KEY` in `.env` matches LiteLLM's configured master key.

## 3. Start services

Sanity-check the rendered compose config first:

```bash
docker compose config >/dev/null
```

Then deploy:

```bash
./deploy.sh --no-logs
```

What this does:
- starts the core OpenWebUI baseline
- starts optional sidecars only if their env flags are enabled
- validates core service health
- keeps OpenWebUI routed to LiteLLM (reference upstream)

## 4. Validate deployment

```bash
./status.sh
docker compose ps
./scripts/testing/audit-no-mock.sh
```

Optional legacy checks (only if `CLIPROXYAPI_ENABLED=true`):

```bash
./scripts/cliproxyapi/check-cliproxyapi.sh
./scripts/cliproxyapi/test-openwebui-cliproxy-routing.sh
```

Baseline validation (fast checks):

```bash
./test-rag.sh --baseline
./test-api.sh --baseline
```

If any of these fail, stop here and use `docs/runbooks/TROUBLESHOOTING.md` before enabling optional flows.

If you enable web search, `test-rag.sh --baseline` probes SearXNG from both the host and the OpenWebUI container.

If `/api/models` is healthy but one chat model is flaky, pin a known-good model for API baseline:

```bash
OPENWEBUI_TEST_MODEL="openrouter/openai/gpt-5-mini" ./test-api.sh --baseline
```

Baseline note:
- `pgvector` is the committed default retrieval backend; `chroma` remains an explicit opt-in override.
- Web search stays disabled by default even though the `searxng` route is preconfigured.
- Jupyter-backed code execution is part of the baseline and is driven by `.env`.
- `CODE_EXECUTION_JUPYTER_AUTH_TOKEN` and `CODE_INTERPRETER_JUPYTER_AUTH_TOKEN` must match `JUPYTER_TOKEN`.
- Use `./test-rag.sh --baseline` after deploy to validate the full default lane.

## 5. After setup

For routine operations after the first deploy:
- use `docs/runbooks/OPERATIONS.md` for start/stop, logs, backups, MCP, RAG tuning, and maintenance
- use `docs/runbooks/TROUBLESHOOTING.md` for auth, model, vector-store, SSL, and sidecar recovery flows
- use `docs/runbooks/EMBEDDING_PROFILES.md` for lane selection and KB/profile lifecycle

Advanced and optional flows intentionally live outside this setup runbook:
- image updates and version review: `./scripts/check-image-versions.sh`, `./update.sh`
- manual document-tuning and multimodal validation: `docs/runbooks/OPERATIONS.md`
- legacy CLIProxyAPI OAuth workflows: `docs/runbooks/OPERATIONS.md`
- archived `vLLM` fallback: `archive/`

If you see an empty model list, a pending activation screen, or repeated `401`/`502` issues:
- take a DB snapshot with `./scripts/admin/backup-openwebui-db.sh`
- then use `./scripts/admin/openwebui-user-admin.sh --email you@example.com` or the troubleshooting runbook, depending on whether the issue is account state or runtime drift
