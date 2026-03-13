# Setup Guide: OpenWebUI + LiteLLM

Last updated: 2026-03-13 (UTC)

This is the canonical setup path for this repository.
Default mode is LiteLLM-first with Docker-managed services.

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
./audit-no-mock.sh
```

Optional legacy checks (only if `CLIPROXYAPI_ENABLED=true`):

```bash
./check-cliproxyapi.sh
./test-openwebui-cliproxy-routing.sh
```

Baseline validation (fast checks):

```bash
./test-rag.sh --baseline
./test-api.sh --baseline
```

If you enable web search, `test-rag.sh --baseline` probes SearXNG from both the host and the OpenWebUI container.

If `/api/models` is healthy but one chat model is flaky, pin a known-good model for API baseline:

```bash
OPENWEBUI_TEST_MODEL="openrouter/openai/gpt-5-mini" ./test-api.sh --baseline
```

If you want to align your pinned image to the latest upstream tag:

```bash
OPENWEBUI_TAG="$(curl -s https://api.github.com/repos/open-webui/open-webui/releases/latest | jq -r '.tag_name')"
sed -i "s|^OPENWEBUI_IMAGE=.*|OPENWEBUI_IMAGE=ghcr.io/open-webui/open-webui:${OPENWEBUI_TAG}|" .env
```

Baseline note:
- `pgvector` is the committed default retrieval backend; `chroma` remains an explicit opt-in override.
- Web search stays disabled by default even though the `searxng` route is preconfigured.
- Jupyter-backed code execution is part of the baseline and is driven by `.env`.
- `CODE_EXECUTION_JUPYTER_AUTH_TOKEN` and `CODE_INTERPRETER_JUPYTER_AUTH_TOKEN` must match `JUPYTER_TOKEN`.
- Use `./test-rag.sh --baseline` after deploy to validate the full default lane.

## Tune "Documents" Settings (Manual-Only)

OpenWebUI stores "Documents" (retrieval) settings persistently in its DB. For reproducible tuning and rollback, use:

```bash
# Snapshot only (dry-run)
OPENWEBUI_API_KEY="<admin-bearer-token>" ./tune-openwebui-documents.sh

# Apply the tuned values
OPENWEBUI_API_KEY="<admin-bearer-token>" ./tune-openwebui-documents.sh --apply

# Restore from a snapshot JSON written under logs/
OPENWEBUI_API_KEY="<admin-bearer-token>" ./tune-openwebui-documents.sh --restore logs/openwebui-documents-snapshot-<timestamp>.json
```

Notes:
- This is intentionally manual-token only (no auto sign-in), to avoid hidden credential assumptions.
- UI path for the same settings: `http://localhost:<WEBUI_PORT>/admin/settings/documents`

### LiteLLM multimodal full validation

When running with LiteLLM (port `4000`), validate multimodal ingestion with a local PDF fixture.

```bash
# 1) Apply multimodal-friendly retrieval settings
OPENWEBUI_API_KEY="<admin-bearer-token>" \
TUNE_CONTENT_EXTRACTION_ENGINE=docling \
TUNE_PDF_EXTRACT_IMAGES=true \
./tune-openwebui-documents.sh --apply

# 2) Run full RAG test suite in LiteLLM mode with a local PDF
OPENWEBUI_URL=http://localhost:${WEBUI_PORT:-3000} \
RAG_UPSTREAM_MODE=litellm \
VLLM_URL=http://localhost:4000 \
RAG_MULTIMODAL_TEST_PDF="data/pdfs/<your-paper>.pdf" \
./test-rag.sh --full
```

If you need strict enforcement (no skip on missing fixture or prereqs), set:

```bash
RAG_MULTIMODAL_STRICT=true
```

## 5. Advanced optional flows

These are not part of the normal baseline setup:

- Legacy CLIProxyAPI OAuth aliases:
```bash
./test-cliproxyapi-oauth.sh
```

- Archived `vLLM` fallback:
```bash
./archive/start-vllm.sh
./archive/check-vllm.sh
```

Use `docs/runbooks/OPERATIONS.md` for the full optional/legacy procedures.

## 6. After setup

For routine operations after the first deploy:
- use `docs/runbooks/OPERATIONS.md` for start/stop, logs, backups, MCP, RAG tuning, and maintenance
- use `docs/runbooks/TROUBLESHOOTING.md` for auth, model, vector-store, SSL, and sidecar recovery flows

If you see an empty model list, a pending activation screen, or repeated `401`/`502` issues:
- take a DB snapshot with `./backup-openwebui-db.sh`
- then use `./openwebui-user-admin.sh --email you@example.com` or the troubleshooting runbook, depending on whether the issue is account state or runtime drift
