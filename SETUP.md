# Setup Guide: OpenWebUI + CLIProxyAPI

Last updated: 2026-02-12 (UTC)

This is the canonical setup path for this repository.
Default mode is CLIProxyAPI-first with Docker-managed services.

## 1. Prepare environment

```bash
cd /LAB/@thesis/openwebui
cp .env.example .env
```

Edit `.env` if needed, then ensure these values:

```bash
OPENAI_API_BASE_URL=http://cliproxyapi:8317/v1
OPENAI_API_BASE_URLS=http://cliproxyapi:8317/v1
CLIPROXYAPI_ENABLED=true
CLIPROXYAPI_DOCKER_MANAGED=true
```

## 2. Configure OAuth credentials

If OAuth credentials are not set yet:

```bash
./configure-cliproxyapi-oauth.sh
```

Manual OAuth in the browser is expected.
Credentials are persisted under `cliproxyapi/auth/`.

## 3. Start services

```bash
./deploy.sh --no-logs
```

What this does:
- starts `cliproxyapi` Docker service
- validates CLIProxyAPI health and models endpoint
- starts OpenWebUI stack
- skips vLLM startup if CLIProxyAPI is healthy

## 4. Validate deployment

```bash
./status.sh
docker compose ps
./check-cliproxyapi.sh
./test-openwebui-cliproxy-routing.sh
./audit-no-mock.sh
```

Baseline validation (fast checks):

```bash
./test-rag.sh --baseline
./test-api.sh --baseline
```

If you enable web search, `test-rag.sh --baseline` probes SearXNG from both the host and the OpenWebUI container.

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
If the host probe works but the container probe fails, the host firewall is likely blocking docker-to-host connections on port 8888.

## 5. Validate OAuth aliases

```bash
./test-cliproxyapi-oauth.sh
```

Expected aliases:
1. `openai-codex`
2. `antigravity-oauth`
3. `qwen-cli`
4. `kimi-cli`

## 6. Optional: enable vLLM fallback

Only do this if you explicitly want fallback behavior.

```bash
./start-vllm.sh
./check-vllm.sh
```

## 7. Useful operations

- Restart CLIProxyAPI:
```bash
./restart-cliproxyapi.sh
```

- Stop CLIProxyAPI:
```bash
./stop-cliproxyapi.sh
```

- See logs:
```bash
./logs.sh
docker compose logs -f cliproxyapi openwebui
```

## 8. Failure recovery

If routing or alias checks fail:

1. Inspect CLIProxyAPI logs:
```bash
docker compose logs --tail=200 cliproxyapi
```

2. Re-run OAuth configuration and retry test:
```bash
./configure-cliproxyapi-oauth.sh
./test-cliproxyapi-oauth.sh
```

3. Re-run OpenWebUI routing test:
```bash
./test-openwebui-cliproxy-routing.sh
```
