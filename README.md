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
2. `antigravity-oauth`
3. `qwen-cli`
4. `kimi-cli`
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

## OpenWebUI baseline enhancements (as of 2026-02-12)

Use this baseline workflow before enabling advanced features:

1. Deploy baseline stack:
```bash
./deploy.sh --no-logs
```

2. Validate runtime baseline:
```bash
./status.sh
```

3. Run fast baseline RAG checks:
```bash
./test-rag.sh --baseline
```

4. Run fast baseline API checks:
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
