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
```

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
