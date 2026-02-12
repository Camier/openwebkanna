# OpenWebUI + CLIProxyAPI Deployment Prerequisites

Last verified: 2026-02-12 (UTC)

This checklist reflects the current default deployment mode in this repo:
- Docker-managed CLIProxyAPI
- Docker-managed OpenWebUI
- vLLM optional fallback only

## Required

1. Docker daemon available
```bash
docker info >/dev/null
```

2. Docker Compose plugin available
```bash
docker compose version
```

3. Core CLI tooling available
```bash
command -v curl
command -v jq
```

4. Repo env file prepared
```bash
cp .env.example .env
```

5. Required `.env` defaults present
```bash
rg -n "^OPENAI_API_BASE_URL=|^OPENAI_API_BASE_URLS=|^CLIPROXYAPI_DOCKER_MANAGED=|^CLIPROXYAPI_ENABLED=" .env
```

Expected values for CLIProxyAPI-first mode:
- `OPENAI_API_BASE_URL=http://cliproxyapi:8317/v1`
- `OPENAI_API_BASE_URLS=http://cliproxyapi:8317/v1`
- `CLIPROXYAPI_DOCKER_MANAGED=true`
- `CLIPROXYAPI_ENABLED=true`

## OAuth prerequisites

For aliases `openai-codex`, `antigravity-oauth`, `qwen-cli`, `kimi-cli`:

1. Complete provider OAuth manually when prompted.
2. Ensure auth artifacts exist under `cliproxyapi/auth/`.
3. Confirm aliases are visible:
```bash
curl -sS -H "Authorization: Bearer ${CLIPROXYAPI_API_KEY:-layra-cliproxyapi-key}" \
  http://127.0.0.1:8317/v1/models | jq .
```

## Optional fallback prerequisites (only if enabling vLLM)

- Python runtime with `vllm` installed
- GPU and CUDA compatibility
- Free port `8000`

Check fallback prerequisites:
```bash
python3 -c "import vllm"
lsof -nP -iTCP:8000 -sTCP:LISTEN || true
```

## Pre-deploy verification

```bash
docker compose config >/dev/null
bash -n deploy.sh start-cliproxyapi.sh check-cliproxyapi.sh test-openwebui-cliproxy-routing.sh
```

## Launch and validate

```bash
./deploy.sh --no-logs
./check-cliproxyapi.sh
./test-openwebui-cliproxy-routing.sh
```

Passing criteria:
- `cliproxyapi` healthy in `docker compose ps`
- OpenWebUI reachable on `http://localhost:${WEBUI_PORT:-3000}`
- Alias list and chat routing checks pass for all required OAuth aliases
