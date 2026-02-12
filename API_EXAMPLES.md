# OpenWebUI + CLIProxyAPI API Examples

Last updated: 2026-02-12 (UTC)

This document provides reproducible HTTP examples for the current stack:
- OpenWebUI in Docker
- CLIProxyAPI in Docker (`http://127.0.0.1:8317` externally)
- OpenWebUI upstream mapped to `http://cliproxyapi:8317/v1` internally

## Runtime environment

```bash
set -a
source ./.env
set +a

export OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
export CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
export CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-layra-cliproxyapi-key}"
export OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-$API_KEY}"
```

## Health checks

CLIProxyAPI root:
```bash
curl -sS -f "$CLIPROXYAPI_BASE_URL/" | jq .
```

CLIProxyAPI models:
```bash
curl -sS -f -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
  "$CLIPROXYAPI_BASE_URL/v1/models" | jq -e '.data | type == "array" and length > 0'
```

OpenWebUI health:
```bash
curl -sS -f "$OPENWEBUI_URL/health" | jq .
```

OpenWebUI models:
```bash
curl -sS -f -H "Authorization: Bearer $OPENWEBUI_API_KEY" \
  "$OPENWEBUI_URL/api/models" | jq .
```

## Required alias presence

```bash
curl -sS -f -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
  "$CLIPROXYAPI_BASE_URL/v1/models" \
| jq -e '
  [.data[].id] as $ids
  | ["openai-codex","antigravity-oauth","qwen-cli","kimi-cli"]
  | all(. as $name | $ids | index($name))
'
```

## Chat completion examples

CLIProxyAPI direct call (`openai-codex`):
```bash
curl -sS -f -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
  -H "Content-Type: application/json" \
  "$CLIPROXYAPI_BASE_URL/v1/chat/completions" \
  -d '{
    "model": "openai-codex",
    "messages": [{"role": "user", "content": "Say ping."}],
    "max_tokens": 24
  }' | jq .
```

OpenWebUI proxied chat (`qwen-cli`):
```bash
curl -sS -f -H "Authorization: Bearer $OPENWEBUI_API_KEY" \
  -H "Content-Type: application/json" \
  "$OPENWEBUI_URL/api/chat/completions" \
  -d '{
    "model": "qwen-cli",
    "messages": [{"role": "user", "content": "Summarize retrieval-augmented generation in one sentence."}],
    "max_tokens": 64
  }' | jq .
```

## Wrapper examples (`cli-proxy-api.sh`)

The wrapper command namespace keeps `vllm` as the OpenAI-compatible endpoint label.
In this repo, `vllm` usually resolves to CLIProxyAPI through `.env` (`OPENAI_API_BASE_URL`).

List models via wrapper against configured OpenAI-compatible endpoint:
```bash
./cli-proxy-api.sh --raw models vllm | jq .
```

Force wrapper to CLIProxyAPI explicitly:
```bash
./cli-proxy-api.sh --raw --url "$CLIPROXYAPI_BASE_URL" models vllm | jq .
```

Health all services:
```bash
./cli-proxy-api.sh health all
```

## Regression test commands

OAuth alias regression:
```bash
./test-cliproxyapi-oauth.sh
```

OpenWebUI routing regression:
```bash
./test-openwebui-cliproxy-routing.sh
```

API smoke suite:
```bash
./test-api.sh
```

RAG smoke suite:
```bash
./test-rag.sh
```

## Notes

- `CLIPROXYAPI_ENABLED=false` is intentionally treated as a failed integration state.
- Manual OAuth is supported; auth material is persisted under `cliproxyapi/auth/`.
- If CLIProxyAPI is healthy, deployment keeps vLLM stopped by default.
