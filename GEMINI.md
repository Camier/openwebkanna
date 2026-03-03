# OpenWebUI + LiteLLM RAG Notes

Last updated: 2026-03-03 (UTC)

## Current architecture

- OpenWebUI runs in Docker.
- LiteLLM is the reference upstream (`http://host.docker.internal:4000/v1`).
- CLIProxyAPI is optional/deprecated as a primary upstream.
- vLLM is optional fallback only.

## Critical commands

```bash
./deploy.sh --no-logs
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

## OAuth aliases

Expected aliases exposed by `/v1/models`:
1. `openai-codex`
2. `qwen-cli`
3. `kimi-cli`

## Guardrails

- Keep deployments LiteLLM-first.
- Keep CLIProxyAPI disabled unless a legacy sidecar workflow is required.
- Use real integration checks (no mocks or bypass modes).
