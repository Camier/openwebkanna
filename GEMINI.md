# OpenWebUI + CLIProxyAPI RAG Notes

Last updated: 2026-02-12 (UTC)

## Current architecture

- OpenWebUI runs in Docker.
- CLIProxyAPI runs in Docker.
- OpenWebUI calls CLIProxyAPI via `http://cliproxyapi:8317/v1`.
- vLLM is optional fallback only.

## Critical commands

```bash
./deploy.sh --no-logs
./status.sh
./check-cliproxyapi.sh
./test-openwebui-cliproxy-routing.sh
./test-cliproxyapi-oauth.sh
```

## OAuth aliases

Expected aliases exposed by `/v1/models`:
1. `openai-codex`
2. `qwen-cli`
3. `kimi-cli`

## Guardrails

- Keep deployments CLIProxyAPI-first.
- Avoid starting vLLM when CLIProxyAPI is healthy.
- Use real integration checks (no mocks or bypass modes).
