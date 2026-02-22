# CLAUDE.md

Last updated: 2026-02-12 (UTC)

## Project mode

Current default mode is **OpenWebUI + Docker-managed CLIProxyAPI**.

- OpenWebUI container: `openwebui`
- CLIProxyAPI container: `cliproxyapi`
- OpenWebUI upstream: `http://cliproxyapi:8317/v1`
- vLLM is optional fallback and should stay stopped when CLIProxyAPI is healthy.

## Primary commands

```bash
./deploy.sh --no-logs
./status.sh
./check-cliproxyapi.sh
./test-openwebui-cliproxy-routing.sh
./test-cliproxyapi-oauth.sh
```

## OAuth aliases expected

- `openai-codex`
- `qwen-cli`
- `kimi-cli`

Manual OAuth setup is acceptable and expected.

## Operational rules

1. Prefer CLIProxyAPI-first path.
2. Do not start vLLM unless explicitly required as fallback.
3. Do not use host-gateway routing for OpenWebUI-to-CLIProxyAPI; use Docker DNS (`cliproxyapi`).
4. Integration checks must use real `/v1/models` and real chat completion calls.
