# CLAUDE.md

Last updated: 2026-03-03 (UTC)

## Project mode

Current default mode is **OpenWebUI + LiteLLM**.

- OpenWebUI container: `openwebui`
- LiteLLM upstream: `http://host.docker.internal:4000/v1`
- CLIProxyAPI container (`cliproxyapi`) is an optional legacy sidecar.
- vLLM remains an optional fallback path only.

## Primary commands

```bash
./deploy.sh --no-logs
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

## OAuth aliases expected

- `openai-codex`
- `qwen-cli`
- `kimi-cli`

Manual OAuth setup is acceptable and expected.

## Operational rules

1. Prefer LiteLLM-first path.
2. Keep CLIProxyAPI disabled unless legacy OAuth aliases are explicitly required.
3. Do not start vLLM unless explicitly required as fallback.
4. Integration checks must use real `/v1/models` and real chat completion calls.
