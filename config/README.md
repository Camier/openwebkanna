# Configuration Map

Canonical runtime configuration lives under `config/`.

This document owns the config edit surface only. Runtime topology and service criticality live in `docs/ssot/stack.md`.

Compatibility entrypoints at the repo root remain intentionally available as regular files/directories:

- `.env.example` -> `config/env/.env.example`
- `docker-compose.yml` -> `config/compose/docker-compose.yml`
- `docker-compose.rg.yml` -> `config/compose/docker-compose.rg.yml`
- `Dockerfile.openwebui-rg` -> `config/docker/Dockerfile.openwebui-rg`
- `jupyter/` -> `config/jupyter/`
- `mcp/` -> `config/mcp/`
- `searxng/` -> `config/searxng/`

If you land in one of the root compatibility directories while browsing locally:

- `jupyter/` means `config/jupyter/`
- `mcp/` means `config/mcp/`
- `searxng/` means `config/searxng/`

Use this tree as the edit target when you need to change runtime behavior:

- `compose/`: Docker Compose topology and service wiring.
- `docker/`: custom image builds for OpenWebUI variants.
- `embeddings/`: OpenWebUI embedding profiles and KB lane bindings.
- `env/`: committed environment template and default values.
- `jupyter/`: Jupyter server configuration for code execution.
- `mcp/`: MCP-to-OpenAPI bridge server registry.
- `searxng/`: local SearXNG configuration.
- `compatibility-copies.txt`: authoritative map of root compatibility copies to canonical `config/` paths.

Editing rule:

- Change files in `config/` first.
- Use `./scripts/sync-compatibility-copies.sh` when you need to refresh the root compatibility copies from the canonical `config/` files.
- Treat root compatibility copies as operator convenience, not the canonical edit surface.
- Re-run the relevant baseline checks after config changes: `./status.sh`, `./test-rag.sh --baseline`, `./test-api.sh --baseline`.
