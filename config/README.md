# Configuration Map

Last updated: 2026-03-13 (UTC)

Canonical runtime configuration lives under `config/`.

This document owns the config edit surface only. Runtime topology and service criticality live in `docs/ssot/stack.md`.

Use this file when you need to know:

- which config path is canonical
- which root files are compatibility copies
- which checks to rerun after changing runtime behavior

Do not use this file as the runtime topology reference or as a troubleshooting manual.
Use `docs/ssot/stack.md` for supported topology and `docs/runbooks/*.md` for procedures.

Common config tasks:
- change Docker topology, ports, profiles, or healthchecks: `compose/docker-compose.yml`
- change committed env defaults: `env/.env.example`
- change Jupyter behavior: `jupyter/jupyter_server_config.py`
- change MCPO server exposure: `mcp/config.json`
- change SearXNG behavior: `searxng/settings.yml`
- change embedding lanes or KB bindings: `embeddings/`

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

After changing config, prefer this validation order:

```bash
docker compose config -q
./scripts/check-doc-consistency.sh
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

Minimal safe workflow:

1. Edit the canonical file under `config/`.
2. Run `./scripts/sync-compatibility-copies.sh` if the matching root compatibility copy exists.
3. Run `docker compose config -q` when compose or env behavior changed.
4. Run `./scripts/check-doc-consistency.sh` if docs or defaults changed.
5. Run `./status.sh`, `./test-rag.sh --baseline`, and `./test-api.sh --baseline` for runtime-impacting changes.
