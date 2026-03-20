# Configuration Map

Last updated: 2026-03-20 (UTC)

Canonical runtime configuration lives under `config/`.

This document owns the config edit surface only. Runtime topology and service criticality live in `docs/ssot/stack.md`.

Use this file when you need to know:

- which config path is canonical
- which checks to rerun after changing runtime behavior

Do not use this file as the runtime topology reference or as a troubleshooting manual.
Use `docs/ssot/stack.md` for supported topology and `docs/runbooks/*.md` for procedures.

Common config tasks:
- change Docker topology, ports, profiles, or healthchecks: `compose/docker-compose.yml`
- change committed env defaults: `env/.env.example`
- change Jupyter behavior: `jupyter/jupyter_server_config.py`
- change MCPO server exposure: `mcp/config.json`
- change SearXNG behavior: `searxng/settings.yml`

Current baseline vs target:
- Current compose-managed baseline RAG path is OpenWebUI retrieval over `VECTOR_DB=pgvector`.
- Canonical future RAG path is `multimodal_retrieval_api` on `POST /api/v1/retrieve`.
- Until migration completes, keep both statements true in docs and config comments. Do not relabel the current baseline as already migrated.

Use this tree as the edit target when you need to change runtime behavior:

- `compose/`: Docker Compose topology and service wiring.
- `docker/`: custom image builds for OpenWebUI variants.
- `env/`: committed environment template and default values.
- `jupyter/`: Jupyter server configuration for code execution.
- `mcp/`: MCP-to-OpenAPI bridge server registry.
- `searxng/`: local SearXNG configuration.
Editing rule:

- Change files in `config/` first.

After changing config, prefer this validation order:

```bash
docker compose --project-directory . -f config/compose/docker-compose.yml config -q
./scripts/check-doc-consistency.sh
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

Minimal safe workflow:

1. Edit the canonical file under `config/`.
2. Run `docker compose --project-directory . -f config/compose/docker-compose.yml config -q` when compose or env behavior changed.
3. Run `./scripts/check-doc-consistency.sh` if docs or defaults changed.
4. Run `./status.sh`, `./test-rag.sh --baseline`, and `./test-api.sh --baseline` for runtime-impacting changes.
5. If you changed canonical future RAG defaults or docs, also verify the host-native retrieval path you documented.
