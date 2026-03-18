# Canonical-Only Layout Migration Note

Last updated: 2026-03-18 (UTC)

This repository no longer uses root compatibility copies for runtime configuration.

## Canonical Paths

- Compose topology: `config/compose/docker-compose.yml`
- Env template: `config/env/.env.example`
- Jupyter config: `config/jupyter/jupyter_server_config.py`
- MCPO registry: `config/mcp/config.json`
- SearXNG config: `config/searxng/settings.yml`
- RG build override: `config/compose/docker-compose.rg.yml`
- RG Dockerfile: `config/docker/Dockerfile.openwebui-rg`

## Operator Changes

- Create the local env file with `cp config/env/.env.example .env`
- Root shell entrypoints still exist: `./deploy.sh`, `./status.sh`, `./cleanup.sh`, `./logs.sh`, `./update.sh`
- Those scripts now target the canonical compose file under `config/compose/`
- Compose execution now uses the repo root as the project directory so `.env` still resolves correctly even though the compose file lives under `config/compose/`

## Removed Root Compatibility Surfaces

The following root copies were removed so the repo has a single committed source of truth:

- `.env.example`
- `docker-compose.yml`
- `docker-compose.rg.yml`
- `Dockerfile.openwebui-rg`
- `jupyter/`
- `mcp/`
- `searxng/`

## Related Cleanup

This migration also removed:

- legacy CLIProxy sidecar surfaces
- archived vLLM fallback surfaces
- compatibility-copy sync machinery
- fake/mock-heavy harness layers tied to the removed legacy paths

## Validation After Migration

Use these checks after layout or runtime changes:

```bash
./scripts/check-doc-consistency.sh
./scripts/testing/verify-scripts.sh
./test-rag.sh --baseline
./test-api.sh --baseline
docker compose --project-directory /LAB/@thesis/openwebui \
  --env-file /LAB/@thesis/openwebui/.env \
  -f /LAB/@thesis/openwebui/config/compose/docker-compose.yml config -q
```

## Where To Look Next

- Runtime truth: `docs/ssot/stack.md`
- Repo layout: `docs/REPO_MAP.md`
- Config ownership: `config/README.md`
- Operator procedures: `docs/runbooks/`
