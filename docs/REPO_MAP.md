# REPO_MAP

Last refreshed: 2026-03-13 (UTC)
Scope: `/LAB/@thesis/openwebui`

This file owns repository shape and operator navigation. It is intentionally not a second runtime spec.

## 1. Fast routing

Use these files in this order so the repo stays legible:

- `README.md`: operator front door, quick start, and the daily command surface.
- `docs/ssot/stack.md`: current runtime topology, critical services, and drift rules.
- `config/README.md`: canonical edit surface under `config/` and the compatibility-copy map.
- `docs/README.md`: documentation routing for runbooks, guides, status notes, plans, and reference snapshots.
- `scripts/README.md`: secondary maintenance and audit tooling.

## 2. Canonical vs compatibility surfaces

This repo intentionally keeps two layers:

- Canonical edit surface:
  - `config/compose/`
  - `config/env/`
  - `config/jupyter/`
  - `config/mcp/`
  - `config/searxng/`
  - `config/embeddings/`
- Root/operator compatibility surface:
  - `.env.example`
  - `docker-compose.yml`
  - `docker-compose.rg.yml`
  - `Dockerfile.openwebui-rg`
  - `jupyter/`
  - `mcp/`
  - `searxng/`

Editing rule:

- Change `config/*` first.
- Keep the root compatibility copies byte-aligned for operator convenience and existing scripts.
- Treat `docs/reference/openwebui/*` as snapshot/reference material, not as local deployment truth.

## 3. Top-level layout

Operationally relevant roots:

- `*.sh`: canonical operator entrypoints.
- `config/`: canonical runtime configuration.
- `docs/`: runbooks, SSOT, guides, status notes, plans, reviews, and reference snapshots.
- `lib/`: shared shell helpers used by root scripts.
- `scripts/`: narrower maintenance, hygiene, backup, and data utilities.
- `cliproxyapi/`: optional legacy sidecar config, local overlay, and auth state.
- `local/`: local-only helper namespace for sidecar binaries, optional tool payloads, and separate workspaces that should not clutter the repo root.
- `certs/`: local TLS material and placeholders; real certificate/key files stay ignored by default.
- `data/`: local corpus, PDFs, extractions, notes, and processing artifacts.
- `artifacts/`: generated audit and evaluation outputs.
- `archive/`: archived fallback flows and retired operator paths.
- `research/`: research-oriented helper assets.

Local helper layout under `local/`:

- `local/bin/`: ignored local sidecar binaries and shims consumed by optional root scripts.
- `local/plugins/`: ignored local tool sources used by optional script-driven integrations, such as the Indigo tool payload.
- `local/entities/`: optional local entity-maintenance workspace; its README is tracked for orientation, while the rest stays ignored by default and should not be treated as baseline runtime topology.

## 4. Operator entrypoints

Canonical operator commands live at the repo root.

Primary daily surface:

- Deploy and status: `deploy.sh`, `status.sh`, `logs.sh`, `cleanup.sh`, `update.sh`
- Baseline validation: `audit-no-mock.sh`, `test-rag.sh --baseline`, `test-api.sh --baseline`, `verify-scripts.sh`
- Admin and repair: `openwebui-user-admin.sh`, `repair-openwebui-tools.sh`, `sync-openwebui-openai-config.sh`, `sync-openwebui-web-search-config.sh`
- Ingest and retrieval tuning: `import-pdfs-to-kb.sh`, `tune-openwebui-documents.sh`, `manage-openwebui-embedding-profiles.sh`

Optional sidecar/tool flows:

- CLIProxyAPI lifecycle/bootstrap/auth: `setup-cliproxyapi.sh`, `start-cliproxyapi.sh`, `stop-cliproxyapi.sh`, `restart-cliproxyapi.sh`, `check-cliproxyapi.sh`, `configure-cliproxyapi-oauth.sh`, `test-cliproxyapi-oauth.sh`, `import-qwen-auth.sh`, `cli-proxy-api.sh`, `local/bin/cliproxyapi`
- Indigo: `start-indigo-service.sh`, `stop-indigo-service.sh`, `restart-indigo-service.sh`, `check-indigo-service.sh`, `enable-indigo-live.sh`, `local/plugins/indigo_chemistry_tool.py`
- Open Terminal: `test-openwebui-open-terminal.sh`
- OpenWebUI tool repair/admin helpers: `audit-openwebui-plugins.sh`, `test-openwebui-tools-endpoints.sh`, `repair-openwebui-tools.sh`, `apply-openwebui-tool-patches.sh`
- Entity-maintenance workspace when present: `local/entities/README.md`, `local/entities/pipelines/`

Use `scripts/` only when you need narrower support tooling or audits.

## 5. Configuration roots

High-signal config paths:

- `config/compose/docker-compose.yml`: compose topology, healthchecks, profiles, and host bindings
- `config/env/.env.example`: committed runtime defaults
- `config/mcp/config.json`: MCPO server registry
- `config/jupyter/jupyter_server_config.py`: Jupyter runtime config
- `config/searxng/settings.yml`: local SearXNG config
- `config/embeddings/profiles.json`: embedding profile catalog
- `config/embeddings/kb-bindings.json`: KB-to-lane binding map
- `config/compatibility-copies.txt`: authoritative root-to-canonical sync map

## 6. Documentation tree

Documentation under `docs/` is split by purpose:

- `docs/ssot/`: runtime source of truth
- `docs/runbooks/`: setup, operations, troubleshooting, and embedding-profile procedures
- `docs/guides/`: focused operator guidance
- `docs/status/`: observed local runtime snapshots and point-in-time deployment notes
- `docs/reviews/`: audits and technical reviews
- `docs/plans/`: historical or pending planning documents
- `docs/reference/openwebui/`: upstream-derived OpenWebUI snapshots
- `docs/legacy/`: compatibility and migration notes

## 7. State and data

Primary persisted state lives in:

- Docker volumes:
  - `openwebui_data`
  - `openwebui_postgres_data`
  - `openwebui_jupyter_data`
- Repo-local state:
  - `data/`
  - `backups/`
  - `logs/`
  - `thesis-exports/`
  - `artifacts/`

## 8. Port inventory

From the compose definition:

- OpenWebUI: `127.0.0.1:${WEBUI_PORT:-3000} -> 8080`
- PostgreSQL: `${POSTGRES_BIND_ADDRESS:-127.0.0.1}:${POSTGRES_PORT:-5432} -> 5432`
- Jupyter: `127.0.0.1:${JUPYTER_PORT:-8890} -> ${JUPYTER_INTERNAL_PORT:-8889}`
- SearXNG: `${SEARXNG_BIND_ADDRESS:-127.0.0.1}:${SEARXNG_PORT:-8888} -> 8080`
- MCPO: `${MCPO_BIND_ADDRESS:-127.0.0.1}:${MCPO_PORT:-8000} -> 8000`
- Docling: `127.0.0.1:5001 -> 5001`
- CLIProxyAPI: `${CLIPROXYAPI_BIND_ADDRESS:-127.0.0.1}:${CLIPROXYAPI_PORT:-8317} -> 8317`
- Open Terminal: `${OPEN_TERMINAL_BIND_ADDRESS:-127.0.0.1}:${OPEN_TERMINAL_PORT:-8320} -> 8000`
- Indigo Service: `${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}:${INDIGO_SERVICE_PORT:-8012} -> 80`

## 9. Orientation path

For a normal repo pass:

1. Read `README.md`.
2. Read `docs/ssot/stack.md`.
3. Check `config/README.md` before editing config.
4. Run `./status.sh`.
5. Run `./test-rag.sh --baseline` and `./test-api.sh --baseline` if runtime behavior changed.
