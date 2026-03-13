# Repository Map

Last refreshed: 2026-03-13 (UTC)
Scope: `/LAB/@thesis/openwebui`

This file owns repository shape and operator navigation. It is intentionally not a second runtime spec.

Use this file when you need to know:
- where a file should live
- whether a path is canonical or only a compatibility surface
- which root commands and local namespaces are operator-facing

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

- root `*.sh`: small daily operator surface for deploy, status, logs, cleanup, update, and baseline validation.
- `config/`: canonical runtime configuration.
- `docs/`: runbooks, SSOT, guides, status notes, plans, reviews, and reference snapshots.
- `lib/`: shared shell helpers used by root scripts.
- `scripts/`: secondary maintenance, admin, testing, ingest, and optional-sidecar workflows.
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

The root now keeps only the daily operator surface. Secondary tooling lives under `scripts/*`.

Primary daily surface:

- Deploy and status: `deploy.sh`, `status.sh`, `logs.sh`, `cleanup.sh`, `update.sh`
- Baseline validation: `test-rag.sh --baseline`, `test-api.sh --baseline`

Secondary scripts under `scripts/`:

- `scripts/testing/`: audit and script-validation helpers such as `scripts/testing/audit-no-mock.sh`, `scripts/testing/audit-openwebui-plugins.sh`, `scripts/testing/test-openwebui-tools-endpoints.sh`, `scripts/testing/test-update-smoke.sh`, `scripts/testing/verify-scripts.sh`
- `scripts/admin/`: OpenWebUI admin, repair, config-sync, and backup helpers such as `scripts/admin/openwebui-user-admin.sh`, `scripts/admin/repair-openwebui-tools.sh`, `scripts/admin/sync-openwebui-openai-config.sh`, `scripts/admin/sync-openwebui-web-search-config.sh`, `scripts/admin/backup-openwebui-db.sh`, `scripts/admin/apply-openwebui-tool-patches.sh`
- `scripts/rag/`: ingest and retrieval-tuning flows such as `scripts/rag/import-pdfs-to-kb.sh`, `scripts/rag/tune-openwebui-documents.sh`, `scripts/rag/manage-openwebui-embedding-profiles.sh`, `scripts/rag/llm-council.sh`

Optional sidecar/tool flows:

- CLIProxyAPI lifecycle/bootstrap/auth: `scripts/cliproxyapi/setup-cliproxyapi.sh`, `scripts/cliproxyapi/start-cliproxyapi.sh`, `scripts/cliproxyapi/stop-cliproxyapi.sh`, `scripts/cliproxyapi/restart-cliproxyapi.sh`, `scripts/cliproxyapi/check-cliproxyapi.sh`, `scripts/cliproxyapi/configure-cliproxyapi-oauth.sh`, `scripts/cliproxyapi/test-cliproxyapi-oauth.sh`, `scripts/cliproxyapi/import-qwen-auth.sh`, `scripts/cliproxyapi/cli-proxy-api.sh`, `local/bin/cliproxyapi`
- Indigo: `scripts/indigo/start-indigo-service.sh`, `scripts/indigo/stop-indigo-service.sh`, `scripts/indigo/restart-indigo-service.sh`, `scripts/indigo/check-indigo-service.sh`, `scripts/indigo/enable-indigo-live.sh`, `local/plugins/indigo_chemistry_tool.py`
- Open Terminal: `scripts/open-terminal/test-openwebui-open-terminal.sh`
- MCP-specific helpers: `scripts/mcp/configure-mcpo-openapi-servers.sh`, `scripts/mcp/test-mcp.sh`
- Entity-maintenance workspace when present: `local/entities/README.md`, `local/entities/pipelines/`

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
