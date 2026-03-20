# Repository Map

Last refreshed: 2026-03-20 (UTC)
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

## 2. Canonical config surfaces

Use `config/*` as the committed runtime edit surface.
Treat `docs/reference/openwebui/*` as snapshot/reference material, not as local deployment truth.

## 3. Top-level layout

Operationally relevant roots:

- root `*.sh`: small daily operator surface for deploy, status, logs, cleanup, update, and baseline validation of the current compose-managed stack.
- `config/`: canonical runtime configuration.
- `docs/`: runbooks, SSOT, guides, status notes, plans, reviews, and reference snapshots.
- `lib/`: shared shell helpers used by root scripts.
- `scripts/`: secondary maintenance, admin, testing, ingest, and optional-sidecar workflows.
- `services/`: tracked runtime APIs and workers that should use framework-native entrypoints instead of shell glue.
- `local/`: local-only helper namespace for sidecar binaries, optional tool payloads, and separate workspaces that should not clutter the repo root.
- `certs/`: local TLS material and placeholders; real certificate/key files stay ignored by default.
- `data/`: local corpus, PDFs, extractions, notes, and processing artifacts.
- `artifacts/`: generated audit and evaluation outputs.
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

- `scripts/testing/`: audit and script-validation helpers such as `scripts/testing/audit-no-mock.sh`, `scripts/testing/audit-openwebui-plugins.sh`, `scripts/testing/test-openwebui-tools-endpoints.sh`, `scripts/testing/verify-scripts.sh`
- `scripts/admin/`: OpenWebUI admin, repair, config-sync, and backup helpers such as `scripts/admin/openwebui-user-admin.sh`, `scripts/admin/repair-openwebui-tools.sh`, `scripts/admin/sync-openwebui-openai-config.sh`, `scripts/admin/sync-openwebui-web-search-config.sh`, `scripts/admin/sync-openwebui-retrieval-config.sh`, `scripts/admin/backup-openwebui-db.sh`, `scripts/admin/apply-openwebui-tool-patches.sh`
- `scripts/eval/`: operator-facing evaluation harnesses such as `scripts/eval/run-retrieval-eval.sh`
- `scripts/rag/`: retrieval helpers such as `scripts/rag/import-pdfs-to-kb.sh`, `scripts/rag/render-multimodal-answer.sh`, `scripts/rag/render-multimodal-answer-v2.sh`, `scripts/rag/migrate_rag_evidence.py` (one-shot text + figure migration driver), and `scripts/rag/build_multimodal_index.py`
- `services/multimodal_retrieval_api/`: canonical future RAG service; it serves `POST /api/v1/retrieve` against the `rag_evidence` one-collection runtime, queries Qdrant directly for text and visual lanes, and only uses a compatibility env file for shared defaults when present

Optional sidecar/tool flows:

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

## 6. Documentation tree

Documentation under `docs/` is split by purpose:

- `docs/ssot/`: runtime source of truth
- `docs/runbooks/`: setup, operations, and troubleshooting procedures
- `docs/guides/`: focused operator guidance
- `docs/status/`: observed local runtime snapshots and point-in-time deployment notes
- `docs/reviews/`: audits and technical reviews
- `docs/plans/`: historical or pending planning documents, including the multimodal redesign target in `docs/plans/MULTIMODAL_RAG_QDRANT_NEMOTRON_REDESIGN.md`
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
- Multimodal retrieval API (host-native target path): `127.0.0.1:${MULTIMODAL_RETRIEVAL_API_PORT:-8510} -> services.multimodal_retrieval_api.app:app`
- PostgreSQL: `${POSTGRES_BIND_ADDRESS:-127.0.0.1}:${POSTGRES_PORT:-5432} -> 5432`
- Jupyter: `127.0.0.1:${JUPYTER_PORT:-8890} -> ${JUPYTER_INTERNAL_PORT:-8889}`
- SearXNG: `${SEARXNG_BIND_ADDRESS:-127.0.0.1}:${SEARXNG_PORT:-8888} -> 8080`
- MCPO: `${MCPO_BIND_ADDRESS:-127.0.0.1}:${MCPO_PORT:-8000} -> 8000`
- Docling: `127.0.0.1:5001 -> 5001`
- Open Terminal: `${OPEN_TERMINAL_BIND_ADDRESS:-127.0.0.1}:${OPEN_TERMINAL_PORT:-8320} -> 8000`
- Indigo Service: `${INDIGO_SERVICE_BIND_ADDRESS:-127.0.0.1}:${INDIGO_SERVICE_PORT:-8012} -> 80`

## 9. Orientation path

For a normal repo pass:

1. Read `README.md`.
2. Read `docs/ssot/stack.md`.
3. Check `config/README.md` before editing config.
4. Run `./status.sh` for the current compose-managed baseline.
5. Run `./test-rag.sh --baseline` and `./test-api.sh --baseline` if baseline runtime behavior changed.
6. For canonical future RAG work, probe `/health`, `/ready`, and `POST /api/v1/retrieve` on `multimodal_retrieval_api`.
