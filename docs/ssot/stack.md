# Stack SSOT

Last updated: 2026-03-20 (UTC)
Scope: `/LAB/@thesis/openwebui`

## North Star + DoD

This repository should remain an operator-facing deployment layer for OpenWebUI-based RAG workflows, with Docker Compose as the supported current baseline, LiteLLM as the primary OpenAI-compatible upstream, the host-native `multimodal_retrieval_api` as the canonical future RAG path, and optional sidecars clearly separated from the core stack.

Definition of done for topology changes in this repo:

- deployment truth is recoverable from `docs/ssot/stack.md`, `docs/ssot/stack.yaml`, `README.md`, `config/env/.env.example`, and `config/compose/docker-compose.yml`
- `README.md` stays a front door, `docs/REPO_MAP.md` stays a layout map, and this document owns runtime topology
- docs distinguish the canonical edit surface (`config/*`) from the daily operator entrypoints
- docs keep the current compose-managed baseline distinct from the target host-native RAG path
- every critical service is named, justified, and backed by evidence
- optional profiles remain explicitly labeled optional
- drift between docs, env defaults, and scripts is recorded instead of hand-waved away

## Run Modes

### Docker primary

Default supported mode. Operators run the root shell entrypoints (`./deploy.sh`, `./status.sh`, `.env`), while compose and committed defaults come from `config/*`.

Evidence:

- `deploy.sh`
- `status.sh`
- `config/compose/docker-compose.yml`
- `config/README.md`
- `scripts/check-doc-consistency.sh`

### Docker RG build variant

Supported build variant for `openwebui` only. `config/compose/docker-compose.rg.yml` overrides the `openwebui` image build to use `config/docker/Dockerfile.openwebui-rg`; it is not a second full-stack topology.

Evidence:

- `config/compose/docker-compose.rg.yml`
- `config/docker/Dockerfile.openwebui-rg`

### Host-assisted upstreams

The Docker stack depends on host-provided or host-mounted resources, especially LiteLLM at `host.docker.internal:4000/v1` and the host Zotero directory mounted into `mcpo`.

Evidence:

- `config/env/.env.example`
- `config/compose/docker-compose.yml`
- `README.md`

### Canonical Future RAG Path

The native `multimodal_retrieval_api` is the repo's canonical future RAG path.
Today it is still run separately on the host rather than as part of the compose
baseline, but RAG architecture and migration work should treat it as the target
runtime over the `rag_evidence` one-collection text and visual lanes.

Evidence:

- `services/README.md`
- `README.md`
- `services/multimodal_retrieval_api/service.py`

### Optional profile sidecars

`open-terminal` and `indigo-service` are opt-in profiles or feature-flagged sidecars. `deploy.sh` starts or stops them only when the corresponding env flags are enabled.

Evidence:

- `config/compose/docker-compose.yml`
- `deploy.sh`
- `README.md`

### Kubernetes / system service modes

Status: unsupported / no evidence found in the tracked operator and config trees during this SSOT refresh. Add a new run mode entry before claiming k8s, Helm, Procfile, or systemd support.

Evidence:

- `command: find config docs scripts archive -type f \( -name 'Procfile' -o -name '*.service' -o -name 'Chart.yaml' -o -name 'values*.yaml' \) -print` returned no results on 2026-03-11

## Service Registry

### Critical services

1. `openwebui`
- Role: primary UI/API and orchestration surface.
- Critical because it is the user entrypoint and depends on `postgres`, `jupyter`, `mcpo`, and `docling`.
- Run mode: `docker-primary`
- Evidence: `config/compose/docker-compose.yml`, `README.md`

2. `litellm`
- Role: primary OpenAI-compatible upstream for model discovery and generation.
- Critical because it is the external adapter on the main E2E query path.
- Run mode: `host-assisted-upstreams`
- Evidence: `config/env/.env.example`, `README.md`, `status.sh`

3. `postgres`
- Role: stateful postgres/pgvector service for application data and the default OpenWebUI vector backend in this repo.
- Critical because it stores persistent state and serves multiple dependents.
- Important nuance: postgres is always deployed and the committed `config/env/.env.example` now defaults OpenWebUI itself to `VECTOR_DB=pgvector`, so pgvector is the baseline retrieval backend for this repo.
- Run mode: `docker-primary`
- Evidence: `config/compose/docker-compose.yml`, `config/env/.env.example`

4. `jupyter`
- Role: code execution and code-interpreter backend.
- Critical because code execution is enabled by default and `openwebui` waits for a healthy Jupyter service.
- Run mode: `docker-primary`
- Evidence: `config/compose/docker-compose.yml`, `config/env/.env.example`, `config/jupyter/jupyter_server_config.py`

5. `mcpo`
- Role: MCP-to-OpenAPI bridge used for OpenWebUI tool exposure.
- Critical because tool integrations depend on it and `openwebui` waits for it during startup.
- Important nuance: the committed env template now pins `MCPO_IMAGE` by digest instead of tracking the floating `main` tag, so fresh deploys stay reproducible across time.
- Run mode: `docker-primary`
- Evidence: `config/compose/docker-compose.yml`, `config/mcp/config.json`

6. `docling`
- Role: multimodal extraction backend.
- Critical because `openwebui` waits for a healthy `docling` service and ingest features depend on it.
- Run mode: `docker-primary`
- Evidence: `config/compose/docker-compose.yml`, `README.md`

### Supporting and optional services

7. `searxng`
- Role: local web-search backend for OpenWebUI web search.
- Non-critical because it is only needed when OpenWebUI web search is enabled.
- Run mode: `optional-profile-sidecars`
- Evidence: `config/compose/docker-compose.yml`, `config/env/.env.example`, `config/searxng/settings.yml`

8. `open-terminal`
- Role: optional terminal and file-server sidecar.
- Non-critical because it is profile-gated and disabled by default.
- Run mode: `optional-profile-sidecars`
- Evidence: `config/compose/docker-compose.yml`, `config/env/.env.example`, `README.md`

9. `indigo-service`
- Role: optional chemistry REST sidecar for Indigo-based tooling.
- Non-critical because it is profile-gated and disabled by default.
- Important nuance: the committed env template now pins `INDIGO_SERVICE_IMAGE` by digest even though the upstream workflow still advertises `latest`.
- Run mode: `optional-profile-sidecars`
- Evidence: `config/compose/docker-compose.yml`, `config/env/.env.example`, `README.md`

10. `multimodal_retrieval_api`
- Role: canonical future retrieval API for one-collection multimodal RAG over `rag_evidence`.
- Non-critical in the current local baseline because it is not yet compose-managed and is still operated separately during the migration window.
- Operator note: docs and config should expose its URL, port, and required host-native env so migration work does not depend on ad hoc local knowledge.
- Important nuance: it now uses native `qdrant-client`, native Transformers/FastEmbed query encoding, lane-based readiness reporting, and a sampled `rag_evidence` completeness probe; it no longer delegates live retrieval through `thesis_graph`.
- Run mode: `host-native-canonical-rag`
- Evidence: `services/README.md`, `README.md`, `services/multimodal_retrieval_api/service.py`

## Critical Paths

### Deploy / run baseline

`operator -> .env -> deploy.sh -> docker compose -> postgres + jupyter + mcpo + docling -> openwebui -> status/test scripts`

Evidence:

- `README.md`
- `deploy.sh`
- `config/compose/docker-compose.yml`
- `status.sh`

### Query / generation (default committed settings)

`browser -> openwebui -> postgres/pgvector -> LiteLLM(host) -> upstream provider`

This is the default OpenWebUI retrieval/generation path implied by `config/env/.env.example` because `VECTOR_DB=pgvector` is now the committed baseline for this repo.

Evidence:

- `config/env/.env.example`
- `config/compose/docker-compose.yml`
- `README.md`

### Native multimodal retrieval API

`client -> multimodal_retrieval_api -> qdrant(rag_evidence) -> text lane + visual lane`

This path is the canonical future RAG architecture for the repo. It is still
operated separately from the compose-managed OpenWebUI pgvector baseline during
the migration window and uses `rag_evidence` as its sole runtime evidence
source. Missing figure payloads are surfaced as runtime diagnostics instead of
triggering local extraction fallbacks.

Operational nuance:

- `/ready` now samples `page` points in `rag_evidence` and reports collection completeness for native `figure_records`.
- A `degraded` readiness state can therefore mean runtime prerequisites are missing, or that sampled native evidence coverage is incomplete.
- `config/env/.env.example` should therefore carry the host-native URL, port, Qdrant URL, and text-model path knobs needed to run and probe this service consistently.

Evidence:

- `services/README.md`
- `README.md`
- `services/multimodal_retrieval_api/service.py`

### Query / generation (chroma fallback mode)

`browser -> openwebui -> chroma in openwebui_data -> LiteLLM(host) -> upstream provider`

This path is still supported when `VECTOR_DB=chroma` is set explicitly.

Evidence:

- `config/env/.env.example`
- `config/compose/docker-compose.yml`

### Ingest

`scripts/rag/import-pdfs-to-kb.sh` -> OpenWebUI APIs -> OpenWebUI data store`

Evidence:

- `scripts/rag/import-pdfs-to-kb.sh`
- `README.md`

### Eval / regression

`local scripts -> compose services -> artifacts`

Evidence:

- `test-rag.sh`
- `test-api.sh`
- `.mise.toml`

## External Dependencies

- Docker daemon and Docker Compose v2
- LiteLLM host endpoint and upstream provider credentials
- optional host Zotero data directory mounted into `mcpo`
Evidence:

- `README.md`
- `config/env/.env.example`
- `config/compose/docker-compose.yml`
- `config/mcp/config.json`

## Drift Radar

- `config/*` is the canonical edit surface for committed runtime behavior.
- Compose environment fallbacks must match the committed env template, especially for baseline retrieval flags such as `VECTOR_DB` and `ENABLE_RAG_HYBRID_SEARCH`.
- Critical and shipped sidecar images should be pinned to immutable digests in `config/env/.env.example` when upstream defaults are floating tags such as `main` or `latest`.
- postgres is always deployed and the committed default now uses `VECTOR_DB=pgvector`; if you switch a local runtime to `chroma`, document that as an intentional override because it changes where vectors persist.
- Optional sidecars remain supported, but the normal baseline should keep them disabled unless a concrete workflow needs them.
- `multimodal_retrieval_api` is the canonical future RAG path; docs must distinguish that target architecture from the current compose-managed OpenWebUI baseline until migration completes.
- Root docs and runbooks should expose direct `/health`, `/ready`, and `POST /api/v1/retrieve` probes for the target service even while baseline validation still centers on OpenWebUI retrieval.
- The native multimodal retrieval path now targets `rag_evidence` directly; do not describe it as delegating to `thesis_graph`.
- Legacy project containers outside the current compose config should be removed instead of being treated as part of the supported topology.

## Observed Local Runtime

Observed from `./status.sh`, `.env`, `docker compose --project-directory . -f config/compose/docker-compose.yml config --services`, and `docker compose --project-directory . -f config/compose/docker-compose.yml ps --format json` on 2026-03-13 after baseline normalization.

- Live UI endpoint: `http://localhost:3000`
- Live OpenWebUI image: `ghcr.io/open-webui/open-webui:v0.8.10`
- Live retrieval mode: `VECTOR_DB=pgvector`
- Live embedding model: `pritamdeka/S-PubMedBert-MS-MARCO`
- Live web-search route: configured via `http://searxng:8080/search?q={query}&format=json`, but disabled in the baseline runtime (`ENABLE_WEB_SEARCH=false`) so the `web-search` profile is intentionally not started
- Live SearXNG host probe: `http://127.0.0.1:8888`
- Live optional sidecars: disabled in the baseline runtime
- Local legacy drift: none expected once orphan containers are removed

## Update Protocol

Update this SSOT when any of the following change:

- `config/compose/docker-compose.yml` or `config/compose/docker-compose.rg.yml`
- `config/env/.env.example` defaults that affect topology, ports, or primary upstreams
- `deploy.sh`, `status.sh`, `update.sh`, or `cleanup.sh`
- profile-gated sidecars and their enablement flags
- CI gates that change the deploy, ingest, or eval path

Minimum refresh checklist:

1. Re-scan `config/compose/docker-compose*.yml`, `config/env/.env.example`, `deploy.sh`, and `status.sh`.
2. Update `docs/ssot/stack.md` and `docs/ssot/stack.yaml`.
3. Update `README.md` if the operator-facing summary changed, and update `docs/REPO_MAP.md` only when repo layout, entrypoint ownership, or config-root navigation changed.
4. Re-run the baseline validation scripts after runtime changes.
5. Run `./scripts/check-doc-consistency.sh`.
6. If runtime behavior changed, run `./status.sh`, `./test-rag.sh --baseline`, and `./test-api.sh --baseline`.

## Changelog

- `2026-03-20`: reclassified `multimodal_retrieval_api` as the canonical future RAG path, removed adjunct wording from the SSOT, and clarified that the current compose-managed baseline still remains the OpenWebUI `pgvector` path during migration.
- `2026-03-20`: promoted the canonical future RAG path across the operator front door, repo map, env template, and runbooks without claiming that the compose-managed baseline migration is already complete.
- `2026-03-15`: added the host-native multimodal retrieval API as a supported adjunct run mode and service, documented the `rag_evidence` one-collection query path, and removed stale references to live `thesis_graph` delegation from the runtime topology.
- `2026-03-13`: aligned `config/compose/docker-compose.yml` environment fallbacks with the committed Docker baseline by restoring `pgvector` and hybrid retrieval as the default fallback values, refreshed the observed-runtime date, and tightened drift guidance so compose defaults, `config/env/.env.example`, and the SSOT stay synchronized.
- `2026-03-13`: pinned `MCPO_IMAGE` and `INDIGO_SERVICE_IMAGE` to verified digests in the committed env template so fresh deploys do not drift with upstream `main` or `latest` tags.
- `2026-03-11`: normalized the repo and local baseline around a regular OpenWebUI RAG topology by making `pgvector` the committed default backend, aligning the committed embedding baseline with the currently indexed OpenWebUI corpus, removing repo-level SMILES retrieval env defaults, returning local runtime ports/image/sidecars to the standard baseline, and treating legacy compose orphans as cleanup targets instead of supported topology.
- `2026-03-10`: refreshed the SSOT against current runtime entrypoints, added LiteLLM and Docling as explicit critical services, separated default `chroma` mode from optional OpenWebUI `pgvector` mode, corrected the live MCPO default back to `ghcr.io/open-webui/mcpo:main`, and removed stale legacy-upstream wording from operator scripts.
- `2026-03-09`: created explicit stack SSOT and aligned repo-navigation docs around canonical config/docs/ops surfaces.
