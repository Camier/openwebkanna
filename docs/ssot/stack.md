# Stack SSOT

Last updated: 2026-03-13 (UTC)
Scope: `/LAB/@thesis/openwebui`

## North Star + DoD

This repository should remain an operator-facing deployment layer for OpenWebUI-based RAG workflows, with Docker Compose as the supported baseline, LiteLLM as the primary OpenAI-compatible upstream, and optional sidecars clearly separated from the core stack.

Definition of done for topology changes in this repo:

- deployment truth is recoverable from `docs/ssot/stack.md`, `docs/ssot/stack.yaml`, `README.md`, `.env.example`, and `docker-compose.yml`
- `README.md` stays a front door, `docs/REPO_MAP.md` stays a layout map, and this document owns runtime topology
- docs distinguish the canonical edit surface (`config/*`) from the live operator entrypoints/root compatibility copies
- every critical service is named, justified, and backed by evidence
- optional profiles remain explicitly labeled optional
- drift between docs, env defaults, and scripts is recorded instead of hand-waved away

## Run Modes

### Docker primary

Default supported mode. Operators run the root compatibility entrypoints (`./deploy.sh`, `./status.sh`, `docker-compose.yml`, `.env`), while `config/*` remains the canonical edit surface and is checked for byte-for-byte alignment.

Evidence:

- `deploy.sh`
- `status.sh`
- `docker-compose.yml`
- `config/README.md`
- `scripts/check-doc-consistency.sh`

### Docker RG build variant

Supported build variant for `openwebui` only. `docker-compose.rg.yml` overrides the `openwebui` image build to use `Dockerfile.openwebui-rg`; it is not a second full-stack topology.

Evidence:

- `docker-compose.rg.yml`
- `Dockerfile.openwebui-rg`

### Host-assisted upstreams

The Docker stack depends on host-provided or host-mounted resources, especially LiteLLM at `host.docker.internal:4000/v1` and the host Zotero directory mounted into `mcpo`.

Evidence:

- `.env.example`
- `docker-compose.yml`
- `README.md`

### Optional profile sidecars

`cliproxyapi`, `open-terminal`, and `indigo-service` are opt-in profiles or feature-flagged sidecars. `deploy.sh` starts or stops them only when the corresponding env flags are enabled.

Evidence:

- `docker-compose.yml`
- `deploy.sh`
- `README.md`

### Archived host fallback

Archived host-side `vLLM` lifecycle scripts remain for fallback/reference only and are not part of the default deploy path.

Evidence:

- `archive/README.md`
- `archive/start-vllm.sh`
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
- Evidence: `docker-compose.yml`, `README.md`

2. `litellm`
- Role: primary OpenAI-compatible upstream for model discovery and generation.
- Critical because it is the external adapter on the main E2E query path.
- Run mode: `host-assisted-upstreams`
- Evidence: `.env.example`, `README.md`, `status.sh`

3. `postgres`
- Role: stateful postgres/pgvector service for application data and the default OpenWebUI vector backend in this repo.
- Critical because it stores persistent state and serves multiple dependents.
- Important nuance: postgres is always deployed and the committed `.env.example` now defaults OpenWebUI itself to `VECTOR_DB=pgvector`, so pgvector is the baseline retrieval backend for this repo.
- Run mode: `docker-primary`
- Evidence: `docker-compose.yml`, `.env.example`

4. `jupyter`
- Role: code execution and code-interpreter backend.
- Critical because code execution is enabled by default and `openwebui` waits for a healthy Jupyter service.
- Run mode: `docker-primary`
- Evidence: `docker-compose.yml`, `.env.example`, `jupyter/jupyter_server_config.py`

5. `mcpo`
- Role: MCP-to-OpenAPI bridge used for OpenWebUI tool exposure.
- Critical because tool integrations depend on it and `openwebui` waits for it during startup.
- Run mode: `docker-primary`
- Evidence: `docker-compose.yml`, `mcp/config.json`

6. `docling`
- Role: multimodal extraction backend.
- Critical because `openwebui` waits for a healthy `docling` service and ingest features depend on it.
- Run mode: `docker-primary`
- Evidence: `docker-compose.yml`, `README.md`

### Supporting and optional services

7. `searxng`
- Role: local web-search backend for OpenWebUI web search.
- Non-critical because it is only needed when OpenWebUI web search is enabled.
- Run mode: `optional-profile-sidecars`
- Evidence: `docker-compose.yml`, `.env.example`, `config/searxng/settings.yml`

8. `cliproxyapi`
- Role: optional legacy OAuth alias sidecar.
- Non-critical because LiteLLM is the primary upstream and this sidecar is profile-gated.
- Run mode: `optional-profile-sidecars`
- Evidence: `docker-compose.yml`, `.env.example`, `cliproxyapi/config.yaml`

9. `open-terminal`
- Role: optional terminal and file-server sidecar.
- Non-critical because it is profile-gated and disabled by default.
- Run mode: `optional-profile-sidecars`
- Evidence: `docker-compose.yml`, `.env.example`, `README.md`

10. `indigo-service`
- Role: optional chemistry REST sidecar for Indigo-based tooling.
- Non-critical because it is profile-gated and disabled by default.
- Run mode: `optional-profile-sidecars`
- Evidence: `docker-compose.yml`, `.env.example`, `README.md`

## Critical Paths

### Deploy / run baseline

`operator -> .env -> deploy.sh -> docker compose -> postgres + jupyter + mcpo + docling -> openwebui -> status/test scripts`

Evidence:

- `README.md`
- `deploy.sh`
- `docker-compose.yml`
- `status.sh`

### Query / generation (default committed settings)

`browser -> openwebui -> postgres/pgvector -> LiteLLM(host) -> upstream provider`

This is the default OpenWebUI retrieval/generation path implied by `.env.example` because `VECTOR_DB=pgvector` is now the committed baseline for this repo.

Evidence:

- `.env.example`
- `docker-compose.yml`
- `README.md`

### Query / generation (chroma fallback mode)

`browser -> openwebui -> chroma in openwebui_data -> LiteLLM(host) -> upstream provider`

This path is still supported when `VECTOR_DB=chroma` is set explicitly.

Evidence:

- `.env.example`
- `docker-compose.yml`

### Ingest

`import-pdfs-to-kb.sh` -> OpenWebUI APIs -> OpenWebUI data store`

Evidence:

- `import-pdfs-to-kb.sh`
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
- `.env.example`
- `docker-compose.yml`
- `mcp/config.json`

## Drift Radar

- `config/*` is the canonical edit surface, but live scripts execute the root compatibility copies; both surfaces must remain byte-aligned.
- Compose environment fallbacks must match the committed env template, especially for baseline retrieval flags such as `VECTOR_DB` and `ENABLE_RAG_HYBRID_SEARCH`.
- postgres is always deployed and the committed default now uses `VECTOR_DB=pgvector`; if you switch a local runtime to `chroma`, document that as an intentional override because it changes where vectors persist.
- Optional sidecars remain supported, but the normal baseline should keep them disabled unless a concrete workflow needs them.
- Legacy project containers outside the current compose config should be removed instead of being treated as part of the supported topology.

## Observed Local Runtime

Observed from `./status.sh`, `.env`, `docker compose config --services`, and `docker compose ps --format json` on 2026-03-13 after baseline normalization.

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

- `docker-compose.yml` or `docker-compose.rg.yml`
- `.env.example` defaults that affect topology, ports, or primary upstreams
- `deploy.sh`, `status.sh`, `update.sh`, or `cleanup.sh`
- profile-gated sidecars and their enablement flags
- CI gates that change the deploy, ingest, or eval path

Minimum refresh checklist:

1. Re-scan `docker-compose*.yml`, `.env.example`, `deploy.sh`, and `status.sh`.
2. Update `docs/ssot/stack.md` and `docs/ssot/stack.yaml`.
3. Update `README.md` if the operator-facing summary changed, and update `docs/REPO_MAP.md` only when repo layout, entrypoint ownership, or config-root navigation changed.
4. Restore any broken root compatibility copies before running consistency checks.
5. Run `./scripts/check-doc-consistency.sh`.
6. If runtime behavior changed, run `./status.sh`, `./test-rag.sh --baseline`, and `./test-api.sh --baseline`.

## Changelog

- `2026-03-13`: aligned `docker-compose.yml` environment fallbacks with the committed Docker baseline by restoring `pgvector` and hybrid retrieval as the default fallback values, refreshed the observed-runtime date, and tightened drift guidance so compose defaults, `.env.example`, and the SSOT stay synchronized.
- `2026-03-11`: normalized the repo and local baseline around a regular OpenWebUI RAG topology by making `pgvector` the committed default backend, aligning the committed embedding baseline with the currently indexed OpenWebUI corpus, removing repo-level SMILES retrieval env defaults, returning local runtime ports/image/sidecars to the standard baseline, and treating legacy compose orphans as cleanup targets instead of supported topology.
- `2026-03-10`: refreshed the SSOT against current runtime entrypoints, added LiteLLM and Docling as explicit critical services, separated default `chroma` mode from optional OpenWebUI `pgvector` mode, restored the root `searxng/` compatibility copy, corrected the live MCPO default back to `ghcr.io/open-webui/mcpo:main`, and removed stale vLLM-primary wording from operator scripts.
- `2026-03-09`: created explicit stack SSOT and aligned repo-navigation docs around canonical config/docs/ops surfaces.
