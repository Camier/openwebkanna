# OpenWebUI Repo Agent Notes

Use this file as a routing aid for agents working in `/LAB/@thesis/openwebui`.
It is intentionally short so it does not compete with the actual repo docs.

## Read Order

Use these files in this order:

- `README.md`: operator front door and daily command surface
- `docs/ssot/stack.md`: canonical supported topology and runtime contract
- `docs/REPO_MAP.md`: repository layout and operator entrypoints
- `config/README.md`: canonical config edit surface and compatibility-copy map
- `docs/README.md`: documentation routing and doc-scope rules
- `docs/runbooks/*.md`: procedure-level operating detail

## Ownership Rules

- Root `*.sh` scripts are the canonical operator surface.
- `config/*` is the canonical edit surface for runtime behavior.
- Root `.env.example`, `docker-compose*.yml`, `mcp/`, `jupyter/`, and `searxng/` are compatibility copies or compatibility entrypoints. Edit `config/*` first.
- `docs/status/*` are observed local-runtime snapshots, not desired-state topology docs.
- `docs/reference/openwebui/*` are upstream/reference snapshots, not local deployment truth.
- `docs/plans/*` and `docs/legacy/*` are historical/context material, not active SSOT.

## Runtime Baseline

Normal local baseline in this repo:

- OpenWebUI on `http://localhost:3000`
- LiteLLM upstream via `http://host.docker.internal:4000/v1`
- PostgreSQL + `pgvector`
- Core services: `openwebui`, `postgres`, `jupyter`, `mcpo`, `docling`
- Optional sidecars disabled by default unless explicitly enabled

Confirm the live baseline in:

- `docs/ssot/stack.md`
- `config/env/.env.example`
- `config/compose/docker-compose.yml`
- `.env`
- `./status.sh`

## Verification

After behavior or config changes, prefer:

```bash
./scripts/check-doc-consistency.sh
docker compose config -q
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

Last updated: 2026-03-12 (UTC)
