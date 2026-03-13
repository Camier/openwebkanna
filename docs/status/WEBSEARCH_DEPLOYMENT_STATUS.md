# Web Search Deployment Status

Observed on 2026-03-11 against the normalized local runtime in `/LAB/@thesis/openwebui`.

This file is an observed local-runtime snapshot, not the canonical supported-topology document for the repo. For desired-state topology and committed defaults, use `docs/ssot/stack.md`, `config/env/.env.example`, and `config/compose/docker-compose.yml`.

## Status

Web search infrastructure is configured, but intentionally not started in the normalized local runtime because OpenWebUI web search is disabled by default.

| Layer | Observed status | Evidence |
|-------|-----------------|----------|
| Compose service | `searxng` not started in baseline runtime | `./status.sh`, `docker compose ps --format json` |
| OpenWebUI integration | disabled in baseline runtime | `.env`, `./status.sh` |
| Container route | `http://searxng:8080/search?q={query}&format=json` | `.env`, `./status.sh` |
| Host probe | `http://127.0.0.1:8888` | `.env`, `./status.sh` |
| OpenWebUI UI | `http://localhost:3000` | `.env`, `./status.sh` |

## Committed Defaults vs Local Runtime

The normalized local runtime now matches the repo baseline for ports, image, and web-search defaults.

| Setting | Committed default (`.env.example`) | Observed local runtime (`.env`) |
|---------|------------------------------------|---------------------------------|
| `WEBUI_PORT` | `3000` | `3000` |
| `SEARXNG_PORT` | `8888` | `8888` reserved for the optional `web-search` profile |
| `ENABLE_WEB_SEARCH` | `false` | `false` |
| `VECTOR_DB` | `pgvector` | `pgvector` |
| `OPENWEBUI_IMAGE` | `ghcr.io/open-webui/open-webui:v0.8.10` | `ghcr.io/open-webui/open-webui:v0.8.10` |

## Current Effective Configuration

Local runtime values observed on 2026-03-11:

```bash
WEBUI_PORT=3000
SEARXNG_PORT=8888
ENABLE_WEB_SEARCH=false
ENABLE_WEBSEARCH=false
WEB_SEARCH_ENGINE=searxng
SEARXNG_QUERY_URL=http://searxng:8080/search?q={query}&format=json
VECTOR_DB=pgvector
RAG_EMBEDDING_MODEL=pritamdeka/S-PubMedBert-MS-MARCO
```

Current URLs:

| Service | URL |
|---------|-----|
| OpenWebUI | `http://localhost:3000` |
| SearXNG host probe | not active in baseline runtime |
| SearXNG container route | `http://searxng:8080` |
| SearXNG query route | `http://searxng:8080/search?q={query}&format=json` |

## Verification Commands

```bash
./status.sh
docker compose --profile web-search ps searxng openwebui
docker compose --profile web-search up -d searxng
curl -fsS "http://127.0.0.1:8888/search?q=openwebui&format=json" | jq '.results | length'
docker exec openwebui env | grep -E '^(ENABLE_WEB_SEARCH|ENABLE_WEBSEARCH|WEB_SEARCH_ENGINE|SEARXNG_QUERY_URL|VECTOR_DB)='
docker exec openwebui curl -fsS "http://searxng:8080/search?q=status+probe&format=json" | jq '.results | length'
```

## Drift Notes

- This document tracks the observed local deployment status after normalization, not a speculative future state.
- The canonical supported topology remains in [../ssot/stack.md](../ssot/stack.md).
- Enable OpenWebUI web search by setting `ENABLE_WEB_SEARCH=true` and `ENABLE_WEBSEARCH=true` in `.env`, then redeploy; `deploy.sh` will start the `web-search` profile automatically.

## Related Documentation

- [Stack SSOT](../ssot/stack.md)
- [Web Search Deep Dive](../reference/openwebui/OPENWEBUI_WEBSEARCH_DEEP_DIVE.md)
- [Environment Variables](../reference/openwebui/openwebui_env_reference.md)
- [RAG Technical Reference](../reference/openwebui/openwebui_rag_technical_reference.md)

Documented against the normalized local runtime on 2026-03-11.
