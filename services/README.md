# Services

Tracked service code lives here when the repository needs a real API or worker,
not another operator wrapper.

Current service roots:

- `multimodal_retrieval_api/`: Starlette retrieval service for the
  one-collection multimodal RAG path. This is the repo's canonical future RAG
  service and it queries `rag_evidence` directly in Qdrant for text and visual
  evidence. Text query encoding is native to this repo through Transformers +
  FastEmbed. A compatibility env file may still be used for shared model or
  Qdrant defaults, but only when explicitly configured.

Direct launch example:

```bash
PYTHONPATH=/LAB/@thesis/openwebui \
MULTIMODAL_RETRIEVAL_API_QDRANT_URL=http://127.0.0.1:6335 \
MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH=/path/to/nemotron-model \
MULTIMODAL_RETRIEVAL_API_TEXT_SPARSE_MODEL_NAME=Qdrant/bm25 \
python3 -m uvicorn \
  services.multimodal_retrieval_api.app:app \
  --host 127.0.0.1 \
  --port 8510
```

Optional shared-defaults file:

```bash
MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE=/absolute/path/to/shared-defaults.env
```

Use this only as a fallback source for `QDRANT_URL`, `QDRANT_API_KEY`,
`NEMOTRON_MODEL_PATH`, `NEMOTRON_DEVICE`, `NEMOTRON_ATTN_IMPLEMENTATION`, and
`SPARSE_MODEL_NAME` when you do not want to export the service-specific env vars
directly.

Readiness:

- `/ready` reports service-level status plus lane-specific runtime readiness for:
  - `text`
  - `visual`
- `/ready.qdrant.collection_completeness` also samples `page` points in
  `rag_evidence` and reports whether native `figure_records` are present.
- Each lane also reports its configured metadata, including vector names and
  model identifiers where applicable.
- `runtime.capabilities` provides a compact operator summary of whether each
  lane is configured and currently usable.
- `status: "degraded"` means at least one runtime lane or config prerequisite is
  missing, or the sampled `rag_evidence` contract is incomplete, even if Qdrant
  connectivity is healthy.

Rules:

- Treat this tree as runtime code, not local scratch space.
- Prefer direct framework interfaces over shell glue.
- Keep contracts explicit and versioned in code.
- Migration-only helpers for extraction backfills belong under `scripts/rag/`, not in this runtime service tree.
