# SMILES + RAG Consolidation

Last updated: 2026-03-03 (UTC)

This document is the operational single source of truth for the current text + SMILES retrieval setup in this repository.

## Active profile

1. Text KB retrieval keeps a dedicated text embedding model via `RAG_EMBEDDING_MODEL`.
2. SMILES retrieval is locked to one stable fingerprint channel: `ecfp4`.
3. Molecular embedding is optional and additive (`SMILES_ENABLE_EMBEDDING_SIGNAL=true`), not a replacement for fingerprints.
4. Multi-channel fusion uses Reciprocal Rank Fusion (RRF).
5. Validation tracks text and SMILES recall independently before and after fusion.

## Configuration baseline

Set these in `.env` (or keep `.env.example` defaults where applicable):

```bash
# LiteLLM is the reference upstream
OPENAI_API_BASE_URL=http://host.docker.internal:4000/v1
OPENAI_API_BASE_URLS=http://host.docker.internal:4000/v1
CLIPROXYAPI_ENABLED=false

# Text retrieval
RAG_EMBEDDING_MODEL=sentence-transformers/all-MiniLM-L6-v2

# Retrieval fusion
RETRIEVAL_FUSION_STRATEGY=rrf
RETRIEVAL_FUSION_RRF_K=60
RETRIEVAL_FUSION_TOP_K=20
# Optional per-channel weighting (name=value CSV)
RETRIEVAL_CHANNEL_WEIGHTS=text_dense=1.0,bm25=1.0,smiles_fingerprint=1.0,smiles_embedding=0.9

# SMILES retrieval
SMILES_RETRIEVAL_ENABLED=true
SMILES_FINGERPRINT_TYPES=ecfp4
SMILES_FINGERPRINT_INDEX_TYPE=ivfflat
SMILES_ENABLE_EMBEDDING_SIGNAL=false
SMILES_EMBEDDING_MODEL=ibm/MoLFormer-XL-both-10pct
```

## Runtime channels

1. `text_dense`: text embedding retrieval from KB chunks.
2. `bm25`: sparse keyword retrieval.
3. `smiles_fingerprint`: structure similarity using ECFP4 vectors and dedicated DB column/index.
4. `smiles_embedding` (optional): molecular encoder signal, enabled explicitly.
5. Fusion output: RRF over available channels.

## API behavior

1. `POST /v1/structure/search` accepts `fingerprint_type=ecfp4` in the stable profile.
2. `POST /v1/structure/batch-search` accepts query param `fingerprint_type=ecfp4` in the stable profile.
3. `POST /v1/retrieval/fuse` fuses `text_dense`, `bm25`, `smiles_fingerprint`, and optional `smiles_embedding`.
4. `GET /v1/structure/stats` exposes combined/per-channel coverage and runtime-enabled fingerprint types.

## Deferred channels

1. `maccs` is intentionally deferred for now to minimize operational drift.
2. Re-enable only after DB column/index parity and eval metrics are validated end-to-end.

## Operations

Activate the working env used for validation:

```bash
source /home/miko/miniforge3/etc/profile.d/conda.sh
conda activate smiles-extraction
```

Validate retrieval metrics with sample data:

```bash
mise run smiles:eval-retrieval
# Eval-only toggle for the mise task flag (--enable-smiles-embedding)
ENABLE_SMILES_EMBEDDING=1 mise run smiles:eval-retrieval
```

Note: `ENABLE_SMILES_EMBEDDING` is only an eval-task switch.
Runtime wiring still uses `.env` key `SMILES_ENABLE_EMBEDDING_SIGNAL`.

Run targeted tests:

```bash
pytest -q smiles-pipeline/tests/test_rrf_fusion.py \
  smiles-pipeline/tests/test_retrieval_eval.py \
  smiles-pipeline/tests/test_rag_fingerprint_storage.py
pytest -q test_rag_smiles_upgrade.py
```

Run stack baseline checks:

```bash
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```

## Verification snapshot (2026-03-03)

1. Env: `smiles-extraction`, Python `3.11.14`.
2. Unit checks: RRF/retrieval/fingerprint tests passed.
3. Retrieval eval:
   1. Text `Recall@5`: pre-fusion `0.5`, post-fusion `1.0`.
   2. SMILES `Recall@5`: pre-fusion `1.0`, post-fusion `1.0`.
4. Baseline integration:
   1. `test-rag --baseline`: passed.
   2. `test-api --baseline`: passed (models listed).

## Known operational note

`check-openwebui-models.sh` is not present in this repository. Use:

1. `./test-api.sh --baseline` to validate model listing and API path.
2. `./refresh-models.sh --wait` to refresh OpenWebUI in-memory model cache.
