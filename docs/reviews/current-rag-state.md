# Current RAG State

Last updated: 2026-03-14 (UTC)

This note captures the current repository state for RAG and multimodal retrieval.
It is a local repo summary, not an upstream reference snapshot and not a live
runtime status record.

Use this file when you want the shortest current-state answer before drilling
into the canonical sources.

## Source Of Truth

Verify deployment-sensitive details in this order:

- `README.md`
- `docs/ssot/stack.md`
- `docs/REPO_MAP.md`
- `config/README.md`
- `config/env/.env.example`
- `config/compose/docker-compose.yml`

## Current Baseline

- The committed OpenWebUI retrieval baseline is still text-first.
- `RAG_EMBEDDING_MODEL=pritamdeka/S-PubMedBert-MS-MARCO` is the main embedding
  control in the committed env template.
- `VECTOR_DB=pgvector` is the committed local baseline.
- Chunking is configured in characters, not tokens, with a committed baseline of
  `CHUNK_SIZE=3000` and `CHUNK_OVERLAP=600`.
- Hybrid search uses the correct local flag name
  `ENABLE_RAG_HYBRID_SEARCH`.

## Multimodal State

- The repo now contains a real multimodal prototype path rather than only OCR or
  extraction helpers.
- A CLIP-backed figure index exists under `scripts/rag/` and can fuse text and
  image similarity for retrieval-oriented answer generation.
- A tracked retrieval service exists under `services/multimodal_retrieval_api/`
  and delegates retrieval to the canonical host `thesis_graph` Qdrant runtime in
  `/LAB/@thesis/wow`.
- The long-term target is still a Qdrant/Nemotron-centered multimodal design as
  documented in `docs/plans/MULTIMODAL_RAG_QDRANT_NEMOTRON_REDESIGN.md`.

## Extraction Backend

- `Docling` is the active canonical extraction backend in this repo's baseline.
- Tika appears in upstream/reference material, but it is not part of the
  default local baseline described by the canonical config and SSOT docs.

## Web Search Naming

- The repo intentionally carries both alias forms for web-search flags and
  engine selection:
  - `ENABLE_WEB_SEARCH` and `ENABLE_WEBSEARCH`
  - `WEB_SEARCH_ENGINE` and `WEBSEARCH_ENGINE`
- This is a compatibility surface, not necessarily a bug, but it is still a
  duplicated naming surface that can confuse readers.

## What Is No Longer Accurate

These older claims should not be repeated as current local truth:

- "The repo has no true vision-language retrieval path."
- "Tika is the active default extraction backend here."
- "The local config still uses `RAG_HYBRID_SEARCH`."

## Related Docs

- Historical upstream-review artifact: `docs/reviews/rag_technical_reference_review.md`
- Upstream snapshots: `docs/reference/openwebui/README.md`
- Target multimodal redesign: `docs/plans/MULTIMODAL_RAG_QDRANT_NEMOTRON_REDESIGN.md`
- Live machine-specific observations: `docs/status/`
