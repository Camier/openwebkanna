# Multimodal RAG Redesign Plan

Last updated: 2026-03-14 (UTC)
Status: proposed target architecture
Scope: replace the current script-orchestrated multimodal path with a native retrieval backend built around the existing host Qdrant collection

## Purpose

This document defines the target architecture for multimodal retrieval in this repository.

The goal is not to improve the current wrapper-based prototype. The goal is to replace it with a proper backend where text, pages, figures, chemical blocks, and `SMILES` are first-class retrieval objects.

## Decision Summary

We will treat the existing host Qdrant collection `pdf_nemotron_hybrid` as the page-level retrieval backbone and build the multimodal system around it.

We will not keep the local CLIP-backed figure index and Python answer orchestration as the long-term architecture.

We will add a dedicated multimodal retrieval service over Qdrant, extend the indexed object model to include figures and chemical blocks, and promote `SMILES` to a native retrieval and evidence field.

OpenWebUI remains the UI and chat client. It should not own the multimodal fusion logic.

## Constraints

- Current local OpenWebUI runtime is still configured for `pgvector` in `.env`.
- OpenWebUI can connect to Qdrant through upstream-supported env vars, but that alone does not provide the multimodal retrieval behavior we need.
- A real host Qdrant deployment already exists and already serves a Nemotron-backed hybrid collection.
- The redesign must prefer official capabilities from Qdrant and Hugging Face models over new local glue layers.

## Non-goals

- Do not keep extending `render-multimodal-answer*.py` as the product path.
- Do not build another local side index as a second source of truth.
- Do not force all retrieval into OpenWebUI's built-in single-vector abstraction.
- Do not assume ColQwen2 is mandatory when a Nemotron/Qdrant collection already exists and works.

## Current State

### Local OpenWebUI facts

- `.env` still sets `VECTOR_DB=pgvector`.
- The current multimodal path is orchestrated by local scripts under `scripts/rag/`.
- The current "multimodal" answer path depends on local orchestration files:
  - `scripts/rag/render_multimodal_answer.py`
  - `scripts/rag/render_multimodal_answer_v2.py`
  - `scripts/rag/multimodal_index.py`
  - `scripts/rag/build-multimodal-index.sh`
  - `scripts/rag/render-multimodal-answer.sh`
  - `scripts/rag/render-multimodal-answer-v2.sh`

This is acceptable as a prototype. It is not the target product architecture.

### Existing host Qdrant backbone

Observed on 2026-03-14 from the existing host deployment:

- collection: `pdf_nemotron_hybrid`
- status: `green`
- points: `4748`
- indexed vectors: `14200`
- named dense vector: `dense_prefetch` (`2560`, cosine)
- named multivector: `late_original` (`2560`, cosine, `MAX_SIM`)
- sparse vector: `sparse_text` with `IDF`
- resolved text field: `document_title`
- resolved document id field: `document_id`
- resolved page field: `page_number`

This means the current host collection already implements the correct broad retrieval pattern:

- dense first-pass retrieval
- sparse companion retrieval
- fusion
- late-interaction reranking

The current gap is not "we have no multimodal search engine". The gap is that the indexed object model is page-level and metadata-level, not figure-level and chemistry-level.

### Existing Nemotron ingestion and retrieval facts

The current external ingestion/retrieval stack already encodes PDF pages into Qdrant with:

- `dense_prefetch` for candidate retrieval
- `late_original` for late-interaction reranking
- payload metadata such as `document_id`, `page_number`, `document_title`, `document_keywords`, `document_entities`, `study_tags`, `compound_mentions`

What it does not currently provide as first-class indexed fields:

- figure ids
- chemical block ids
- image region provenance
- `SMILES`
- canonical `SMILES`
- OCSR provenance and confidence

## Why The Current Prototype Must Be Retired

The current prototype breaks the system into separate local concerns:

1. OpenWebUI retrieval
2. a local figure index
3. local score fusion
4. local answer rendering

That creates the wrong ownership boundaries:

- retrieval truth is split across systems
- figures are attached after retrieval instead of being indexed objects
- `SMILES` can appear in answers without being native retrieval signals
- OpenWebUI is forced to depend on local orchestration instead of a retrieval service contract

The redesign must remove that split.

## Reference Sources

The target architecture is based on official docs and model references.

### Qdrant

- Hybrid queries, nested prefetch, and fusion:
  - https://qdrant.tech/documentation/concepts/hybrid-queries/
- Multivectors and late-interaction retrieval:
  - https://qdrant.tech/documentation/tutorials-search-engineering/using-multivector-representations/

These official docs are the key reason to reuse Qdrant as the retrieval backbone instead of building another local index manager.

### Hugging Face

- ColQwen2 model card:
  - https://huggingface.co/vidore/colqwen2-v1.0
  - https://huggingface.co/vidore/colqwen2-v1.0-hf
- Hugging Face task page for multimodal document/image-text pipelines:
  - https://huggingface.co/tasks/image-text-to-text
- CLIP model card:
  - https://huggingface.co/openai/clip-vit-base-patch32

Interpretation:

- ColQwen2 remains a strong reference design for document-visual retrieval and reranking.
- CLIP is acceptable for prototyping but should not define the long-term retrieval backbone here.
- Because a working Nemotron/Qdrant collection already exists, ColQwen2 is a later option, not a prerequisite for the redesign.

### Hugging Face MCP note

Hugging Face MCP access was attempted during this design pass but was unavailable because the MCP endpoints required authentication in this session. The Hugging Face evidence used here therefore comes from public official Hugging Face pages and model cards rather than direct MCP document fetches.

## Target Architecture

### 1. Canonical Retrieval Substrate

Qdrant becomes the retrieval source of truth.

Use:

- the existing `pdf_nemotron_hybrid` collection for page-level retrieval
- one additional figure/chemical collection for image-region retrieval and chemistry enrichment

Recommended split:

- `pdf_nemotron_hybrid`
  - unit: page
  - role: document/page candidate generation and late rerank
- `pdf_nemotron_figures`
  - unit: figure or chemical block crop
  - role: figure retrieval, chemistry evidence, and image-grounded evidence objects

If Qdrant payload size or collection ownership makes a second collection undesirable, the alternative is a single mixed collection with strong type fields and named vectors. The default recommendation remains two collections because the object granularity is materially different.

### 2. Canonical Object Model

The retrieval system must promote the following object types to first-class records:

- `document`
- `page`
- `figure`
- `chemical_block`
- `molecule_evidence`

Minimum payload for figure and chemical records:

- `document_id`
- `source_pdf_sha256`
- `page_number`
- `block_id`
- `object_type`
- `caption_text`
- `surrounding_text`
- `figure_kind`
- `image_uri` or `image_storage_key`
- `has_smiles`
- `smiles`
- `canonical_smiles`
- `smiles_source_backend`
- `smiles_confidence`
- `smiles_review_status`

### 3. Retrieval API

Add a dedicated multimodal retrieval service in front of Qdrant.

Inputs:

- text query
- image query
- `SMILES` query
- mixed query

Responsibilities:

- retrieve page candidates from `pdf_nemotron_hybrid`
- retrieve figure or chemical candidates from the figure collection
- apply metadata filters and joins by `document_id`, `page_number`, and `block_id`
- apply exact and canonical `SMILES` matching when present
- fuse page, figure, and chemistry evidence into a typed evidence set

Output contract:

- typed evidence objects, not free-form merged blobs
- per-object provenance
- per-object modality
- per-object retrieval scores
- optional answer-context serialization for the downstream generator

### 4. Answer API

The answer stage should consume typed evidence objects and produce:

- citations
- figure references
- image attachments
- `SMILES` where available
- explicit provenance back to page and block ids

This answer stage should live behind a service contract, not behind `render_multimodal_answer*.py`.

### 5. OpenWebUI Role

OpenWebUI becomes a client of the multimodal retrieval and answer APIs.

OpenWebUI can still use upstream-supported Qdrant configuration when useful, but it should not be the place where:

- page retrieval
- figure retrieval
- `SMILES` retrieval
- score fusion
- provenance stitching

are all recombined ad hoc.

## Retrieval Behavior By Query Type

### Text query

1. Query `pdf_nemotron_hybrid` using its existing hybrid pattern.
2. Use document and page hits to constrain figure/chemical lookup.
3. Retrieve figure or chemical evidence from the figure collection.
4. Fuse results into a single evidence set.

### Image query

1. Query the figure collection directly with the visual retrieval path.
2. Join matched figure records back to pages and documents.
3. Return figures as primary evidence and pages as support evidence.

### `SMILES` query

1. Canonicalize the query `SMILES`.
2. Exact-match `canonical_smiles` in the figure/chemical collection.
3. Fall back to raw `smiles` and metadata aliases if needed.
4. Join back to page and document evidence.

`SMILES` must be treated as both:

- a precise structured filter
- a retrieval and answer evidence field

## Migration Plan

### Phase 0: Freeze the prototype boundary

Treat these files as transitional only:

- `scripts/rag/render_multimodal_answer.py`
- `scripts/rag/render_multimodal_answer_v2.py`
- `scripts/rag/multimodal_index.py`
- `scripts/rag/build-multimodal-index.sh`
- `scripts/rag/render-multimodal-answer.sh`
- `scripts/rag/render-multimodal-answer-v2.sh`

No new product behavior should depend on them.

### Phase 1: Qdrant introspection and contract

Build a first-class local contract for the existing host collection:

- collection profile
- vector names
- sparse vector names
- payload field contract
- join keys

This contract should be versioned and testable.

### Phase 2: Figure and chemical indexing

Add figure and chemical objects as first-class Qdrant records.

Requirements:

- deterministic `block_id`
- image region persistence
- caption and local text payload
- optional figure embedding path aligned with the host retrieval stack

### Phase 3: Chemistry enrichment

Persist OCSR output as indexed payload, not as sidecar-only JSON.

Required chemistry fields:

- raw `smiles`
- canonical `smiles`
- backend name
- confidence
- review status
- source crop provenance

### Phase 4: Retrieval service

Introduce a dedicated retrieval API that:

- queries Qdrant directly
- understands page and figure records
- understands `SMILES`
- returns typed evidence objects

### Phase 5: Answer service and OpenWebUI integration

Move answer assembly behind a stable service contract and make OpenWebUI consume it.

At this point the current local multimodal render scripts become migration tools only.

## Exit Criteria

The redesign is done when all of the following are true:

- page retrieval uses the existing host Qdrant collection directly
- figure and chemical evidence are indexed in Qdrant as first-class records
- `SMILES` are retrievable as structured evidence
- OpenWebUI no longer depends on `render_multimodal_answer*.py` for product behavior
- the answer path consumes typed multimodal evidence instead of stitching local artifacts

## Risks

- The current Nemotron collection is page-level, so figure-level indexing is a real schema and ingestion change, not a config toggle.
- If the host Nemotron stack cannot encode figure crops in a way that matches retrieval quality needs, a second visual encoder family may still be needed later.
- Moving `SMILES` into indexed payload without a review status would contaminate retrieval precision.

## Immediate Next Step

Do not add more wrapper logic.

The next implementation artifact should be a contract document and service skeleton for:

- Qdrant-backed page retrieval
- Qdrant-backed figure retrieval
- `SMILES`-aware evidence objects
- OpenWebUI client integration against that service
