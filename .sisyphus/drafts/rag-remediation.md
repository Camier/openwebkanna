# Draft: RAG Stack Remediation Plan

## Status: IN PROGRESS — awaiting external research agents

## What Exists Today

### services/rag/ — Honest Inventory

| File | Status | Reason |
|---|---|---|
| text_query_runtime.py | ✅ ALIVE | Official: transformers + fastembed |
| text_lane_qdrant.py | ✅ ALIVE | Official: qdrant_client SDK |
| vision_lane_qdrant.py | ✅ ALIVE | Official: qdrant_client SDK |
| vision_colmodernvbert.py | ✅ ALIVE | Official: transformers |
| colmodernvbert_profile.py | ✅ ALIVE | Pure config dataclass |
| fusion_policy.py | ✅ ALIVE (chem removed) | 50-line DIY RRF math |
| qdrant_schema.py | ✅ ALIVE | Pure config/setup |
| manifest.py | ⚠️ BROKEN | imports MaterializedRecord from DELETED materialize_evidence |
| embed_and_upsert.py | ⚠️ BROKEN | imports MaterializedRecord from DELETED materialize_evidence |
| query_router.py | ⚠️ BROKEN | imports normalize_smiles from DELETED file; ZERO production callers |
| materialize_evidence.py | DELETED | Was deleted in Phase 1 — causes breakage above |

### Broken Dependency Chain

```
manifest.py line 10:     from .materialize_evidence import MaterializedRecord
embed_and_upsert.py line 8: from .materialize_evidence import MaterializedRecord, make_qdrant_point_id
query_router.py line 19:    from .normalize_smiles import NormalizedSmiles, normalize_smiles

services/__init__.py (parent pkg): re-exports QueryRouter, ParsedQuery, QueryRouterConfig → broken
```

### Live Callers of services/rag/

| File | Imports | Status |
|---|---|---|
| `services/multimodal_retrieval_api/service.py` | ColModernVBERTConfig, probe_colmodernvbert_runtime | ✅ WORKS |
| `services/multimodal_retrieval_api/adapter.py` | text_lane + vision_lane classes | ✅ WORKS |
| `scripts/rag/materialize_rag_evidence.py` | 9 symbols from services.rag | ⚠️ BREAKS |
| `scripts/rag/index_colmodernvbert_figures.py` | materialize_evidence (DELETED) | ⚠️ BREAKS |
| `scripts/rag/migrate_legacy_text_points.py` | materialize_evidence (DELETED) | ⚠️ BREAKS |
| `scripts/rag/migrate_rag_evidence.py` | ZERO imports from services | ✅ SELF-CONTAINED |

### Scripts with Independent make_point_id

- `scripts/rag/migrate_rag_evidence.py`: defines local `make_point_id` (line 35) — NO dependency on deleted file
- `scripts/rag/migrate_legacy_text_points.py`: imports from DELETED `services.rag.materialize_evidence` — BREAKS
- `scripts/rag/index_colmodernvbert_figures.py`: imports from DELETED `services.rag.materialize_evidence` — BREAKS

### Chem Fields Still Present (not yet removed)

- `contracts.py`: `smiles`, `canonical_smiles`, `molecule_id`, `smiles_review_status` fields in `RetrievalEvidenceObject`
- `qdrant_schema.py`: `enable_chem_dense=True` by default, `chem_dense_dim=768`
- `qdrant_schema.py`: `has_smiles` payload index still created
- `service.py` lines 404-412: extracts SMILES fields from visual hits
- `text_lane_qdrant.py` line 63: outputs `molecule_id` in normalized hit

---

## Official Library Usage (confirmed)

### ✅ Using Official Libraries (No Change Needed)

1. **Text query encoding** (`text_query_runtime.py`):
   - `transformers.AutoModel.from_pretrained()` — correct
   - `fastembed.SparseTextEmbedding` — correct
   - Pure PyTorch tensor math — correct

2. **Qdrant hybrid retrieval** (`text_lane_qdrant.py`):
   - `qdrant_client.query_points()` with `Prefetch` + `FusionQuery.RRF` — correct Qdrant pattern
   - Qdrant named multivector API for vision — correct

3. **ColModernVBERT** (`vision_colmodernvbert.py`):
   - Official HuggingFace Transformers interface — correct
   - `ColModernVBertProcessor` + `ColModernVBertForRetrieval` — correct

### ❌ DIY Code (Needs Replacement)

1. **`query_router.py`** — 248 lines, broken import, ZERO callers
2. **`manifest.py`** — broken import chain, dataclass boilerplate
3. **`embed_and_upsert.py`** — broken import chain, batching logic
4. **`fusion_policy.py`** — 50 lines of RRF math, CHEM REMOVED but file still exists
5. **`materialize_evidence.py`** — DELETED, causing breakage in 2 files + manifest

---

## Phase 1: Fix Import Breakage (Immediate — No Research Needed)

### Step 1.1: Recreate `services/rag/materialize_evidence.py`

**What to recreate** (from script usage analysis):
- `PageInput`, `FigureInput`, `MoleculeInput` — dataclasses for input validation
- `SourceDocument` — dataclass for document metadata
- `MaterializedRecord` — dataclass wrapping payload + point_id + text_for_embedding
- `materialize_pages(pages, document)` → list[MaterializedRecord]
- `materialize_figures(figures, document)` → list[MaterializedRecord]
- `materialize_molecules(molecules, document)` → list[MaterializedRecord]
- `make_point_id(prefix, value)` → `f"{prefix}#{value}"`
- `make_qdrant_point_id(point_id)` → numeric hash or uuid

**What NOT to recreate** (chemistry-related, deleted):
- `normalize_smiles` — REMOVED
- `exact_chemistry_lane` — REMOVED
- `live_one_collection_router` — REMOVED

**File location**: `services/rag/materialize_evidence.py`

### Step 1.2: Fix `services/rag/manifest.py`

- Remove `from .materialize_evidence import MaterializedRecord`
- Import `MaterializedRecord` from `.materialize_evidence` (now recreated)
- Or: define a minimal `ManifestRecord` dataclass inline

### Step 1.3: Fix `services/rag/embed_and_upsert.py`

- Remove `from .materialize_evidence import MaterializedRecord, make_qdrant_point_id`
- Import `MaterializedRecord, make_qdrant_point_id` from `.materialize_evidence` (now recreated)

### Step 1.4: Delete `services/rag/query_router.py`

- Zero production callers — confirmed by grep
- Has broken `normalize_smiles` import
- `services/__init__.py` re-exports it — must remove from `__all__`

### Step 1.5: Rewrite `services/rag/__init__.py`

Remove broken/dead exports:
- Remove `QueryRouter`, `ParsedQuery`, `QueryRouterConfig` (query_router — DELETED)
- Remove `EmbeddingModelSpec`, `build_manifest`, `write_jsonl` (manifest — broken import chain, but rebuild manifest.py not delete)
- Actually: fix manifest.py first, then update __init__.py exports

### Step 1.6: Fix `scripts/rag/materialize_rag_evidence.py`

Update `from services.rag import` to match new exports after Step 1.1-1.3.

### Step 1.7: Fix `scripts/rag/index_colmodernvbert_figures.py`

Replace `from services.rag.materialize_evidence import make_point_id, make_qdrant_point_id` with local implementations — OR create a standalone utils module for these ID generation functions.

### Step 1.8: Fix `scripts/rag/migrate_legacy_text_points.py`

Replace `from services.rag.materialize_evidence import make_point_id, make_qdrant_point_id` with local implementations.

---

## Phase 2: Remove Remaining Chemistry Fields

### Step 2.1: Strip chemistry from `contracts.py`

`RetrievalEvidenceObject` still has:
```python
smiles: list[str]
canonical_smiles: list[str]
smiles_review_status: list[str]
smiles_backends: list[str]
smiles_confidences: list[float]
molecule_id: str | None
```

Remove all molecule/chemistry fields from `RetrievalEvidenceObject`.

### Step 2.2: Strip chemistry from `qdrant_schema.py`

- `enable_chem_dense=True` → `False`
- `chem_dense_dim=768` → remove
- `has_smiles` payload index → remove
- `canonical_smiles` payload index → remove (if molecule-specific)
- `inchikey` payload index → remove (if molecule-specific)

### Step 2.3: Strip SMILES from `service.py`

Lines 404-412: Remove SMILES field extraction from `_figure_evidence_from_visual_hit`.
`_select_qdrant_figures` SMILES field extractions.

### Step 2.4: Strip molecule_id from `text_lane_qdrant.py`

Line 63: Remove `molecule_id` from normalized hit output dict.

### Step 2.5: Update `migrate_rag_evidence.py`

Remove `--chem-dense-dim` arg, `chem_dense` vector handling, molecule record materialization.

---

## Phase 3: Simplify / Replace DIY Components

### Step 3.1: Evaluate `fusion_policy.py`

Current: 77 lines, weighted RRF for 2 lanes (text=1.0, visual=1.1).
After chem removal: 54 lines.

Options:
- **Keep minimal DIY**: ~10 lines of pure Python math — acceptable, no library needed
- **Replace with LlamaIndex**: Would require `llama-index-core` dependency for 10 lines of math

Decision: Keep minimal DIY. LlamaIndex would be heavier than the problem.

### Step 3.2: Evaluate `embed_and_upsert.py`

Current: 213 lines, batching pipeline around Qdrant client.
Protocol-based embedder interfaces (DenseTextEmbedder, etc.).

Options:
- **Keep as-is**: Works, batching logic is legitimate
- **Replace with LlamaIndex**: Would require adapting the pipeline to LlamaIndex's `IngestionPipeline`

Decision: Keep as-is. The batching logic is working code.

### Step 3.3: Evaluate `manifest.py`

Current: 120 lines, builds a manifest dataclass from materialized records.
Imports `MaterializedRecord` from deleted file.

After fixing: Just dataclass boilerplate — could be replaced with a plain dict or Pydantic model.

Decision: Rewrite to fix import, keep structure. The manifest format is a repo convention.

---

## Phase 4: Gap Implementations (from Phase 2 Research)

Pending external research agent results.

### Gap 1: Cross-encoder Reranking
- LlamaIndex: `CrossEncoderReranker(model_name="ms-marco-MiniLM-L-6-v2", top_n=30)`
- Insert between RRF output and final top-k
- ~30 lines integration

### Gap 2: HyDE Selective Activation
- LlamaIndex: `HyDEQueryTransform(llm=..., include_original=True)`
- Trigger heuristic: query length < 4 tokens OR has rare chemical name
- ~15 lines conditional

### Gap 3: RAGAS Evaluation
- Install `ragas`
- `Faithfulness` metric: requires `user_input`, `response`, `retrieved_contexts`
- Pipeline produces all three — bridge to RAGAS evaluate()
- ~50 lines

### Gap 4: Chemistry Lane (Future)
- NOT implementing until exact-match path is stable in notebook
- Plan: exact canonical SMILES via RDKit + Qdrant payload filter

---

## Research Findings (Pending)

### From LlamaIndex docs (bg_573c7ac1) — PENDING

### From Qdrant docs (bg_542ece7a) — PENDING

### From RAGAS docs (bg_f02259fa) — PENDING

### From GitHub examples (bg_a7ffdbdd) — PENDING

### From LlamaIndex HyDE docs (bg_7698fd11) — PENDING

### From GitHub ColModernVBERT (bg_d703a26a) — PENDING

---

## Files to Modify / Delete

### DELETE
- `services/rag/query_router.py` (dead, broken)

### REWRITE (new content)
- `services/rag/materialize_evidence.py` (recreate deleted file)
- `services/rag/__init__.py` (remove dead exports)

### REWRITE (fix broken imports)
- `services/rag/manifest.py`
- `services/rag/embed_and_upsert.py`

### EDIT (strip chem fields)
- `services/multimodal_retrieval_api/contracts.py`
- `services/rag/qdrant_schema.py`
- `services/multimodal_retrieval_api/service.py`
- `services/rag/text_lane_qdrant.py`

### EDIT (fix script imports)
- `scripts/rag/materialize_rag_evidence.py`
- `scripts/rag/index_colmodernvbert_figures.py`
- `scripts/rag/migrate_legacy_text_points.py`
- `scripts/rag/migrate_rag_evidence.py`

---

## Not Using LlamaIndex For

| Component | Reason |
|---|---|
| Qdrant hybrid search | Raw `qdrant_client` with `Prefetch` + `FusionQuery.RRF` is the correct, fully-documented pattern. LlamaIndex's `QdrantVectorStore` adds indirection. |
| ColModernVBERT late interaction | Custom model (not in LlamaIndex model registry), must use raw transformers + Qdrant multivector API |
| Fusion policy | 10 lines of math. LlamaIndex would be heavier than the problem. |
| Data model schemas | `MaterializedRecord`, `PageInput`, etc. are repo-specific document shapes, not library concerns |

## Using LlamaIndex For

| Gap | LlamaIndex Component | Integration Points |
|---|---|---|
| Reranking | `CrossEncoderReranker` | service.py retrieve() — insert between text_lane() and evidence building |
| HyDE | `HyDEQueryTransform` | Conditional on query length / type heuristic |
| Evaluation | `ragas` + `llama_index.evaluation` | scripts/eval/run_retrieval_eval.py extension |
| BM25 sparse | `BM25Retriever` | Could replace fastembed dependency (optional) |
