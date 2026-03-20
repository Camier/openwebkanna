# RAG Stack Remediation Plan

## TL;DR

> Fix broken import chain from deleted chemistry files. Strip remaining chemistry fields from contracts and schema. Remove dead `query_router.py`. Recreate deleted `materialize_evidence.py` dataclasses. Then add cross-encoder reranking via LlamaIndex `CrossEncoderReranker`.

**Deliverables:**
- `services/rag/` imports cleanly (no broken references)
- Chemistry payload stripped from contracts, schema, and service
- Dead code (`query_router.py`) removed
- `scripts/rag/` imports fixed
- Cross-encoder reranking gap bridged via LlamaIndex

**Estimated Effort:** Short (import fixes) + Short (chem strip) + Medium (reranking gap)
**Parallel Execution:** YES — 3 waves
**Critical Path:** Fix imports → Strip chem → Add reranking

---

## Context

### Problem 1: Import Chain Is Broken

Three files in `services/rag/` import from the DELETED `materialize_evidence.py`:

```
manifest.py line 10:    from .materialize_evidence import MaterializedRecord
embed_and_upsert.py line 8: from .materialize_evidence import MaterializedRecord, make_qdrant_point_id
```

`query_router.py` line 19: `from .normalize_smiles import NormalizedSmiles` — normalize_smiles was deleted.

`services/__init__.py` re-exports `QueryRouter`, `ParsedQuery`, `QueryRouterConfig` from the broken `query_router.py`.

This cascades: `from services.rag import X` (for ANY X) triggers `services/rag/__init__.py` → imports `QueryRouter` from broken `query_router.py` → **ImportError**.

### Problem 2: `QueryRouter` Has ZERO Production Callers

Grep across entire codebase: `QueryRouter` appears only in `services/__init__.py` (re-export), `services/rag/__init__.py` (re-export), and `services/rag/query_router.py` (definition). Not used by `service.py`, `adapter.py`, or any script. Dead code with a broken import.

### Problem 3: Chemistry Fields Still Present After Deletion

`contracts.py`: `smiles`, `canonical_smiles`, `smiles_review_status`, `smiles_backends`, `smiles_confidences`, `molecule_id` — all still in `RetrievalEvidenceObject`.

`qdrant_schema.py`: `enable_chem_dense=True` by default, `chem_dense_dim=768`, `has_smiles` payload index.

`service.py` lines 404-412: extracts SMILES fields from visual hit payloads.

`text_lane_qdrant.py` line 63: outputs `molecule_id` orphan field.

### What Is Working Correctly

| Component | Verdict |
|---|---|
| `text_query_runtime.py` | ✅ Official: `transformers.AutoModel` + `fastembed.SparseTextEmbedding` |
| `text_lane_qdrant.py` | ✅ Official: `qdrant_client` with `Prefetch` + `FusionQuery.RRF` — correct Qdrant pattern |
| `vision_lane_qdrant.py` | ✅ Official: `qdrant_client` multivector query API |
| `vision_colmodernvbert.py` | ✅ Official: HuggingFace Transformers interface |
| `fusion_policy.py` | ✅ Alive (chem removed) — 10-line DIY RRF is fine |

### LlamaIndex Usage Strategy

**Use LlamaIndex for:**
- Cross-encoder reranking (`CrossEncoderReranker` — replaces DIY reranking gap)
- HyDE query transform (optional gap, conditional)
- RAGAS evaluation bridge (optional gap)

**Do NOT use LlamaIndex for:**
- Qdrant hybrid search — raw `qdrant_client` SDK with `Prefetch` + `FusionQuery.RRF` is the canonical Qdrant pattern, fully documented, battle-tested. LlamaIndex's `QdrantVectorStore` adds indirection.
- ColModernVBERT late interaction — custom model, must use raw transformers + Qdrant multivector API.
- Fusion policy — 10 lines of math, no library needed.

---

## Work Objectives

### Concrete Deliverables

1. `services/rag/materialize_evidence.py` — recreated with surviving dataclasses
2. `services/rag/__init__.py` — clean exports, no broken imports
3. `services/rag/manifest.py` — fixed import
4. `services/rag/embed_and_upsert.py` — fixed import
5. `services/rag/query_router.py` — DELETED (dead code)
6. `contracts.py` — chemistry fields stripped from `RetrievalEvidenceObject`
7. `qdrant_schema.py` — `enable_chem_dense=False`, payload indexes cleaned
8. `service.py` — SMILES field extraction removed from `_figure_evidence_from_visual_hit`
9. `text_lane_qdrant.py` — `molecule_id` orphan field removed
10. Scripts fixed: `materialize_rag_evidence.py`, `index_colmodernvbert_figures.py`, `migrate_legacy_text_points.py`, `migrate_rag_evidence.py`
11. Cross-encoder reranking integration via LlamaIndex

### Definition of Done

- [ ] `python -c "from services.rag import QdrantTextHybridLane"` → succeeds
- [ ] `python -c "from services.rag import FusionPolicy"` → succeeds
- [ ] `python -c "from services.rag import manifest"` → succeeds
- [ ] `python -c "from services.rag import embed_and_upsert"` → succeeds
- [ ] `check-doc-consistency.sh` → all OK
- [ ] `contracts.py` — `RetrievalEvidenceObject` has zero chemistry fields
- [ ] `qdrant_schema.py` — `chem_dense` vector removed, `enable_chem_dense=False`
- [ ] `QueryRouter` absent from codebase

### Must Have

- All imports resolve without errors
- No chemistry fields in evidence objects
- Chem vector disabled in Qdrant schema
- `QueryRouter` dead code removed

### Must NOT Have

- No LlamaIndex as a dependency for core retrieval
- No new chemistry code
- No mock/fallback paths

---

## Verification Strategy

### Test Decision

- **Infrastructure exists**: YES — Python import chain, Qdrant client, existing test scripts
- **Automated tests**: NO — existing test harness uses `--baseline` scripts
- **Agent-Executed QA Scenarios**: YES — every task verified by running the affected import/command

### QA Policy

Every task includes agent-executed QA — no human confirmation needed.

**Import smoke test**: `python -c "from services.rag import X"` per modified module
**Doc consistency**: `./scripts/check-doc-consistency.sh`
**Qdrant schema**: verify `enable_chem_dense=False` in source

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Foundation — can all run in parallel):
├── T1: Recreate services/rag/materialize_evidence.py
├── T2: Rewrite manifest.py (fix import)
├── T3: Rewrite embed_and_upsert.py (fix import)
├── T4: Rewrite services/rag/__init__.py (clean exports)
└── T5: DELETE services/rag/query_router.py

Wave 2 (Chem strip + scripts — after Wave 1 completes):
├── T6: Strip chemistry from contracts.py
├── T7: Strip chemistry from qdrant_schema.py
├── T8: Strip SMILES from service.py
├── T9: Strip molecule_id from text_lane_qdrant.py
├── T10: Fix scripts (materialize_rag_evidence.py, index_*.py, migrate_*.py)
└── T11: Strip migrate_rag_evidence.py chem handling

Wave 3 (Gap implementation — after Wave 1 completes):
├── T12: Add LlamaIndex CrossEncoderReranker to service.py retrieve()
└── T13: Update eval script with RAGAS faithfulness metric

Wave FINAL (Verification — after ALL waves):
├── T14: Full import smoke test
├── T15: check-doc-consistency.sh
└── T16: Scope fidelity check
```

**Critical Path:** Wave 1 (T1-T5) unblocks everything. T6-T11 can start after T1 completes.
**Parallel Speedup:** ~60% faster than sequential

---

## TODOs

---

## Final Verification Wave

- [ ] F1. **Import smoke test** — `oracle`: Run `python -c "from services.rag import QdrantTextHybridLane, FusionPolicy, manifest, embed_and_upsert"` — must succeed with zero ImportError.
- [ ] F2. **Check-doc-consistency** — `./scripts/check-doc-consistency.sh` — all checks pass.
- [ ] F3. **Chem fields audit** — `grep -n "smiles\|molecule_id\|chem_dense" services/multimodal_retrieval_api/contracts.py services/rag/qdrant_schema.py services/multimodal_retrieval_api/service.py` — zero matches.
- [ ] F4. **QueryRouter audit** — `grep -rn "QueryRouter\|query_router" services/` — zero matches.
- [ ] F5. **Reranking integration** — verify `CrossEncoderReranker` is wired into `service.py retrieve()` flow.

---

## Commit Strategy

- 1 commit for all Wave 1 changes
- 1 commit for all Wave 2 changes
- 1 commit for Wave 3 (reranking)
- Message: `refactor(rag): fix import chain, strip chemistry, add reranking`

---

## Success Criteria

- `from services.rag import QdrantTextHybridLane` → no ImportError
- `RetrievalEvidenceObject` → no chemistry fields (`smiles`, `canonical_smiles`, etc.)
- `qdrant_schema.py` → `enable_chem_dense=False`
- `QueryRouter` → deleted from codebase
- Cross-encoder reranking → integrated into `service.py`
