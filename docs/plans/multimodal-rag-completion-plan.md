# Multimodal RAG Completion Plan

## Current Status (Phase 0 Complete)

| Component | Status |
|-----------|--------|
| Collection `rag_evidence` | ✅ Recreated with full schema |
| `text_dense` (2560-dim) | 🟡 Empty |
| `vision_li` (128-dim MV) | 🟡 Empty |
| `chem_dense` (768-dim) | 🟡 Empty |
| `text_sparse` (BM25) | 🟡 Empty |

---

## Actionable Plan

### Phase 1: Text Migration (30 min)
- [ ] 1.1 Migrate 4,704 pages from `pdf_nemotron_hybrid` → `rag_evidence`
- [ ] 1.2 Verify text vectors indexed

### Phase 2: Figure Indexing (45 min)
- [ ] 2.1 Download ColModernVBERT model
- [ ] 2.2 Run figure indexing (1,505 figures → vision_li vectors)
- [ ] 2.3 Verify figure vectors indexed

### Phase 3: Molecule Pipeline (45 min)
- [ ] 3.1 Decide: RDKit fingerprints vs neural (ChemBERTa)
- [ ] 3.2 Generate molecule embeddings
- [ ] 3.3 Upsert to `chem_dense` vector space

### Phase 4: Payload & Validation (15 min)
- [ ] 4.1 Backfill document metadata
- [ ] 4.2 Validate collection completeness
- [ ] 4.3 Test retrieval queries

---

## Commands

### Phase 1: Text Migration
```bash
cd /LAB/@thesis/openwebui
python scripts/rag/migrate_legacy_text_points.py \
    --source-collection pdf_nemotron_hybrid \
    --target-collection rag_evidence \
    --source-dense-vector dense_prefetch \
    --source-sparse-vector sparse_text \
    --target-dense-vector text_dense \
    --target-sparse-vector text_sparse
```

### Phase 2: Figure Indexing
```bash
# Download model first
# Then run indexing
python scripts/rag/index_colmodernvbert_figures.py \
    --figures-jsonl artifacts/rag/materialized/generated.figures.jsonl \
    --collection-name rag_evidence \
    --text-dense-dim 2560 \
    --chem-dense-dim 768 \
    --device cuda
```

### Phase 3: Molecule Upsert
```bash
python scripts/rag/upsert_resolved_molecule_points.py \
    --resolved-records artifacts/rag/smiles_identity/ambiguous_resolved_external_records.parquet \
    --collection-name rag_evidence \
    --apply
```

---

## Dependencies

- [x] Qdrant (running on localhost:6333)
- [x] GPU (RTX 5000 16GB)
- [ ] ColModernVBERT model
- [ ] Python env with dependencies

---

## Target State

After completion:
- 4,704 page points with text_dense + text_sparse
- 1,505 figure points with vision_li multivectors
- 55+ molecule points with chem_dense
- Full payload: document_id, page_id, figure_id, molecule_id, SMILES, InChIKey
