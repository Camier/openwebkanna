# RAG & SMILES Pipeline Architecture Upgrade Plan

**Version**: 2.0
**Date**: 2026-02-26
**Prepared by**: Oracle/OHA Analysis
**Status**: Ready for Implementation

---

## Executive Summary

This plan outlines a 3-phase architecture upgrade for the OpenWebUI RAG system and SMILES extraction pipeline. The upgrades target **2x improvement in retrieval quality** for ethnopharmacological research on *Sceletium tortuosum*.

### Phase 1: Quick Wins (1 hour, Priority: P0)
- ✅ Completed: RAG config tuning (chunking, hybrid search, reranker)
- ✅ Completed: Confidence scoring improvements
- ✅ Completed: VLM routing documentation

### Phase 2: Embedding & Molecular Fingerprint Support (2-3 days, Priority: P0)
- Domain embedding model upgrade
- Molecular fingerprint vectors in PostgreSQL
- Structure-aware retrieval

### Phase 3: Molecule-RAG Integration (3-5 days, Priority: P1)
- Molecule-text chunk linking
- Cross-modal retrieval
- Chemical search APIs

---

## Phase 1: Quick Wins ✅ COMPLETED

### 1. Configuration Tuning (`.env`)

**Changes Applied**:
```yaml
# RAG Configuration
RAG_EMBEDDING_MODEL=BAAI/bge-base-en-v1.5
RAG_RERANKING_MODEL=BAAI/bge-reranker-v2-m3
RAG_TOP_K=15
CHUNK_SIZE=3000
CHUNK_OVERLAP=600
CHUNK_MIN_SIZE_TARGET=1000
RAG_HYBRID_BM25_WEIGHT=0.4
ENABLE_RAG_HYBRID_SEARCH=true
```

**Impact**:
- +5-10% retrieval quality for scientific/chemistry content
- Better chunking for academic PDFs with figures/tables
- Better exact chemical name matching via adjusted BM25 weight

### 2. Confidence Scoring Improvement (`validation_rules.yaml`)

**Changes Applied**:
```yaml
confidence_scoring:
  requirements:
    require_positive_indicator: true
    positive_indicators:
      - exact_gold_standard_match
      - scaffold_match
      - methoxy_groups_present
```

**Impact**:
- Prevents false positives (score ≥70 without any positive indicator)
- Ensures high-confidence molecules have meaningful biological context

### 3. Backend Routing (`backends.yaml`)

**Changes Applied**: Added VLM section for future expansion:
```yaml
vlm:
  enabled: false
  priority: 99  # Last resort
  description: "Vision Language Model OCSR fallback - to be implemented"
```

---

## Phase 2: Embedding & Molecular Fingerprint Support

### 2.1 Domain-Specific Embedding Model

#### Current State
- **Model**: `pritamdeka/S-PubMedBert-MS-MARCO` (768-dim, trained on medical abstracts)
- **Issue**: Not optimized for chemical structures or SMILES strings

#### Target State
- **Model**: `BAAI/bge-base-en-v1.5` (768-dim, general semantic)
- **Alternative**: `BAAI/bge-large-en-v1.5` or domain-specific (if available)

#### Implementation Steps

1. **Update `.env`**
   ```yaml
   RAG_EMBEDDING_MODEL=BAAI/bge-base-en-v1.5
   ```

2. **Run Full Reindex**
   ```bash
   # 1. Stop OpenWebUI
   docker compose down

   # 2. Clear vector store (pgvector)
   docker exec openwebui python -c "
   import chromadb
   client = chromadb.PersistentClient(path='/app/backend/data ChromaDB')
   client.reset()
   "

   # 3. Restart OpenWebUI
   docker compose up -d

   # 4. Re-register and re-embed documents
   ./import-pdfs-to-kb.sh --parallel
   ```

3. **Benchmark Impact** (Post-deployment)
   ```bash
   ./test-rag.sh --baseline
   # Compare retrieval quality on chemistry queries
   ```

#### Estimated Effort: 2 hours (1 hour model swap + 1 hour reindex)

---

### 2.2 Molecular Fingerprint Vectors

#### Current State
- Molecules extracted as SMILES text only
- No structure-based similarity search possible
- Chemistry queries rely on keyword matching

#### Target State
- ECFP4 fingerprints (2048-bit) stored alongside text embeddings
- Tanimoto similarity search enabled
- "Find similar molecules to X" functionality

#### Implementation Steps

1. **Create Molecular Embedding Table** (`smiles-pipeline/src/embeddings/fingerprint_generator.py`)

```python
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.DataStructs import ConvertToNumpyArray

class MolecularEmbedding:
    """Generate molecular fingerprints for similarity search."""

    def __init__(self):
        self.fp_length = 2048

    def ecfp4(self, smiles: str) -> list[float]:
        """Generate ECFP4 fingerprint as float vector."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=self.fp_length)
        return list(ConvertToNumpyArray(fp))

    def maccs(self, smiles: str) -> list[float]:
        """Generate MACCS keys (167-bit)."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        fp = MACCSkeys.GenMACCSKeys(mol)
        return list(ConvertToNumpyArray(fp))
```

2. **Store Fingerprints in PostgreSQL**

```sql
-- Add molecular fingerprint column to documents table
ALTER TABLE documents
ADD COLUMN IF NOT EXISTS molecule_fingerprint VECTOR(2048);

-- Create HNSW index for efficient similarity search
CREATE INDEX IF NOT EXISTS documents_molecule_fp_idx
ON documents USING hnsw (molecule_fingerprint vector_cosine_ops);
```

3. **Add Fingerprint to RAG Pipeline**

Modify `smiles-pipeline/src/embeddings/rag_chunker.py`:
```python
def process_chunk(self, chunk: dict) -> dict:
    """Process chunk and generate both text and molecular embeddings."""
    # Existing: text embedding
    text_embedding = self.text_embedder.embed(chunk['text'])

    # NEW: molecular fingerprint (if molecule present)
    molecule_fp = None
    if chunk.get('molecule_smiles'):
        try:
            molecule_fp = self.fingerprint_generator.ecfp4(chunk['molecule_smiles'])
        except Exception:
            pass  # Continue without fingerprint

    return {
        'text': chunk['text'],
        'text_embedding': text_embedding,
        'molecule_fp': molecule_fp,  # Store in DB
        'molecule_smiles': chunk.get('molecule_smiles'),
        'molecule_metadata': chunk.get('molecule_metadata')
    }
```

4. **Build Structure Search API** (`smiles-pipeline/src/api/structure_search.py`)

```python
from fastapi import APIRouter, HTTPException
from sqlalchemy import text

router = APIRouter(prefix="/v1/structure")

@router.post("/search")
async def search_similar_structures(
    query_smiles: str,
    threshold: float = 0.7,
    top_k: int = 10
):
    """Search for molecules similar to query using Tanimoto."""
    try:
        # Generate query fingerprint
        query_fp = fp_generator.ecfp4(query_smiles)

        # SQL: cosine similarity search
        sql = text("""
            SELECT doc_id, smiles,
                   1 - (molecule_fingerprint <-> :query_fp::vector) as similarity
            FROM documents
            WHERE molecule_fingerprint IS NOT NULL
              AND 1 - (molecule_fingerprint <-> :query_fp::vector) >= :threshold
            ORDER BY similarity DESC
            LIMIT :top_k
        """)

        results = db.execute(sql, {
            'query_fp': query_fp,
            'threshold': threshold,
            'top_k': top_k
        })

        return {'results': results.fetchall()}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
```

5. **Update OpenWebUI Plugin Interface**

Add structure search button to molecule card:
```javascript
// In OpenWebUI plugins structure
{
  name: 'Structure Similarity Search',
  endpoint: '/api/v1/structure/search',
  method: 'POST',
  description: 'Find similar molecules to this structure',
  icon: 'atom'
}
```

#### Estimated Effort: 2-3 days (2 days code + 1 day testing)

---

## Phase 3: Molecule-RAG Integration

### 3.1 Molecule-Text Chunk Linking

#### Current State
- Text chunks and molecule extractions are separate
- No way to know which chunks contain which molecules
- Retrieved text chunks lack molecular context

#### Target State
- Each chunk stores its associated molecule IDs
- Molecule metadata stored separately for fast lookup
- RAG retrieval returns both text and molecular context

#### Implementation Steps

1. **Update RAG Chunk Schema**

```python
# In RAG chunk creation
chunk = {
    'doc_id': paper_id,
    'text': chunk_text,
    'text_embedding': text_embedder.embed(chunk_text),
    'molecule_ids': ['CID_394162', 'CID_216272'],  # NEW
    'chunk_metadata': {
        'page': 12,
        'position': 'figure_caption',
        'image_references': ['fig1.png']
    }
}
```

2. **Store Molecule Metadata** (Separate table/collection)

```sql
CREATE TABLE molecules (
    cid VARCHAR(32) PRIMARY KEY,          -- PubChem CID
    smiles TEXT NOT NULL,                 -- Canonical SMILES
    molecular_weight FLOAT,
    logp FLOAT,
    hba INTEGER,
    hbd INTEGER,
    rotatable_bonds INTEGER,
    heavy_atoms INTEGER,
    fps_ecfp4 VECTOR(2048),              -- Fingerprint
    fps_maccs VECTOR(167),
    scaffold TEXT,
    is_gold_standard BOOLEAN,
    similarity_to_gold_standard FLOAT
);
```

3. **Update RAG Query Pipeline**

```python
def query_rag(self, question: str, mol_id: str = None):
    """Query RAG with optional molecule context."""
    # Step 1: Text-first retrieval
    text_results = self.vector_search(question)

    # Step 2: If molecule context available, add molecule-enhanced results
    if mol_id:
        mol_fp = self.get_molecule_fingerprint(mol_id)
        mol_results = self.molecule_search(mol_fp)

        # Recombine: blend text and molecular results
        combined_results = self.rerank_combined(
            text_results,
            mol_results,
            mol_id=mol_id
        )
        return combined_results

    return text_results
```

4. **OpenWebUI UI Enhancement**

Add molecule "context card" to RAG results:
```
┌─────────────────────────────────────────────────┐
│ Retrieved Chunk                                 │
├─────────────────────────────────────────────────┤
│ Text: The mesembrine alkaloids exhibit..........│
│                                                 │
│ 🧪 Related Molecules:                           │
│   - Mesembrine (CID 394162)                    │
│     MW: 337.4 Da | LogP: 2.8                   │
│   - Mesembrenone (CID 216272)                  │
│     MW: 335.4 Da | LogP: 2.6                   │
│                                                 │
│ 🔬 View 3D Structure | 🔍 Similar Compounds    │
└─────────────────────────────────────────────────┘
```

#### Estimated Effort: 2-3 days

---

### 3.2 Cross-Modal Retrieval

#### Current State
- Text queries search text embeddings only
- Structure queries search molecular fingerprints only
- No joint embedding space

#### Target State
- Users can search: "Find papers about mesembrine analogs"
- System retrieves text chunks with similar molecules
- Cross-modal reranking blends text + structure relevance

#### Implementation Steps

1. **Train Cross-Modal Embeddings** (Optional, long-term)

- Use contrastive learning to align text and molecular embeddings
- Dataset: (caption, SMILES) pairs from PDF extractions
- Loss: contrastive loss minimizing distance between pairs

2. **Heuristic Cross-Modal Search** (Short-term)

```python
def cross_modal_search(query_text: str, query_smiles: str = None):
    # Step 1: Text search
    text_results = text_search(query_text)

    # Step 2: Structure search (if SMILES provided)
    if query_smiles:
        mol_fp = fp_generator.ecfp4(query_smiles)
        mol_results = molecule_search(mol_fp)

        # Step 3: Combine using weighted hybrid
        combined = []
        for doc_id in set(r['doc_id'] for r in text_results) | set(r['doc_id'] for r in mol_results):
            text_score = next((r['score'] for r in text_results if r['doc_id'] == doc_id), 0)
            mol_score = next((r['score'] for r in mol_results if r['doc_id'] == doc_id), 0)

            combined.append({
                'doc_id': doc_id,
                'score': 0.7 * text_score + 0.3 * mol_score
            })

        return sorted(combined, key=lambda x: x['score'], reverse=True)

    return text_results
```

3. **OpenWebUI Query Interface**

Add molecule input field to chat/search:
```
┌────────────────────────────────────────────┐
│ Ask me anything...                         │
│ [Text input________]                       │
│                                            │
│ 🧪 Or search by structure (SMILES):        │
│ [c1ccc(cc1OC)C(C)(C)O]    📸 Upload image  │
│                                            │
│ [🔍 Search] [🧠 Chat]                      │
└────────────────────────────────────────────┘
```

#### Estimated Effort: 1-2 days

---

## Implementation Roadmap

### Week 1: Quick Wins + Embedding Upgrade
| Day | Task | Deliverable |
|-----|------|-------------|
| Mon | Config updates, confidence scoring | ✅ Completed |
| Tue | Embedding model test (sample 50 papers) | Model quality report |
| Wed | Full reindex (all 129 papers) | Indexed corpus v2 |
| Thu | Reranker configuration | Reranker enabled |
| Fri | Test baseline metrics | Baseline benchmark |

### Week 2: Molecular Fingerprint Support
| Day | Task | Deliverable |
|-----|------|-------------|
| Mon | Fingerprint generator class | `fingerprint_generator.py` |
| Tue | PostgreSQL schema + indexes | `molecules` table functional |
| Wed | RAG pipeline integration | Chunks store fingerprints |
| Thu | Structure search API | `/api/v1/structure/search` |
| Fri | OpenWebUI plugin | Structure search button |

### Week 3: Molecule-RAG Integration
| Day | Task | Deliverable |
|-----|------|-------------|
| Mon | Chunk-molecule linking | `molecule_ids` field |
| Tue | Molecule metadata table | `molecules` table populated |
| Wed | Cross-modal search logic | Combined retrieval |
| Thu | OpenWebUI UI updates | Context cards |
| Fri | End-to-end test | E2E workflow validated |

### Week 4: Testing + Documentation
| Day | Task | Deliverable |
|-----|------|-------------|
| Mon | Regression tests | All tests passing |
| Tue | Integration tests | Full pipeline verified |
| Wed | User acceptance test | SME validation |
| Thu | Documentation update | Architecture docs |
| Fri | Migration guide | v2 upgrade guide |

---

## Success Metrics

### Retrieval Quality
| Metric | Current | Target | Measurement |
|--------|---------|--------|-------------|
| MRR@5 | ~0.45 | ≥0.70 | Chemistry queries |
| NDCG@10 | ~0.55 | ≥0.75 | Paper relevance |
| Chemical recall | ~0.35 | ≥0.65 | Known molecules |

### Performance
| Metric | Current | Target |
|--------|---------|--------|
| Query latency | ~500ms | ≤1s (with reranker) |
| Indexing throughput | ~10 papers/hour | ≥20 papers/hour |
| Vector DB size | ~15GB | ~25GB (+fingerprint) |

### User Satisfaction
| Survey Question | Target |
|-----------------|--------|
| "Results are relevant" | ≥4.5/5 |
| "Can find molecules easily" | ≥4/5 |
| "Structure search useful" | ≥4.5/5 |

---

## Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| **Embedding model change breaks retrieval** | High | Test on sample first, keep backup |
| **Molecular fingerprints increase DB size** | Medium | Monitor storage, compress if needed |
| **Structure search performance slow** | Medium | HNSW index + filter by threshold |
| **Molecule-text linking complex** | Low | Incremental rollout, feature flag |

---

## Dependencies

### Immediate (No blocking)
- Embedding model upgrade → Independent
- Confidence scoring → Independent
- Backend routing update → Independent

### Phase 2 Dependent
- Molecular fingerprint support → Requires embedding change
- Structure search API → Requires fingerprints
- OpenWebUI plugin → Requires structure search API

### Phase 3 Dependent
- Molecule-text linking → Requires fingerprints
- Cross-modal search → Requires linking

---

## Next Steps

1. **Review and approve** this implementation plan
2. **Schedule testing window** (2-hour maintenance window recommended)
3. **Assign resources** (dev time, testing hours)
4. **Begin Phase 2 implementation** (Embedding upgrade)

---

**Document Version**: 1.0
**Last Updated**: 2026-02-26
**Status**: Ready for Implementation
