# RAG Remediation Learnings

## Task T1: Recreate materialize_evidence.py

### What was done
- Created `/LAB/@thesis/openwebui/services/rag/materialize_evidence.py` with minimal surviving exports
- Exports: MaterializedRecord, make_point_id, make_qdrant_point_id, SourceDocument, PageInput, FigureInput, MoleculeInput, materialize_pages, materialize_figures, materialize_molecules

### Key observations
1. **Pre-existing import failure**: The `services.rag` package has a broken import chain due to `query_router.py` trying to import from deleted `normalize_smiles.py`. This prevents standard `from services.rag.materialize_evidence import ...` from working.
2. **Verification workaround**: Direct module loading via `importlib.util.spec_from_file_location` bypasses the broken package init and confirms the module is syntactically correct and exports all required symbols.
3. **Chemistry stripped**: `materialize_molecules()` returns empty list as expected - chemistry is being removed in this wave but the function signature is preserved for script compatibility.
4. **smiles_for_embedding kept**: The `MaterializedRecord` dataclass still has `smiles_for_embedding: str | None = None` because `embed_and_upsert.py` line 182 uses it - will be stripped in Wave 2.

### Dependencies unblocked
- T2 (manifest.py) - now can import MaterializedRecord
- T3 (embed_and_upsert.py) - now can import MaterializedRecord and make_qdrant_point_id
- T10 (scripts) - now can import make_point_id and make_qdrant_point_id

### Remaining blocker
- T5 (query_router.py deletion) needs to run to fix the broken package import chain before standard imports will work
