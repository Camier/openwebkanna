#!/bin/bash
# Recreate rag_evidence collection with full multimodal schema (text + vision + molecule)
set -e

# Load Qdrant credentials
source /home/miko/.config/qdrant/qdrant.env

QDRANT_URL="${QDRANT_URL:-http://127.0.0.1:6333}"
QDRANT_API_KEY="${QDRANT__SERVICE__API_KEY}"

echo "=== Step 1: Delete existing rag_evidence collection ==="
curl -s -X DELETE "${QDRANT_URL}/collections/rag_evidence" \
    -H "Api-Key: ${QDRANT_API_KEY}" &&
    echo "Deleted" || echo "Already gone or didn't exist"

echo ""
echo "=== Step 2: Create rag_evidence with full schema ==="

# Create collection with text_dense, vision_li, and chem_dense vectors
curl -s -X PUT "${QDRANT_URL}/collections/rag_evidence" \
    -H "Api-Key: ${QDRANT_API_KEY}" \
    -H "Content-Type: application/json" \
    -d '{
        "vectors": {
            "text_dense": {
                "size": 2560,
                "distance": "Cosine"
            },
            "vision_li": {
                "size": 128,
                "distance": "Cosine",
                "multivector_config": {
                    "comparator": "max_sim"
                }
            },
            "chem_dense": {
                "size": 768,
                "distance": "Cosine"
            }
        },
        "sparse_vectors": {
            "text_sparse": {
                "modifier": "idf"
            }
        },
        "on_disk_payload": true
    }'

echo ""
echo "=== Step 3: Create payload indexes ==="

# Create payload indexes for filtering
for field in object_type document_id page_id page_number figure_id molecule_id canonical_smiles inchikey review_status figure_kind has_smiles; do
    curl -s -X PUT "${QDRANT_URL}/collections/rag_evidence/index" \
        -H "Api-Key: ${QDRANT_API_KEY}" \
        -H "Content-Type: application/json" \
        -d "{
            \"field_name\": \"${field}\",
            \"field_schema\": \"keyword\"
        }" 2>/dev/null || true
done

# Integer field
curl -s -X PUT "${QDRANT_URL}/collections/rag_evidence/index" \
    -H "Api-Key: ${QDRANT_API_KEY}" \
    -H "Content-Type: application/json" \
    -d '{
        "field_name": "page_number",
        "field_schema": "integer"
    }' 2>/dev/null || true

# Boolean field
curl -s -X PUT "${QDRANT_URL}/collections/rag_evidence/index" \
    -H "Api-Key: ${QDRANT_API_KEY}" \
    -H "Content-Type: application/json" \
    -d '{
        "field_name": "has_smiles",
        "field_schema": "bool"
    }' 2>/dev/null || true

echo ""
echo "=== Collection created successfully ==="
echo "Schema: text_dense (2560) + vision_li (128 multivector) + chem_dense (768)"
echo "Sparse: text_sparse (BM25)"
echo ""
echo "Next steps:"
echo "1. Run text migration: python scripts/rag/migrate_legacy_text_points.py --target-collection rag_evidence"
echo "2. Run figure indexing: python scripts/rag/index_colmodernvbert_figures.py ..."
echo "3. Run molecule upsert: python scripts/rag/upsert_resolved_molecule_points.py ..."
