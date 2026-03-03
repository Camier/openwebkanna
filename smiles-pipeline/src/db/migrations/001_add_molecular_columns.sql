-- Migration: Add molecular fingerprint support to document_chunk table
-- Phase 2, Task 2.1
-- Date: 2026-02-26
--
-- This migration adds columns to store molecular fingerprints (ECFP4)
-- alongside text embeddings for structure similarity search.
--
-- TECHNICAL NOTE: pgvector v0.8.1 supports max 2000 dimensions for 'vector' type.
-- ECFP4 fingerprints are 2048-dim, so we use halfvec (half-precision) which supports up to 4000 dims.
-- Reference: https://github.com/pgvector/pgvector#indexing

BEGIN;

-- Add molecular columns to document_chunk table (embedded schema)
-- These columns store:
--   molecule_fingerprint: ECFP4 2048-bit as half-precision vector
--   molecule_smiles: Canonical SMILES string
--   molecule_metadata: Additional molecule properties (MW, LogP, etc.)

ALTER TABLE document_chunk
ADD COLUMN IF NOT EXISTS molecule_fingerprint HALFVEC(2048),
ADD COLUMN IF NOT EXISTS molecule_smiles TEXT,
ADD COLUMN IF NOT EXISTS molecule_metadata JSONB;

-- Create IVFFlat index using half-precision casting for 2048-dim ECFP4
-- Note: Direct indexing of vector(2048) is not supported (max 2000 dims)
-- Solution: Cast to halfvec during indexing (supports up to 4000 dims)
--
-- Performance tradeoff: half-precision has slightly lower accuracy but much smaller index size
-- For chemical fingerprints, this is typically acceptable (binary-like nature of ECFP4)

CREATE INDEX IF NOT EXISTS document_chunk_mol_fp_idx
ON document_chunk USING ivfflat (molecule_fingerprint halfvec_cosine_ops)
WITH (lists = 100);

-- Add comments for documentation
COMMENT ON COLUMN document_chunk.molecule_fingerprint IS 'ECFP4 molecular fingerprint (2048-bit) for structure similarity search';
COMMENT ON COLUMN document_chunk.molecule_smiles IS 'Canonical SMILES string for the molecule';
COMMENT ON COLUMN document_chunk.molecule_metadata IS 'JSONB metadata: {source, generated_at, properties, ...}';
COMMENT ON INDEX document_chunk_mol_fp_idx IS 'IVFFlat index on halfvec-cast molecular fingerprint (cosine distance)';

-- Grant permissions (adjust role name if needed)
GRANT SELECT ON document_chunk TO openwebui;

-- Verify the migration
-- Uncomment for debugging:
-- \d document_chunk
-- SELECT indexname, indexdef FROM pg_indexes WHERE tablename = 'document_chunk';

COMMIT;
