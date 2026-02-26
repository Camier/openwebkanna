-- Migration: Add molecular fingerprint support to documents table
-- Phase 2, Task 2.1
-- Date: 2026-02-26
--
-- This migration adds columns to store molecular fingerprints (ECFP4)
-- alongside text embeddings for structure similarity search.

BEGIN;

-- Add molecular columns to documents table (embedded schema)
-- These columns store:
--   molecule_fingerprint: ECFP4 2048-bit vector for similarity search
--   molecule_smiles: Canonical SMILES string
--   molecule_metadata: Additional molecule properties (MW, LogP, etc.)

ALTER TABLE documents
ADD COLUMN IF NOT EXISTS molecule_fingerprint VECTOR(2048),
ADD COLUMN IF NOT EXISTS molecule_smiles TEXT,
ADD COLUMN IF NOT EXISTS molecule_metadata JSONB;

-- Create HNSW index for efficient molecular similarity search
-- Parameters:
--   m = 16: Number of connections per node (higher = more accurate, slower)
--   ef_construction = 64: Size of dynamic candidate list during index build
--
-- This enables sub-second similarity searches even with millions of molecules

CREATE INDEX IF NOT EXISTS documents_mol_fp_idx
ON documents USING hnsw (molecule_fingerprint vector_cosine_ops)
WITH (m = 16, ef_construction = 64);

-- Add comments for documentation
COMMENT ON COLUMN documents.molecule_fingerprint IS 'ECFP4 molecular fingerprint (2048-bit) for structure similarity search';
COMMENT ON COLUMN documents.molecule_smiles IS 'Canonical SMILES string for the molecule';
COMMENT ON COLUMN documents.molecule_metadata IS 'JSONB metadata: {source, generated_at, properties, ...}';
COMMENT ON INDEX documents_mol_fp_idx IS 'HNSW index for molecular fingerprint similarity search (cosine distance)';

-- Grant permissions (adjust role name if needed)
GRANT SELECT ON documents TO openwebui;

-- Verify the migration
-- Uncomment for debugging:
-- \d documents
-- SELECT indexname, indexdef FROM pg_indexes WHERE tablename = 'documents';

COMMIT;
