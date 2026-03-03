BEGIN;

DO $$
BEGIN
    IF EXISTS (
        SELECT 1
        FROM information_schema.columns
        WHERE table_name = 'document_chunk'
          AND column_name = 'molecule_fingerprint'
          AND udt_name = 'vector'
    ) THEN
        ALTER TABLE document_chunk
        ALTER COLUMN molecule_fingerprint TYPE halfvec(2048)
        USING (molecule_fingerprint::halfvec(2048));
    END IF;
END $$;

DO $$
BEGIN
    IF EXISTS (
        SELECT 1
        FROM information_schema.columns
        WHERE table_name = 'document_chunk'
          AND column_name = 'molecule_fingerprint_maccs'
          AND udt_name = 'vector'
    ) THEN
        ALTER TABLE document_chunk
        ALTER COLUMN molecule_fingerprint_maccs TYPE halfvec(167)
        USING (molecule_fingerprint_maccs::halfvec(167));
    END IF;
END $$;

DROP INDEX IF EXISTS document_chunk_mol_fp_idx;
DROP INDEX IF EXISTS document_chunk_mol_maccs_idx;

CREATE INDEX IF NOT EXISTS document_chunk_mol_fp_idx
ON document_chunk USING ivfflat (molecule_fingerprint halfvec_cosine_ops)
WITH (lists = 100);

CREATE INDEX IF NOT EXISTS document_chunk_mol_maccs_idx
ON document_chunk USING ivfflat (molecule_fingerprint_maccs halfvec_cosine_ops)
WITH (lists = 100);

COMMENT ON INDEX document_chunk_mol_fp_idx IS 'IVFFlat index on halfvec molecular fingerprint (cosine distance)';
COMMENT ON INDEX document_chunk_mol_maccs_idx IS 'IVFFlat index on halfvec MACCS fingerprint (cosine distance)';

COMMIT;
