# SMILES Operator Must-Read (Top 5)

Updated: 2026-02-27

This is the shortest possible official reading list for running and maintaining the SMILES pipeline in production.

1. RDKit MolStandardize API
   - https://mail.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html
   - Why operators care: defines the exact cleanup/fragment/charge/tautomer standardization behavior used for deterministic outputs.

2. Indigo Documentation
   - https://lifescience.opensource.epam.com/indigo/
   - Why operators care: syntax parsing and canonicalization behavior for Level-1 gate, including runtime options and error modes.

3. InChI Trust Technical FAQ
   - https://www.inchi-trust.org/technical-faq/
   - Why operators care: identity policy for cross-system interoperability (`standard_inchi`, `standard_inchikey`).

4. DECIMER-Image_Transformer (Official Repository)
   - https://github.com/Kohulan/DECIMER-Image_Transformer
   - Why operators care: fallback OCSR backend behavior, release notes, and version-sensitive runtime expectations.

5. pgvector (Official Documentation/Repository)
   - https://github.com/pgvector/pgvector
   - Why operators care: vector/halfvec indexing and cosine-distance query semantics used for molecular retrieval.

## How To Use This List

- Before changing extraction logic: read 1 + 2 + 4.
- Before changing identity/storage fields: read 3 + 5.
- Before major rollout: verify pinned runtime versions match these docs and current pipeline manifests.
