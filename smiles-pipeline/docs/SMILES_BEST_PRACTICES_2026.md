# SMILES & Cheminformatics Best Practices (2026)

As of **2026-02-26**, this document consolidates current technical best practices for SMILES-centric pipelines, with priority on production extraction, validation, and RAG interoperability.

## Scope

- SMILES extraction and normalization for cheminformatics workflows
- OCSR (image/PDF -> SMILES) production architecture
- Molecular identity, reproducibility, and ML/evaluation hygiene
- Toolkit operations for RDKit, Open Babel, Indigo, and CDK

## Gold Standard Practices

1. Run a deterministic 3-stage curation pipeline before downstream use.
Why: reduces silent structural noise and keeps outputs reproducible.
Recommended order: syntax/parse validation -> structure standardization -> parent extraction (salt/solvent handling).

2. Do not treat canonical SMILES as a universal global identifier.
Why: canonicalization is toolkit/version-specific.
Use Standard InChI and InChIKey as interoperability identity, while retaining toolkit-specific canonical SMILES for local processing.

3. Pin toolkit versions and parser options explicitly.
Why: parser, aromaticity perception, tautomer behavior, and canonicalization can change across releases.
Record exact toolkit versions and critical parse/standardize flags in run metadata.

4. Keep interchange representations in OpenSMILES-compatible core.
Why: extensions (for example CXSMILES/vendor-specific annotations) are not guaranteed to round-trip across toolchains.
Use extension fields only when producer and consumer both explicitly support them.

5. Make tautomer/protonation/charge policy explicit and stable.
Why: deduplication, retrieval, and labels can drift if normalization policy changes.
Define whether canonical tautomerization is used, and document charge/fragment parent logic.

6. Require post-extraction chemistry validation for OCSR outputs.
Why: OCR/vision outputs can be syntactically parseable but chemically implausible.
Use strict parse + sanitization + plausibility checks before acceptance.

7. Use multi-engine OCSR with composite confidence.
Why: [Inference] consensus from independent engines plus chemistry validation is more robust than a single score.
Composite confidence should include:
- model-native confidence
- cross-engine structural agreement
- chemistry validator pass/fail and warning severity

8. Maintain multi-representation storage for ML and retrieval.
Recommended minimum fields:
- raw_smiles
- canonical_smiles (toolkit/version scoped)
- standard_inchi
- standard_inchikey
- provenance (backend, model version, thresholds, extraction time)

9. Evaluate models with multiple split regimes.
Why: single split modes can misrepresent generalization.
Report at least random split and scaffold split; add one harder shift split (cluster or temporal) when available.

10. Treat reproducibility metadata as mandatory.
Minimum metadata:
- toolkit/backend versions
- standardization policy revision
- split definitions and random seeds
- dropped-invalid counts with reasons
- exact train/validation/test IDs for benchmarkable tasks

## OCSR Architecture Pattern (Recommended)

1. Ingest: PDF/image intake and page rasterization
2. Detect: structure region segmentation
3. Recognize: parallel OCSR engines (primary + fallback)
4. Validate: syntax + chemistry + domain rules
5. Standardize: canonical/parent/tautomer policy
6. Score: composite confidence and triage labels
7. Review: manual queue for low-confidence/failed structures
8. Persist: identity fields + provenance + quality indicators

## Validated Gold Standard Path (2026-02-26)

This path consolidates architecture decisions validated against recent OCSR/cheminformatics documentation and peer-reviewed references.

1. Segment and classify first.
- Use DECIMER segmentation/classification before recognition for multi-panel and noisy figures.

2. Run multi-engine recognition with deterministic reconciliation.
- Primary: MolScribe.
- Secondary: DECIMER.
- Optional modality-specific fallback (for example OSRA for some multi-structure layouts).

3. Canonicalize and validate in strict stages.
- Stage A: Indigo syntax parsing + wildcard/placeholder policy.
- Stage B: RDKit standardization (cleanup -> fragment parent -> charge/tautomer policy).
- Stage C: RDKit chemical plausibility.
- Stage D: domain constraints + gold standard comparison.

4. Use composite confidence with triage labels.
- Combine backend confidence, domain/chemistry indicators, and optional cross-backend agreement.
- Route into three states:
  - `READY`
  - `MANUAL_REVIEW`
  - `REJECT_OR_DEEP_REVIEW`

5. Persist interoperability identities and full provenance.
- Always store `standardized_smiles`, `standard_inchi`, and `standard_inchikey`.
- Store toolchain versions and policy IDs (RDKit/Indigo/OCSR versions, config hash, thresholds, run timestamp).

6. Close the loop with calibration and review feedback.
- Calibrate confidence on held-out in-domain sets (reliability curves / post-hoc calibration).
- Feed reviewed corrections back into training/evaluation sets.

### Short Action Plan

1. Keep conservative consensus check policy in `config/backends.yaml` and enable it in stricter profiles for low-confidence outputs.
Status: `Implemented (config-gated)` via `routing.consensus_check`.
2. Keep runtime/toolkit manifest in run outputs (`extraction_stats.json` + per-record provenance).
Status: `Implemented` (run id, timestamps, dependency versions, config hashes, backend runtime metadata).
3. Add structured rejection taxonomy fields (`stage`, `reason_code`, `reason_detail`) to all rejected paths.
Status: `Implemented` (runtime rejection events + QA aggregation by stage/reason/backend).
4. Add evaluation slices (single/multi/hand-drawn/chirality) and track exact-match + graph similarity + stereo retention.
Status: `Next`.
5. Introduce periodic confidence calibration job and drift checks.
Status: `Next`.
6. Add explicit Indigo runtime preflight and option pinning (timeout, aromaticity/stereo options, per-run option reset policy).
Status: `Implemented (baseline)` via `--preflight` and `level_1_syntax.runtime_options` / `set_options`.

## Current Ecosystem Snapshot (Checked 2026-02-26)

- RDKit latest public release line: `2025.9.5` (PyPI, released 2026-02-16)
- Open Babel latest release: `3.1.1`
- Indigo latest listed release: `1.40.0` (released 2026-02-23)
- DECIMER-Image_Transformer latest package release: `2.8.0` (PyPI, released 2025-12-02)
- DECIMER-Image_Transformer latest GitHub tag/release: `v2.7.2` (2025-12-02)
- DECIMER-Image-Segmentation latest PyPI wheel release: `1.5.0` (2025-08-12)

Note: Production reproducibility should use runtime manifests produced by the pipeline, not this snapshot.

## Source Links

- ChEMBL Structure Pipeline: <https://github.com/chembl/ChEMBL_Structure_Pipeline>
- ChEMBL NAR 2024 update: <https://academic.oup.com/nar/article/52/D1/D1180/7337608>
- PubChem compounds standardization note: <https://pubchem.ncbi.nlm.nih.gov/docs/compounds>
- OpenSMILES specification: <https://opensmiles.org/opensmiles.html>
- InChI standard overview: <https://www.inchi-trust.org/about-the-inchi-standard/>
- InChI technical FAQ: <https://www.inchi-trust.org/technical-faq/>
- IUPAC InChI 1.07 release note: <https://iupac.org/inchi-1-07-available-on-github/>
- RDKit documentation index: <https://mail.rdkit.org/docs/>
- RDKit releases: <https://github.com/rdkit/rdkit/releases>
- RDKit PyPI releases: <https://pypi.org/project/rdkit/#history>
- RDKit backward incompatible changes: <https://mail.rdkit.org/docs/BackwardsIncompatibleChanges.html>
- RDKit parser params: <https://www.rdkit.org/docs/cppapi/structRDKit_1_1v1_1_1SmilesParserParams.html>
- RDKit MolStandardize API: <https://mail.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html>
- Open Babel docs: <https://open-babel.readthedocs.io/>
- Open Babel canonical notes: <https://openbabel.org/api/3.0/canonical_code_algorithm.shtml>
- Open Babel strict Smiley parser: <https://open-babel.readthedocs.io/en/latest/FileFormats/SMILES_format_using_Smiley_parser.html>
- Indigo docs: <https://lifescience.opensource.epam.com/indigo/>
- Indigo release notes index: <https://lifescience.opensource.epam.com/indigo/release-notes/>
- Indigo 1.40.0 release note: <https://lifescience.opensource.epam.com/indigo/release-notes/indigo-1.40.0.html>
- CDK website: <https://cdk.github.io/>
- TDC split docs: <https://tdcommons.ai/functions/data_split/>
- DeepChem splitters: <https://deepchem.readthedocs.io/en/2.8.0/api_reference/splitters.html>
- Scaffold split analysis (2024): <https://arxiv.org/abs/2406.00873>
- SELFIES docs: <https://selfiesv2.readthedocs.io/en/latest/index.html>
- MolScribe repo: <https://github.com/thomas0809/MolScribe>
- MolScribe PyPI: <https://pypi.org/project/molscribe/>
- DECIMER-Image_Transformer: <https://github.com/Kohulan/DECIMER-Image_Transformer>
- DECIMER-Image_Transformer releases: <https://github.com/Kohulan/DECIMER-Image_Transformer/releases>
- DECIMER package releases: <https://pypi.org/project/decimer/#history>
- DECIMER-Image-Segmentation: <https://github.com/Kohulan/DECIMER-Image-Segmentation>
- DECIMER-Segmentation package releases: <https://www.piwheels.org/project/decimer-segmentation>
- DECIMER.ai paper: <https://www.nature.com/articles/s41467-023-40782-0>
- OSRA project files: <https://sourceforge.net/projects/osra/files/osra/>

## Notes

- Items marked as [Inference] are architectural recommendations inferred from multiple official sources, not direct single-source prescriptions.
- This document is intentionally implementation-agnostic. For repository-specific status, see `SMILES_PIPELINE_COMPLIANCE_CHECKLIST.md`.
