# SMILES Pipeline Compliance Checklist (Mapped to Current Code)

As of **2026-02-26**, this checklist maps best-practice expectations to the current implementation under `smiles-pipeline/`.

Status legend:
- `PASS`: implemented and observable in code path
- `PARTIAL`: partially implemented or implemented with caveats
- `GAP`: not implemented in active runtime path

## Summary Scorecard

| Area | Status | Notes |
|------|--------|-------|
| 3-layer validation gate | PASS | Indigo -> RDKit plausibility -> Domain validator in pipeline |
| Gold standards path SSOT | PASS | CLI/env/default precedence implemented in `DomainValidator` |
| Multi-backend OCSR fallback | PASS | Ordered fallback `molscribe,decimer` in pipeline |
| Strict sanitization taxonomy | PARTIAL | Parse checks exist, but no explicit `SanitizeMol(..., catchErrors=True)` error taxonomy |
| Standardization policy (tautomer/parent/charge) | PASS | Deterministic RDKit MolStandardize stage in `validate_molecule()` |
| Cross-system identity (`InChI`, `InChIKey`) | PASS | Persisted as `standard_inchi` and `standard_inchikey` on accepted records |
| Composite confidence (model + agreement + validation) | PASS | Score now includes backend confidence penalties, validation requirements, and optional cross-backend agreement bonus |
| OCSR confidence gating | PASS | `fallback.min_confidence` is enforced; low-confidence outputs only pass if cross-backend agreement is `true` |
| Config-file-driven behavior (`validation_rules.yaml`, `backends.yaml`) | PASS | Runtime now loads YAML config and maps thresholds/backend policy into validators/extractors |
| Reproducibility metadata capture | PASS | Run-level manifest now records toolkit versions, config checksums, run id/timestamps, and effective config |
| Structured rejection taxonomy | PASS | Rejections now persist `stage`, `reason_code`, `reason_detail` + aggregated counters |

## Component Mapping

| Best-practice control | Current implementation | Status | Evidence |
|----------------------|------------------------|--------|----------|
| Layered extraction + validation flow | `ExtractionPipeline.process_image()` runs extract -> validate -> enrich | PASS | `scripts/extract_smiles_pipeline.py` |
| Syntax validation gate | Indigo parse/canonical/placeholder/wildcard checks | PASS | `src/validators/indigo_validator.py` |
| Chemical plausibility gate | RDKit property/rules filtering | PASS | `src/validators/chemical_validator.py` |
| Domain + gold standard comparison | ECFP4/Tanimoto against configured references | PASS | `src/validators/domain_validator.py` |
| Gold standards path precedence | `--gold-standards-file` -> env -> canonical path | PASS | `src/validators/domain_validator.py` + CLI in `scripts/extract_smiles_pipeline.py` |
| OCSR engine fallback | sequential backend order; fallback on exceptions | PASS | `scripts/extract_smiles_pipeline.py` |
| Multi-engine consensus scoring | Optional second-backend consensus check can be enabled in `backends.yaml` (`routing.consensus_check`) | PARTIAL | `extract_from_image()` compares canonical outputs when consensus is enabled |
| Deterministic standardization | Cleanup -> FragmentParent -> Uncharge -> CanonicalTautomer | PASS | `src/validators/standardization_validator.py` + `scripts/extract_smiles_pipeline.py` |
| InChI/InChIKey persistence | identity fields persisted in output records | PASS | output record schema in `scripts/extract_smiles_pipeline.py` |
| Chemistry error taxonomy | rejection reasons exist; sanitize failure flags not explicitly captured | PARTIAL | `src/validators/chemical_validator.py` |
| Confidence policy alignment | Composite score enforces low-confidence penalties and positive-indicator requirement for high-confidence routing | PASS | confidence scoring and requirement enforcement in `scripts/extract_smiles_pipeline.py` |
| DECIMER confidence semantics | heuristic confidence estimation, not model-native score | PARTIAL | `src/extractors/decimer_extractor.py` |
| Config as runtime SSOT | Pipeline loads config files and applies policy defaults/overrides | PASS | YAML loader and config mapping in `scripts/extract_smiles_pipeline.py` |
| Extraction acceptance gate | Rejects below `fallback.min_confidence` unless consensus check explicitly agrees | PASS | `process_image()` confidence gate in `scripts/extract_smiles_pipeline.py` |
| Runtime/toolkit manifest | `extraction_stats.json` includes dependency versions + config hashes + run ids/timestamps | PASS | stats payload assembly in `scripts/extract_smiles_pipeline.py` |
| QA rejection analytics | `generate_qa_report.py` computes rejection counts by stage/reason/backend | PASS | `compute_rejection_analysis()` + report output |
| Indigo runtime dependency readiness | Pipeline now supports runtime preflight checks and actionable import failures | PASS | `--preflight` in `scripts/extract_smiles_pipeline.py` + lazy Indigo import guard |

## Prioritized Remediation Plan

1. Add structured rejection taxonomy for every failed path.
Status: `Done`
Delivered:
- persists `stage`, `reason_code`, `reason_detail` uniformly
- includes backend + confidence context for triage analytics

2. Add Indigo runtime preflight + option pinning.
Target:
- dependency readiness check with actionable install guidance
- explicit Indigo option policy (`timeout`, stereo/aromaticity options) recorded in manifest
Status: `Delivered (baseline)` (preflight checks + runtime option pinning for `timeout`/reset policy + manifest capture)

3. Add native DECIMER confidence integration.
Target:
- use model-native confidence output when available, keep heuristic fallback as backup

4. Add evaluation slices and calibration reporting.
Target:
- evaluate by modality (`single`, `multi`, `hand-drawn`, stereochemistry-heavy)
- track exact match + graph similarity + stereo retention + calibration error

## Definition of Done for "Best-Practice Compliant"

- All accepted molecules include:
`smiles`, `canonical_smiles`, `standardized_smiles`, `standard_inchi`, `standard_inchikey`, provenance metadata.
- Runtime behavior is driven by `config/backends.yaml` and `config/validation_rules.yaml` (no hidden constants).
- Confidence routing has 3 states:
`READY`, `MANUAL_REVIEW`, `REJECT`.
- Pipeline run output includes effective configuration, run id/timestamps, and toolkit version manifest.
- QA report includes rejection taxonomy with counts by reason and backend.

## Verification Commands

```bash
# Confirm CLI exposes gold standard override
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py --help | rg "gold-standards-file"

# Confirm active references in code
rg -n "SMILES_GOLD_STANDARDS_FILE|gold_standards|MolStandardize|InChI|InChIKey|validation_rules|backends.yaml" smiles-pipeline/src smiles-pipeline/scripts

# Syntax check key modules
python -m py_compile \
  smiles-pipeline/scripts/extract_smiles_pipeline.py \
  smiles-pipeline/src/validators/indigo_validator.py \
  smiles-pipeline/src/validators/chemical_validator.py \
  smiles-pipeline/src/validators/domain_validator.py
```
