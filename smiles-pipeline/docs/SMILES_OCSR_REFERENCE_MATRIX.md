# OCSR Reference Matrix (2026-02-24)

This matrix captures source-backed implementation references for image-to-SMILES extraction.

## Shortlist

| Project | Role | Strength | Caveat | Recency signal |
| --- | --- | --- | --- | --- |
| `thomas0809/MolScribe` | Primary extraction backend | Strong practical extraction surface, checkpoints, and confidence output support | Heavy training/eval profile if retraining | Active repo; latest commit observed in 2025 |
| `audivir/molscribe-server` | Service wrapper pattern | Fast path to API-style deployment around MolScribe | Older dependency pinning risk | Updated in 2025 |
| `Kohulan/DECIMER-Image_Transformer` | Complementary extraction backend | Canonical DECIMER successor docs and maintained transformer line | Different stack behavior vs MolScribe; benchmark before promotion | Repo/docs current in 2026 |
| `Kohulan/DECIMER-Image-to-SMILES` | Legacy/historical baseline | Useful for backward context and older workflows | Explicit migration pointer to newer DECIMER repo | Older baseline |
| `hustvl/MolSight` | Modern multimodal OCSR candidate | Strong modern approach signal (multimodal + RL + stereo focus) | Research-grade integration cost; needs local benchmark proof | Recent releases in late 2025 to early 2026 |
| `chemvision/chemvision` | Multimodal benchmark/tooling context | Useful for multimodal chemistry workflow ideas | Not a drop-in OCR extraction service | Updated in 2025 |
| `epam/Indigo` | Cheminformatics interoperability toolkit | Canonical SMILES tooling (`indigo-cano`), broad language bindings, search/render utilities | Not an OCSR extractor; best as interoperability/validation adjunct | Latest release observed: 1.40.0 (2026-02-23) |
| `hsiaoyi0504/awesome-cheminformatics` | Discovery index | Curated map to candidate libraries/tools for future evaluation | Not authoritative for performance; entries require independent validation | Active curation; no release cadence |

## Recommended default posture

1. Keep `molscribe,decimer,vlm` backend order for production extraction.
2. Continue RDKit canonicalization and strict high-confidence filtering before KB ingestion.
3. Keep Indigo as optional interoperability support (format conversion/search adjunct), not a replacement for OCSR extraction backends.
4. Run MolSight as an offline benchmark track only, then promote if it improves valid/high-confidence yield on local papers.

## Validation references

- RDKit API: `https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html`
- RDKit cookbook: `https://www.rdkit.org/docs/Cookbook.html`

## Note

This matrix summarizes documentation and repository evidence gathered via search-mode background agents on 2026-02-24.
