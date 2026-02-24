# SMILES Pipeline for OpenWebUI Retrieval

This repository now includes a reproducible molecular extraction pipeline that can publish high-confidence molecule artifacts into an OpenWebUI Knowledge Base.

## What it produces

Per paper directory under `prod_max/`:

- `smiles_extracted.json`
- `molecules.jsonl`
- `molecules_high_confidence.jsonl`

Run-level output under `prod_max/`:

- `smiles_qa_summary.json`

## Main workflow

Run the full pipeline (extract -> embed -> validate -> high-confidence filter -> QA summary):

```bash
LIMIT=3 EXTRACT_ARGS='--backend-order molscribe,decimer,vlm --max-new-tokens 128' mise run smiles:pipeline
```

Publish high-confidence artifacts to OpenWebUI KB:

```bash
OPENWEBUI_URL=http://localhost:3010 mise run smiles:publish-kb
```

If duplicate KB entries exist for the same name:

```bash
OPENWEBUI_URL=http://localhost:3010 ARGS='--keep-id <kb-id> --apply' mise run smiles:kb-dedupe
```

## Tasks

- `smiles:extract`
- `smiles:extract-all`
- `smiles:embed`
- `smiles:validate`
- `smiles:high-confidence`
- `smiles:qa-summary`
- `smiles:publish-kb`
- `smiles:kb-dedupe`
- `smiles:pipeline`

## Environments

- `scripts/environment-smiles-extraction.yml` for extraction + RDKit + DECIMER
- `scripts/environment-molscribe-legacy.yml` for isolated MolScribe subprocess execution

The extractor calls MolScribe through `scripts/molscribe_predict.py` in the legacy environment to avoid dependency conflicts in the main runtime.
