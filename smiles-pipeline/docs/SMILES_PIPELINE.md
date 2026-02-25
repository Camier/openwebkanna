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

Publish markdown-rendered high-confidence artifacts (A/B retrieval path):

```bash
OPENWEBUI_URL=http://localhost:3010 mise run smiles:publish-kb-markdown
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
- `smiles:publish-kb-markdown`
- `smiles:kb-dedupe`
- `smiles:pipeline`

## Environments

- `scripts/environment-smiles-extraction.yml` for extraction + RDKit + DECIMER
- `scripts/environment-molscribe-legacy.yml` for isolated MolScribe subprocess execution

The extractor calls MolScribe through `scripts/molscribe_predict.py` in the legacy environment to avoid dependency conflicts in the main runtime.

## Source-backed OCSR guidance (2026-02-24)

This section captures the latest canonical references gathered in search mode for DECIMER, MolScribe, RDKit validation, and modern multimodal OCSR.

Primary references:

- DECIMER Image Transformer docs: `https://decimer-image-transformer.readthedocs.io/en/latest/`
- DECIMER legacy repo (migration context): `https://github.com/Kohulan/DECIMER-Image-to-SMILES`
- MolScribe canonical repo: `https://github.com/thomas0809/MolScribe`
- Indigo toolkit repo: `https://github.com/epam/Indigo`
- Awesome cheminformatics index: `https://github.com/hsiaoyi0504/awesome-cheminformatics`
- RDKit rdmolfiles API docs: `https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html`
- RDKit cookbook: `https://www.rdkit.org/docs/Cookbook.html`
- Modern multimodal OCSR candidate (MolSight): `https://github.com/hustvl/MolSight`

### Practical stack choice

- Default production path in this repository: MolScribe + RDKit canonicalization and filtering.
- Keep DECIMER as a complementary backend in `SMILES_BACKEND_ORDER` for recovery and coverage.
- Treat legacy DECIMER V1 (`DECIMER-Image-to-SMILES`) as historical context, not the primary target for new deployment work.
- Keep modern multimodal OCSR (for example MolSight) as an evaluation track, then promote only after benchmark wins on local paper images.

### Validation profile

- Keep RDKit validation enabled by default unless explicitly disabled.
- Keep Indigo as an optional interoperability and toolkit fallback for format/search workflows outside primary RDKit validation path.
- Canonicalize with `MolFromSmiles` and `MolToSmiles` before writing molecule rows.
- Maintain strict placeholder rejection and post-filtering (`molecules_high_confidence.jsonl`) before KB ingestion.
- Use awesome-cheminformatics as a discovery index only, then validate each candidate tool against local extraction and retrieval benchmarks before adoption.

### Recommended operating modes

Baseline extraction:

```bash
LIMIT=3 EXTRACT_ARGS='--backend-order molscribe,decimer,vlm --write-molecules-jsonl' mise run smiles:pipeline
```

Confidence-aware extraction (captures backend confidence metadata where available):

```bash
LIMIT=3 EXTRACT_ARGS='--backend-order molscribe,decimer,vlm --write-molecules-jsonl --molscribe-return-confidence --decimer-return-confidence' mise run smiles:pipeline
```

Candidate modernization track (research only, separate from baseline):

- Benchmark modern multimodal OCSR alternatives against the same `prod_max` subset.
- Compare valid SMILES yield, high-confidence yield, and downstream KB retrieval utility.

## OpenWebUI ingestion behavior notes

The current importer (`import-smiles-to-kb.sh`) already follows the documented OpenWebUI sequence:

1. upload file (`POST /api/v1/files/?process=true&process_in_background=false`)
2. check processing status (`GET /api/v1/files/{id}/process/status`)
3. attach file to KB (`POST /api/v1/knowledge/{id}/file/add`)

Operational guidance:

- Keep `process_in_background=false` for deterministic batch automation.
- Keep duplicate-content handling idempotent: if the API reports duplicate content while attaching, treat it as already indexed for retry-safe imports.
- Raw JSON/JSONL ingestion works, but retrieval can improve when records are transformed into stable, sectioned markdown with consistent headers.
- If testing markdown conversion, run side-by-side benchmarks on the same paper subset and compare high-confidence retrieval utility.

### A/B retrieval snapshot (2026-02-24)

Bounded live ingestion was run with `--limit 1` on the same source paper to compare retrieval behavior:

- Markdown KB: `Molecular SMILES High Confidence Markdown AB 2026-02-24` (`90972141-1cae-44b2-ad6a-7913c74a9261`)
- Raw JSONL KB: `Molecular SMILES High Confidence JSONL AB 2026-02-24` (`2ce11ad4-7035-4277-9b79-2c987fbf75ae`)

Observed retrieval shape from `POST /api/v1/retrieval/query/doc` for four queries (`k=3`):

- Both paths returned 3 documents for each query.
- Markdown path returned sectioned, molecule-first snippets (`## Molecules`, `### Molecule N`, `SMILES: ...`).
- Raw path returned JSON object strings (`{"paper": ..., "paper_dir": ..., ...}`), which are less readable in retrieved chunks.
- Top distance values were close overall; markdown was slightly better in 3/4 probes and slightly worse in 1/4 probe.

Operator recommendation:

- Keep `smiles:publish-kb` (raw JSONL) as the default for compatibility and minimal transform risk.
- Use `smiles:publish-kb-markdown` when retrieval readability and stable section-level grounding are priorities.
- For production promotion, rerun this comparison on a larger bounded set (for example `--limit 10`) and track query-level hit quality alongside distance.

### A/B retrieval snapshot expansion (2026-02-24, `--limit 10`)

Expanded bounded ingestion was run with `--limit 10` for both paths. Under the current `SOURCE_ROOT` dataset this yielded 3 available artifacts per run (`Artifacts seen: 3`), so the comparison remained like-for-like across the same 3 papers:

- Markdown KB: `SMILES A-B Markdown Limit10 2026-02-24` (`8c62ec87-cb90-4fbf-8781-ebb500ecfcf0`)
- Raw JSONL KB: `SMILES A-B Raw Limit10 2026-02-24` (`ffbb27df-3202-4436-8d54-d346446dccbc`)

Retrieval probes (`k=3`) used these queries:

- `mesembrine alkaloid structure`
- `Sceletium tortuosum mesembrenone`
- `SMILES notation for mesembrine`

Observed behavior from `POST /api/v1/retrieval/query/doc`:

- Both paths returned 3 documents for each query (no coverage delta in this bounded run).
- Markdown results were consistently sectioned and molecule-first (`## Molecules`, `### Molecule N`, `SMILES: ...`).
- Raw JSONL results were consistently single-line JSON object strings, less readable for direct analyst grounding.
- Distance trends were mixed across queries. Markdown had tighter clustering (~0.926 to ~0.931 average per query), while raw had higher variance (~0.915 to ~0.933 average per query).

### RAG chunks ingestion benchmark (2026-02-24, `--find-max-depth 3`)

To enable larger bounded runs, the importer now supports depth control for artifact discovery:

- `--find-min-depth N`: Min depth for artifact discovery (default: 2)
- `--find-max-depth N`: Max depth for artifact discovery (default: 2)

Using `rag/chunks.jsonl` (PDF extraction chunks at depth 3) with `--find-max-depth 3 --limit 10`:

```bash
# Raw RAG chunks ingestion
./import-smiles-to-kb.sh \
  --source-root /LAB/@thesis/openwebui/prod_max \
  --file-name "chunks.jsonl" \
  --find-max-depth 3 \
  --kb-name "SMILES Raw RAG-Chunks Limit10 2026-02-24" \
  --limit 10 \
  --url http://localhost:3010

# Markdown RAG chunks ingestion
./import-smiles-to-kb.sh \
  --source-root /LAB/@thesis/openwebui/prod_max \
  --file-name "chunks.jsonl" \
  --find-max-depth 3 \
  --convert-markdown \
  --markdown-title "SMILES RAG Chunks" \
  --kb-name "SMILES Markdown RAG-Chunks Limit10 2026-02-24" \
  --limit 10 \
  --url http://localhost:3010
```

Results:

- Raw KB: `5c95c381-4465-4135-b568-2a1b05d12386` (10 seen, 7 attached, 3 duplicate-skipped)
- Markdown KB: `2a54163d-b888-4abd-82bb-2dd1b958e945` (10 seen, 7 attached, 3 duplicate-skipped)

Retrieval comparison (`k=3`, 4 queries):

| Query | Raw distances | Raw snippet type | MD distances | MD snippet type |
|-------|---------------|------------------|--------------|-----------------|
| mesembrine alkaloid structure | 0.960, 0.957, 0.956 | JSON objects | 0.959, 0.959, 0.959 | Molecule sections |
| Sceletium tortuosum mesembrenone | 0.950, 0.939, 0.938 | JSON objects | 0.952, 0.941, 0.939 | Molecule sections |
| SMILES notation for mesembrine | 0.939, 0.933, 0.930 | JSON objects | 0.944, 0.943, 0.936 | Molecule sections |
| joubertiamine class alkaloids | 0.950, 0.948, 0.947 | JSON objects | 0.950, 0.949, 0.947 | Molecule sections |

**Important caveat**: The default markdown converter (`scripts/jsonl_to_markdown_smiles.py`) is designed for molecule extraction JSONL (`molecules_high_confidence.jsonl`) with `smiles`, `backend`, `confidence` fields. When applied to `rag/chunks.jsonl` (PDF extraction chunks with `id`, `block_type`, `html` fields), molecule-specific fields render as `N/A`/`unknown`. For general RAG chunks, either:

1. Use raw ingestion (no `--convert-markdown`), or
2. Create a dedicated chunk-to-markdown converter for PDF extraction output.

Updated operator recommendation:

- Keep `smiles:publish-kb` (raw JSONL) as the default compatibility path.
- Prefer `smiles:publish-kb-markdown` for retrieval-facing KBs where human-readable grounding in returned chunks matters.
- Before changing the default, rerun on a truly larger corpus slice (more than 3 artifacts available in source) and evaluate query-level relevance judgments, not distance alone.

Markdown conversion controls in importer:

- `--convert-markdown`: convert each JSONL artifact to markdown before upload.
- `--markdown-file`: output markdown filename per paper.
- `--markdown-title`: top-level markdown heading.
- `--markdown-script`: custom converter path (default: `scripts/jsonl_to_markdown_smiles.py`).
