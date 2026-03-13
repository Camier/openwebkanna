# Local Entities Workspace

This directory isolates entity-oriented extraction maintenance from the main RAG runtime.
It now lives under `local/` so the repo root stays focused on baseline operator surfaces.

## Scope
- Entity data source (canonical): `data/extractions/<paper>/entities.json` and `data/extractions/<paper>/.errors.jsonl`
- Entity maintenance pipelines: `local/entities/pipelines/`
- Entity run outputs and manifests: `local/entities/artifacts/`

## Layout
```text
local/entities/
├── README.md
├── docs/
│   └── ARCHITECTURE.md
├── data/
│   └── README.md
├── pipelines/
│   ├── remediate-extraction-error-logs.py
│   ├── retry-entity-extraction.py
│   └── run-entities-quality-audit.sh
├── artifacts/
│   ├── error-log-remediation/
│   ├── retry-entity-extraction/
│   └── quality-assurance/
│       ├── extractions-base-audit/
│       ├── extractions-outlier-qa/
│       ├── post-error-log-remediation-audit/
│       └── post-error-log-remediation-gate/
```

## Naming Conventions
- Pipelines: verb-first snake-case with clear intent, e.g. `retry-entity-extraction.py`.
- Artifact folders: kebab-case by workflow, then timestamp for runs, e.g. `error-log-remediation/20260305T110434Z`.
- QA snapshots: grouped under `local/entities/artifacts/quality-assurance/` by audit family.

## Pipeline Order
1. Build unresolved manifests from extraction failures:
```bash
python local/entities/pipelines/remediate-extraction-error-logs.py \
  --root data/extractions \
  --prune-covered
```

2. Retry unresolved chunk extraction from latest manifest:
```bash
python local/entities/pipelines/retry-entity-extraction.py \
  --model cliproxy/gpt-5-codex-mini \
  --max-chunks 40
```

3. Apply resolved retries to `entities.json` and prune `.errors.jsonl`:
```bash
python local/entities/pipelines/retry-entity-extraction.py \
  --model cliproxy/gpt-5-codex-mini \
  --max-chunks 40 \
  --apply
```

4. Run entity QA snapshot into dedicated artifact path:
```bash
local/entities/pipelines/run-entities-quality-audit.sh \
  --name extractions-base-audit
```

## Notes
- This workspace is prepared for future graph-oriented work.
- Current non-graph RAG remains unchanged.
