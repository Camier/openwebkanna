# Generated Artifacts

`artifacts/` stores generated outputs from audits, extraction runs, repo scans, and evaluation harnesses.

Current producers include:

- `data-quality/`: curation and corpus quality reports.
- `extraction-pipeline-operator/`: operator logs and delta reports for extraction waves.
- `repo-cartographer/`: repository scan logs used to refresh architecture docs.
- `retrieval-evaluation-harness*/`: retrieval benchmark outputs for CI, live runs, and local harness executions.

Usage rules:

- Treat this directory as generated evidence, not as source code or runtime configuration.
- Keep stable report formats when other scripts consume them.
- Avoid hard-coding artifact filenames in docs unless the path is part of an explicit workflow.
- Git policy: generated contents under `artifacts/` are ignored; keep only this README tracked.
