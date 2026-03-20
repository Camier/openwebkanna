# Generated Artifacts

`artifacts/` stores generated outputs from audits, extraction runs, repo scans, and evaluation harnesses.

Usage rules:

- Treat this directory as generated evidence, not as source code or runtime configuration.
- Keep stable report formats when other scripts consume them.
- Avoid hard-coding artifact filenames in docs unless the path is part of an explicit workflow.
- Git policy: generated contents under `artifacts/` are ignored; keep only this README tracked.

## Producers

### Evaluation

- `retrieval-evaluation-harness/`: Local retrieval benchmark harness.
- `retrieval-evaluation-harness-ci/`: CI retrieval benchmark outputs.
- `retrieval-evaluation-harness-live/`: Live production retrieval benchmark outputs.

### RAG and Multimodal

- `multimodal-answer/`: Rendered multimodal RAG answers (v1).
- `multimodal-answer-v2/`: Rendered multimodal RAG answers (v2, CLIP-backed).
- `multimodal-index/`: CLIP-backed multimodal figure index build outputs.
- `rag/`: RAG materialization outputs, evidence manifests, and `migrate_rag_evidence` reports.

### Quality and Extraction

- `data-quality/`: Corpus quality reports, curation audits, and resolution records.
- `data-quality-guardian/`: Live quality monitoring dashboards and issue logs.
- `extraction-pipeline-operator/`: Extraction wave operator logs, delta reports, and manual review records.

### Repository Maintenance

- `dependency-env-doctor/`: Python dependency version snapshots for environment auditing.
- `repo-cartographer/`: Repository scan logs used to refresh architecture docs.
- `repo-hygiene/`: Hygiene run manifests and quarantined stale/rotten artifacts.
