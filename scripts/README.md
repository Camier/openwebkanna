# Utility Scripts

`scripts/` contains secondary utilities, audits, and maintenance helpers. It is not the primary operator entrypoint surface.

Use these first for daily operations:

- Root `*.sh` commands for the primary operator workflows.

Use `scripts/` when you need narrower support tooling:

- Quality and hygiene: `check-doc-consistency.sh`, `repo-hygiene.sh`, `cleanup-empty-dirs.sh`
- Evaluation: `eval/run-retrieval-eval.sh`, `eval/run-chemical-ocsr-eval.sh` (`run.json`, `disagreements.json`, `report.md`), `eval/run-chemical-ocsr-disagreement-analysis.sh` (`disagreement_analysis.json`, `disagreement_analysis.md`), `eval/run-chemical-ocsr-fusion-analysis.sh` (`fusion_candidates.json`, `fusion_report.md`, `review_queue.json`, `review_queue.md`, `review_queue.csv`, `review_queue.jsonl`, `review_queue/*.json`, `review_queue/*.md`, `review_queue/*.csv`, `review_queue/*.jsonl`), `eval/run-chemical-ocsr-adjudication-analysis.sh` (`adjudication_manual_*` and `adjudication_defaulted_*` summary/report/accepted/rejected/pending exports plus `accepted_catalog` JSONL/CSV with compact provenance)
- Data and corpus utilities: `audit-data-quality.sh`, `curate-corpus-chunks.sh`, `data-quality-gate.sh`
- Thesis export and backup helpers: `export-thesis-chats.sh`, `thesis-backup.sh`
- Image utilities: `embed_images_in_chunks.py`, `blip_caption.py`, `rag/build-multimodal-index.sh`, `rag/render-multimodal-answer.sh`, `rag/render-multimodal-answer-v2.sh`, `rag/extract-chemical-smiles.sh` (`MolDetv2` detector by default, `--backend molgrapher|molscribe`)
- Maintenance and monitoring: `audit-dependencies.sh`, `check-image-versions.sh`
- RAG/admin sync helpers: `scripts/admin/sync-openwebui-web-search-config.sh`, `scripts/admin/sync-openwebui-retrieval-config.sh`

Conventions:

- Prefer shell scripts here for repo-local maintenance and validation.
- Prefer Python helpers here for deterministic data processing utilities.
- Keep these tools callable directly from the repo root without additional path assumptions.
