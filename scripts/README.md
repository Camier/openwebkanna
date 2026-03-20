# Utility Scripts

`scripts/` contains secondary utilities, audits, and maintenance helpers. It is not the primary operator entrypoint surface.

Use these first for daily operations:

- Root `*.sh` commands for the primary operator workflows.

Use `scripts/` when you need narrower support tooling:

- Quality and hygiene: `check-doc-consistency.sh`, `repo-hygiene.sh`, `cleanup-empty-dirs.sh`
- Evaluation: `eval/run-retrieval-eval.sh`
- Data and corpus utilities: `audit-data-quality.sh`, `curate-corpus-chunks.sh`, `data-quality-gate.sh`
- Thesis export and backup helpers: `export-thesis-chats.sh`, `thesis-backup.sh`
- Image utilities: `embed_images_in_chunks.py`, `blip_caption.py`, `rag/build_multimodal_index.py`, `rag/render-multimodal-answer.sh`, `rag/render-multimodal-answer-v2.sh`
- Maintenance and monitoring: `check-image-versions.sh`
- RAG/admin sync helpers: `scripts/admin/sync-openwebui-web-search-config.sh`, `scripts/admin/sync-openwebui-retrieval-config.sh`

Conventions:

- Prefer shell scripts here for repo-local maintenance and validation.
- Prefer Python helpers here for deterministic data processing utilities.
- Keep these tools callable directly from the repo root without additional path assumptions.
