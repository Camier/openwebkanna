# Embedding Profiles Architecture

This runbook defines a reusable, production-ready embedding architecture for multiple KB use cases.

## Objective

Keep retrieval quality stable across heterogeneous knowledge bases by treating embeddings as a managed deployment surface (not a per-query ad hoc toggle).

## Core Principle

OpenWebUI uses one runtime query embedding model at a time.
All KB collections queried by that runtime must be indexed with the same model for reliable recall.

## Architecture Layers

1. Control plane (versioned profile registry)
- File: `config/embeddings/profiles.json`
- Stores named embedding profiles (model, engine, batch size, async settings)
- One `active_profile` tracked locally for operator visibility

2. Policy plane (KB bindings)
- File: `config/embeddings/kb-bindings.json`
- Maps each KB ID to a target profile
- Optional `kb_key` adds a short lane alias (`smiles`, `sceletium`, etc.)
- Used to prevent accidental cross-model retrieval drift

3. Runtime plane (OpenWebUI)
- Endpoint: `POST /api/v1/retrieval/embedding/update`
- Active query embedding model used by retrieval

4. Data plane (pgvector)
- Table: `document_chunk`
- Chunks carry `vmetadata.embedding_config` showing indexing model

## Operator CLI

Use:
```bash
./manage-openwebui-embedding-profiles.sh <command>
```

Commands:
- `list`: show profiles + runtime model
- `lanes`: show KB lane shortcuts from bindings (`kb_key -> kb/profile/model`)
- `prewarm --profile ID`: pre-download model in container
- `apply --profile ID [--prewarm]`: switch runtime model safely
- `bindings`: list KB->profile bindings
- `bind-kb --kb-id UUID --profile ID [--kb-name NAME] [--lane KEY] [--replace-lane] [--use-case TEXT]`
- `use-kb (--kb-id UUID | --lane KEY) [--prewarm] [--profile ID]`: one-shot KB activation (resolve->apply->diagnose)
- `diagnose [--kb-id UUID]`: detect runtime/index/binding mismatches

Notes:
- `--lane` is the recommended operator-facing lane key.
- `--kb-key` remains supported as a compatibility alias.
- `--replace-lane` intentionally transfers an existing lane binding to the new KB during reingest/cutover.
- `--use-case` is optional legacy metadata, not part of the daily switch flow.

## Single-Instance Workflow (SSOT)

This repo standard is now:
- one OpenWebUI instance,
- one runtime embedding model at a time,
- controlled switch per KB lane with `use-kb`.

1. Check runtime and lane map
```bash
./manage-openwebui-embedding-profiles.sh list
./manage-openwebui-embedding-profiles.sh lanes
```

2. Activate the Sceletium lane
```bash
./manage-openwebui-embedding-profiles.sh use-kb --lane sceletium --prewarm
```

3. If a KB has no lane key yet, bind it once
```bash
./manage-openwebui-embedding-profiles.sh bind-kb \
  --kb-id <KB_UUID> \
  --profile biomed_pubmed_msmarco \
  --lane sceletium \
  --kb-name "Sceletium Corpus"
```

4. Diagnose compatibility before querying
```bash
./manage-openwebui-embedding-profiles.sh diagnose
```

This pattern avoids manual UUID copy/paste and keeps the profile switch explicit.

- `biomed_pubmed_msmarco`
  - Model: `pritamdeka/S-PubMedBert-MS-MARCO`
  - Use for: pharmacology papers, mechanisms, biomedical prose queries

## Serving Patterns

Pattern A: Single OpenWebUI instance (serial profile switching)
- Recommended for this personal setup.
- Before querying a KB domain, switch runtime profile to matching index profile using `use-kb --lane <KEY>`.

Pattern B: Multiple OpenWebUI instances (recommended for concurrent domains)
- Run one instance per embedding profile.
- Route users/use-cases to the instance pinned to the matching profile.
- Avoids runtime switching contention.

## Migration Rules

- If changing profile for an existing KB, re-ingest/re-index documents for that KB.
- Do not assume old chunks remain semantically compatible after profile switch.
- Validate with:
```bash
./manage-openwebui-embedding-profiles.sh diagnose --kb-id <KB_UUID>
```

## Notes

- This architecture is intentionally model-agnostic and reusable across biomedical and chemistry KBs.
- Keep `profiles.json` and `kb-bindings.json` under version control for traceability.
