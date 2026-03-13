# Legacy Documentation Notes

This directory tracks migration notes and historical documentation decisions.

## Canonical Documentation Roots

- Runbooks: `docs/runbooks/`
- Status docs: `docs/status/`
- Guides: `docs/guides/`
- Plans: `docs/plans/`
- Reviews: `docs/reviews/`
- OpenWebUI reference snapshots: `docs/reference/openwebui/`

## Current Rule

- New canonical docs live under `docs/`, not at the repo root.
- The old root markdown compatibility stubs were removed during the consolidation pass.
- Repo-root markdown should be limited to true top-level docs such as `README.md`, `AGENTS.md`, `SECURITY.md`, `API_EXAMPLES.md`, and agent-facing notes.

## Migration Note

If you find an old path in local notes or shell history, replace it with the matching canonical file under `docs/`.

Historical planning docs under `docs/plans/` may still be useful for context, but they are not active source-of-truth documents.
