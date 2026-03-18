# local/

Last updated: 2026-03-13 (UTC)

This subtree groups helper assets that are intentionally local to one machine or clone.
It keeps the repo root focused on the canonical operator scripts and versioned config/docs.

## Layout

- `local/bin/`: local sidecar binaries and wrappers
- `local/plugins/`: local tool payloads used by optional integrations
- `local/entities/`: separate entity-maintenance workspace and its run outputs

## Rules

- Treat everything under `local/` as local-only by default.
- Root scripts may use files from `local/`, but the baseline deployment does not require this subtree.
- `local/README.md` and `local/entities/README.md` are the tracked orientation points here; the rest of `local/` stays ignored by default.
