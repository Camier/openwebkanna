# Documentation Map

This repository uses structured documentation tiers:

- `runbooks`: operational guides and workflows used to run/maintain the stack.
- `reference`: upstream-derived technical references and deep-dive docs. Treat these as snapshot material, not the runtime source of truth for this repository.
- `guides`: topic-focused guidance and curated how-to material.
- `plans`: implementation plans and migration strategies.
- `status`: observed local runtime snapshots and rollout records for this machine.
- `reviews`: audit and technical review reports.
- `legacy`: compatibility and migration notes for old paths/layouts.

## Runtime Source Of Truth

Read these in this order so the docs do not compete with each other:

- `README.md`: operator front door, quick-start path, and high-signal commands.
- `docs/ssot/stack.md`: canonical runtime topology, service registry, critical paths, and drift notes.
- `docs/ssot/stack.yaml`: machine-readable stack registry.
- `docs/REPO_MAP.md`: repository layout, entrypoints, config roots, and port inventory.
- `config/README.md`: canonical config edit surface and compatibility-copy ownership.
- `docs/runbooks/*.md`: procedure-level operating detail.
- `config/env/.env.example` and `config/compose/docker-compose.yml`: committed runtime defaults and exact service wiring.

Do not use these as canonical desired-state truth:

- `docs/reference/openwebui/*`: upstream/reference snapshots for background only.
- `docs/plans/*`: planning artifacts that may describe intended or historical changes rather than the current state.
- `docs/status/*`: observed local runtime snapshots. Use them to understand what this specific machine was running at a specific date, not to redefine supported topology or committed defaults.

## Local Directory Maps

Use these when you need subtree ownership and navigation without scanning the whole repo:

- `config/README.md`
- `scripts/README.md`
- `artifacts/README.md`
- `cliproxyapi/README.md`
- `research/README.md`

## Standalone Root Docs

These live at the repo root on purpose, but they are not runtime SSOT files:

- `API_EXAMPLES.md`: repo-specific HTTP/API examples
- `MCP_INTEGRATION_GUIDE.md`: MCP setup and integration guide
- `SECURITY.md`: security policy and scan posture
- `CLAUDE.md`, `GEMINI.md`: agent-facing quick entrypoints

## Current Doc Surfaces

- `docs/reference/openwebui/`
  Use [docs/reference/openwebui/README.md](reference/openwebui/README.md) as the entrypoint for upstream/reference snapshots instead of maintaining a second index here.
- `docs/runbooks/`
  Active operator procedures live here: setup, operations, troubleshooting, prerequisites, and embedding-profile workflows.
- `docs/guides/`
  Curated topic guides that are narrower than the runbooks.
- `docs/plans/`
  Planning artifacts and migration notes that may be historical or aspirational.
- `docs/status/`
  Point-in-time local runtime observations, not desired-state topology.
- `docs/reviews/`
  Audit and technical review outputs.
- `docs/legacy/`
  Compatibility and migration notes for retired layouts or paths.
