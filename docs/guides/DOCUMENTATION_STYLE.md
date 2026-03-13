# Documentation Style Guide

Last updated: 2026-03-13 (UTC)

Use this guide when you are writing, reviewing, or restructuring active technical documentation in this repository.
It defines the minimum editorial contract for operator docs, not a general writing theory.

## Core Rules

- Keep one source of truth per fact. If code, config, and docs drift, update the docs to point at the real owner instead of copying stale details.
- Keep one primary job per document. A README is a front door, a runbook is a procedure, a map is a locator, and the SSOT defines current-state runtime truth.
- Write task-first. The reader should know within the first screen whether the document is for setup, day-2 operations, troubleshooting, config edits, or architecture truth.
- Prefer commands and observable outcomes over abstract prose.
- Keep optional, legacy, or archived flows out of baseline runbooks unless the section is explicitly marked as optional.

## Required Structure

For active docs, prefer this shape:

1. Title
2. Freshness marker
3. Scope line or usage block
4. Task-oriented body
5. Handoff, exit criteria, or validation loop when applicable

Timestamp rules:

- Use `Last updated:` for guides, maps, runbooks, and SSOT pages.
- Use `Last verified:` for checklists that primarily assert current readiness or compatibility.

Scope rules:

- Start with `Use this file when...` or `Use this runbook when...` for front doors, maps, and procedures.
- Add `Do not use this file...` only when it prevents a common routing mistake.
- Do not repeat multiple near-identical routing blocks in the same document.

## Writing Rules

- Use direct imperative wording for procedures.
- Keep paragraphs short and information-dense.
- Prefer bullets for operator choices, prerequisites, and expected outcomes.
- Use exact filenames, commands, ports, and env vars.
- Use absolute facts with concrete paths rather than vague references like "the config file" or "the service script".
- Use exact dates such as `2026-03-13`, not relative wording like "recently" or "today", when freshness matters.

Avoid:

- restating the same routing guidance in multiple sections of the same file
- copying full env baselines into multiple docs when `config/env/.env.example` is the canonical owner
- mixing first-deploy, day-2 operations, troubleshooting, and legacy flows in one procedure unless the split is explicit
- long narrative explanations before the first actionable step

## Command And Snippet Rules

- Prefer the canonical root script surface before raw low-level commands.
- Use `docker compose ...` instead of hard-coded container names when low-level Docker commands are still needed.
- Mark raw `curl`, `psql`, `docker exec`, or SQL snippets as advanced debugging when they are not the default operator path.
- Keep snippets copy-paste safe. If a placeholder is required, name it explicitly in the snippet.
- Keep validation commands close to the task they validate.

## Source-Of-Truth Boundaries

- `README.md` owns the operator front door.
- `docs/ssot/stack.md` owns current-state runtime truth.
- `docs/REPO_MAP.md` owns repository shape and path ownership.
- `config/README.md` owns the canonical config edit surface.
- `docs/runbooks/*.md` own procedures.
- `docs/status/*`, `docs/plans/*`, and `docs/reference/openwebui/*` are not canonical desired-state runtime truth.

## Change Checklist

After documentation edits:

```bash
./scripts/check-doc-consistency.sh
```

After docs that change or describe runtime behavior materially:

```bash
docker compose config -q
./status.sh
./test-rag.sh --baseline
./test-api.sh --baseline
```
