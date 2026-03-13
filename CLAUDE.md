# CLAUDE.md

Last updated: 2026-03-12 (UTC)

This file is a quick agent-facing entrypoint, not the repository source of truth.

## Read These First

- [README.md](README.md) for current deployment posture, baseline workflow, and operator-facing commands
- [AGENTS.md](AGENTS.md) for the short agent-routing note and local guardrails
- [docs/REPO_MAP.md](docs/REPO_MAP.md) for current topology, services, ports, and ownership
- [docs/runbooks/OPERATIONS.md](docs/runbooks/OPERATIONS.md) for day-to-day commands

## Minimal Repo Posture

- LiteLLM-first deployment
- OpenWebUI runs in Docker
- CLIProxyAPI is optional legacy sidecar only
- Use root scripts for operator workflows

## Guardrails

- Keep integration checks real: use actual `/v1/models` and real chat-completion calls
- Keep CLIProxyAPI disabled unless a legacy OAuth alias workflow explicitly needs it
- Do not treat this file as canonical for versions, ports, or compose topology
