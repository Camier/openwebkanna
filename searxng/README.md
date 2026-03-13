# SearXNG Compatibility Copy

This directory is a root compatibility surface for operators.

Canonical edit paths:

- `config/searxng/settings.yml`
- `config/searxng/limiter.toml`

Current compatibility copies in this directory:

- `settings.yml`
- `limiter.toml`

Rule:

- Edit `config/searxng/` first.
- Keep this root copy aligned because compose and operator workflows still expose `searxng/` at the repo root.
