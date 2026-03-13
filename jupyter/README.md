# Jupyter Compatibility Copy

This directory is a root compatibility surface for operators.

Canonical edit path:

- `config/jupyter/jupyter_server_config.py`

Current compatibility copy in this directory:

- `jupyter_server_config.py`

Rule:

- Edit `config/jupyter/` first.
- Keep this root copy aligned because compose and operator workflows still expose `jupyter/` at the repo root.
