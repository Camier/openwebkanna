# CLIProxyAPI Sidecar

`cliproxyapi/` contains configuration and local runtime state for the optional legacy CLIProxyAPI sidecar.

Current status:

- LiteLLM is the primary OpenAI-compatible upstream for this repository.
- CLIProxyAPI remains supported only for legacy OAuth alias workflows.
- Do not treat this directory as the default path unless `CLIPROXYAPI_ENABLED=true`.

Directory layout:

- `config.yaml`: tracked base configuration and alias map.
- `config.local.yaml`: local overlay for real credentials and machine-specific values.
- `auth/`: OAuth material and local auth state.
- `logs/`: sidecar runtime logs.
- `vendor/`: vendored sidecar assets when needed.
- `upstream-release.txt`: tracked upstream version marker.

Handling rules:

- Never commit secrets or rotate credentials by editing tracked files directly.
- Use the root provider scripts directly (`start-cliproxyapi.sh`, `stop-cliproxyapi.sh`, `restart-cliproxyapi.sh`, `check-cliproxyapi.sh`, `configure-cliproxyapi-*.sh`).
- Keep documentation and `.env.example` aligned with its legacy/optional role.
