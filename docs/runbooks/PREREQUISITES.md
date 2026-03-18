# OpenWebUI + LiteLLM Deployment Prerequisites

Last verified: 2026-03-13 (UTC)

This checklist reflects the current default deployment mode in this repo:
- LiteLLM as primary OpenAI-compatible upstream
- Docker-managed OpenWebUI
- optional sidecars disabled by default

Use this checklist before the first deploy or before rebuilding a machine from scratch.
If a required item fails, fix it before running `./deploy.sh`.

## Required

1. Docker daemon available
```bash
docker info >/dev/null
```

1.a Host cgroup-runtime sanity check
```bash
docker run --rm --network none busybox:latest true
```
This validates the host/runtime can start a container before full deploy. If this fails with a cgroup-scoped error, fix Docker cgroup delegation first (systemd + cgroup driver mismatch is the typical root cause in this environment).

2. Docker Compose plugin available
```bash
docker compose version
```

3. Core CLI tooling available
```bash
command -v curl
command -v jq
```

4. Repo env file prepared
```bash
cp config/env/.env.example .env
```

5. Required `.env` defaults present
```bash
for key in \
  WEBUI_SECRET_KEY \
  JUPYTER_TOKEN \
  CODE_EXECUTION_JUPYTER_AUTH_TOKEN \
  CODE_INTERPRETER_JUPYTER_AUTH_TOKEN \
  POSTGRES_PASSWORD \
  OPENAI_API_BASE_URL \
  OPENAI_API_BASE_URLS \
  OPENAI_API_KEY
do
  grep -Eq "^${key}=.+$" .env || echo "missing: ${key}"
done
```

6. Jupyter auth tokens aligned
```bash
test "$(grep '^JUPYTER_TOKEN=' .env | cut -d= -f2-)" = "$(grep '^CODE_EXECUTION_JUPYTER_AUTH_TOKEN=' .env | cut -d= -f2-)" &&
test "$(grep '^JUPYTER_TOKEN=' .env | cut -d= -f2-)" = "$(grep '^CODE_INTERPRETER_JUPYTER_AUTH_TOKEN=' .env | cut -d= -f2-)"
```

7. Secret and credential placeholders replaced
```bash
test "$(grep '^WEBUI_SECRET_KEY=' .env | cut -d= -f2- | wc -c)" -ge 33 &&
! grep -Eq '^POSTGRES_PASSWORD=<|^OPENAI_API_KEY=<' .env
```

Expected values for LiteLLM-first mode:
- `WEBUI_SECRET_KEY=<stable-random-secret-at-least-32-chars>`
- `JUPYTER_TOKEN=<random-jupyter-token>`
- `CODE_EXECUTION_JUPYTER_AUTH_TOKEN=<same-as-JUPYTER_TOKEN>`
- `CODE_INTERPRETER_JUPYTER_AUTH_TOKEN=<same-as-JUPYTER_TOKEN>`
- `POSTGRES_PASSWORD=<strong-local-password>`
- `OPENAI_API_BASE_URL=http://host.docker.internal:4000/v1`
- `OPENAI_API_BASE_URLS=http://host.docker.internal:4000/v1`
- `OPENAI_API_KEY=<litellm-master-key>`
- `OPEN_TERMINAL_ENABLED=false`
- `INDIGO_SERVICE_ENABLED=false`

For the complete committed baseline and image defaults, use `config/env/.env.example`.

Advanced optional prerequisites are not required for the baseline runtime.

## Pre-deploy verification

```bash
docker compose config >/dev/null
bash -n deploy.sh status.sh test-rag.sh test-api.sh
```

Expected outcome:
- no missing required env keys
- Jupyter auth tokens match
- placeholder secrets have been replaced with real values
- `docker compose config` renders successfully
- the operator scripts parse cleanly

## Handoff to setup

If the checks above pass, continue with `docs/runbooks/SETUP.md`.
