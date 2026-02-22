# Troubleshooting: OpenWebUI + CLIProxyAPI (Auth, Models, RAG/Vector DB)

This repository runs:
- OpenWebUI in Docker (`openwebui`)
- CLIProxyAPI in Docker (`cliproxyapi`) as the OpenAI-compatible upstream
- Optional host services (vLLM, SearXNG, PostgreSQL/pgvector)

When something looks "gone" in the UI (models, settings), the most common root cause is auth/session drift, not data loss.

## Quick triage checklist (2 minutes)

1. Service health:
```bash
./status.sh
docker compose ps
curl -fsS http://127.0.0.1:${WEBUI_PORT:-3000}/health >/dev/null && echo "openwebui health ok"
```

2. CLIProxyAPI models reachable from host:
```bash
./check-cliproxyapi.sh
```

3. OpenWebUI can route to CLIProxyAPI (end-to-end):
```bash
./test-openwebui-cliproxy-routing.sh
```

If those pass, your "models disappeared" problem is almost always browser auth state (cookies/localStorage) or user role state in OpenWebUI DB.

## Symptom: "Models disappeared" / model dropdown is empty

### What it usually means
- OpenWebUI UI is not authorized to call `GET /api/models`.
- The browser gets `401` and the UI renders an empty list.

### Confirm quickly
- Open browser DevTools, Network tab.
- Reload OpenWebUI.
- Look for `GET /api/models`.
If it returns `401`, you are not correctly authenticated.

### Fix (preferred)
1. Clear site data for OpenWebUI.
- URL: `http://127.0.0.1:${WEBUI_PORT:-3000}`
- Clear cookies and localStorage for that origin.

2. Hard refresh and sign in again.

### Fix (server-side, when you cannot sign in)
Promote/activate the user and optionally reset password (this edits SQLite DB inside the container and restarts OpenWebUI):
```bash
./openwebui-user-admin.sh --email you@example.com --role admin
./openwebui-user-admin.sh --email you@example.com --reset-password
```

Before edits, it auto-creates a DB backup under `./backups/`.

## Symptom: "Activation du compte en attente" / pending activation

### What it means
OpenWebUI thinks the account role is not active (often `pending`), even if the auth row exists.

### Fix
Promote the account role explicitly:
```bash
./openwebui-user-admin.sh --email you@example.com --role user
```

If you want admin UI access:
```bash
./openwebui-user-admin.sh --email you@example.com --role admin
```

Then clear site data and sign in again.

## Symptom: Auto-reconnects as "User" / cannot stay logged out

### What it means
- Your browser still has a valid session token/cookie for another account.
- Or `WEBUI_AUTH=false` is set and OpenWebUI operates in "single-user/no-auth" mode.

### Fix
1. Ensure `.env` has `WEBUI_AUTH=true`.
2. Redeploy:
```bash
docker compose up -d
```
3. Clear browser site data for the OpenWebUI origin and sign in again.

If you need to invalidate all sessions:
- Change `WEBUI_SECRET_KEY` in `.env` to a new strong random value (>= 32 bytes).
- Redeploy.
This forces all clients to re-authenticate (expected).

## Symptom: 502 from reverse proxy / OpenWebUI not reachable

### What it means
Usually the container is crash-looping or not healthy.

### Triage
```bash
docker compose ps
docker compose logs --tail 200 openwebui
```

### Common causes in this repo context
- Corrupted or inconsistent SQLite DB values from manual edits.
- OpenWebUI image tag mismatch and a migration edge case.

### Recovery steps
1. Backup DB:
```bash
./backup-openwebui-db.sh
```
2. Restart:
```bash
docker compose restart openwebui
```
3. If it keeps failing, restore from a known-good backup:
- Stop OpenWebUI.
- Copy a backup sqlite file back into the volume location.

Note: This repo defaults to a named Docker volume `openwebui_data` mapped to `/app/backend/data` inside the container.
The DB file is `/app/backend/data/webui.db`.

## Symptom: RAG "Knowledge Base" exists but vectors not where you expect

### Current expectation in this repo
Vector backend is chosen by `.env`:
- `VECTOR_DB=chroma` stores vectors inside OpenWebUI container data.
- `VECTOR_DB=pgvector` stores vectors in external PostgreSQL (host or separate container).

### Confirm which backend you are using
```bash
./status.sh
```

### When using `pgvector`
Required env vars:
- `VECTOR_DB=pgvector`
- `PGVECTOR_DB_URL=postgresql://...`

Common connectivity rule:
- If PostgreSQL runs on the host, from the container you must use `host.docker.internal`, not `localhost`.
- This repo config includes `extra_hosts: host.docker.internal:host-gateway` for that reason.

## Symptom: Web search / web_scraper fails with SSL certificate errors

### What it means
- The OpenWebUI container can't validate the HTTPS certificate chain for the target host.
- Common causes: corporate proxy / TLS interception, private CA, or a missing CA bundle in the container.

Typical error strings:
- `SSLError`
- `certificate verify failed`
- `unable to get local issuer certificate`

### Fix (recommended): mount a CA bundle + set `REQUESTS_CA_BUNDLE`
This repo mounts `./certs` into the OpenWebUI container at `/certs` (see `docker-compose.yml`).

1. Place a PEM bundle at `certs/ca-bundle.pem` (gitignored by default).
2. In `.env`, set:
```bash
REQUESTS_CA_BUNDLE=/certs/ca-bundle.pem
```
3. Restart OpenWebUI:
```bash
docker compose up -d --force-recreate openwebui
```

If a specific library requires it, you can also set:
```bash
SSL_CERT_FILE=/certs/ca-bundle.pem
```

### Temporary workaround (not recommended)
For the custom `web_scraper` tool, you can relax verification via its valves:
- `SSL_VERIFY=false`, or
- `ALLOW_INSECURE_SSL_FALLBACK=true`

This weakens TLS security and is only meant for short-lived debugging.

## Symptom: SearXNG web search fails from container

### What it means
- SearXNG is reachable from the host but not from the OpenWebUI container.
- Typically a host firewall blocking docker-to-host connections.

### Triage
The baseline test (`./test-rag.sh --baseline`) probes SearXNG from both the host and the container.
If the host probe works but the container probe fails, the host firewall is likely blocking docker-to-host connections on the SearXNG port (default 8888).

### Fix
1. Allow docker-to-host traffic on the SearXNG port:
```bash
# Example for ufw
sudo ufw allow from 172.16.0.0/12 to any port 8888
# Or for firewalld
sudo firewall-cmd --permanent --zone=trusted --add-port=8888/tcp
sudo firewall-cmd --reload
```

2. Or bind SearXNG to all interfaces and ensure Docker can reach it via `host.docker.internal`.

## Symptom: Tools not working

### What it means
- Tool plugins may not be registered, have broken endpoints, or have configuration drift.

### Triage
1. Check if tools are registered:
```bash
./audit-openwebui-plugins.sh
```

2. Repair tool registrations:
```bash
./repair-openwebui-tools.sh
```

3. Verify tool endpoints are reachable:
```bash
./test-openwebui-tools-endpoints.sh
```

### Common causes
- Tools were added manually but not registered in the OpenWebUI database.
- Tool container (e.g., Jupyter) is not healthy.
- Tool endpoint URLs changed after initial registration.

## Symptom: API key rotation needed

### What it means
- The CLIProxyAPI local key may be compromised, expired, or needs rotation for security compliance.

### Fix
1. Rotate the CLIProxyAPI local key:
```bash
./rotate-cliproxyapi-local-key.sh
```

2. Update any services or scripts that use the old key:
   - Check `.env` for references to the old key.
   - Update external integrations that call CLIProxyAPI.
   - Redeploy affected containers:
```bash
docker compose up -d --force-recreate
```

### After rotation
- Verify the new key works:
```bash
./check-cliproxyapi.sh
```

## Symptom: Code execution not working

### What it means
- The Jupyter code execution backend is unhealthy or misconfigured.

### Triage
1. Check Jupyter container health:
```bash
docker compose ps jupyter
docker compose logs --tail 50 jupyter
```

2. Verify `JUPYTER_TOKEN` is set correctly in `.env`.

3. Confirm OpenWebUI's `CODE_EXECUTION_JUPYTER_AUTH_TOKEN` matches `JUPYTER_TOKEN`:
   - Both must have the same value for authentication to succeed.

### Common fixes
- Restart the Jupyter container:
```bash
docker compose restart jupyter
```

- Regenerate and sync tokens if needed:
  1. Generate a new token.
  2. Set `JUPYTER_TOKEN` in `.env`.
  3. Set `CODE_EXECUTION_JUPYTER_AUTH_TOKEN` to the same value.
  4. Restart both containers:
```bash
docker compose up -d --force-recreate jupyter openwebui
```

## Symptom: Post-update validation failures

### What it means
- After an OpenWebUI update, something may be broken or misconfigured.

### Triage
1. Run the smoke test suite:
```bash
./test-update-smoke.sh
```

2. Check container logs for errors:
```bash
./logs.sh
docker compose logs --tail 200 openwebui
```

### Common post-update issues
- Database migration failures (check logs for migration errors).
- Configuration changes requiring `.env` updates.
- Plugin/tool compatibility with new OpenWebUI version.

### Recovery
- If smoke tests fail, check the specific test output for guidance.
- For DB issues, restore from backup:
```bash
./backup-openwebui-db.sh  # Create fresh backup first
# Then restore from a known-good backup
```

## Safety notes (do not skip)

- Never commit `.env`.
- Do not paste API keys, tokens, or DB URLs containing credentials into logs/issues.
- Before touching OpenWebUI SQLite DB, take a backup:
```bash
./backup-openwebui-db.sh
```
