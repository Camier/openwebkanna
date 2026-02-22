# OpenWebUI + CLIProxyAPI API Examples

Last updated: 2026-02-12 (UTC)

This document provides reproducible HTTP examples for the current stack:
- OpenWebUI in Docker
- CLIProxyAPI in Docker (`http://127.0.0.1:8317` externally)
- OpenWebUI upstream mapped to `http://cliproxyapi:8317/v1` internally

## Runtime environment

```bash
set -a
source ./.env
set +a

export OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
export CLIPROXYAPI_BASE_URL="${CLIPROXYAPI_BASE_URL:-http://127.0.0.1:8317}"
export CLIPROXYAPI_API_KEY="${CLIPROXYAPI_API_KEY:-layra-cliproxyapi-key}"
export OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-$API_KEY}"
```

## Authentication

### Sign In (Get JWT Token)

```bash
curl -sS -X POST "${OPENWEBUI_URL}/api/v1/auths/signin" \
  -H "Content-Type: application/json" \
  -d '{
    "email": "user@example.com",
    "password": "your-password"
  }' | jq .
```

**Response:**
```json
{
  "token": "eyJhbGciOiJIUzI1NiIs...",
  "id": "user-uuid",
  "email": "user@example.com",
  "name": "User Name",
  "role": "user",
  "profile_image_url": "/user.png"
}
```

### Sign Up (Create Account)

```bash
curl -sS -X POST "${OPENWEBUI_URL}/api/v1/auths/signup" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "New User",
    "email": "newuser@example.com",
    "password": "secure-password"
  }' | jq .
```

**Response:** Same as sign in (201 Created).

### Get Current User Info

```bash
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/users/user/info" | jq .
```

**Response:**
```json
{
  "id": "user-uuid",
  "email": "user@example.com",
  "name": "User Name",
  "role": "user",
  "profile_image_url": "/user.png"
}
```

### API Key Authentication

API keys can be used instead of JWT tokens. Generate API keys in Settings > Account.

```bash
curl -sS -H "Authorization: Bearer sk-xxxxxxxxxxxxxxxxxxxxxxxx" \
  "${OPENWEBUI_URL}/api/models" | jq .
```

**Requirements:**
1. Admin must enable "Enable API Keys" in Admin Panel > Settings
2. User must have `features.api_keys` permission

### Token Expiration

JWT tokens expire based on `JWT_EXPIRES_IN` environment variable:
- `"-1"`: Never expire (not recommended for production)
- `"4w"`: 4 weeks
- `"24h"`: 24 hours
- `"30m"`: 30 minutes

When a token expires, you'll receive:
```json
{
  "detail": "Invalid token or expired token."
}
```
With HTTP 401 status code.

## Health checks

CLIProxyAPI root:
```bash
curl -sS -f "$CLIPROXYAPI_BASE_URL/" | jq .
```

CLIProxyAPI models:
```bash
curl -sS -f -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
  "$CLIPROXYAPI_BASE_URL/v1/models" | jq -e '.data | type == "array" and length > 0'
```

OpenWebUI health:
```bash
curl -sS -f "$OPENWEBUI_URL/health" | jq .
```

OpenWebUI models:
```bash
curl -sS -f -H "Authorization: Bearer $OPENWEBUI_API_KEY" \
  "$OPENWEBUI_URL/api/models" | jq .
```

## Required alias presence

```bash
curl -sS -f -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
  "$CLIPROXYAPI_BASE_URL/v1/models" \
| jq -e '
  [.data[].id] as $ids
  | ["openai-codex","qwen-cli","kimi-cli"]
  | all(. as $name | $ids | index($name))
'
```

## Chat completion examples

CLIProxyAPI direct call (`openai-codex`):
```bash
curl -sS -f -H "Authorization: Bearer $CLIPROXYAPI_API_KEY" \
  -H "Content-Type: application/json" \
  "$CLIPROXYAPI_BASE_URL/v1/chat/completions" \
  -d '{
    "model": "openai-codex",
    "messages": [{"role": "user", "content": "Say ping."}],
    "max_tokens": 24
  }' | jq .
```

OpenWebUI proxied chat (`qwen-cli`):
```bash
curl -sS -f -H "Authorization: Bearer $OPENWEBUI_API_KEY" \
  -H "Content-Type: application/json" \
  "$OPENWEBUI_URL/api/chat/completions" \
  -d '{
    "model": "qwen-cli",
    "messages": [{"role": "user", "content": "Summarize retrieval-augmented generation in one sentence."}],
    "max_tokens": 64
  }' | jq .
```

## Wrapper examples (`cli-proxy-api.sh`)

The wrapper command namespace keeps `vllm` as the OpenAI-compatible endpoint label.
In this repo, `vllm` usually resolves to CLIProxyAPI through `.env` (`OPENAI_API_BASE_URL`).

List models via wrapper against configured OpenAI-compatible endpoint:
```bash
./cli-proxy-api.sh --raw models vllm | jq .
```

Force wrapper to CLIProxyAPI explicitly:
```bash
./cli-proxy-api.sh --raw --url "$CLIPROXYAPI_BASE_URL" models vllm | jq .
```

Health all services:
```bash
./cli-proxy-api.sh health all
```

## Regression test commands

OAuth alias regression:
```bash
./test-cliproxyapi-oauth.sh
```

OpenWebUI routing regression:
```bash
./test-openwebui-cliproxy-routing.sh
```

LLM council evaluation (candidate + judge pass):
```bash
./llm-council.sh \
  --models "glm-5 minimax/chat-elite" \
  --judges "glm-5 minimax/chat-elite" \
  --prompt "Summarize retrieval quality risks in one paragraph."
```

LLM council with stricter judge retries (useful for reasoning-heavy models):
```bash
./llm-council.sh \
  --models "glm-5 minimax/chat-elite" \
  --judges "glm-5 minimax/chat-elite" \
  --prompt "Reply with exactly council-ok" \
  --judge-max-tokens 256 --judge-retries 3 --judge-force-max-tokens 96
```

LLM council anti-position-bias mode (A/B then B/A):
```bash
./llm-council.sh \
  --models "glm-5 minimax/chat-elite" \
  --judges "glm-5 minimax/chat-elite" \
  --prompt "Pick the stronger incident response plan." \
  --position-swap
```

LLM council with strict JSON judge output:
```bash
./llm-council.sh \
  --models "glm-5 minimax/chat-elite" \
  --judges "glm-5 minimax/chat-elite" \
  --prompt "Reply with exactly council-ok" \
  --judge-output-format json --judge-retries 3
```

API smoke suite:
```bash
./test-api.sh
```

RAG smoke suite:
```bash
./test-rag.sh
```

## Additional API Examples

### File Download

```bash
# Download file content
curl -sS -f -H "Authorization: Bearer $OPENWEBUI_API_KEY" \
  "${OPENWEBUI_URL}/api/v1/files/{file_id}/content" \
  -o downloaded_file.pdf
```

### Chat Management

```bash
# Create new chat
curl -sS -X POST "${OPENWEBUI_URL}/api/v1/chats/new" \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  -H "Content-Type: application/json" \
  -d '{
    "chat": {
      "title": "New Research Chat",
      "models": ["llama3.2"],
      "messages": [],
      "history": {}
    }
  }' | jq .

# List user chats
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/chats" | jq .

# Get specific chat
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/chats/{chat_id}" | jq .

# Update chat
curl -sS -X POST "${OPENWEBUI_URL}/api/v1/chats/{chat_id}" \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  -H "Content-Type: application/json" \
  -d '{
    "title": "Updated Chat Title",
    "messages": [{"role": "user", "content": "Hello"}],
    "history": {"messages": {"message-id": {"role": "user", "content": "Hello"}}}
  }' | jq .

# Delete chat
curl -sS -X DELETE -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/chats/{chat_id}"
```

### Tools and Functions

```bash
# List tools
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/tools/" | jq .

# Get tool details
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/tools/id/{tool_id}" | jq .

# List functions
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/functions/" | jq .
```

### Ollama Proxy

```bash
# Ollama generate
curl -sS -X POST "${OPENWEBUI_URL}/ollama/api/generate" \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.2",
    "prompt": "Why is the sky blue?",
    "stream": false
  }' | jq .

# Ollama chat
curl -sS -X POST "${OPENWEBUI_URL}/ollama/api/chat" \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.2",
    "messages": [
      {"role": "user", "content": "Hello!"}
    ],
    "stream": false
  }' | jq .

# List Ollama models
curl -sS -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/ollama/api/tags" | jq .

# Generate embeddings via Ollama
curl -sS -X POST "${OPENWEBUI_URL}/ollama/api/embed" \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.2",
    "input": ["Text to embed", "Another text"]
  }' | jq .
```

## Error Handling Examples

```bash
# Example: Handle 404 error for non-existent file
curl -sS -w "\nHTTP Status: %{http_code}\n" \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  "${OPENWEBUI_URL}/api/v1/files/nonexistent-id/content" \
  -o /dev/null 2>&1 | tail -2

# Example: Handle 401 unauthorized (no token)
curl -sS -w "\nHTTP Status: %{http_code}\n" \
  "${OPENWEBUI_URL}/api/v1/chats" \
  2>&1 | tail -2

# Example: Handle 400 bad request (malformed JSON)
curl -sS -X POST "${OPENWEBUI_URL}/api/v1/chats/new" \
  -H "Authorization: Bearer ${OPENWEBUI_API_KEY}" \
  -H "Content-Type: application/json" \
  -d 'invalid json' | jq .
```

## Notes

- `CLIPROXYAPI_ENABLED=false` is intentionally treated as a failed integration state.
- Manual OAuth is supported; auth material is persisted under `cliproxyapi/auth/`.
- If CLIProxyAPI is healthy, deployment keeps vLLM stopped by default.
