# OpenWebUI + LiteLLM API Examples

Last updated: 2026-03-12 (UTC)

This document provides reproducible HTTP examples for the current stack:
- OpenWebUI in Docker
- LiteLLM as the primary OpenAI-compatible upstream (`http://127.0.0.1:4000/v1`)

Scope note:
- This file is a repo-specific examples guide, not the runtime source of truth.
- Validate ports, image tags, and default env values against `README.md`, `docs/ssot/stack.md`, `config/env/.env.example`, and `config/compose/docker-compose.yml`.

## Runtime environment

```bash
set -a
source ./.env
set +a

export OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
export LITELLM_BASE_URL="${LITELLM_BASE_URL:-http://127.0.0.1:4000/v1}"
export LITELLM_API_KEY="${LITELLM_API_KEY:-${OPENAI_API_KEY:-}}"
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

LiteLLM models:
```bash
curl -sS -f -H "Authorization: Bearer $LITELLM_API_KEY" \
  "$LITELLM_BASE_URL/models" | jq -e '.data | type == "array" and length > 0'
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

## Chat completion examples

LiteLLM direct call (replace `MODEL_ID` with one from `/v1/models`):
```bash
curl -sS -f -H "Authorization: Bearer $LITELLM_API_KEY" \
  -H "Content-Type: application/json" \
  "$LITELLM_BASE_URL/chat/completions" \
  -d '{
    "model": "MODEL_ID",
    "messages": [{"role": "user", "content": "Say ping."}],
    "max_tokens": 24
  }' | jq .
```

OpenWebUI proxied chat:
```bash
curl -sS -f -H "Authorization: Bearer $OPENWEBUI_API_KEY" \
  -H "Content-Type: application/json" \
  "$OPENWEBUI_URL/api/chat/completions" \
  -d '{
    "model": "MODEL_ID",
    "messages": [{"role": "user", "content": "Summarize retrieval-augmented generation in one sentence."}],
    "max_tokens": 64
  }' | jq .
```

## Regression test commands

LiteLLM-first baseline:
```bash
./test-rag.sh --baseline
./test-api.sh --baseline
```

Full API smoke suite:
```bash
./test-api.sh
```

Full RAG smoke suite:
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

- LiteLLM is the reference upstream for this repository.
