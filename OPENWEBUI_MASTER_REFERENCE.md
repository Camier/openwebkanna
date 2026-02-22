# OpenWebUI Master Technical Reference

> **OpenWebUI Version:** v0.8.3
> **Repository:** CLIProxyAPI-integrated deployment
> **Last Updated:** 2026-02-18

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Architecture Overview](#2-architecture-overview)
3. [Environment Variables Complete Reference](#3-environment-variables-complete-reference)
4. [Tools & Functions Development](#4-tools--functions-development)
5. [RAG System Deep Dive](#5-rag-system-deep-dive)
6. [API Endpoints Reference](#6-api-endpoints-reference)
7. [Pipelines & Extensibility](#7-pipelines--extensibility)
8. [Quick Start Guides](#8-quick-start-guides)

---

## Related Documentation

| Document | Purpose | Size |
|----------|---------|------|
| [Environment Variables](openwebui_env_reference.md) | Complete env var reference (150+ variables) | 1120 lines |
| [Tools & Functions Guide](openwebui-tools-functions-guide.md) | Development guide with examples | 1365 lines |
| [RAG Technical Reference](openwebui_rag_technical_reference.md) | RAG system deep dive | 923 lines |
| [Pipelines Guide](OPENWEBUI_PIPELINES_GUIDE.md) | Pipelines deployment and development | 877 lines |
| [API Examples](API_EXAMPLES.md) | Repository-specific API examples | 166 lines |
| [README.md](README.md) | Main project documentation | 448 lines |
| [AGENTS.md](AGENTS.md) | Repository map for AI agents | 546 lines |
| [TROUBLESHOOTING.md](TROUBLESHOOTING.md) | Common issues and solutions | 326 lines |

> **Note:** This Master Reference provides an overview. See the specialized
> documents above for detailed technical information.

---

## 1. Executive Summary

OpenWebUI is a feature-rich, self-hosted AI interface supporting multiple LLM backends. This reference covers:

| Area | Key Capabilities |
|------|------------------|
| **Configuration** | 400+ environment variables across 20 categories |
| **Extensibility** | Tools, Functions, Pipelines, Filters |
| **RAG** | Multi-vector DB support, hybrid search, web integration |
| **API** | OpenAI-compatible + native endpoints |
| **Integrations** | Ollama, OpenAI, Gemini, 20+ web search engines |

---

## 2. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         OpenWebUI Architecture                          │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌──────────────┐     ┌──────────────┐     ┌──────────────────────┐    │
│  │   Web UI     │────▶│   Backend    │────▶│   LLM Providers      │    │
│  │   (Svelte)   │     │   (FastAPI)  │     │   - Ollama           │    │
│  └──────────────┘     └──────┬───────┘     │   - OpenAI           │    │
│                              │             │   - Gemini           │    │
│                              ▼             │   - Azure            │    │
│  ┌──────────────┐     ┌──────────────┐     └──────────────────────┘    │
│  │   Vector DB  │◄────│   RAG Engine │                               │
│  │   (Chroma/   │     │   (Embedding │     ┌──────────────────────┐    │
│  │    pgvector) │     │    + Search) │────▶│   Tools/Functions    │    │
│  └──────────────┘     └──────────────┘     │   - Web Search       │    │
│                                            │   - Code Execution   │    │
│  ┌──────────────┐     ┌──────────────┐     │   - Image Gen        │    │
│  │   Database   │◄────│   Auth/      │     │   - Custom Tools     │    │
│  │(SQLite/PG/MySQL)   │   Permissions│     └──────────────────────┘    │
│  └──────────────┘     └──────────────┘                               │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Environment Variables Complete Reference

### Critical Security Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `WEBUI_SECRET_KEY` | ✅ | Session signing (≥32 bytes) |
| `WEBUI_AUTH` | ✅ | Enable authentication |
| `JWT_EXPIRES_IN` | - | Token lifetime (default: `4w` = 4 weeks). Security best practice: Use explicit expiration, not `-1` |
| `ENABLE_SIGNUP` | - | Allow new registrations |

### Core Application

| Variable | Default | Description |
|----------|---------|-------------|
| `PORT` | 8080 | Internal container port |
| `WEBUI_URL` | - | External URL for links |
| `ENV` | production | dev/prod mode |
| `CUSTOM_NAME` | Open WebUI | Instance name |

### Model Configuration

| Variable | Description |
|----------|-------------|
| `OLLAMA_BASE_URL` | Ollama server URL |
| `OPENAI_API_BASE_URL` | OpenAI-compatible endpoint |
| `OPENAI_API_KEY` | API key for OpenAI providers |
| `ENABLE_OLLAMA_API` | Enable Ollama integration |
| `ENABLE_OPENAI_API` | Enable OpenAI API compatibility |

### RAG / Vector Database

| Variable | Default | Options |
|----------|---------|---------|
| `VECTOR_DB` | chroma | chroma, pgvector, faiss, milvus, qdrant, weaviate, opensearch, elasticsearch, pinecone, oracle23ai |
| `RAG_EMBEDDING_MODEL` | sentence-transformers/all-MiniLM-L6-v2 | Any sentence-transformers model |
| `RAG_TOP_K` | 5 | Number of chunks to retrieve |
| `CHUNK_SIZE` | 1500 | Characters per chunk |
| `CHUNK_OVERLAP` | 150 | Overlap between chunks |
| `RAG_SYSTEM_CONTEXT` | true | Inject RAG context into system prompt for KV cache optimization |

### Web Search (20+ Engines)

| Variable | Description |
|----------|-------------|
| `ENABLE_WEB_SEARCH` | Enable web search feature |
| `WEB_SEARCH_ENGINE` | searxng, google, brave, bing, duckduckgo, etc. |
| `SEARXNG_QUERY_URL` | Custom SearXNG endpoint |
| `GOOGLE_PSE_API_KEY` | Google Programmable Search |
| `BRAVE_SEARCH_API_KEY` | Brave Search API |

### Code Execution

| Variable | Default | Description |
|----------|---------|-------------|
| `ENABLE_CODE_EXECUTION` | false | Enable code interpreter |
| `CODE_EXECUTION_ENGINE` | pyodide | pyodide or jupyter |
| `CODE_EXECUTION_JUPYTER_URL` | - | Jupyter kernel URL |
| `CODE_EXECUTION_JUPYTER_AUTH` | token | Authentication method |

### Image Generation

| Variable | Description |
|----------|-------------|
| `ENABLE_IMAGE_GENERATION` | Enable image generation |
| `IMAGE_GENERATION_ENGINE` | openai, comfyui, automatic1111 |
| `OPENAI_API_BASE_URL` | For DALL-E integration |

### Audio / Speech

| Variable | Description |
|----------|-------------|
| `WHISPER_MODEL` | STT model (base, small, medium, large) |
| `AUDIO_TTS_ENGINE` | openai, azure, elevenlabs |
| `AUDIO_TTS_VOICE` | Voice selection |

### Database

| Variable | Default | Description |
|----------|---------|-------------|
| `DATABASE_URL` | - | PostgreSQL/MySQL connection |
| `DATABASE_POOL_SIZE` | 10 | Connection pool size |
| `DATABASE_POOL_MAX_OVERFLOW` | 20 | Overflow connections |

### Redis (Multi-Instance)

| Variable | Description |
|----------|-------------|
| `REDIS_URL` | Single instance URL |
| `REDIS_SENTINEL_HOSTS` | Sentinel host list |
| `REDIS_SENTINEL_MASTER_NAME` | Master name |

### Performance Tuning

| Variable | Default | Description |
|----------|---------|-------------|
| `AIOHTTP_CLIENT_TIMEOUT` | 300 | HTTP client timeout |
| `AIOHTTP_CLIENT_TIMEOUT_MODEL_LIST` | 10 | Model list fetch timeout |
| `MAX_BODY_SIZE` | - | Max request body |
| `MAX_UPLOAD_SIZE` | - | Max file upload |

### User Permissions (50+ Granular Controls)

| Variable Pattern | Description |
|------------------|-------------|
| `USER_PERMISSIONS_CHAT_{FILE/ARCHIVE/DELETE}` | Chat permissions |
| `USER_PERMISSIONS_FEATURES_{WEB_SEARCH/CODE/TOOLS}` | Feature access |
| `USER_PERMISSIONS_WORKSPACE_{MODELS/PROMPTS/KNOWLEDGE}` | Workspace access |

---

## 4. Tools & Functions Development

### Tools vs Functions Comparison

| Aspect | Tools | Functions |
|--------|-------|-----------|
| **Purpose** | LLM function calling | Lifecycle hooks |
| **Triggered by** | LLM decision | Request/response flow |
| **Use cases** | External APIs, search | Modify requests, custom models |
| **Types** | Single type | Pipe, Filter, Action |

### Tools Development

**Basic Structure:**
```python
from pydantic import BaseModel, Field

class Tools:
    class Valves(BaseModel):
        api_key: str = Field(default="", description="API key")
        endpoint: str = Field(default="https://api.example.com")

    def __init__(self):
        self.valves = self.Valves()

    async def search_web(self, query: str) -> str:
        """Search the web for information"""
        # Implementation
        return results
```

### Functions Development

**Pipe Function (Custom Model):**
```python
class Pipe:
    class Valves(BaseModel):
        model: str = Field(default="gpt-4o")
        temperature: float = Field(default=0.7)

    def pipe(self, body: dict) -> Union[str, Generator, Iterator]:
        # Process the request
        messages = body.get("messages", [])
        # Call your LLM
        return response
```

**Filter Function (Request/Response Interception):**
```python
class Filter:
    def inlet(self, body: dict) -> dict:
        # Modify request before sending to LLM
        return body

    def outlet(self, body: dict) -> dict:
        # Modify response before returning to user
        return body
```

**Action Function (UI Button):**
```python
class Action:
    def action(self, body: dict) -> None:
        # Triggered by UI button click
        pass
```

### Injected Parameters

| Parameter | Description |
|-----------|-------------|
| `__user__` | User info (id, name, email, role) |
| `__event_emitter__` | Status updates and citations |
| `__metadata__` | Chat metadata |
| `__files__` | Uploaded files |
| `__tools__` | Available tools (in Pipes) |
| `__model__` | Current model info |

### Event Emitter Usage

```python
async def my_tool(self, query: str, __event_emitter__: callable) -> str:
    # Send status update
    await __event_emitter__({
        "type": "status",
        "data": {"description": "Searching...", "done": False}
    })

    # Do work
    result = await search(query)

    # Send citation
    await __event_emitter__({
        "type": "citation",
        "data": {"source": {"name": "Web Search", "url": result.url}}
    })

    # Mark complete
    await __event_emitter__({
        "type": "status",
        "data": {"description": "Done", "done": True}
    })

    return result.content
```

---

## 5. RAG System Deep Dive

### RAG Pipeline Architecture

```
Document Upload → Content Extraction → Chunking → Embedding → Vector Storage
                                                        ↓
User Query → Query Embedding → Vector Search → Context Assembly → LLM
                        ↓
            (Optional: BM25 + Reranking = Hybrid Search)
```

### Vector Database Options

| Backend | Best For | Notes |
|---------|----------|-------|
| **Chroma** | Simple setups, default | Embedded, no external deps |
| **pgvector** | Production, scale | PostgreSQL extension |
| **FAISS** | Performance, local | Facebook AI library |
| **Milvus** | Enterprise scale | Cloud-native |
| **Qdrant** | Rust-based, fast | Good for high throughput |
| **Weaviate** | GraphQL interface | Module ecosystem |

### Embedding Models

| Model | Dimensions | Speed | Quality |
|-------|------------|-------|---------|
| all-MiniLM-L6-v2 | 384 | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| all-mpnet-base-v2 | 768 | ⭐⭐⭐ | ⭐⭐⭐⭐ |
| BAAI/bge-small-en-v1.5 | 384 | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| BAAI/bge-base-en-v1.5 | 768 | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |

### Chunking Strategies

| Strategy | Use Case |
|----------|----------|
| **Token Splitter** | Code, structured text |
| **Character Splitter** | General documents |
| **Markdown Splitter** | Markdown files (preserves headers) |

### Hybrid Search Configuration

```python
# RAG System Context enables KV cache prefix caching
RAG_SYSTEM_CONTEXT=true

# Hybrid search combines:
# 1. Dense (vector similarity)
# 2. Sparse (BM25 keyword)
# 3. CrossEncoder reranking
```

### Knowledge Base Permissions

| Access Level | Description |
|--------------|-------------|
| **Private** | Owner only |
| **Public** | All users |
| **Group** | Specific user groups |

### Web Search + RAG

```python
# Combine web search results with knowledge base
# Web pages are processed same as uploaded documents
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=searxng
SEARXNG_QUERY_URL=http://searxng:8080/search?q=<query>&format=json
```

---

## 6. API Endpoints Reference

### Base URLs

| Prefix | Purpose |
|--------|---------|
| `/api` | Primary management APIs |
| `/api/v1` | File & RAG management |
| `/v1` | OpenAI-compatible |
| `/ollama` | Ollama proxy |

### Authentication

**Bearer Token:**
```bash
Authorization: Bearer <token>
```

**Get JWT Token:**
```bash
POST /api/v1/auths/signin
{"email": "user@example.com", "password": "password"}
```

### OpenAI-Compatible Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/v1/models` | GET | List available models |
| `/v1/chat/completions` | POST | Chat completion |
| `/v1/embeddings` | POST | Generate embeddings |

### Chat Completions Example

```bash
curl -X POST http://localhost:3000/v1/chat/completions \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.1",
    "messages": [{"role": "user", "content": "Hello!"}],
    "stream": false,
    "temperature": 0.7
  }'
```

### File Management

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/files/` | POST | Upload file |
| `/api/v1/files/` | GET | List files |
| `/api/v1/files/{id}` | DELETE | Delete file |
| `/api/v1/files/{id}/process/status` | GET | Check processing status |

### Knowledge Base

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/knowledge/create` | POST | Create KB |
| `/api/v1/knowledge` | GET | List KBs |
| `/api/v1/knowledge/{id}` | GET | Get KB details |
| `/api/v1/knowledge/{id}/file/add` | POST | Add file to KB |
| `/api/v1/knowledge/{id}/delete` | DELETE | Delete KB |

### RAG Query

```bash
POST /api/chat/completions
{
  "model": "llama3.1",
  "messages": [{"role": "user", "content": "Query"}],
  "files": [{
    "id": "knowledge-id",
    "type": "collection",
    "status": "processed"
  }]
}
```

---

## 7. Pipelines & Extensibility

### When to Use What

| Use Case | Solution |
|----------|----------|
| Simple request/response modification | Functions (Filter) |
| LLM function calling | Functions (Tools) |
| Custom model provider | Functions (Pipe) |
| Heavy computation, external libraries | Pipelines |
| Multi-stage processing chains | Pipelines |

### Pipeline Architecture

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   OpenWebUI     │────▶│   Pipelines     │────▶│   LLM Providers │
│   (Port 3000)   │     │   (Port 9099)   │     │                 │
└─────────────────┘     └─────────────────┘     └─────────────────┘
```

### Pipeline Types

| Type | Purpose | Key Method |
|------|---------|------------|
| **Filter** | Modify requests/responses | `inlet()`, `outlet()` |
| **Pipe** | Custom model/provider | `pipe()` |
| **Manifold** | Multi-model provider | `pipe()` with model list |
| **Tools** | Function calling | Inherits `FunctionCallingBlueprint` |

### Pipeline Valves Configuration

```python
class Valves(BaseModel):
    pipelines: List[str] = ["*"]  # Apply to all pipelines
    priority: int = 0  # Execution order (0 = first)
    api_key: str = Field(default="", description="API Key")
```

### Pipeline Security Warning

⚠️ **NEVER install pipelines from untrusted sources**
- Pipelines execute arbitrary Python code
- Can access filesystem, network, environment
- Can exfiltrate data or mine cryptocurrency

### Pipeline Development Example

```python
from pydantic import BaseModel
from typing import Optional

class Filter:
    class Valves(BaseModel):
        target_user_roles: List[str] = ["user"]

    def __init__(self):
        self.valves = self.Valves()

    def inlet(self, body: dict) -> dict:
        # Log or modify request
        print(f"Request: {body}")
        return body

    def outlet(self, body: dict) -> dict:
        # Log or modify response
        print(f"Response: {body}")
        return body
```

---

## 8. Quick Start Guides

### 8.1 Deploy OpenWebUI (Repository-Specific)

This repository uses CLIProxyAPI for LLM provider integration with OAuth authentication.

```bash
# 1. Prepare environment
cp .env.example .env
# Edit .env to set required values:
# - WEBUI_SECRET_KEY (>=32 bytes)
# - OPENAI_API_KEY (for CLIProxyAPI)
# - JUPYTER_TOKEN
# - POSTGRES_PASSWORD

# 2. Deploy complete stack
./deploy.sh --no-logs

# 3. Verify services
./status.sh
./check-cliproxyapi.sh

# 4. Run baseline tests
./test-rag.sh --baseline
./test-api.sh --baseline
```

**Note:** For OAuth setup, see `./configure-cliproxyapi-oauth.sh`

### 8.2 Enable RAG with Custom Settings

```bash
# docker-compose.yml environment
environment:
  - RAG_EMBEDDING_MODEL=sentence-transformers/all-MiniLM-L6-v2
  - RAG_TOP_K=5
  - CHUNK_SIZE=1500
  - CHUNK_OVERLAP=150
  - VECTOR_DB=chroma
  - RAG_SYSTEM_CONTEXT=true
```

### 8.3 Create a Simple Tool

```python
# tools/calculator.py
import ast
import operator

class Tools:
    def __init__(self):
        # Supported operators for safe evaluation
        self.operators = {
            ast.Add: operator.add,
            ast.Sub: operator.sub,
            ast.Mult: operator.mul,
            ast.Div: operator.truediv,
            ast.Pow: operator.pow,
            ast.USub: operator.neg,
        }

    def calculate(self, expression: str) -> str:
        """Safely calculate mathematical expressions."""
        try:
            tree = ast.parse(expression.strip(), mode='eval')
            result = self._safe_eval(tree.body)
            return f"Result: {result}"
        except Exception as e:
            return f"Error: {str(e)}"

    def _safe_eval(self, node):
        if isinstance(node, ast.Constant):  # Python 3.8+
            if isinstance(node.value, (int, float)):
                return node.value
            raise ValueError("Only numbers allowed")
        elif hasattr(ast, 'Num') and isinstance(node, ast.Num):  # Python < 3.8
            return node.n
        elif isinstance(node, ast.BinOp):
            op_type = type(node.op)
            if op_type not in self.operators:
                raise ValueError(f"Operator not allowed: {op_type.__name__}")
            left = self._safe_eval(node.left)
            right = self._safe_eval(node.right)
            return self.operators[op_type](left, right)
        elif isinstance(node, ast.UnaryOp):
            if type(node.op) not in self.operators:
                raise ValueError("Unary operator not allowed")
            operand = self._safe_eval(node.operand)
            return self.operators[type(node.op)](operand)
        else:
            raise ValueError(f"Node type not allowed: {type(node).__name__}")
```

> ⚠️ **Security Note**: This example uses AST-based safe evaluation instead of
> `eval()` to prevent code injection attacks. Never use `eval()` on user input.

### 8.4 RAG Workflow via API

```bash
# 1. Upload file
FILE_ID=$(curl -s -X POST \
  -H "Authorization: Bearer $TOKEN" \
  -F "file=@document.pdf" \
  http://localhost:3000/api/v1/files/ | jq -r '.id')

# 2. Wait for processing
while [ "$(curl -s -H "Authorization: Bearer $TOKEN" \
  http://localhost:3000/api/v1/files/$FILE_ID/process/status | jq -r '.status')" != "processed" ]; do
  sleep 1
done

# 3. Create KB and add file
KB_ID=$(curl -s -X POST \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name": "My KB"}' \
  http://localhost:3000/api/v1/knowledge/create | jq -r '.id')

curl -s -X POST \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{\"file_id\": \"$FILE_ID\"}" \
  http://localhost:3000/api/v1/knowledge/$KB_ID/file/add

# 4. Query with RAG
curl -X POST \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{
    \"model\": \"llama3.1\",
    \"messages\": [{\"role\": \"user\", \"content\": \"Summarize the document\"}],
    \"files\": [{\"id\": \"$KB_ID\", \"type\": \"collection\"}]
  }" \
  http://localhost:3000/api/chat/completions
```

---

---

**OpenWebUI Version:** v0.8.3
**Last Updated:** 2026-02-18
*Sources: Official OpenWebUI docs, GitHub discussions, community guides*
