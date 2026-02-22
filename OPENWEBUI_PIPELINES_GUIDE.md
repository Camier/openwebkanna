# Comprehensive Guide to OpenWebUI Pipelines and Extensibility

> **OpenWebUI Version:** v0.8.3 | **Last Updated:** 2026-02-18

**Related Documents:**
- [Master Reference](OPENWEBUI_MASTER_REFERENCE.md) - Overview and navigation
- [Environment Variables](openwebui_env_reference.md) - Complete configuration reference
- [Tools & Functions Guide](openwebui-tools-functions-guide.md) - Development patterns
- [RAG Technical Reference](openwebui_rag_technical_reference.md) - RAG system details
- [API Examples](API_EXAMPLES.md) - Repository-specific API usage

---

> **Research Summary**: This guide compiles findings from official documentation, community resources, and practical examples to provide a complete understanding of OpenWebUI's extensibility mechanisms.

---

## Table of Contents

1. [Introduction to OpenWebUI Extensibility](#introduction)
2. [Pipelines vs Tools vs Functions](#comparison)
3. [What Are Pipelines?](#what-are-pipelines)
4. [Pipeline Architecture & Components](#architecture)
5. [Pipeline Types Explained](#pipeline-types)
6. [Creating Custom Pipelines](#creating-pipelines)
7. [Valves: Configuration System](#valves)
8. [Integration with External Services](#integration)
9. [Deployment & Configuration](#deployment)
10. [Production Deployment](#production-deployment)
11. [Community Pipelines & Examples](#community-examples)
12. [Best Practices](#best-practices)
13. [Advanced Patterns](#advanced-patterns)
14. [Decision Matrix: When to Use What](#decision-matrix)

---

> **See Also:** [Tools & Functions Guide](openwebui-tools-functions-guide.md) for simpler extensibility options, [Master Reference](OPENWEBUI_MASTER_REFERENCE.md#7-pipelines--extensibility) for overview, and [RAG Technical Reference](openwebui_rag_technical_reference.md) for RAG-specific pipelines.

---

## 1. Introduction to OpenWebUI Extensibility {#introduction}

OpenWebUI is more than just a chat interface—it's a **full-fledged, extensible AI platform** that allows customization of virtually every component:

- **Models & API integration**
- **User interface & behavior**
- **Prompts & system messages**
- **RAG (Retrieval-Augmented Generation)**
- **Custom workflows & logic**

The platform provides three main extensibility mechanisms:
- **Functions** (built-in, easiest)
- **Tools** (function calling for LLMs)
- **Pipelines** (separate service, most powerful)

---

## 2. Pipelines vs Tools vs Functions {#comparison}

### Quick Comparison Table

| Feature | Functions | Tools | Pipelines |
|---------|-----------|-------|-----------|
| **Location** | Built into OpenWebUI | Built into OpenWebUI | Separate service/container |
| **Complexity** | Low | Medium | High |
| **Performance Impact** | On main instance | On main instance | Offloaded to separate service |
| **Use Case** | Simple filters/integrations | Function calling | Heavy computation, complex workflows |
| **Setup** | Copy-paste in UI | Copy-paste in UI | Requires separate container |
| **Python Libraries** | Limited | Limited | Any Python library |
| **Scalability** | Limited by main instance | Limited by main instance | Independently scalable |

### When to Use Each

**Use FUNCTIONS when:**
- Adding basic filters (e.g., profanity filters, simple transformations)
- Simple API integrations
- Modifying prompts before they reach the LLM
- Quick prototypes and simple customizations

**Use TOOLS when:**
- You want the LLM to call external functions
- Implementing calculators, search, APIs that LLM decides to use
- Function calling scenarios
- The LLM needs to perform actions based on user queries

**Use PIPELINES when:**
- Running computationally heavy tasks (large models, complex logic)
- Need to offload processing from main OpenWebUI instance
- Implementing custom RAG systems
- Building multi-step workflows
- Need any Python library
- Message monitoring (Langfuse, Opik)
- Rate limiting across multiple requests
- Real-time translation
- Custom LLM provider integration

---

## 3. What Are Pipelines? {#what-are-pipelines}

**Pipelines** is a **UI-agnostic OpenAI API Plugin Framework** that brings modular, customizable workflows to any client supporting OpenAI API specs.

### Key Characteristics

- **Standalone Service**: Runs in its own Docker container
- **OpenAI API Compatible**: Works with any UI/client supporting OpenAI API specs
- **Arbitrary Code Execution**: Full Python environment with any library
- **Modular**: Add/remove pipelines without restarting OpenWebUI

### ⚠️ Important Warning from Official Docs

> **"DO NOT USE PIPELINES IF!"**
>
> If your goal is simply to add support for additional providers like Anthropic or basic filters, you likely **don't need Pipelines**. For those cases, Open WebUI **Functions** are a better fit—it's built-in, much more convenient, and easier to configure.
>
> Pipelines comes into play when you're dealing with **computationally heavy tasks** (e.g., running large models or complex logic) that you want to **offload** from your main Open WebUI instance for better performance and scalability.

### What You Can Build with Pipelines

- **Function Calling Pipeline**: Handle function calls with custom logic
- **Custom RAG Pipeline**: Sophisticated retrieval-augmented generation
- **Message Monitoring**: Langfuse, Opik integration for analytics
- **Rate Limit Filter**: Control request flow
- **Real-Time Translation**: LibreTranslate integration
- **Toxic Message Filter**: Content moderation
- **Custom LLM Providers**: Integrate new model providers
- **Multi-Model Routing**: Load balancing between models

---

## 4. Pipeline Architecture & Components {#architecture}

### Architecture Overview

```
┌─────────────────┐         ┌──────────────────┐         ┌─────────────────┐
│   OpenWebUI     │◄───────►│   Pipelines      │◄───────►│  External APIs  │
│   (Container)   │  HTTP   │   (Container)    │  HTTP   │  (LLMs, etc.)   │
│                 │         │   Port: 9099     │         │                 │
└─────────────────┘         └──────────────────┘         └─────────────────┘
                                    │
                                    ▼
                           ┌──────────────────┐
                           │  /pipelines dir  │
                           │  (Python files)  │
                           └──────────────────┘
```

### Core Components

1. **Pipeline Server**: FastAPI-based service running on port 9099
2. **Pipeline Files**: Python modules in `/pipelines` directory
3. **Valves**: Configuration system for user-customizable settings
4. **Hooks**: Lifecycle methods (`on_startup`, `on_shutdown`, `inlet`, `outlet`, `pipe`)

### Request Flow

```
User Request → OpenWebUI → Pipelines (inlet) → External LLM → Pipelines (outlet) → User
                  │              │                                  │
                  │              ▼                                  ▼
                  │         Filter/Transform                  Filter/Transform
                  │                                              │
                  └──────────────────────────────────────────────┘
```

---

## 5. Pipeline Types Explained {#pipeline-types}

### 5.1 Filter Pipeline

**Purpose**: Intercept and modify requests before/after LLM processing

**When to use**:
- Pre-process user messages (add context, modify prompts)
- Post-process LLM responses (format, filter, translate)
- RAG context injection
- Safety/content moderation
- Prompt injection detection

**Key Methods**:
- `inlet(body, user)`: Called BEFORE request goes to LLM
- `outlet(body, user)`: Called AFTER response comes from LLM

**Example Structure**:
```python
class Pipeline:
    class Valves(BaseModel):
        pipelines: List[str] = []  # Target pipeline IDs
        priority: int = 0          # Execution order (lower = earlier)

    def __init__(self):
        self.type = "filter"
        self.name = "My Filter"
        self.valves = self.Valves(**{"pipelines": ["*"]})

    async def inlet(self, body: dict, user: Optional[dict] = None) -> dict:
        # Modify request before LLM
        return body

    async def outlet(self, body: dict, user: Optional[dict] = None) -> dict:
        # Modify response after LLM
        return body
```

### 5.2 Pipe Pipeline

**Purpose**: Take complete control over the chat flow

**When to use**:
- Custom LLM provider integration
- Complex multi-step workflows
- Custom RAG implementation
- When you want to completely handle the request

**Key Methods**:
- `pipe(user_message, model_id, messages, body)`: Main processing function

**Return Types**:
- `str`: Complete response as string
- `Generator`: Streaming response
- `Iterator`: Streaming chunks

**Example Structure**:
```python
class Pipeline:
    def __init__(self):
        self.valves = self.Valves()

    def pipe(
        self,
        user_message: str,
        model_id: str,
        messages: List[dict],
        body: dict
    ) -> Union[str, Generator, Iterator]:
        # Implement your logic
        return "Response"
```

### 5.3 Manifold Pipeline (Special Pipe)

**Purpose**: Integrate new LLM providers with model selection

**When to use**:
- Adding new LLM provider (Anthropic, Groq, etc.)
- Exposing multiple models from a single provider
- Dynamic model discovery

**Key Methods**:
- `pipelines()`: Return list of available models
- `pipe()`: Handle the actual inference

**Example Structure**:
```python
class Pipeline:
    def __init__(self):
        self.type = "manifold"

    def pipelines(self) -> List[dict]:
        return ["model-1", "model-2", "model-3"]

    def pipe(
        self,
        user_message: str,
        model_id: str,
        messages: List[dict],
        body: dict
    ):
        # Use model_id to route to specific model
        pass
```

### 5.4 Tools (Special Filter)

**Purpose**: Enable LLM function calling with automatic tool selection

**When to use**:
- Calculator tools
- Search tools
- API integrations where LLM decides when to call

**Implementation**: Inherit from `FunctionCallingBlueprint`

**Example Structure**:
```python
from blueprints.function_calling_blueprint import Pipeline as FunctionCallingBlueprint
import ast
import operator

class Pipeline(FunctionCallingBlueprint):
    class Tools:
        def __init__(self, pipeline):
            self.pipeline = pipeline
            # Supported operators for safe evaluation
            self.operators = {
                ast.Add: operator.add,
                ast.Sub: operator.sub,
                ast.Mult: operator.mul,
                ast.Div: operator.truediv,
                ast.Pow: operator.pow,
                ast.USub: operator.neg,
            }

        def calculator(self, equation: str) -> str:
            """Safely calculate mathematical expressions."""
            try:
                tree = ast.parse(equation.strip(), mode='eval')
                result = self._safe_eval(tree.body)
                return f"Result: {result}"
            except Exception as e:
                return f"Error: Invalid expression - {str(e)}"

        def _safe_eval(self, node):
            if isinstance(node, ast.Constant):  # Python 3.8+ (check first)
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
                    raise ValueError(f"Unary operator not allowed")
                operand = self._safe_eval(node.operand)
                return self.operators[type(node.op)](operand)
            else:
                raise ValueError(f"Node type not allowed: {type(node).__name__}")

    def __init__(self):
        super().__init__()
        self.tools = self.Tools(self)
```

> ⚠️ **Security Warning**: Never use `eval()` on user input. Use safe evaluation
> methods like `ast.literal_eval()` or the pattern shown above.

### 5.5 Action (UI Feature)

**Purpose**: Add interactive buttons to chat interface

**When to use**:
- Suggested follow-up questions
- Quick action buttons
- UI shortcuts

**Note**: Actions are typically implemented as Functions within OpenWebUI rather than in Pipelines.

---

## 6. Creating Custom Pipelines {#creating-pipelines}

### Step-by-Step Development Guide

#### Step 1: Set Up Development Environment

**Option A: Docker (Recommended for Production)**
```bash
# Clone repository
git clone https://github.com/open-webui/pipelines.git
cd pipelines

# Run with Docker
docker run -d -p 9099:9099 \
  --add-host=host.docker.internal:host-gateway \
  -v pipelines:/app/pipelines \
  --name pipelines \
  --restart always \
  ghcr.io/open-webui/pipelines:main
```

**Option B: Local Python**
```bash
# Requirements
# - Python 3.11 (only officially supported version)

git clone https://github.com/open-webui/pipelines.git
cd pipelines
pip install -r requirements.txt
sh ./start.sh
```

#### Step 2: Create Your Pipeline File

Create a Python file in the `/pipelines` directory:

```python
"""
title: My Custom Pipeline
author: Your Name
date: 2024-01-01
version: 1.0
license: MIT
description: A description of what this pipeline does
requirements: requests, pydantic
"""

from typing import List, Optional
from pydantic import BaseModel

class Pipeline:
    # Configuration valves
    class Valves(BaseModel):
        api_key: str = ""
        model: str = "gpt-4"
        temperature: float = 0.7

    def __init__(self):
        self.type = "filter"  # or "pipe", "manifold"
        self.name = "My Pipeline"
        self.valves = self.Valves()

    async def on_startup(self):
        """Called when server starts"""
        print(f"on_startup:{__name__}")

    async def on_shutdown(self):
        """Called when server stops"""
        print(f"on_shutdown:{__name__}")

    # For Filter pipelines
    async def inlet(self, body: dict, user: Optional[dict] = None) -> dict:
        """Process before LLM call"""
        return body

    async def outlet(self, body: dict, user: Optional[dict] = None) -> dict:
        """Process after LLM call"""
        return body

    # For Pipe pipelines
    def pipe(self, user_message: str, model_id: str,
             messages: List[dict], body: dict):
        """Handle the complete request"""
        return "Response"
```

#### Step 3: Add Dependencies (Front-matter)

Specify requirements in the file header:

```python
"""
requirements: requests, beautifulsoup4, numpy
"""
```

Dependencies are automatically installed when the pipeline loads.

#### Step 4: Test Your Pipeline

1. Save file to `/pipelines` directory
2. Restart pipelines container (or it may auto-reload)
3. Check logs: `docker logs -f pipelines`
4. Connect OpenWebUI to your pipeline URL
5. Test through OpenWebUI interface

---

## 7. Valves: Configuration System {#valves}

### What are Valves?

**Valves** provide a mechanism for configuring pipelines through the OpenWebUI interface. They allow administrators to customize pipeline behavior without modifying code.

### Valve Types

1. **Pipeline Valves**: Configure which pipelines a filter applies to
   - `pipelines`: List of target pipeline IDs (use `["*"]` for all)
   - `priority`: Execution order (lower number = higher priority)

2. **Custom Valves**: Your own configuration parameters
   - API keys
   - Model names
   - Threshold values
   - Feature flags

### Valve Definition Example

```python
from pydantic import BaseModel, Field

class Valves(BaseModel):
    # Pipeline targeting
    pipelines: List[str] = Field(
        default=["*"],
        description="Target pipelines for this filter"
    )
    priority: int = Field(
        default=0,
        description="Execution priority (lower = earlier)"
    )

    # Custom configuration
    api_key: str = Field(
        default="",
        description="API key for external service"
    )
    temperature: float = Field(
        default=0.7,
        description="Sampling temperature"
    )
    enabled: bool = Field(
        default=True,
        description="Enable this pipeline"
    )
```

### Accessing Valves in Your Pipeline

```python
# In your pipeline methods
api_key = self.valves.api_key
temperature = self.valves.temperature
```

### Configuring Valves in OpenWebUI

1. Go to **Admin Panel > Settings > Pipelines**
2. Select your pipeline
3. Modify valve values in the UI
4. Changes apply immediately (no restart needed)

---

## 8. Integration with External Services {#integration}

### Pattern for External API Integration

```python
import requests
from typing import List, Optional
from pydantic import BaseModel

class Pipeline:
    class Valves(BaseModel):
        api_key: str = ""
        base_url: str = "https://api.example.com"

    def __init__(self):
        self.valves = self.Valves()

    def call_external_api(self, message: str):
        headers = {
            "Authorization": f"Bearer {self.valves.api_key}",
            "Content-Type": "application/json"
        }
        response = requests.post(
            f"{self.valves.base_url}/endpoint",
            json={"message": message},
            headers=headers
        )
        return response.json()
```

### Common Integration Patterns

#### RAG with External Vector DB
```python
# Inject context before LLM call
async def inlet(self, body: dict, user: Optional[dict] = None) -> dict:
    messages = body.get("messages", [])
    user_message = self.get_last_user_message(messages)

    # Retrieve relevant documents
    context = self.vector_search(user_message)

    # Prepend context to system message
    messages.insert(0, {
        "role": "system",
        "content": f"Context: {context}"
    })

    body["messages"] = messages
    return body
```

#### Multi-Model Routing
```python
def pipe(self, user_message: str, model_id: str,
         messages: List[dict], body: dict):
    # Route based on message complexity
    if len(user_message) > 1000:
        return self.call_expensive_model(body)
    else:
        return self.call_cheap_model(body)
```

#### Message Monitoring
```python
async def outlet(self, body: dict, user: Optional[dict] = None) -> dict:
    # Send to monitoring service
    self.log_to_langfuse(body)
    return body
```

---

## 9. Deployment & Configuration {#deployment}

### Docker Deployment

#### Basic Deployment
```bash
docker run -d -p 9099:9099 \
  --add-host=host.docker.internal:host-gateway \
  -v pipelines:/app/pipelines \
  --name pipelines \
  --restart always \
  ghcr.io/open-webui/pipelines:main
```

#### With Custom Pipelines
```bash
docker run -d -p 9099:9099 \
  --add-host=host.docker.internal:host-gateway \
  -e PIPELINES_URLS="https://github.com/user/repo/blob/main/pipeline.py" \
  -v pipelines:/app/pipelines \
  --name pipelines \
  --restart always \
  ghcr.io/open-webui/pipelines:main
```

#### Docker Compose
```yaml
services:
  openwebui:
    image: ghcr.io/open-webui/open-webui:main
    ports:
      - "3000:8080"
    volumes:
      - open-webui:/app/backend/data
    environment:
      - OPENAI_API_BASE_URL=http://pipelines:9099
      - OPENAI_API_KEY=0p3n-w3bu!

  pipelines:
    image: ghcr.io/open-webui/pipelines:main
    volumes:
      - pipelines:/app/pipelines
    restart: always
    environment:
      - PIPELINES_API_KEY=0p3n-w3bu!

volumes:
  open-webui: {}
  pipelines: {}
```

### Connecting OpenWebUI to Pipelines

1. Go to **Admin Panel > Settings > Connections**
2. Click **+** to add a connection
3. Set:
   - **API URL**: `http://localhost:9099` (or `http://pipelines:9099` if both in Docker)
   - **API Key**: `0p3n-w3bu!`
4. Verify connection (pipelines icon appears)

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `PIPELINES_DIR` | Directory for pipeline files | `/app/pipelines` |
| `PIPELINES_URLS` | Semicolon-separated URLs to install pipelines | "" |
| `PIPELINES_API_KEY` | API key for authentication | `0p3n-w3bu!` |

---

## 10. Production Deployment {#production-deployment}

### Docker Compose Configuration

```yaml
version: '3.8'
services:
  pipelines:
    image: ghcr.io/open-webui/pipelines:main
    container_name: openwebui_pipelines
    restart: unless-stopped
    ports:
      - "127.0.0.1:9099:9099"
    environment:
      - PIPELINES_API_KEY=${PIPELINES_API_KEY:-change-me-in-production}
      - PIPELINES_DIR=/app/pipelines
    volumes:
      - pipelines_data:/app/pipelines
    # Resource limits
    cpus: "${PIPELINES_CPUS:-2.0}"
    mem_limit: "${PIPELINES_MEM_LIMIT:-4g}"
    pids_limit: 512
    # Health check
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:9099/health"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 40s
    # Security
    security_opt:
      - no-new-privileges:true
    read_only: true
    tmpfs:
      - /tmp:noexec,nosuid,size=100m
    networks:
      - openwebui_network

volumes:
  pipelines_data:
    driver: local

networks:
  openwebui_network:
    external: true  # Use existing OpenWebUI network
```

### Security Hardening Checklist

- [ ] Change default API key (`0p3n-w3bu!`)
- [ ] Enable TLS/HTTPS (use reverse proxy like Nginx or Traefik)
- [ ] Set resource limits (CPU, memory, PIDs)
- [ ] Run with read-only filesystem
- [ ] Use non-root user (if supported by image)
- [ ] Network isolation (dedicated Docker network)
- [ ] Regular security updates for base image
- [ ] Review all pipeline code before deployment
- [ ] Enable container health checks
- [ ] Set up log aggregation and monitoring

### Network Configuration

Pipelines must be on the same Docker network as OpenWebUI:

```bash
# Create network (if not exists)
docker network create openwebui_network

# Connect OpenWebUI (if not already connected)
docker network connect openwebui_network openwebui

# Connect Pipelines
docker network connect openwebui_network openwebui_pipelines
```

### Health Checks

Pipelines provides a health endpoint:

```bash
curl http://localhost:9099/health
# Returns: {"status": "healthy"}
```

### Monitoring

**Prometheus Metrics (if enabled):**
```yaml
# Add to docker-compose
environment:
  - ENABLE_PROMETHEUS=true
```

**Structured Logging:**
```yaml
environment:
  - LOG_LEVEL=info
  - LOG_FORMAT=json
```

### Backup and Restore

**Backup pipeline configurations:**
```bash
docker exec openwebui_pipelines tar czf - /app/pipelines > pipelines_backup_$(date +%Y%m%d).tar.gz
```

**Restore:**
```bash
docker exec -i openwebui_pipelines tar xzf - -C / < pipelines_backup_YYYYMMDD.tar.gz
```

### High Availability (Advanced)

For production deployments requiring high availability:

```yaml
services:
  pipelines-lb:
    image: nginx:alpine
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf:ro
    depends_on:
      - pipelines-1
      - pipelines-2

  pipelines-1:
    image: ghcr.io/open-webui/pipelines:main
    # ... config

  pipelines-2:
    image: ghcr.io/open-webui/pipelines:main
    # ... config
```

### Troubleshooting

**Container won't start:**
```bash
docker logs openwebui_pipelines
docker inspect openwebui_pipelines
```

**OpenWebUI can't connect to Pipelines:**
```bash
# Test from OpenWebUI container
docker exec openwebui curl http://pipelines:9099/health

# Check network
docker network inspect openwebui_network
```

**Pipeline not loading:**
```bash
# Check pipeline files exist
docker exec openwebui_pipelines ls -la /app/pipelines/

# Check pipeline logs
docker logs openwebui_pipelines | grep -i error
```

---

## 11. Community Pipelines & Examples {#community-examples}

### Official Examples Repository

The official pipelines repo contains examples at:
https://github.com/open-webui/pipelines/tree/main/examples

### Common Example Categories

#### Filters
- **Detoxify Filter**: Toxic content detection
- **LibreTranslate Filter**: Real-time translation
- **Rate Limit Filter**: Request throttling
- **Prompt Injection Filter**: Security filtering

#### Pipes
- **Haystack RAG**: Custom RAG implementation
- **Function Calling**: Tool use examples
- **Langfuse Monitoring**: Message tracking
- **Opik Monitoring**: Alternative tracking

#### Providers (Manifolds)
- **Anthropic**: Claude integration
- **Groq**: Groq API integration
- **Custom**: Template for new providers

### Example: LibreTranslate Filter

```python
"""
title: LibreTranslate Filter
author: open-webui
date: 2024-05-30
version: 1.0
license: MIT
description: Real-time translation using LibreTranslate
requirements: requests
"""

import requests
from typing import List, Optional
from pydantic import BaseModel

class Pipeline:
    class Valves(BaseModel):
        pipelines: List[str] = ["*"]
        priority: int = 0
        libretranslate_url: str = "http://localhost:5000"
        target_language: str = "en"

    def __init__(self):
        self.type = "filter"
        self.name = "LibreTranslate Filter"
        self.valves = self.Valves()

    async def inlet(self, body: dict, user: Optional[dict] = None) -> dict:
        # Translate user message to target language
        messages = body.get("messages", [])
        # ... translation logic ...
        return body

    async def outlet(self, body: dict, user: Optional[dict] = None) -> dict:
        # Translate response back to user's language
        return body
```

---

## 11. Best Practices {#best-practices}

### Security

⚠️ **CRITICAL WARNING**:
> Pipelines execute arbitrary Python code. **Never install pipelines from untrusted sources.** A malicious pipeline could:
> - Access your file system
> - Exfiltrate data
> - Mine cryptocurrency
> - Compromise your system

**Best Practices:**
- Always review pipeline source code before installing
- Use official examples as templates
- Run pipelines in isolated containers
- Limit network access where possible
- Use API keys stored in valves, not hardcoded

### Performance

- Offload heavy computation to pipelines (that's their purpose)
- Use streaming responses (`Generator`, `Iterator`) for long outputs
- Set appropriate timeouts for external API calls
- Consider caching frequently accessed data

### Development

- Use front-matter for dependencies
- Implement proper error handling
- Add logging for debugging
- Test incrementally
- Check Docker logs when debugging: `docker logs -f pipelines`

### Configuration

- Use valves for all configurable parameters
- Provide sensible defaults
- Add descriptions to valve fields
- Use proper types (List, Optional, etc.)

### Error Handling

```python
async def inlet(self, body: dict, user: Optional[dict] = None) -> dict:
    try:
        # Your logic
        return body
    except Exception as e:
        print(f"Error in inlet: {e}")
        # Return unchanged body on error
        return body
```

---

## 12. Advanced Patterns {#advanced-patterns}

### Multi-Model Routing

Route requests to different models based on content:

```python
def pipe(self, user_message: str, model_id: str,
         messages: List[dict], body: dict):

    # Analyze complexity
    complexity = self.assess_complexity(user_message)

    if complexity == "high":
        # Route to powerful model
        return self.call_model("gpt-4", body)
    elif complexity == "medium":
        # Route to balanced model
        return self.call_model("gpt-3.5-turbo", body)
    else:
        # Route to fast/cheap model
        return self.call_model("local-model", body)
```

### Request Modification

Add system prompts dynamically:

```python
async def inlet(self, body: dict, user: Optional[dict] = None) -> dict:
    messages = body.get("messages", [])

    # Add system context based on user
    if user and user.get("role") == "admin":
        context = "You are speaking to an administrator."
    else:
        context = "You are speaking to a regular user."

    messages.insert(0, {"role": "system", "content": context})
    body["messages"] = messages

    return body
```

### Chaining Multiple Pipelines

Use priority to control execution order:

```python
# Pipeline 1: Translation (priority: 0)
class TranslationFilter:
    class Valves(BaseModel):
        priority: int = 0

# Pipeline 2: Content Filter (priority: 1)
class ContentFilter:
    class Valves(BaseModel):
        priority: int = 1

# Pipeline 3: Logging (priority: 2)
class LoggingFilter:
    class Valves(BaseModel):
        priority: int = 2
```

### Custom RAG Pipeline

```python
def pipe(self, user_message: str, model_id: str,
         messages: List[dict], body: dict):

    # 1. Retrieve relevant documents
    docs = self.vector_store.search(user_message, k=5)

    # 2. Build context
    context = "\n".join([doc.content for doc in docs])

    # 3. Augment prompt
    augmented_prompt = f"""Context:
{context}

User Question: {user_message}
"""

    # 4. Call LLM with augmented prompt
    messages[-1]["content"] = augmented_prompt
    return self.call_llm(messages)
```

---

## 13. Decision Matrix: When to Use What {#decision-matrix}

| Use Case | Solution | Complexity | Example |
|----------|----------|------------|---------|
| Simple prompt modification | Function | Low | Add timestamp to prompts |
| Profanity filter | Function | Low | Basic word filtering |
| LLM function calling | Tool | Medium | Calculator, search |
| External API integration (simple) | Function | Medium | Weather lookup |
| External API integration (complex) | Pipeline | High | Multi-step workflows |
| Custom LLM provider | Pipeline | High | New AI service |
| Heavy ML models | Pipeline | High | Custom embeddings |
| Message monitoring | Pipeline | Medium | Langfuse integration |
| Multi-model routing | Pipeline | High | Load balancing |
| Complex RAG | Pipeline | High | Custom retrieval |
| Rate limiting | Pipeline | Medium | Request throttling |
| Real-time translation | Pipeline | Medium | LibreTranslate |
| Safety filtering | Pipeline | Medium | LlamaGuard |

---

## Summary

OpenWebUI provides a powerful, layered extensibility system:

1. **Functions** → Built-in, simple, quick modifications
2. **Tools** → Function calling for LLM-driven actions
3. **Pipelines** → Separate service for heavy computation and complex workflows

**The key insight**: Start with Functions, move to Pipelines only when you need to offload heavy processing or require capabilities beyond what the built-in system can provide.

**Remember**: Pipelines = Power + Flexibility - Convenience. Use them wisely!

---

## References

1. [Official Pipelines Documentation](https://docs.openwebui.com/features/pipelines/)
2. [Pipelines GitHub Repository](https://github.com/open-webui/pipelines)
3. [OpenWebUI Main Repository](https://github.com/open-webui/open-webui)
4. [Zohaib's Pipeline Guide](https://zohaib.me/extending-openwebui-using-pipelines/)
5. [OpenWebUI Community](https://openwebui.com/)

---

**OpenWebUI Version:** v0.8.3
**Last Updated:** 2026-02-18

*Generated: February 2026*
*Research compiled from official docs, GitHub repositories, and community resources*
