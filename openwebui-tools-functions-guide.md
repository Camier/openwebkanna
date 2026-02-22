# OpenWebUI Tools and Functions Development Guide

> **OpenWebUI Version:** v0.8.3 | **Last Updated:** 2026-02-18

**Related Documents:**
- [Master Reference](OPENWEBUI_MASTER_REFERENCE.md) - Overview and navigation
- [Environment Variables](openwebui_env_reference.md) - Complete configuration reference
- [RAG Technical Reference](openwebui_rag_technical_reference.md) - RAG system details
- [Pipelines Guide](OPENWEBUI_PIPELINES_GUIDE.md) - Pipeline deployment
- [API Examples](API_EXAMPLES.md) - Repository-specific API usage

---

A comprehensive technical guide for developing OpenWebUI plugins - covering Tools, Functions (Pipe, Filter, Action), Valves system, and best practices.

---

## Table of Contents

1. [Overview: Tools vs Functions](#1-overview-tools-vs-functions)
2. [The Valves System](#2-the-valves-system)
3. [Tools Development](#3-tools-development)
4. [Functions Development](#4-functions-development)
5. [Common Patterns and Examples](#5-common-patterns-and-examples)
6. [Testing and Debugging](#6-testing-and-debugging)
7. [Best Practices](#7-best-practices)
8. [Community Resources](#8-community-resources)

---

> **See Also:** [Pipelines Guide](OPENWEBUI_PIPELINES_GUIDE.md) for advanced extensibility patterns, [Master Reference](OPENWEBUI_MASTER_REFERENCE.md#4-tools--functions-development) for overview, and [RAG Technical Reference](openwebui_rag_technical_reference.md) for RAG integration patterns.

---

## 1. Overview: Tools vs Functions

### Key Differences

| Aspect | Tools | Functions |
|--------|-------|-----------|
| **Purpose** | Extend LLM capabilities with external actions | Modify/extend OpenWebUI behavior |
| **Execution** | LLM decides when to call (function calling) | Hook into request/response lifecycle |
| **Types** | Single type (callable tools) | Pipe, Filter, Action |
| **UI** | Appear as available capabilities | Can add buttons, filters, custom models |
| **Scope** | Per-request tool invocation | Global or per-model behavior |

### When to Use What

**Use Tools when:**
- You need the LLM to decide when to perform an action
- The action requires external data/API calls (weather, search, calculations)
- You want function-calling behavior
- Examples: Web search, calculator, database lookup, API integrations

**Use Functions when:**
- You need to intercept/modify requests or responses
- You want to create custom "models" (Pipe)
- You need to add UI buttons (Action)
- You want to filter/transform all messages (Filter)
- Examples: Custom agent workflows, message preprocessing, UI enhancements

### Architecture Overview

```
User Request
     |
     v
[Filter Function] ----> Modifies request before LLM
     |
     v
[LLM Processing] <----> [Tools called as needed]
     |
     v
[Filter Function] ----> Modifies response after LLM
     |
     v
  Response

[Pipe Function] -------> Replaces entire LLM call
[Action Function] -----> Adds UI buttons
```

---

## 2. The Valves System

Valves provide a declarative configuration mechanism for all OpenWebUI plugins. They automatically generate UI controls based on Python type hints and Pydantic field definitions.

### Types of Valves

| Valve Type | Access Level | Set In | Use Case |
|------------|--------------|--------|----------|
| `Valves` | Admin | Tools/Functions settings | API keys, endpoints, global settings |
| `UserValves` | User | Chat interface | Personal preferences, per-session config |

### Basic Valves Structure

```python
from pydantic import BaseModel, Field
from typing import Optional, Literal

class Tools:
    """Tools class with Valves configuration."""

    class Valves(BaseModel):
        """Admin-configurable settings."""
        api_key: str = Field(
            default="",
            description="API key for the service"
        )
        endpoint: str = Field(
            default="https://api.example.com",
            description="Service endpoint URL"
        )
        timeout: int = Field(
            default=30,
            description="Request timeout in seconds"
        )
        debug_mode: bool = Field(
            default=False,
            description="Enable debug logging"
        )
        pass  # Always end with pass for parsing reliability

    class UserValves(BaseModel):
        """User-configurable settings."""
        result_count: int = Field(
            default=5,
            description="Number of results to return"
        )
        language: Literal["en", "es", "fr"] = Field(
            default="en",
            description="Response language"
        )
        pass

    def __init__(self):
        self.valves = self.Valves()
        # Access with: self.valves.api_key
```

### UI Control Mapping by Type

| Python Type | UI Control | Example |
|-------------|------------|---------|
| `str` | Text input | `name: str = Field(default="")` |
| `int` | Number input | `count: int = Field(default=5)` |
| `float` | Decimal input | `temperature: float = Field(default=0.7)` |
| `bool` | Toggle switch | `enabled: bool = Field(default=True)` |
| `Literal[...]` | Dropdown | `mode: Literal["fast", "slow"]` |
| `Optional[T]` | Optional field | `api_key: Optional[str] = None` |

### Runtime Access Patterns

```python
class Tools:
    class Valves(BaseModel):
        api_key: str = Field(default="")
        pass

    class UserValves(BaseModel):
        max_results: int = Field(default=5)
        pass

    def __init__(self):
        self.valves = self.Valves()  # Admin valves accessible here

    async def search(
        self,
        query: str,
        __user__: dict = None  # User dict contains UserValves
    ) -> str:
        # Access admin valves
        api_key = self.valves.api_key

        # Access user valves (if UserValves defined)
        if __user__ and "valves" in __user__:
            max_results = __user__["valves"]["max_results"]

        # Your logic here
        return result
```

### Special Valves for Filter Functions

```python
class Filter:
    class Valves(BaseModel):
        # Control which models this filter applies to
        pipelines: List[str] = Field(
            default=["*"],
            description="Target pipeline IDs (* = all)"
        )
        # Execution priority (lower = earlier)
        priority: int = Field(
            default=0,
            description="Filter priority (lower executes first)"
        )
        pass
```

### Environment Variable Integration

```python
import os

class Tools:
    class Valves(BaseModel):
        api_key: str = Field(
            default=os.getenv("MY_API_KEY", ""),
            description="API key (can set via env var)"
        )
        pass
```

---

## 3. Tools Development

### Complete Tools Class Structure

```python
"""
title: Weather Search Tool
description: Get current weather information for any location
author: Your Name
version: 1.0.0
"""

import requests
import json
from pydantic import BaseModel, Field
from typing import Optional, Callable, Any

class Tools:
    """
    Tools class containing all callable tools.
    Each public method (not starting with _) becomes a tool.
    """

    # =========================================================================
    # VALVES CONFIGURATION
    # =========================================================================
    class Valves(BaseModel):
        """Admin-configurable settings."""
        api_key: str = Field(
            default="",
            description="OpenWeatherMap API key"
        )
        units: str = Field(
            default="metric",
            description="Temperature units (metric/imperial)"
        )
        pass

    class UserValves(BaseModel):
        """User-configurable settings."""
        include_forecast: bool = Field(
            default=False,
            description="Include 5-day forecast"
        )
        pass

    # =========================================================================
    # INITIALIZATION
    # =========================================================================
    def __init__(self):
        """Initialize the Tool instance."""
        self.valves = self.Valves()
        # Disable automatic citations if doing custom citation handling
        # self.citation = False

    # =========================================================================
    # TOOL METHODS
    # =========================================================================
    async def get_weather(
        self,
        location: str,
        __event_emitter__: Callable[[dict], Any] = None,
        __user__: dict = None
    ) -> str:
        """
        Get current weather for a location.

        :param location: City name or coordinates (e.g., "London, UK")
        :param __event_emitter__: Injected by Open WebUI for status updates
        :param __user__: User information including UserValves
        :return: Weather information as formatted string
        """
        # Emit status update
        if __event_emitter__:
            await __event_emitter__({
                "type": "status",
                "data": {"description": "Fetching weather data...", "done": False}
            })

        try:
            # Access valves
            api_key = self.valves.api_key
            if not api_key:
                return "Error: API key not configured"

            # Access user valves
            user_valves = __user__.get("valves", {}) if __user__ else {}
            include_forecast = user_valves.get("include_forecast", False)

            # Make API call
            url = f"https://api.openweathermap.org/data/2.5/weather"
            params = {
                "q": location,
                "appid": api_key,
                "units": self.valves.units
            }

            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            # Format result
            weather = data["weather"][0]["description"]
            temp = data["main"]["temp"]
            humidity = data["main"]["humidity"]

            result = f"Weather in {location}: {weather}, {temp}°C, {humidity}% humidity"

            # Emit citation
            if __event_emitter__:
                await __event_emitter__({
                    "type": "citation",
                    "data": {
                        "source": {"name": "OpenWeatherMap", "url": url},
                        "document": [result],
                        "metadata": {"location": location}
                    }
                })

            # Mark status as done
            if __event_emitter__:
                await __event_emitter__({
                    "type": "status",
                    "data": {"description": "Weather data retrieved", "done": True}
                })

            return result

        except requests.RequestException as e:
            error_msg = f"Error fetching weather: {str(e)}"
            if __event_emitter__:
                await __event_emitter__({
                    "type": "status",
                    "data": {"description": error_msg, "done": True}
                })
            return error_msg

    def _helper_method(self, data: dict) -> dict:
        """
        Private helper method (not exposed as tool).
        Prefix with _ to hide from LLM function calling.
        """
        return data
```

### Available Injected Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `__user__` | `dict` | User info: `id`, `name`, `email`, `role`, `valves` |
| `__event_emitter__` | `Callable` | Emit status updates and citations |
| `__event_call__` | `Callable` | Same as event emitter (legacy) |
| `__metadata__` | `dict` | Chat metadata |
| `__messages__` | `list` | Current conversation messages |
| `__files__` | `list` | Uploaded files in current request |
| `__oauth_token__` | `str` | OAuth token if using OAuth |

### Event Emitter Types

```python
# Status update
await __event_emitter__({
    "type": "status",
    "data": {
        "description": "Processing...",
        "done": False  # Set True when complete
    }
})

# Citation
await __event_emitter__({
    "type": "citation",
    "data": {
        "source": {
            "name": "Source Name",
            "url": "https://example.com"
        },
        "document": ["Relevant text snippet"],
        "metadata": {"key": "value"}
    }
})

# Message (stream to chat)
await __event_emitter__({
    "type": "message",
    "data": {"content": "Partial response..."}
})
```

### Native Function Calling Mode

```python
class Tools:
    """Native function calling with proper response format."""

    def __init__(self):
        self.valves = self.Valves()
        # Disable automatic citations for manual control
        self.citation = False

    async def search_web(
        self,
        query: str,
        __event_emitter__: Callable = None
    ) -> str:
        """Search the web and return results."""
        # Perform search
        results = await self._do_search(query)

        # Emit citations for each result
        for idx, result in enumerate(results, 1):
            await __event_emitter__({
                "type": "citation",
                "data": {
                    "source": {
                        "name": result["title"],
                        "url": result["url"]
                    },
                    "document": [result["snippet"]],
                    "metadata": {"index": idx}
                }
            })

        return json.dumps(results)
```

---

## 4. Functions Development

Functions come in three types: **Pipe**, **Filter**, and **Action**.

### 4.1 Pipe Functions

Pipe Functions create custom "models/agents" that intercept the entire request/response cycle.

```python
"""
title: Custom LLM Pipe
author: Your Name
version: 1.0.0
"""

from pydantic import BaseModel, Field
from typing import Optional, Callable, Any, Union, Generator, Iterator
import requests

class Pipe:
    """
    Pipe function - creates a custom model in OpenWebUI.
    The 'pipe' method handles all requests to this model.
    """

    class Valves(BaseModel):
        """Configuration options."""
        target_model: str = Field(
            default="llama3.2:latest",
            description="Model to route to"
        )
        api_endpoint: str = Field(
            default="http://localhost:11434",
            description="Ollama API endpoint"
        )
        temperature: float = Field(
            default=0.7,
            description="Generation temperature"
        )
        pass

    def __init__(self):
        self.type = "pipe"  # Required
        self.valves = self.Valves()

    async def pipe(
        self,
        body: dict,
        __user__: Optional[dict] = None,
        __metadata__: Optional[dict] = None,
        __event_emitter__: Optional[Callable] = None,
        __tools__: Optional[dict] = None,
        __files__: Optional[list] = None
    ) -> Union[str, Generator, Iterator]:
        """
        Process the request and return a response.

        :param body: Full request body with messages, model, etc.
        :param __user__: User information
        :param __metadata__: Chat metadata
        :param __event_emitter__: For status updates
        :param __tools__: Available tools for function calling
        :param __files__: Uploaded files
        :return: Response string or generator for streaming
        """
        # Access user info
        user_id = __user__.get("id") if __user__ else None
        user_valves = __user__.get("valves", {}) if __user__ else {}

        # Modify body as needed
        body["model"] = self.valves.target_model

        # Add custom parameters
        if "options" not in body:
            body["options"] = {}
        body["options"]["temperature"] = self.valves.temperature

        # Call your custom logic or API
        response = await self._call_llm(body)

        return response

    async def _call_llm(self, body: dict) -> str:
        """Internal method to call LLM."""
        url = f"{self.valves.api_endpoint}/api/chat"
        response = requests.post(url, json=body)
        return response.json()["message"]["content"]
```

### Pipe Method Body Structure

```python
body = {
    "model": "pipe-model-name",
    "messages": [
        {"role": "system", "content": "..."},
        {"role": "user", "content": "..."}
    ],
    "stream": True,
    "options": {
        "temperature": 0.7,
        # ... other options
    }
}
```

### 4.2 Filter Functions

Filter Functions intercept and modify requests (inlet) and responses (outlet).

```python
"""
title: Message Filter
author: Your Name
version: 1.0.0
"""

from pydantic import BaseModel, Field
from typing import Optional, List

class Filter:
    """
    Filter function - intercepts requests and responses.
    inlet() modifies before LLM, outlet() modifies after LLM.
    """

    class Valves(BaseModel):
        """Configuration."""
        pipelines: List[str] = Field(
            default=["*"],
            description="Apply to all pipelines"
        )
        priority: int = Field(
            default=0,
            description="Execution priority (lower = earlier)"
        )
        add_system_prompt: bool = Field(
            default=True,
            description="Add system prompt to all requests"
        )
        system_prompt: str = Field(
            default="You are a helpful assistant.",
            description="System prompt to add"
        )
        pass

    def __init__(self):
        self.type = "filter"  # Required
        self.name = "MessageFilter"
        self.valves = self.Valves()

    async def inlet(
        self,
        body: dict,
        __user__: Optional[dict] = None,
        __metadata__: Optional[dict] = None
    ) -> dict:
        """
        Process request BEFORE it goes to the LLM.

        :param body: Request body with messages
        :param __user__: User information
        :param __metadata__: Chat metadata
        :return: Modified body
        """
        messages = body.get("messages", [])

        # Add system prompt if enabled
        if self.valves.add_system_prompt:
            # Check if system message exists
            has_system = any(m.get("role") == "system" for m in messages)
            if not has_system:
                messages.insert(0, {
                    "role": "system",
                    "content": self.valves.system_prompt
                })

        # Modify user message
        for message in reversed(messages):
            if message.get("role") == "user":
                # Add prefix to user message
                content = message.get("content", "")
                message["content"] = f"[Filtered] {content}"
                break

        body["messages"] = messages
        return body

    async def outlet(
        self,
        body: dict,
        __user__: Optional[dict] = None,
        __metadata__: Optional[dict] = None
    ) -> dict:
        """
        Process response AFTER it comes from the LLM.

        :param body: Response body with messages
        :param __user__: User information
        :param __metadata__: Chat metadata
        :return: Modified body
        """
        messages = body.get("messages", [])

        # Modify assistant's response
        for message in reversed(messages):
            if message.get("role") == "assistant":
                content = message.get("content", "")
                # Add suffix to response
                message["content"] = f"{content}\n\n[Response filtered]"
                break

        body["messages"] = messages
        return body
```

### Filter Priority System

When multiple filters are installed, they execute in priority order:

```python
# Priority 0 runs first, then priority 1, etc.
class HighPriorityFilter:
    class Valves(BaseModel):
        priority: int = 0  # Runs first

class LowPriorityFilter:
    class Valves(BaseModel):
        priority: int = 10  # Runs later
```

### 4.3 Action Functions

Action Functions add custom buttons to the chat interface.

```python
"""
title: Explain Like I'm Five
author: Your Name
version: 1.0.0
"""

from pydantic import BaseModel, Field
from typing import Optional, Callable

class Action:
    """
    Action function - adds a button to the message toolbar.
    Clicking the button triggers the action() method.
    """

    class Valves(BaseModel):
        """Configuration."""
        target_model: str = Field(
            default="",
            description="Model to use for the action"
        )
        pass

    def __init__(self):
        self.valves = self.Valves()

    async def action(
        self,
        body: dict,
        __user__: Optional[dict] = None,
        __metadata__: Optional[dict] = None,
        __event_emitter__: Optional[Callable] = None
    ) -> str:
        """
        Execute when user clicks the action button.

        :param body: Contains the selected message
        :param __user__: User information
        :param __metadata__: Chat metadata
        :param __event_emitter__: For UI updates
        :return: Result to display
        """
        # Get the message that was clicked
        messages = body.get("messages", [])
        if not messages:
            return "No message selected"

        # Get the last user or assistant message
        target_message = None
        for msg in reversed(messages):
            if msg.get("role") in ["user", "assistant"]:
                target_message = msg.get("content", "")
                break

        if not target_message:
            return "No content to process"

        # Process the message
        explanation = await self._explain_like_im_five(target_message)

        return explanation

    async def _explain_like_im_five(self, text: str) -> str:
        """Generate simple explanation."""
        # Your implementation here
        return f"Simple explanation of: {text[:100]}..."
```

### Action Button Configuration

The button appears in the message toolbar with:
- Title from the docstring title
- Description from the docstring description

---

## 5. Common Patterns and Examples

### Pattern 1: Web Search Tool

```python
"""
title: Web Search Tool
description: Search the web using DuckDuckGo
"""

from duckduckgo_search import DDGS
from pydantic import BaseModel, Field
from typing import Optional, Callable, Any

class Tools:
    class Valves(BaseModel):
        max_results: int = Field(default=5, description="Max search results")
        region: str = Field(default="wt-wt", description="Search region")
        safesearch: str = Field(default="moderate", description="Safe search level")
        pass

    def __init__(self):
        self.valves = self.Valves()

    async def web_search(
        self,
        query: str,
        __event_emitter__: Callable[[dict], Any] = None
    ) -> str:
        """Search the web for information."""
        if __event_emitter__:
            await __event_emitter__({
                "type": "status",
                "data": {"description": "Searching the web...", "done": False}
            })

        try:
            with DDGS() as ddgs:
                results = list(ddgs.text(
                    query,
                    max_results=self.valves.max_results,
                    region=self.valves.region
                ))

            # Format results
            formatted = []
            for idx, r in enumerate(results, 1):
                formatted.append(f"{idx}. {r['title']}\n{r['href']}\n{r['body']}")

                # Emit citation
                if __event_emitter__:
                    await __event_emitter__({
                        "type": "citation",
                        "data": {
                            "source": {"name": r['title'], "url": r['href']},
                            "document": [r['body']]
                        }
                    })

            if __event_emitter__:
                await __event_emitter__({
                    "type": "status",
                    "data": {"description": f"Found {len(results)} results", "done": True}
                })

            return "\n\n".join(formatted)

        except Exception as e:
            return f"Search error: {str(e)}"
```

### Pattern 2: Calculator Tool

```python
"""
title: Calculator
description: Perform mathematical calculations safely
"""

import ast
import operator
from pydantic import BaseModel, Field

class Tools:
    class Valves(BaseModel):
        max_precision: int = Field(default=10, description="Decimal precision")
        pass

    def __init__(self):
        self.valves = self.Valves()

    async def calculate(
        self,
        expression: str
    ) -> str:
        """
        Safely evaluate a mathematical expression.

        :param expression: Math expression (e.g., "2 + 2 * 5")
        :return: Calculation result
        """
        try:
            result = self._safe_eval(expression)
            return f"{expression} = {result:.{self.valves.max_precision}f}"
        except Exception as e:
            return f"Error: {str(e)}"

    def _safe_eval(self, expression: str) -> float:
        """Safely evaluate math expression using AST."""
        allowed_ops = {
            ast.Add: operator.add,
            ast.Sub: operator.sub,
            ast.Mult: operator.mul,
            ast.Div: operator.truediv,
            ast.Pow: operator.pow,
            ast.USub: operator.neg,
        }

        def eval_node(node):
            if isinstance(node, ast.Constant):
                if isinstance(node.value, (int, float)):
                    return node.value
                raise ValueError("Unsupported constant type")
            # For Python 3.7 compatibility, also handle ast.Num:
            elif hasattr(ast, 'Num') and isinstance(node, ast.Num):
                return node.n
            elif isinstance(node, ast.BinOp):
                return allowed_ops[type(node.op)](eval_node(node.left), eval_node(node.right))
            elif isinstance(node, ast.UnaryOp):
                return allowed_ops[type(node.op)](eval_node(node.operand))
            else:
                raise ValueError(f"Unsupported operation: {type(node)}")

        tree = ast.parse(expression.strip(), mode='eval')
        return eval_node(tree.body)
```

### Pattern 3: Custom Agent Pipe

```python
"""
title: Research Agent
author: Your Name
description: Multi-step research agent with planning
"""

from pydantic import BaseModel, Field
from typing import Optional, Callable, Any
import json

class Pipe:
    class Valves(BaseModel):
        model: str = Field(default="gpt-4", description="Model for agent")
        max_iterations: int = Field(default=5, description="Max research steps")
        pass

    def __init__(self):
        self.type = "pipe"
        self.valves = self.Valves()

    async def pipe(
        self,
        body: dict,
        __user__: Optional[dict] = None,
        __event_emitter__: Optional[Callable] = None
    ) -> str:
        """Execute multi-step research."""
        messages = body.get("messages", [])
        query = self._get_last_user_message(messages)

        if __event_emitter__:
            await __event_emitter__({
                "type": "status",
                "data": {"description": "Planning research steps...", "done": False}
            })

        # Step 1: Create plan
        plan = await self._create_plan(query)

        # Step 2: Execute each step
        results = []
        for idx, step in enumerate(plan, 1):
            if __event_emitter__:
                await __event_emitter__({
                    "type": "status",
                    "data": {"description": f"Step {idx}/{len(plan)}: {step}", "done": False}
                })

            result = await self._execute_step(step)
            results.append(result)

        # Step 3: Synthesize
        if __event_emitter__:
            await __event_emitter__({
                "type": "status",
                "data": {"description": "Synthesizing results...", "done": False}
            })

        final = await self._synthesize(query, results)

        if __event_emitter__:
            await __event_emitter__({
                "type": "status",
                "data": {"description": "Research complete", "done": True}
            })

        return final

    def _get_last_user_message(self, messages: list) -> str:
        """Extract last user message."""
        for msg in reversed(messages):
            if msg.get("role") == "user":
                return msg.get("content", "")
        return ""

    async def _create_plan(self, query: str) -> list:
        """Create research plan."""
        # Implementation
        return ["Search for background", "Find recent developments", "Analyze trends"]

    async def _execute_step(self, step: str) -> str:
        """Execute single research step."""
        # Implementation
        return f"Result for: {step}"

    async def _synthesize(self, query: str, results: list) -> str:
        """Synthesize final answer."""
        # Implementation
        return f"Synthesized answer for: {query}\n" + "\n".join(results)
```

### Pattern 4: Rate Limiting Filter

```python
"""
title: Rate Limit Filter
description: Limit requests per user
"""

from pydantic import BaseModel, Field
from typing import Optional
import time

class Filter:
    class Valves(BaseModel):
        pipelines: list = Field(default=["*"])
        priority: int = Field(default=0)
        max_requests: int = Field(default=10, description="Requests per window")
        window_seconds: int = Field(default=60, description="Time window")
        pass

    def __init__(self):
        self.type = "filter"
        self.name = "RateLimitFilter"
        self.valves = self.Valves()
        self.request_log = {}  # user_id -> [timestamps]

    async def inlet(
        self,
        body: dict,
        __user__: Optional[dict] = None
    ) -> dict:
        """Check rate limit before processing."""
        user_id = __user__.get("id") if __user__ else "anonymous"
        now = time.time()

        # Get user's request history
        if user_id not in self.request_log:
            self.request_log[user_id] = []

        # Clean old entries
        window = self.valves.window_seconds
        self.request_log[user_id] = [
            ts for ts in self.request_log[user_id]
            if now - ts < window
        ]

        # Check limit
        if len(self.request_log[user_id]) >= self.valves.max_requests:
            raise Exception(f"Rate limit exceeded. Max {self.valves.max_requests} requests per {window}s")

        # Log this request
        self.request_log[user_id].append(now)

        return body

    async def outlet(self, body: dict, __user__: Optional[dict] = None) -> dict:
        """Pass through unchanged."""
        return body
```

---

## 6. Testing and Debugging

### Debugging Techniques

#### 1. Use Print Statements

```python
class Tools:
    def __init__(self):
        self.valves = self.Valves()
        print(f"[DEBUG] Tools initialized with valves: {self.valves}")

    async def my_tool(self, query: str, __user__: dict = None):
        print(f"[DEBUG] Tool called with query: {query}")
        print(f"[DEBUG] User: {__user__}")
        # ... logic ...
        print(f"[DEBUG] Returning result: {result}")
        return result
```

Print output appears in the OpenWebUI container logs:
```bash
docker logs open-webui 2>&1 | grep "\[DEBUG\]"
```

#### 2. Event Emitter Status Updates

```python
async def my_tool(self, query: str, __event_emitter__: Callable = None):
    if __event_emitter__:
        await __event_emitter__({
            "type": "status",
            "data": {"description": f"Processing: {query}", "done": False}
        })

    # Step 1
    if __event_emitter__:
        await __event_emitter__({
            "type": "status",
            "data": {"description": "Step 1 complete", "done": False}
        })

    # ... more steps ...
```

#### 3. Return Debug Info

```python
async def my_tool(self, query: str, __user__: dict = None) -> str:
    try:
        result = await self._process(query)
        return result
    except Exception as e:
        # Return full error info for debugging
        import traceback
        return f"Error: {str(e)}\n\n{traceback.format_exc()}"
```

### Accessing Logs

```bash
# View OpenWebUI container logs
docker logs -f open-webui

# Filter for specific output
docker logs open-webui 2>&1 | grep "your-tool-name"

# Save logs to file
docker logs open-webui > openwebui.log 2>&1
```

### Testing Tools

1. **Workspace Editor Testing:**
   - Go to Workspace > Tools
   - Create/edit your tool
   - Use the built-in test interface

2. **Direct Chat Testing:**
   - Enable the tool in a conversation
   - Ask questions that should trigger the tool
   - Check the tool call indicator

3. **API Testing:**
   ```bash
   curl -X POST http://localhost:3000/api/chat/completions \
     -H "Content-Type: application/json" \
     -H "Authorization: Bearer YOUR_API_KEY" \
     -d '{
       "model": "your-model",
       "messages": [{"role": "user", "content": "test query"}],
       "tools": ["your-tool-name"]
     }'
   ```

### Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| Tool not appearing | Check docstring has `title` and `description` |
| Valves not saving | Ensure `pass` at end of Valves class |
| Import errors | Use only standard library + pre-installed packages |
| Async errors | Mark tool methods with `async def` |
| No citations shown | Set `self.citation = False` for manual control |
| Filter not applying | Check `pipelines` valve includes target model |

---

## 7. Best Practices

### Code Organization

1. **Document Metadata:**
```python
"""
title: Clear Tool Name
author: Your Name
version: 1.0.0
description: Clear description of what this does
license: MIT
"""
```

2. **Use Type Hints:**
```python
async def search(
    self,
    query: str,
    max_results: int = 5,
    __user__: dict = None
) -> str:
```

3. **Validate Inputs:**
```python
async def search(self, query: str) -> str:
    if not query or not query.strip():
        return "Error: Empty query"
    if len(query) > 1000:
        return "Error: Query too long (max 1000 chars)"
    # ... continue
```

4. **Handle Errors Gracefully:**
```python
async def api_call(self, url: str) -> str:
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.Timeout:
        return "Error: Request timed out"
    except requests.RequestException as e:
        return f"Error: {str(e)}"
```

### Security Considerations

```python
class Tools:
    class Valves(BaseModel):
        api_key: str = Field(
            default="",
            description="API key (stored securely)"
        )
        pass

    def _sanitize_input(self, text: str) -> str:
        """Remove potentially dangerous characters."""
        import re
        # Remove script tags, etc.
        return re.sub(r'<script.*?>.*?</script>', '', text, flags=re.DOTALL)

    def _validate_url(self, url: str) -> bool:
        """Check URL is in allowed domains."""
        allowed = ["api.example.com", "api2.example.com"]
        from urllib.parse import urlparse
        parsed = urlparse(url)
        return parsed.netloc in allowed
```

### Performance Tips

1. **Use Async I/O:**
```python
import aiohttp

async def fetch(self, url: str) -> str:
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            return await response.text()
```

2. **Cache Results:**
```python
from functools import lru_cache

class Tools:
    def __init__(self):
        self._cache = {}

    async def expensive_operation(self, key: str) -> str:
        if key in self._cache:
            return self._cache[key]

        result = await self._compute(key)
        self._cache[key] = result
        return result
```

3. **Set Timeouts:**
```python
response = requests.get(url, timeout=(5, 30))  # (connect, read)
```

---

## 8. Community Resources

### Official Resources

- **Documentation:** https://docs.openwebui.com/
- **GitHub:** https://github.com/open-webui/open-webui
- **Community Hub:** https://openwebui.com/

### Community Tool Collections

- **Haervwe's Tools:** https://github.com/Haervwe/open-webui-tools
  - arXiv Search, Perplexica Search, Image Generation
  - Planner Agent, Multi-Model Conversations
  - ComfyUI Integration

- **Mplogas Collection:** https://github.com/mplogas/open-webui
  - Various tools, functions, and prompts

### Useful Packages

Tools can use these pre-installed packages:
- `requests` - HTTP requests
- `aiohttp` - Async HTTP
- `pydantic` - Data validation
- `duckduckgo-search` - Web search
- `beautifulsoup4` - HTML parsing

### Version Compatibility

| OpenWebUI | Features Available |
|-----------|-------------------|
| 0.4.x | Basic Tools, Valves |
| 0.5.x | + __request__ parameter, improved Valves |
| 0.6.0+ | Full Functions (Pipe/Filter/Action), UserValves |

### Migration Notes (0.4 to 0.5+)

```python
# OLD (0.4.x)
def pipe(self, user_message: str, model_id: str, messages: list, body: dict):
    pass

# NEW (0.5.x+)
async def pipe(self, body: dict, __user__: dict = None):
    messages = body.get("messages", [])
    # ...
```

---

## Quick Reference

### Tools Template
```python
"""title: My Tool"""
from pydantic import BaseModel, Field
from typing import Callable, Any

class Tools:
    class Valves(BaseModel):
        setting: str = Field(default="")
        pass

    def __init__(self):
        self.valves = self.Valves()

    async def tool_name(
        self,
        param: str,
        __event_emitter__: Callable[[dict], Any] = None
    ) -> str:
        """Tool description."""
        return "result"
```

### Filter Template
```python
"""title: My Filter"""
from pydantic import BaseModel, Field

class Filter:
    class Valves(BaseModel):
        pipelines: list = Field(default=["*"])
        priority: int = Field(default=0)
        pass

    def __init__(self):
        self.type = "filter"
        self.valves = self.Valves()

    async def inlet(self, body: dict, __user__: dict = None):
        return body

    async def outlet(self, body: dict, __user__: dict = None):
        return body
```

### Pipe Template
```python
"""title: My Pipe"""
from pydantic import BaseModel, Field

class Pipe:
    class Valves(BaseModel):
        model: str = Field(default="gpt-4")
        pass

    def __init__(self):
        self.type = "pipe"
        self.valves = self.Valves()

    async def pipe(self, body: dict, __user__: dict = None):
        return "response"
```

---

**OpenWebUI Version:** v0.8.3
**Last Updated:** 2026-02-18

*This guide was compiled from official OpenWebUI documentation, community discussions, and practical examples from the OpenWebUI ecosystem.*
