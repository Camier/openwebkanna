#!/bin/bash
#
# Patch OpenWebUI tools in-place via the admin API.
#
# Why:
# - OpenWebUI tool "specs" are generated from *all* callables on the Tools object.
#   Helper methods like `_foo()` accidentally become model-callable tool functions.
# - This script updates selected tools so only intended functions are exposed.
#
# Safety:
# - Never prints secrets (tokens, API keys).
# - Uses admin signin (or OPENWEBUI_API_KEY if provided).
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/init.sh"
load_env_defaults

require_cmds() {
    local missing=0
    for cmd in curl jq; do
        if ! command_exists "$cmd"; then
            print_error "Missing dependency: $cmd"
            missing=1
        fi
    done
    if [ "$missing" -ne 0 ]; then
        exit 1
    fi
}

uri_encode() {
    local value="$1"
    jq -nr --arg v "$value" '$v|@uri'
}

OPENWEBUI_URL="${OPENWEBUI_URL:-http://localhost:${WEBUI_PORT:-3000}}"
OPENWEBUI_SIGNIN_PATH="${OPENWEBUI_SIGNIN_PATH:-/api/v1/auths/signin}"
OPENWEBUI_SIGNIN_EMAIL="${OPENWEBUI_SIGNIN_EMAIL:-admin@localhost}"
OPENWEBUI_SIGNIN_PASSWORD="${OPENWEBUI_SIGNIN_PASSWORD:-admin}"
OPENWEBUI_API_KEY="${OPENWEBUI_API_KEY:-}"
OPENWEBUI_AUTO_AUTH="${OPENWEBUI_AUTO_AUTH:-true}"
CURL_TIMEOUT="${CURL_TIMEOUT:-45}"

API_TOKEN=""
TMP_DIR=""

cleanup() {
    if [ -n "$TMP_DIR" ] && [ -d "$TMP_DIR" ]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT

request_api() {
    local method="$1"
    local path="$2"
    local payload_file="${3:-}"
    local output_file="$4"
    local url="${OPENWEBUI_URL%/}${path}"

    local -a args=(
        -sS
        -m "$CURL_TIMEOUT"
        -X "$method"
        "$url"
        -H "Accept: application/json"
        -o "$output_file"
        -w "%{http_code}"
    )

    if [ -n "$API_TOKEN" ]; then
        args+=(-H "Authorization: Bearer $API_TOKEN")
    fi
    if [ -n "$payload_file" ]; then
        args+=(-H "Content-Type: application/json" -d @"$payload_file")
    fi

    curl "${args[@]}" 2>/dev/null || true
}

signin() {
    if [ -n "$OPENWEBUI_API_KEY" ]; then
        API_TOKEN="$OPENWEBUI_API_KEY"
        print_success "Using OPENWEBUI_API_KEY for auth"
        return 0
    fi

    if [ "${OPENWEBUI_AUTO_AUTH}" != "true" ]; then
        print_warning "OPENWEBUI_AUTO_AUTH=false and OPENWEBUI_API_KEY is empty; continuing unauthenticated"
        return 0
    fi

    print_step "Signing in to OpenWebUI to get bearer token"
    local payload_file="$TMP_DIR/signin.payload.json"
    local output_file="$TMP_DIR/signin.response.json"
    local code=""

    jq -n --arg email "$OPENWEBUI_SIGNIN_EMAIL" --arg password "$OPENWEBUI_SIGNIN_PASSWORD" \
        '{email:$email, password:$password}' >"$payload_file"

    code="$(request_api "POST" "$OPENWEBUI_SIGNIN_PATH" "$payload_file" "$output_file")"
    if [ "$code" != "200" ]; then
        local detail
        detail="$(jq -r '.detail // .error.message // .message // "signin failed"' "$output_file" 2>/dev/null || true)"
        print_warning "Signin failed (HTTP $code): $detail"
        return 0
    fi

    API_TOKEN="$(jq -r '.token // empty' "$output_file" 2>/dev/null || true)"
    if [ -n "$API_TOKEN" ]; then
        print_success "OpenWebUI bearer token acquired"
        return 0
    fi

    print_warning "Signin response did not include token; continuing unauthenticated"
    return 0
}

patch_tool_by_id() {
    local tool_id="$1"
    local content_file="$2"

    local tool_id_encoded
    tool_id_encoded="$(uri_encode "$tool_id")"

    print_step "Patching tool: $tool_id"

    local get_file="$TMP_DIR/tool_${tool_id}.get.json"
    local get_code
    get_code="$(request_api "GET" "/api/v1/tools/id/${tool_id_encoded}" "" "$get_file")"
    if [ "$get_code" != "200" ]; then
        local detail
        detail="$(jq -r '.detail // .error.message // .message // "failed to fetch tool"' "$get_file" 2>/dev/null || true)"
        print_error "Failed to fetch tool ${tool_id} (HTTP $get_code): $detail"
        return 1
    fi

    local name meta_json update_payload update_file update_code
    name="$(jq -r '.name // empty' "$get_file" 2>/dev/null || true)"
    meta_json="$(jq -c '.meta // {}' "$get_file" 2>/dev/null || echo '{}')"
    if [ -z "$name" ]; then
        print_error "Tool ${tool_id} missing name in response"
        return 1
    fi

    update_payload="$TMP_DIR/tool_${tool_id}.update.payload.json"
    update_file="$TMP_DIR/tool_${tool_id}.update.response.json"

    jq -n \
        --arg id "$tool_id" \
        --arg name "$name" \
        --rawfile content "$content_file" \
        --argjson meta "$meta_json" \
        '{id:$id, name:$name, content:$content, meta:$meta}' >"$update_payload"

    update_code="$(request_api "POST" "/api/v1/tools/id/${tool_id_encoded}/update" "$update_payload" "$update_file")"
    if [ "$update_code" != "200" ]; then
        local detail
        detail="$(jq -r '.detail // .error.message // .message // "update failed"' "$update_file" 2>/dev/null || true)"
        print_error "Tool ${tool_id} update failed (HTTP $update_code): $detail"
        return 1
    fi

    # Basic validation: no specs should start with "_" after patching.
    if jq -r '.specs[]?.name // empty' "$update_file" 2>/dev/null | grep -q '^_'; then
        print_error "Tool ${tool_id} still exposes private helper functions in specs"
        return 1
    fi

    print_success "Tool patched: $tool_id"
    return 0
}

write_tool_sources() {
    local llm_file="$TMP_DIR/llm_council.py"
    local scraper_file="$TMP_DIR/web_scraper.py"

    cat >"$llm_file" <<'PY_LLM_COUNCIL'
#!/usr/bin/env python3
import asyncio
import json
import os
import urllib.parse
from typing import Any, Dict, List, Optional, Tuple

import httpx
from pydantic import BaseModel


class Valves(BaseModel):
    # Comma-separated model IDs that must match the IDs exposed by /v1/models.
    MODELS: str = "gpt-5.3-codex,minimax/chat-elite,minimax/chat-thinking,minimax/chat-quality,minmax,glm-5,qwen3-coder-flash"

    # Semicolon-separated OpenAI-compatible base URLs (with or without /v1 suffix).
    # Default is the local CLIProxyAPI container.
    API_BASE_URLS: str = "http://cliproxyapi:8317/v1"

    # Semicolon-separated API keys matching API_BASE_URLS.
    # NOTE: For local CLIProxyAPI base URLs, you can leave entries empty and this tool
    # will fall back to OPENAI_API_KEY (or CLIPROXYAPI_API_KEY) from the OpenWebUI env.
    API_KEYS: str = ""

    # Optional: force routing per model to avoid ambiguous discovery when multiple
    # providers expose the same model id.
    #
    # JSON object: { "<model_id>": "<base_url>|<index>" }
    # Examples:
    #   {"glm-5": "http://cliproxyapi:8317/v1", "minimax/chat-elite": 0}
    MODEL_ROUTE_MAP: str = ""

    TIMEOUT: int = 120


def _split_csv(value: str, sep: str = ",") -> List[str]:
    return [item.strip() for item in (value or "").split(sep) if item and item.strip()]


def _default_local_api_key() -> str:
    # IMPORTANT: only use this fallback for known-local base URLs to avoid leaking a
    # local key to an external host.
    for env_name in ("CLIPROXYAPI_API_KEY", "OPENAI_API_KEY"):
        value = os.getenv(env_name, "").strip()
        if value:
            return value
    return ""


def _is_local_cliproxy_base_url(base_url: str) -> bool:
    try:
        parsed = urllib.parse.urlsplit(base_url)
    except Exception:
        return False

    host = (parsed.hostname or "").strip().lower()
    port = parsed.port
    if not host or port is None:
        return False

    return host in {"cliproxyapi", "localhost", "127.0.0.1", "host.docker.internal"} and port == 8317


def _provider_configs(api_base_urls: str, api_keys: str) -> List[Tuple[str, str]]:
    base_urls = _split_csv(api_base_urls, sep=";") or ["http://cliproxyapi:8317/v1"]
    keys = _split_csv(api_keys, sep=";")

    local_key = _default_local_api_key()

    providers: List[Tuple[str, str]] = []
    for idx, raw_base_url in enumerate(base_urls):
        base_url = raw_base_url.rstrip("/")
        key = keys[idx] if idx < len(keys) else ""
        if not key and local_key and _is_local_cliproxy_base_url(base_url):
            key = local_key
        providers.append((base_url, key))
    return providers


def _models_endpoint(base_url: str) -> str:
    if base_url.endswith("/v1"):
        return f"{base_url}/models"
    return f"{base_url}/v1/models"


def _chat_endpoint(base_url: str) -> str:
    if base_url.endswith("/v1"):
        return f"{base_url}/chat/completions"
    return f"{base_url}/v1/chat/completions"


def _parse_model_route_map(value: str) -> Tuple[Dict[str, object], str]:
    if not value or not value.strip():
        return {}, ""

    try:
        payload = json.loads(value)
    except Exception as exc:
        return {}, f"Invalid MODEL_ROUTE_MAP JSON: {exc}"

    if not isinstance(payload, dict):
        return {}, "MODEL_ROUTE_MAP must be a JSON object (model_id -> base_url or index)"

    return payload, ""


def _match_provider_target(target: str, providers: List[Tuple[str, str]]) -> Optional[Tuple[str, str]]:
    desired = (target or "").strip().rstrip("/")
    if not desired:
        return None

    # Be lenient about /v1 suffix in the mapping.
    candidates = {desired}
    if desired.endswith("/v1"):
        candidates.add(desired[:-3].rstrip("/"))
    else:
        candidates.add(f"{desired}/v1")

    for base_url, api_key in providers:
        normalized = base_url.rstrip("/")
        if normalized in candidates:
            return base_url, api_key
        if normalized.endswith("/v1") and normalized[:-3].rstrip("/") in candidates:
            return base_url, api_key
    return None


def _forced_routes_from_map(
    route_map: Dict[str, object], providers: List[Tuple[str, str]]
) -> Tuple[Dict[str, Tuple[str, str]], List[str]]:
    forced: Dict[str, Tuple[str, str]] = {}
    errors: List[str] = []

    for model_id, target in route_map.items():
        if not isinstance(model_id, str) or not model_id.strip():
            errors.append("MODEL_ROUTE_MAP contains a non-string model id key")
            continue

        if isinstance(target, int):
            if 0 <= target < len(providers):
                forced[model_id] = providers[target]
            else:
                errors.append(f"MODEL_ROUTE_MAP[{model_id}] index out of range: {target}")
            continue

        if isinstance(target, str):
            match = _match_provider_target(target, providers)
            if match:
                forced[model_id] = match
            else:
                errors.append(f"MODEL_ROUTE_MAP[{model_id}] base_url not found in API_BASE_URLS: {target}")
            continue

        errors.append(f"MODEL_ROUTE_MAP[{model_id}] must be a string base_url or an integer index")

    return forced, errors


async def _discover_model_routes(client: Any, providers: List[Tuple[str, str]]):
    routes = {}
    errors: List[str] = []

    async def probe(base_url: str, api_key: str):
        headers = {}
        if api_key:
            headers["Authorization"] = f"Bearer {api_key}"
        try:
            response = await client.get(_models_endpoint(base_url), headers=headers)
            response.raise_for_status()
            payload = response.json()
            models = payload.get("data", [])

            ids: List[str] = []
            for model in models:
                if isinstance(model, dict):
                    model_id = model.get("id")
                    if isinstance(model_id, str) and model_id:
                        ids.append(model_id)
            return base_url, api_key, ids, None
        except Exception as exc:
            return base_url, api_key, [], f"{base_url}: {exc}"

    results = await asyncio.gather(*(probe(url, key) for url, key in providers))
    for base_url, api_key, model_ids, error in results:
        if error:
            errors.append(error)
            continue
        for model_id in model_ids:
            routes.setdefault(model_id, (base_url, api_key))

    return routes, errors


async def _query_model(
    client: Any, model_id: str, base_url: str, api_key: str, question: str
) -> str:
    headers = {}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    try:
        response = await client.post(
            _chat_endpoint(base_url),
            headers=headers,
            json={
                "model": model_id,
                "messages": [{"role": "user", "content": question}],
                "max_tokens": 800,
            },
        )
        response.raise_for_status()
        payload = response.json()
        choices = payload.get("choices", [])
        if not choices:
            return f"## {model_id}\nError: empty choices in response"

        message = choices[0].get("message", {})
        content = message.get("content")
        if not content:
            return f"## {model_id}\nError: empty content in response"

        return f"## {model_id}\n{content}"
    except Exception as exc:
        return f"## {model_id}\nError: {exc}"


class Tools:
    def __init__(self):
        self.valves = Valves()
        # OpenWebUI loads Tools() instances (not modules), so expose Valves on the instance
        # to make /api/v1/tools/id/{id}/valves/spec return a schema.
        self.Valves = Valves

    async def council(self, question: str) -> str:
        """Query configured models and aggregate responses with per-model error handling."""
        model_ids = _split_csv(self.valves.MODELS, sep=",")
        if not model_ids:
            return "No model IDs configured in MODELS"

        providers = _provider_configs(self.valves.API_BASE_URLS, self.valves.API_KEYS)
        timeout = float(self.valves.TIMEOUT if self.valves.TIMEOUT else 120)

        route_map, route_map_error = _parse_model_route_map(self.valves.MODEL_ROUTE_MAP)
        forced_routes, forced_errors = _forced_routes_from_map(route_map, providers) if route_map else ({}, [])

        async with httpx.AsyncClient(timeout=timeout) as client:
            routes, discovery_errors = await _discover_model_routes(client, providers)
            if forced_routes:
                routes.update(forced_routes)
            tasks = []
            results = []

            for model_id in model_ids:
                route = routes.get(model_id)
                if not route:
                    results.append(
                        f"## {model_id}\nError: model not found on configured providers"
                    )
                    continue

                base_url, api_key = route
                tasks.append(_query_model(client, model_id, base_url, api_key, question))

            if tasks:
                results.extend(await asyncio.gather(*tasks))

        if route_map_error:
            results.append(f"## Route map\n- {route_map_error}")
        if forced_errors:
            errors_block = "\n".join(f"- {err}" for err in forced_errors)
            results.append(f"## Route map\n{errors_block}")

        if discovery_errors:
            errors_block = "\n".join(f"- {err}" for err in discovery_errors)
            results.append(f"## Provider discovery\n{errors_block}")

        if not results:
            return "No responses received from models"

        return "### Multi-Model Council Response\n\n" + "\n\n---\n\n".join(results)
PY_LLM_COUNCIL

    cat >"$scraper_file" <<'PY_WEB_SCRAPER'
#!/usr/bin/env python3
import json
import os
import urllib.parse

import requests
from bs4 import BeautifulSoup
from pydantic import BaseModel
from requests.exceptions import RequestException, SSLError


class Valves(BaseModel):
    MAX_LENGTH: int = 10000
    TIMEOUT: int = 30
    SSL_VERIFY: bool = True
    # Safer default: do not silently downgrade to verify=False on SSL errors.
    ALLOW_INSECURE_SSL_FALLBACK: bool = False


def _request_with_ssl_fallback(valves: Valves, method: str, url: str, **kwargs):
    timeout = kwargs.pop("timeout", valves.TIMEOUT)
    # IMPORTANT:
    # - Use verify=True by default so requests can honor REQUESTS_CA_BUNDLE / CURL_CA_BUNDLE.
    # - Do not force certifi.where() which breaks custom CA setups.
    verify = kwargs.pop("verify", None)
    if verify is None:
        verify = bool(valves.SSL_VERIFY)

    try:
        return requests.request(method, url, timeout=timeout, verify=verify, **kwargs)
    except SSLError:
        if verify is not False and valves.ALLOW_INSECURE_SSL_FALLBACK:
            return requests.request(method, url, timeout=timeout, verify=False, **kwargs)
        raise


def _ensure_searxng_json_url(url: str) -> str:
    parsed = urllib.parse.urlsplit(url)
    query = urllib.parse.parse_qs(parsed.query, keep_blank_values=True)

    # SearXNG returns HTML by default unless format=json is set.
    if "format" not in query or not query["format"] or not any(v for v in query["format"]):
        query["format"] = ["json"]

    encoded_query = urllib.parse.urlencode(query, doseq=True)
    return urllib.parse.urlunsplit((parsed.scheme, parsed.netloc, parsed.path, encoded_query, parsed.fragment))


def _searxng_search(valves: Valves, query: str, num_results: int):
    template = os.getenv("SEARXNG_QUERY_URL", "").strip()
    if not template:
        raise ValueError("SEARXNG_QUERY_URL is not set")

    encoded_query = urllib.parse.quote_plus(query)
    if "{query}" in template:
        searx_url = template.replace("{query}", encoded_query)
    else:
        separator = "&" if "?" in template else "?"
        searx_url = f"{template}{separator}q={encoded_query}"

    searx_url = _ensure_searxng_json_url(searx_url)

    headers = {"Accept": "application/json", "User-Agent": "OpenWebUI-WebScraper/1.0"}
    response = _request_with_ssl_fallback(valves, "GET", searx_url, headers=headers)
    response.raise_for_status()

    payload = response.json()
    raw_results = payload.get("results", [])
    results = []
    for item in raw_results:
        if len(results) >= num_results:
            break
        title = item.get("title") or "Untitled"
        link = item.get("url") or item.get("link") or ""
        snippet = item.get("content") or item.get("snippet") or ""
        if link:
            results.append({"title": title, "url": link, "snippet": snippet})

    return {
        "query": query,
        "source": "searxng",
        "results": results,
    }


def _duckduckgo_fallback(valves: Valves, query: str, num_results: int):
    headers = {"Accept": "application/json", "User-Agent": "OpenWebUI-WebScraper/1.0"}
    params = {
        "q": query,
        "format": "json",
        "no_redirect": "1",
        "no_html": "1",
        "skip_disambig": "1",
    }
    response = _request_with_ssl_fallback(
        valves,
        "GET",
        "https://api.duckduckgo.com/",
        headers=headers,
        params=params,
    )
    response.raise_for_status()

    payload = response.json()
    results = []

    abstract_text = payload.get("AbstractText") or ""
    abstract_url = payload.get("AbstractURL") or ""
    heading = payload.get("Heading") or "DuckDuckGo Instant Answer"
    if abstract_text or abstract_url:
        results.append({
            "title": heading,
            "url": abstract_url,
            "snippet": abstract_text,
        })

    related_topics = payload.get("RelatedTopics", [])
    for entry in related_topics:
        if len(results) >= num_results:
            break
        topics = entry.get("Topics") if isinstance(entry, dict) and "Topics" in entry else [entry]
        for topic in topics:
            if len(results) >= num_results:
                break
            if not isinstance(topic, dict):
                continue
            text = topic.get("Text") or ""
            url = topic.get("FirstURL") or ""
            if not text and not url:
                continue
            title = text.split(" - ")[0].strip() if " - " in text else (text[:120] or "DuckDuckGo Result")
            results.append({"title": title, "url": url, "snippet": text})

    return {
        "query": query,
        "source": "duckduckgo_instant_answer",
        "results": results[:num_results],
    }


class Tools:
    def __init__(self):
        self.valves = Valves()
        # OpenWebUI loads Tools() instances (not modules), so expose Valves on the instance
        # to make /api/v1/tools/id/{id}/valves/spec return a schema.
        self.Valves = Valves

    async def search_web(self, query: str, num_results: int = 5) -> str:
        searx_error = None

        try:
            payload = _searxng_search(self.valves, query, num_results)
            return json.dumps(payload, ensure_ascii=False, indent=2)
        except SSLError as error:
            searx_error = f"SSL error: {error}"
        except Exception as error:
            searx_error = str(error)

        try:
            fallback_payload = _duckduckgo_fallback(self.valves, query, num_results)
            fallback_payload["fallback_reason"] = (
                f"SearXNG unavailable: {searx_error}" if searx_error else "SearXNG unavailable"
            )
            return json.dumps(fallback_payload, ensure_ascii=False, indent=2)
        except SSLError as fallback_error:
            ssl_help = (
                "SSL verification failed. If you're using a corporate proxy or a private CA, set "
                "REQUESTS_CA_BUNDLE or SSL_CERT_FILE in the OpenWebUI container to point at your CA bundle. "
                "As a temporary workaround you can set ALLOW_INSECURE_SSL_FALLBACK=true (not recommended)."
            )
            return json.dumps(
                {
                    "query": query,
                    "source": "none",
                    "results": [],
                    "error": "Search failed due to SSL errors (SearXNG + fallback).",
                    "searxng_error": searx_error,
                    "fallback_error": f"SSL error: {fallback_error}",
                    "ssl_help": ssl_help,
                },
                ensure_ascii=False,
                indent=2,
            )
        except Exception as fallback_error:
            return json.dumps(
                {
                    "query": query,
                    "source": "none",
                    "results": [],
                    "error": "Search failed for both SearXNG and DuckDuckGo fallback",
                    "searxng_error": searx_error,
                    "fallback_error": str(fallback_error),
                },
                ensure_ascii=False,
                indent=2,
            )

    async def scrape_url(self, url: str, max_length: int = None) -> str:
        limit = max_length if max_length is not None else self.valves.MAX_LENGTH

        try:
            headers = {"User-Agent": "Mozilla/5.0 (compatible; OpenWebUI-WebScraper/1.0)"}
            response = _request_with_ssl_fallback(self.valves, "GET", url, headers=headers)
            response.raise_for_status()

            content_type = (response.headers.get("content-type") or "").lower()
            if "application/json" in content_type:
                try:
                    payload = response.json()
                    text = json.dumps(payload, ensure_ascii=False, indent=2)
                except Exception:
                    text = response.text or ""
            else:
                soup = BeautifulSoup(response.text, "html.parser")
                for tag in soup(["script", "style", "noscript"]):
                    tag.decompose()
                text = soup.get_text(separator="\n", strip=True)

            if len(text) > limit:
                text = text[:limit] + "..."

            return f"Content from {url}:\n\n{text}"
        except SSLError as error:
            return (
                "Error: SSL verification failed while fetching the URL.\n"
                "Fix: install/attach the correct CA bundle in the OpenWebUI container and set REQUESTS_CA_BUNDLE "
                "or SSL_CERT_FILE; or set SSL_VERIFY=false / ALLOW_INSECURE_SSL_FALLBACK=true as a temporary workaround.\n"
                f"Details: {error}"
            )
        except RequestException as error:
            return f"Error: {str(error)}"
        except Exception as error:
            return f"Error: {str(error)}"
PY_WEB_SCRAPER

    echo "$llm_file"
    echo "$scraper_file"
}

main() {
    require_cmds
    print_header

    TMP_DIR="$(mktemp -d)"

    signin

    local llm_file scraper_file
    llm_file="$(write_tool_sources | sed -n '1p')"
    scraper_file="$(write_tool_sources | sed -n '2p')"

    patch_tool_by_id "llm_council" "$llm_file"
    patch_tool_by_id "web_scraper" "$scraper_file"

    print_step "Done"
    print_success "Tool patches applied"
}

main "$@"
