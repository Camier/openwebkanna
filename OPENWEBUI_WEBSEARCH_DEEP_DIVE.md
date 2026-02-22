# OpenWebUI Web Search - Complete Deep Dive

> Comprehensive guide to web search integration, configuration, and optimization for OpenWebUI deployments.

---

## Table of Contents

1. [Overview & Architecture](#1-overview--architecture)
2. [Supported Search Engines](#2-supported-search-engines)
3. [Environment Variables](#3-environment-variables)
4. [Search Engine Setup Guides](#4-search-engine-setup-guides)
5. [Web Search + RAG Integration](#5-web-search--rag-integration)
6. [API Endpoints](#6-api-endpoints)
7. [Troubleshooting](#7-troubleshooting)
8. [Security Considerations](#8-security-considerations)
9. [Recommendations](#9-recommendations)

---

## 1. Overview & Architecture

### What is Web Search in OpenWebUI?

Web search enables LLMs to access real-time information from the internet, extending their knowledge beyond training data cutoffs. When enabled, OpenWebUI can:

- Search the web for current information
- Fetch and process web pages
- Extract content and add to RAG context
- Provide citations for web sources

### Architecture Flow

```
User Query
    │
    ├──▶ Local Knowledge Base (RAG)
    │         │
    │         ▼
    │    Vector Search Results
    │
    └──▶ Web Search Engine API
              │
              ▼
         Search Results
              │
              ▼
         Web Page Fetching
              │
              ▼
         Content Extraction
              │
              ▼
         Text Chunking
              │
              ▼
         Combined Local + Web Context
              │
              ▼
         LLM Response with Citations
```

### Key Components

| Component | Purpose |
|-----------|---------|
| **Search Engine** | Performs web queries (SearXNG, Brave, Google, etc.) |
| **Web Loader** | Fetches web page content |
| **Content Extractor** | Parses HTML to extract text |
| **Chunker** | Splits content for embedding |
| **Embedder** | Creates vector representations |
| **Reranker** | Re-ranks results by relevance |

---

## 2. Supported Search Engines

OpenWebUI supports **20+ search engines**:

### By Category

| Category | Engines |
|----------|---------|
| **Self-Hosted** | SearXNG, YaCy |
| **Privacy-Focused** | Brave, DuckDuckGo, Mojeek, Kagi |
| **Major Commercial** | Google PSE, Bing, SerpAPI |
| **AI-Optimized** | Tavily, Perplexity, Exa, Jina |
| **Other** | Serper, Serply, SearchAPI, Serpstack, Bocha |

### Quick Comparison

| Engine | Cost | Privacy | Setup | Best For |
|--------|------|---------|-------|----------|
| **SearXNG** | Free | ⭐⭐⭐⭐⭐ | Medium | Self-hosted, maximum privacy |
| **Brave** | Freemium | ⭐⭐⭐⭐⭐ | Easy | Privacy + reliability |
| **DuckDuckGo** | Free | ⭐⭐⭐⭐⭐ | Easy | Free option (rate limited) |
| **Google PSE** | Paid | ⭐⭐ | Easy | Google results |
| **Tavily** | Freemium | ⭐⭐⭐ | Easy | AI/RAG applications |
| **Perplexity** | Paid | ⭐⭐⭐ | Easy | AI answers with citations |

---

## 3. Environment Variables

### Core Web Search Variables

```bash
# Master toggle
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=searxng  # or: brave, google_pse, tavily, etc.

# Search behavior
WEB_SEARCH_RESULT_COUNT=3          # Number of results to retrieve
WEB_SEARCH_CONCURRENT_REQUESTS=10  # Parallel requests (0 = unlimited)

# Legacy alias (also supported)
ENABLE_WEBSEARCH=true
WEBSEARCH_ENGINE=searxng
```

### SearXNG Configuration (Recommended)

```bash
# Query URL - use host.docker.internal when running in Docker
SEARXNG_QUERY_URL=http://host.docker.internal:8888/search?q={query}&format=json
SEARXNG_LANGUAGE=en
```

### Commercial API Keys

```bash
# Brave Search
BRAVE_SEARCH_API_KEY=your_api_key

# Google PSE
GOOGLE_PSE_API_KEY=your_api_key
GOOGLE_PSE_ENGINE_ID=your_engine_id

# Bing
BING_SEARCH_V7_SUBSCRIPTION_KEY=your_key

# Tavily
TAVILY_API_KEY=your_api_key
TAVILY_EXTRACT_DEPTH=basic  # or 'comprehensive'

# Perplexity
PERPLEXITY_API_KEY=your_api_key
PERPLEXITY_MODEL=sonar

# SerpAPI
SERPAPI_API_KEY=your_api_key

# Add other keys as needed...
```

### Web Loader Configuration

```bash
# Engine selection
WEB_LOADER_ENGINE=playwright  # or: firecrawl, external

# Performance
WEB_LOADER_CONCURRENT_REQUESTS=10
WEB_LOADER_TIMEOUT=30

# Security
ENABLE_WEB_LOADER_SSL_VERIFICATION=true
```

### Proxy & Network

```bash
# HTTP proxy support
HTTP_PROXY=http://proxy.company.com:8080
HTTPS_PROXY=http://proxy.company.com:8080
NO_PROXY=localhost,127.0.0.1

# Trust proxy environment
WEB_SEARCH_TRUST_ENV=true
```

---

## 4. Search Engine Setup Guides

### 4.1 SearXNG (Self-Hosted) - RECOMMENDED

**Why SearXNG?**
- ✅ Completely free
- ✅ Maximum privacy
- ✅ Aggregates 70+ search engines
- ✅ No API rate limits
- ✅ Already configured in this repository

**Docker Compose Setup:**

```yaml
services:
  searxng:
    image: searxng/searxng:latest
    container_name: searxng
    ports:
      - "8888:8080"
    volumes:
      - ./searxng:/etc/searxng
    environment:
      - BASE_URL=http://localhost:8888/
    restart: unless-stopped
    networks:
      - openwebui_network
```

**OpenWebUI Configuration:**

```bash
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=searxng
SEARXNG_QUERY_URL=http://searxng:8080/search?q={query}&format=json
```

**Test:**
```bash
curl "http://localhost:8888/search?q=openwebui&format=json"
```

---

### 4.2 Brave Search API

**Why Brave?**
- ✅ Privacy-focused
- ✅ Independent web index
- ✅ 2,000 free queries/month
- ✅ Fast API response

**Setup:**
1. Visit [api.search.brave.com](https://api.search.brave.com)
2. Create account and get API key
3. Subscribe to "Data for AI" (free tier available)

**Configuration:**

```bash
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=brave
BRAVE_SEARCH_API_KEY=your_api_key
WEB_SEARCH_RESULT_COUNT=3
```

**Pricing:**
- Free: 2,000 queries/month
- Paid: $5 per 1,000 requests

---

### 4.3 DuckDuckGo (No API Key)

**Why DuckDuckGo?**
- ✅ No API key required
- ✅ Free
- ✅ Privacy-focused

**Setup:**
```bash
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=ddgs
DDGS_BACKEND=auto  # auto, html, or lite
```

**Limitations:**
- ⚠️ Unofficial API (~100 req/day limit)
- ⚠️ Rate limiting errors common
- ⚠️ May break without warning

---

### 4.4 Google Programmable Search Engine (PSE)

**Why Google PSE?**
- ✅ Google search results
- ✅ 100 free queries/day
- ✅ Custom search engines

**Setup:**
1. Go to [Google Cloud Console](https://console.cloud.google.com/)
2. Enable "Custom Search API"
3. Create API key
4. Go to [Programmable Search Engine](https://programmablesearchengine.google.com/)
5. Create search engine and get Search Engine ID (cx)

**Configuration:**

```bash
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=google_pse
GOOGLE_PSE_API_KEY=your_api_key
GOOGLE_PSE_ENGINE_ID=your_cx_id
```

**Pricing:**
- Free: 100 queries/day
- Paid: $5 per 1,000 queries

---

### 4.5 Tavily (AI-Optimized)

**Why Tavily?**
- ✅ Built for LLM applications
- ✅ Clean, AI-ready content
- ✅ 1,000 free credits/month

**Setup:**
1. Visit [tavily.com](https://tavily.com)
2. Sign up and get API key

**Configuration:**

```bash
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=tavily
TAVILY_API_KEY=your_api_key
TAVILY_EXTRACT_DEPTH=basic  # or 'comprehensive'
```

**Pricing:**
- Free: 1,000 credits/month
- Paid: $0.008/credit

---

## 5. Web Search + RAG Integration

### How It Works

When web search is enabled with RAG:

1. **Dual Retrieval**: Query triggers both local RAG and web search
2. **Web Processing**: Pages fetched, extracted, chunked like documents
3. **Unified Search**: Local + web embeddings searched together
4. **Hybrid Re-ranking**: Results combined and re-ranked
5. **Context Assembly**: Top-K results fed to LLM

### Configuration

```bash
# Enable both
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=searxng

# RAG settings
RAG_TOP_K=5                    # Total chunks (local + web)
RAG_TOP_K_RERANKER=3           # Final after re-ranking
ENABLE_RAG_HYBRID_SEARCH=true  # BM25 + vector search
RAG_SYSTEM_CONTEXT=true        # KV cache optimization

# Chunking (applies to web pages too)
CHUNK_SIZE=1500
CHUNK_OVERLAP=150
```

### Context Window Requirements

| Use Case | Minimum | Recommended |
|----------|---------|-------------|
| Local RAG only | 4096 | 8192 |
| Web search + RAG | 8192 | 16384+ |
| Large web pages | 16384 | 32768+ |

### Web Search as Additional Source

Web search results appear as additional documents in the RAG context:

```json
{
  "source": "web_search",
  "url": "https://example.com/page",
  "title": "Page Title",
  "engine": "searxng"
}
```

---

## 6. API Endpoints

### Get Web Search Configuration

```bash
GET /api/v1/retrieval/config
Authorization: Bearer <token>
```

**Response includes:**
- `ENABLE_WEB_SEARCH`
- `WEB_SEARCH_ENGINE`
- `WEB_SEARCH_RESULT_COUNT`
- SSL verification settings

### Execute Web Search

```bash
POST /api/v1/retrieval/process/web/search
Authorization: Bearer <token>
Content-Type: application/json

{
  "query": "search query string"
}
```

### Chat with Web Search

Enable via UI toggle or use `#` prefix:

```bash
POST /api/chat/completions
{
  "model": "llama3.1",
  "messages": [
    {"role": "user", "content": "#Search for recent AI developments"}
  ],
  "stream": true
}
```

---

## 7. Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| Web search not working | Check `ENABLE_WEB_SEARCH=true` |
| SSL certificate errors | Mount CA bundle, set `REQUESTS_CA_BUNDLE` |
| Container can't reach SearXNG | Use `host.docker.internal` not `localhost` |
| "No results found" | Check API key, increase timeout |
| Rate limit errors | Reduce `WEB_SEARCH_RESULT_COUNT` |
| Proxy errors | Set `WEB_SEARCH_TRUST_ENV=true` |

### Verification Commands

```bash
# Test SearXNG
curl "http://localhost:8888/search?q=test&format=json"

# Test from OpenWebUI container
docker exec openwebui curl "http://searxng:8080/health"

# Check network
docker network inspect openwebui_network
```

### Debug Mode

```bash
# Enable debug logging
docker logs -f openwebui

# Look for web search related messages
docker logs openwebui 2>&1 | grep -i "web\|search"
```

---

## 8. Security Considerations

### API Key Security

- Store in `.env` file (gitignored)
- Never commit to version control
- Rotate keys periodically
- Use environment-specific keys

### Rate Limiting

```bash
# Prevent abuse
WEB_SEARCH_CONCURRENT_REQUESTS=10
WEB_SEARCH_RESULT_COUNT=3  # Not too high
```

### Domain Filtering

```bash
# Restrict to trusted domains
WEB_SEARCH_DOMAIN_FILTER_LIST="wikipedia.org,arxiv.org,github.com"
```

### SSL Verification

```bash
# Keep enabled in production
ENABLE_WEB_LOADER_SSL_VERIFICATION=true
```

---

## 9. Recommendations

### For This Deployment

**Current Setup:** SearXNG is already configured in `docker-compose.yml`

**Recommended Configuration:**

```bash
# .env
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=searxng
SEARXNG_QUERY_URL=http://searxng:8080/search?q={query}&format=json
SEARXNG_LANGUAGE=en

WEB_SEARCH_RESULT_COUNT=3
WEB_SEARCH_CONCURRENT_REQUESTS=10

# RAG optimization
RAG_SYSTEM_CONTEXT=true
CHUNK_SIZE=1500
CHUNK_OVERLAP=150
```

### By Use Case

| Scenario | Recommended Engine | Why |
|----------|-------------------|-----|
| **Maximum Privacy** | SearXNG | Self-hosted, no data sharing |
| **Production/Reliability** | Brave | Good free tier, reliable API |
| **AI/RAG Applications** | Tavily | Built for LLMs |
| **Free Option** | DuckDuckGo | No API key (but rate limited) |
| **Google Results** | Google PSE | Direct Google access |

### Multi-Engine Strategy

For production deployments:

1. **Primary**: SearXNG (self-hosted, no limits)
2. **Backup**: Brave Search (if SearXNG fails)
3. **Development**: DuckDuckGo (no setup)

Switch via environment variables:
```bash
# Quick failover
WEB_SEARCH_ENGINE=brave
BRAVE_SEARCH_API_KEY=backup_key
```

---

## Related Documentation

- [Environment Variables Reference](openwebui_env_reference.md)
- [RAG Technical Reference](openwebui_rag_technical_reference.md)
- [Master Reference](OPENWEBUI_MASTER_REFERENCE.md)
- [API Examples](API_EXAMPLES.md)

---

*Compiled: 2026-02-18*
*OpenWebUI Version: v0.8.3*
