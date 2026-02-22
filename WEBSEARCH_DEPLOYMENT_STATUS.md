# Web Search Deployment Status

> SearXNG web search integration for OpenWebUI - DEPLOYED ✅

**Date:** 2026-02-18
**Status:** ✅ OPERATIONAL

---

## 🎯 Deployment Summary

| Component | Status | Details |
|-----------|--------|---------|
| **SearXNG Instance** | ✅ Running | /LAB/@ai_hub/searxng |
| **OpenWebUI Integration** | ✅ Configured | Using host.docker.internal:8888 |
| **API Connectivity** | ✅ Working | 1,100+ results for test query |
| **Recent Activity** | ✅ Active | 3 web searches in last logs |

---

## 🔧 Configuration

### Current Setup

**SearXNG (Local Instance):**
- **Location:** `/LAB/@ai_hub/searxng`
- **Port:** 8888 (localhost)
- **Status:** Running as systemd service
- **Health:** OK

**OpenWebUI (.env):**
```bash
ENABLE_WEB_SEARCH=true
ENABLE_WEBSEARCH=true
WEB_SEARCH_ENGINE=searxng
WEB_SEARCH_RESULT_COUNT=5
WEB_SEARCH_CONCURRENT_REQUESTS=10
SEARXNG_QUERY_URL=http://host.docker.internal:8888/search?q={query}&format=json
```

### SearXNG Settings

**File:** `searxng/settings.yml`
```yaml
use_default_settings: true
server:
  secret_key: "7bba336a0d250fbfd17fb2fb4869124165906ceab892b530e8ab99819dd60e3d"
  limiter: false
  image_proxy: true
  public_instance: false
  bind_address: "0.0.0.0"
  port: 8080
  base_url: "http://searxng:8080/"
search:
  formats:
    - html
    - json
```

---

## 🧪 Test Results

### API Test
```bash
Query: Sceletium tortuosum
Results: 1,100
First result: Sceletium tortuosum — Wikipédia
```

### Health Check
```bash
SearXNG Health: OK
OpenWebUI Integration: Active
```

### Recent Activity
OpenWebUI logs show web search activity:
- Web page embeddings generated
- Content extracted and chunked
- Results added to vector database
- RAG queries returning web results

---

## 🚀 How to Use Web Search

### Method 1: Web UI Toggle

1. Open OpenWebUI: http://localhost:3010
2. Start a new chat
3. Click the **globe icon** 🌐 (Web Search toggle)
4. Ask your question

### Method 2: # Prefix

Type `#` before your query:
```
#What are the latest research findings on Sceletium tortuosum?
```

### Method 3: API

```bash
curl -X POST http://localhost:3010/api/chat/completions \
  -H "Authorization: Bearer YOUR_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "your-model",
    "messages": [{"role": "user", "content": "#Search for recent AI developments"}]
  }'
```

---

## 📊 Performance Characteristics

| Metric | Value |
|--------|-------|
| Search Results | 5 per query (configurable) |
| Concurrent Requests | 10 |
| Response Time | ~1-3 seconds |
| Result Sources | 70+ aggregated engines |
| Cost | Free (self-hosted) |

---

## 🔗 Access URLs

| Service | URL |
|---------|-----|
| OpenWebUI | http://localhost:3010 |
| SearXNG API | http://localhost:8888 |
| SearXNG Health | http://localhost:8888/healthz |
| SearXNG Search | http://localhost:8888/search?q=test&format=json |

---

## 🛠️ Management Commands

### Check SearXNG Status
```bash
# Check process
ps aux | grep searxng

# Test API
curl "http://localhost:8888/search?q=openwebui&format=json"

# Check health
curl http://localhost:8888/healthz
```

### Check OpenWebUI Integration
```bash
# View web search logs
docker logs openwebui 2>&1 | grep -i "web\|search"

# Verify environment
docker exec openwebui env | grep -E "(WEB_SEARCH|SEARXNG)"
```

---

## ⚠️ Troubleshooting

### Issue: Web search not working
**Solution:**
1. Check SearXNG: `curl http://localhost:8888/healthz`
2. Check OpenWebUI can reach SearXNG:
   ```bash
   docker exec openwebui curl http://host.docker.internal:8888/healthz
   ```
3. Verify `.env` has correct `SEARXNG_QUERY_URL`

### Issue: Slow responses
**Solution:**
- Reduce `WEB_SEARCH_RESULT_COUNT` to 3
- Check SearXNG instance resources
- Verify network connectivity

### Issue: No results
**Solution:**
- Check SearXNG logs: `journalctl -u searxng -f`
- Verify SearXNG has working search backends
- Test directly: `curl "http://localhost:8888/search?q=test&format=json"`

---

## 🎓 Features Enabled

- ✅ Web search toggle in chat UI
- ✅ Automatic web page fetching
- ✅ Content extraction and chunking
- ✅ Embedding generation for web content
- ✅ Integration with RAG pipeline
- ✅ Source citations in responses
- ✅ Parallel local + web retrieval

---

## 📚 Related Documentation

- [Web Search Deep Dive](OPENWEBUI_WEBSEARCH_DEEP_DIVE.md)
- [Environment Variables](openwebui_env_reference.md)
- [RAG Technical Reference](openwebui_rag_technical_reference.md)
- [SearXNG Documentation](https://docs.searxng.org/)

---

## 🎉 Deployment Complete!

Web search is **fully operational** and ready for use in your ethnopharmacological research.

**Next Steps:**
1. Open http://localhost:3010
2. Start a chat with web search enabled
3. Ask questions about Sceletium tortuosum and other research topics
4. Web results will be automatically retrieved and integrated into responses

---

*Deployed: 2026-02-18*
*OpenWebUI Version: v0.8.3*
*SearXNG: Latest*
