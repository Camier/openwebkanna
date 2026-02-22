# OpenWebUI RAG (Retrieval-Augmented Generation) Technical Reference

> **OpenWebUI Version:** v0.8.3 | **Last Updated:** 2026-02-18

**Related Documents:**
- [Master Reference](OPENWEBUI_MASTER_REFERENCE.md) - Overview and navigation
- [Environment Variables](openwebui_env_reference.md) - Complete configuration reference
- [Tools & Functions Guide](openwebui-tools-functions-guide.md) - Development patterns
- [Pipelines Guide](OPENWEBUI_PIPELINES_GUIDE.md) - Pipeline deployment
- [API Examples](API_EXAMPLES.md) - Repository-specific API usage

---

## Table of Contents
1. [Architecture Overview](#architecture-overview)
2. [Document Processing Pipeline](#document-processing-pipeline)
3. [Vector Database Options](#vector-database-options)
4. [Embedding Models](#embedding-models)
5. [Chunking Strategies](#chunking-strategies)
6. [Knowledge Base Management](#knowledge-base-management)
7. [RAG Query Flow](#rag-query-flow)
8. [Hybrid Search](#hybrid-search)
9. [Web Search Integration](#web-search-integration)
10. [Performance Optimization](#performance-optimization)
11. [Troubleshooting Guide](#troubleshooting-guide)
12. [Environment Variables Reference](#environment-variables-reference)

---

> **See Also:** [Environment Variables Reference](openwebui_env_reference.md#rag--knowledge-retrieval) for complete RAG configuration options, [Master Reference](OPENWEBUI_MASTER_REFERENCE.md#5-rag-system-deep-dive) for overview, and [Pipelines Guide](OPENWEBUI_PIPELINES_GUIDE.md) for custom RAG pipelines.

---

## Architecture Overview

Retrieval-Augmented Generation (RAG) in OpenWebUI combines LLMs with retrieved knowledge from external sources. The system retrieves relevant data from uploaded documents or knowledge bases, enhancing the quality and accuracy of responses.

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           OpenWebUI RAG Pipeline                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────────┐    ┌──────────────────────────────────────────────────┐   │
│  │   Document   │───▶│         Content Extraction Engine               │   │
│  │   Upload     │    │  (Tika / Docling / Default Text Extractor)      │   │
│  └──────────────┘    └──────────────────┬───────────────────────────────┘   │
│                                         │                                   │
│                                         ▼                                   │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                     Text Chunking & Splitting                        │   │
│  │  (Token-based / Character-based / Markdown Header Splitter)         │   │
│  └────────────────────────────────┬────────────────────────────────────┘   │
│                                   │                                         │
│                                   ▼                                         │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                     Embedding Generation                             │   │
│  │  (SentenceTransformers / OpenAI / External API)                     │   │
│  └────────────────────────────────┬────────────────────────────────────┘   │
│                                   │                                         │
│                                   ▼                                         │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                     Vector Database Storage                          │   │
│  │  (Chroma / FAISS / pgvector)                                        │   │
│  └────────────────────────────────┬────────────────────────────────────┘   │
│                                   │                                         │
│                                   ▼                                         │
│  ┌──────────────┐    ┌──────────────────────────────────────────────────┐   │
│  │ User Query   │───▶│              Retrieval Engine                    │   │
│  │  + Context   │    │  (Vector Search + BM25 Hybrid + Re-ranking)     │   │
│  └──────────────┘    └──────────────────┬───────────────────────────────┘   │
│                                         │                                   │
│                                         ▼                                   │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                Context Injection into LLM Prompt                     │   │
│  └────────────────────────────────┬────────────────────────────────────┘   │
│                                   │                                         │
│                                   ▼                                         │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                      LLM Response Generation                         │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Key Components

| Component | Purpose | Configuration Location |
|-----------|---------|------------------------|
| Content Extraction | Extract text from PDFs, DOCX, etc. | Admin Settings > Documents |
| Text Splitter | Divide documents into chunks | Admin Settings > Documents |
| Embedding Model | Convert text to vectors | Admin Settings > Documents |
| Vector Database | Store and search embeddings | Environment variables |
| Retrieval Engine | Find relevant chunks | Admin Settings > Documents |
| Re-ranker | Score and sort results | Admin Settings > Documents |

---

## Document Processing Pipeline

### 1. Document Upload

**Supported File Types:**
- PDF (`.pdf`)
- Word Documents (`.docx`, `.doc`)
- PowerPoint (`.pptx`, `.ppt`)
- Text files (`.txt`, `.md`, `.mdx`)
- HTML files (`.html`, `.htm`)
- And more via content extraction engines

**Upload Methods:**
1. **Chat Uploads**: Direct file upload in chat interface
2. **Knowledge Base Uploads**: Organized document collections
3. **Folder Uploads**: Batch document uploads
4. **API Uploads**: Programmatic file uploads via REST API

**Upload Limits (Configurable):**
| Limit | Environment Variable | Default |
|-------|---------------------|---------|
| Max file size | `RAG_FILE_MAX_SIZE` | Unlimited |
| Max file count (chat) | `RAG_FILE_MAX_COUNT` | Unlimited |
| Max file count (folder) | `FOLDER_MAX_FILE_COUNT` | 100 |
| Allowed extensions | `RAG_ALLOWED_FILE_EXTENSIONS` | All |

### 2. Content Extraction

**Available Extraction Engines:**

| Engine | Best For | Pros | Cons |
|--------|----------|------|------|
| **Apache Tika** | Complex documents (PDF, DOCX, PPT) | Handles 1000+ formats, metadata extraction | Requires separate container |
| **Docling** | Academic/scientific papers | Good table extraction, structured output | May need additional setup |
| **Default** | Simple text files | Built-in, no dependencies | Limited format support |

**Apache Tika Setup:**
```yaml
# docker-compose.yml addition
services:
  tika:
    image: apache/tika:latest-full
    container_name: tika
    restart: unless-stopped
    ports:
      - "9998:9998"
```

**Configuration:**
- Navigate to: `Admin Settings > Documents > Content Extraction Engine`
- Set Tika URL: `http://tika:9998/tika` (if using Docker)

### 3. Text Chunking

**Available Text Splitters:**

| Splitter | Description | Use Case |
|----------|-------------|----------|
| **Token (Tiktoken)** | Splits by token count | General purpose, OpenAI-compatible |
| **Character** | Splits by character count | Simple documents |
| **Markdown Header** | Splits by markdown headers | Structured markdown documents |

**Chunking Parameters:**

| Parameter | Description | Default | Recommended |
|-----------|-------------|---------|-------------|
| `CHUNK_SIZE` | Size of each chunk | 1500 chars | 500-2000 tokens |
| `CHUNK_OVERLAP` | Overlap between chunks | 150 chars | 10-20% of chunk size |
| `CHUNK_MIN_SIZE_TARGET` | Minimum chunk size | 0 | 50-60% of chunk size |

**Chunking Best Practices:**
- Use **Token-based splitting** for most documents
- Set overlap to preserve context between chunks
- Increase `CHUNK_MIN_SIZE_TARGET` to avoid tiny fragments (e.g., table of contents entries)
- Recommended: `CHUNK_MIN_SIZE_TARGET = 1000` (for 1500-2000 token chunks)

### 4. Embedding Generation

**Process Flow:**
1. Text chunks are passed to the embedding model
2. Model converts text to high-dimensional vectors
3. Vectors capture semantic meaning
4. Embeddings are stored in vector database

**Configuration:**
- Navigate to: `Admin Settings > Documents > Embedding`
- Select embedding engine and model
- After changing models, **reindex all documents**

---

## Vector Database Options

OpenWebUI supports three vector database backends:

### 1. Chroma (Default)

**Best for:** Local deployments, simple setups, quick prototyping

```bash
VECTOR_DB=chroma
```

| Pros | Cons |
|------|------|
| Zero configuration | Single-node only |
| Embedded (no separate service) | Limited scalability |
| Good performance for small-medium datasets | Not suitable for high-availability |
| Persistent storage | |

**Storage Location:**
- Docker: Inside container volume (persisted via `openwebui_data` volume)
- Local: SQLite database file

### 2. FAISS

**Best for:** Research, experimentation, CPU-optimized similarity search

```bash
VECTOR_DB=faiss
```

| Pros | Cons |
|------|------|
| Facebook's efficient similarity search | In-memory only (data lost on restart) |
| Optimized for CPU | No persistence by default |
| Fast approximate search | Requires manual index management |

### 3. pgvector (PostgreSQL)

**Best for:** Production deployments, multi-node setups, existing PostgreSQL infrastructure

```bash
VECTOR_DB=pgvector
PGVECTOR_DB_URL=postgresql://user:pass@host:5432/dbname
PGVECTOR_CREATE_EXTENSION=true
PGVECTOR_INITIALIZE_MAX_VECTOR_LENGTH=1536
```

| Pros | Cons |
|------|------|
| Production-grade reliability | Requires PostgreSQL setup |
| ACID compliance | Additional infrastructure |
| Scalable to multiple nodes | Network latency overhead |
| Full SQL support | |
| Concurrent access | |

**PostgreSQL Setup:**
```sql
-- Enable pgvector extension
CREATE EXTENSION IF NOT EXISTS vector;

-- Create table for embeddings (OpenWebUI manages this automatically)
-- Example dimension sizes based on embedding models:
-- - all-MiniLM-L6-v2: 384 dimensions
-- - text-embedding-3-small: 1536 dimensions
-- - text-embedding-3-large: 3072 dimensions
```

### Vector Database Comparison

| Feature | Chroma | FAISS | pgvector |
|---------|--------|-------|----------|
| Setup Complexity | Low | Low | Medium |
| Persistence | ✅ | ❌ | ✅ |
| Scalability | Limited | Limited | High |
| ACID Compliance | ✅ | ❌ | ✅ |
| Network Access | Local only | Local only | Remote |
| Best Use Case | Development | Research | Production |

---

## Embedding Models

### Supported Embedding Engines

| Engine | Description | Best For |
|--------|-------------|----------|
| **SentenceTransformers** | Local Hugging Face models | Privacy, offline use, cost control |
| **OpenAI** | OpenAI API embeddings | Quality, multilingual |
| **External** | Custom embedding API | Specialized models, existing infrastructure |

### Recommended Local Models

| Model | Dimensions | Size | Quality | Speed |
|-------|------------|------|---------|-------|
| `sentence-transformers/all-MiniLM-L6-v2` | 384 | ~80MB | Good | Fast |
| `sentence-transformers/all-mpnet-base-v2` | 768 | ~420MB | Better | Medium |
| `BAAI/bge-small-en-v1.5` | 384 | ~130MB | Good | Fast |
| `BAAI/bge-base-en-v1.5` | 768 | ~440MB | Better | Medium |
| `BAAI/bge-large-en-v1.5` | 1024 | ~1.3GB | Best | Slower |

### Cloud Embedding Models

| Model | Provider | Dimensions | Context |
|-------|----------|------------|---------|
| `text-embedding-3-small` | OpenAI | 1536 | 8191 tokens |
| `text-embedding-3-large` | OpenAI | 3072 | 8191 tokens |
| `text-embedding-ada-002` | OpenAI | 1536 | 8191 tokens |

### Configuration Example

```bash
# Local SentenceTransformers
RAG_EMBEDDING_MODEL=sentence-transformers/all-MiniLM-L6-v2

# OpenAI (requires API key)
RAG_EMBEDDING_MODEL=text-embedding-3-small
OPENAI_API_KEY=sk-...
```

### Custom Embedding Models

OpenWebUI supports loading custom embedding models in SentenceTransformers format:
1. Place model in accessible directory
2. Reference by path or Hugging Face model ID
3. Ensure model is in SBERT-compatible format

---

## Chunking Strategies

### Strategy Selection Guide

| Document Type | Recommended Splitter | Chunk Size | Overlap |
|---------------|---------------------|------------|---------|
| General text | Token (Tiktoken) | 500-1000 | 100-150 |
| Code | Token (Tiktoken) | 800-1200 | 100-200 |
| Markdown docs | Markdown Header | 1000-2000 | 150-200 |
| Legal/Contracts | Character | 1500-2000 | 200-300 |
| Academic papers | Token (Tiktoken) | 500-800 | 100 |

### Advanced Chunking Options

**Markdown Header Splitter:**
- Preserves document structure
- Splits on H1, H2, H3 headers
- Good for documentation and knowledge bases

**Token-based Splitter:**
- Uses Tiktoken (same as OpenAI)
- Ensures chunks fit within token limits
- Best for LLM context management

**Handling Small Fragments:**
- Set `CHUNK_MIN_SIZE_TARGET` to merge tiny chunks
- Prevents table of contents, headers from becoming separate chunks
- Recommended: `CHUNK_MIN_SIZE_TARGET = ~50-60% of CHUNK_SIZE`

---

## Knowledge Base Management

### Creating Knowledge Bases

1. Navigate to: `Workspace > Knowledge`
2. Click `+ Create a Knowledge Base`
3. Configure:
   - **Name**: Descriptive title
   - **Description**: Purpose and contents
   - **Visibility**: Private, Public, or Group-specific
   - **Access Groups**: Control who can use it

### Knowledge Base Features

| Feature | Description |
|---------|-------------|
| **Document Upload** | Drag & drop multiple files |
| **Folder Upload** | Bulk upload entire directories |
| **Auto-indexing** | Automatic embedding on upload |
| **Reindex** | Reprocess all documents (after model changes) |
| **Sync** | Update from source files |

### Permissions and Access Control

**Visibility Levels:**
- **Private**: Only creator can access
- **Public**: All users can access
- **Group**: Specific groups only

**RBAC Integration:**
- Admin can control knowledge base creation
- Group-based access management
- Model-level knowledge source assignment

### Using Knowledge Bases

**Method 1: Knowledge Command**
```
#Reference knowledge base in chat
#<knowledge_base_name>

Example:
#Open WebUI Documentation
How do I configure environment variables?
```

**Method 2: Model Assignment**
- Create custom model
- Assign knowledge base as source
- All queries automatically use knowledge base

**Method 3: URL/Document Reference**
```
#URL processing
#https://example.com/document

#File reference
#file_name.pdf
```

---

## RAG Query Flow

### Standard Query Pipeline

```
1. User Query Input
         │
         ▼
2. Query Embedding Generation
         │
         ▼
3. Vector Similarity Search (Top K)
         │
         ▼
4. [Optional] BM25 Keyword Search
         │
         ▼
5. [Optional] Re-ranking (CrossEncoder)
         │
         ▼
6. Top K Selection
         │
         ▼
7. Context Injection into Prompt
         │
         ▼
8. LLM Generation
         │
         ▼
9. Response with Source Citations
```

### Context Injection Modes

**Standard Mode (Default):**
- RAG context injected into user message
- Good for single-turn queries
- May cause KV cache invalidation in multi-turn chats

**System Context Mode:**
```bash
RAG_SYSTEM_CONTEXT=true
```
- RAG context injected into system message
- Fixed position throughout conversation
- Enables KV cache prefix caching
- Faster follow-up responses
- Better for multi-turn conversations

**Full Context Mode:**
- Bypasses embedding/retrieval
- Sends full document content
- Use for small documents only
- Requires large context window models

### Query Configuration

| Parameter | Description | Default |
|-----------|-------------|---------|
| `RAG_TOP_K` | Number of chunks to retrieve | 5 |
| `RAG_TOP_K_RERANKER` | Chunks after re-ranking | 5 |
| `RAG_RELEVANCE_THRESHOLD` | Minimum relevance score | 0.0 |

---

## Hybrid Search

OpenWebUI implements hybrid search combining:
1. **Dense Retrieval** (Vector similarity)
2. **Sparse Retrieval** (BM25 keyword matching)
3. **Re-ranking** (Cross-encoder scoring)

### Dense Search (Vector Similarity)

- Uses cosine similarity between query and document embeddings
- Captures semantic meaning
- Works across languages and synonyms
- May miss exact keyword matches

### Sparse Search (BM25)

- Traditional keyword-based ranking
- Ensures exact term matches
- Complements semantic search
- Good for technical terms, names, codes

### Re-ranking

**Purpose:** Score and sort retrieved chunks for maximum relevance

**Options:**

| Engine | Model | Requirements |
|--------|-------|--------------|
| Local (SentenceTransformers) | `BAI/bge-reranker-v2-m3` | GPU recommended |
| External | Jina AI Reranker | API key required |

**Configuration:**
```
Reranking Engine: Default (SentenceTransformers) or External
Reranking Model: BAAI/bge-reranker-v2-m3
Top K Reranker: 10 (candidates to re-rank)
Top K: 5 (final chunks to use)
```

**External Re-ranking (Jina AI):**
- URL: `https://api.jina.ai/v1/rerank`
- Model: `jina-reranker-v2-base-multilingual`
- No GPU required
- Trade-off: Network latency vs. local compute

### Hybrid Search Benefits

1. **Better Recall**: Dense + sparse covers more relevant documents
2. **Precision**: Re-ranking surfaces most relevant chunks
3. **Robustness**: Works across document types and query styles

---

## Web Search Integration

### Overview

OpenWebUI can combine local RAG with web search for comprehensive answers:

```
Query ──▶ Local Knowledge Base
     ──▶ Web Search (Brave, SearXNG, etc.)
              │
              ▼
      Retrieved Web Content
              │
              ▼
      Combined Context ──▶ LLM
```

### Supported Web Search Engines

| Engine | Description | Configuration |
|--------|-------------|---------------|
| **SearXNG** | Self-hosted, privacy-focused | `WEB_SEARCH_ENGINE=searxng` |
| **Brave** | Privacy-focused search API | Requires API key |
| **Google PSE** | Programmable Search Engine | Requires API key |
| **DuckDuckGo** | Privacy search | No API key |
| **SerpAPI** | Unified search API | Requires API key |
| **Serper** | Google Search API | Requires API key |

### SearXNG Configuration

```bash
# Enable web search
ENABLE_WEB_SEARCH=true
WEB_SEARCH_ENGINE=searxng

# SearXNG endpoint
SEARXNG_QUERY_URL=http://host.docker.internal:8888/search?q={query}&format=json
SEARXNG_LANGUAGE=en

# Search parameters
WEB_SEARCH_RESULT_COUNT=3
WEB_SEARCH_CONCURRENT_REQUESTS=10
```

### Web Search RAG Pipeline

1. User query triggers both local RAG and web search
2. Web pages are fetched and processed
3. Content extracted and chunked
4. Combined with local knowledge base results
5. Re-ranked together
6. Top results sent to LLM

### Web Search Considerations

| Challenge | Solution |
|-----------|----------|
| Large web pages | Use large context models (8192+ tokens) |
| Irrelevant content | Adjust relevance threshold |
| Rate limits | Use SearXNG self-hosted |
| Slow responses | Enable concurrent requests |

---

## Performance Optimization

### 1. Context Window Management

**Problem:** Default Ollama models limited to 2048 tokens

**Solutions:**

| Approach | Configuration |
|----------|--------------|
| Increase context length | Model Settings > Advanced Parameters > Context Length |
| Use Full Context Mode | For small documents only |
| Enable RAG System Context | `RAG_SYSTEM_CONTEXT=true` |
| Use cloud models | GPT-4, Claude 3 (128k context) |

**Recommended Context Sizes:**
- Local development: 4096 tokens
- Production: 8192+ tokens
- Web search: 16384+ tokens
- Complex analysis: 32768+ tokens

### 2. Embedding Optimization

| Strategy | Impact |
|----------|--------|
| Use smaller models | Faster indexing, less memory |
| Batch embeddings | `EMBEDDING_BATCH_SIZE` |
| GPU acceleration | CUDA for SentenceTransformers |
| External embeddings | Offload compute |

### 3. Retrieval Optimization

| Strategy | Configuration |
|----------|--------------|
| Adjust Top K | Balance precision vs. recall |
| Enable hybrid search | Better coverage |
| Use re-ranking | Improved relevance |
| Tune relevance threshold | Filter noise |

### 4. Storage Optimization

| Database | Optimization |
|----------|--------------|
| Chroma | Regular backups, volume management |
| FAISS | Not recommended for production |
| pgvector | Index tuning, connection pooling |

### 5. Caching Strategies

| Cache Type | Benefit |
|------------|---------|
| KV Cache Prefix | Faster follow-up responses |
| Embedding Cache | Avoid re-computation |
| Document Cache | Skip re-processing |

---

## Troubleshooting Guide

### Issue 1: Model "Can't See" Content

**Symptoms:** Hallucinations, generic responses, missing document info

**Causes & Solutions:**

| Cause | Solution |
|-------|----------|
| Poor content extraction | Switch to Apache Tika or Docling |
| Extraction preview shows blanks | Check extraction engine settings |
| Wrong file format support | Verify extraction engine capabilities |

**Steps:**
1. Go to `Admin Settings > Documents`
2. Check extraction engine
3. Upload test document and preview extracted content
4. If blank/missing, change extraction engine

### Issue 2: Only Small Part of Document Used

**Symptoms:** Incomplete answers, missing sections

**Causes:**
- Token limit too short (default 2048 for Ollama)
- Aggressive context trimming

**Solutions:**

| Option | Configuration |
|--------|--------------|
| Enable Full Context Mode | Admin Settings > Documents > Full Context Mode |
| Bypass Embedding | Admin Settings > Documents > Bypass Embedding |
| Increase context length | Model Settings > Advanced Parameters > Context Length |

**Warning:** Ensure model supports larger context before increasing

### Issue 3: Token Limit Too Short

**Symptoms:** Responses cut off, incomplete processing

**Context Limit Guidelines:**

| Use Case | Minimum | Recommended |
|----------|---------|-------------|
| Simple docs | 2048 | 4096 |
| Web search | 4096 | 8192-16384 |
| Complex analysis | 8192 | 16384-32768 |
| Large documents | 16384 | 32768-128000 |

**Ollama Context Extension:**
```
Admin Panel > Models > Settings (model) > Advanced Parameters
Modify context length to 8192+ (if model supports it)
```

**Alternative:** Use external LLMs with larger context:
- GPT-4o (128k tokens)
- Claude 3 (200k tokens)
- Gemini 1.5 (1M+ tokens)

### Issue 4: Poor Quality Embeddings

**Symptoms:** Irrelevant retrieval, bad matches

**Solutions:**

| Action | Steps |
|--------|-------|
| Change embedding model | Admin Settings > Documents > Embedding Model |
| Use higher quality model | e.g., `BAAI/bge-large-en-v1.5` |
| Reindex documents | After changing embedding model |
| Check model download | Ensure model downloaded completely |

### Issue 5: 'NoneType' object has no attribute 'encode'

**Cause:** Misconfigured or missing embedding model

**Solution:**
1. Go to: `Admin Settings > Documents > Embedding Model`
2. Save the embedding model again (even if already selected)
3. For external models, verify API accessibility
4. Check logs for download errors

### Issue 6: Fragmented/Tiny Chunks

**Symptoms:** Poor retrieval of table of contents, headers as separate chunks

**Solution:**
1. Go to: `Admin Settings > Documents`
2. Increase `Chunk Min Size Target`
3. Recommended: 1000 (or 50-60% of chunk size)
4. Reindex documents

### Issue 7: Slow Follow-up Responses

**Symptoms:** First response fast, subsequent queries slow

**Cause:** KV Cache invalidation due to context position shifting

**Solution:**
```bash
RAG_SYSTEM_CONTEXT=true
```

This injects RAG context into system message (fixed position), enabling:
- KV prefix caching
- Prompt caching (OpenAI/Vertex AI)
- Nearly instant follow-up responses

### Issue 8: API Upload "Content is Empty" Error

**Symptoms:** `400: The content provided is empty`

**Cause:** Race condition - file processing is asynchronous

**Solution:** Poll processing status before adding to knowledge base:

```python
import requests
import time

def wait_for_processing(token, file_id, timeout=300):
    url = f'http://localhost:3000/api/v1/files/{file_id}/process/status'
    headers = {'Authorization': f'Bearer {token}'}
    start_time = time.time()

    while time.time() - start_time < timeout:
        response = requests.get(url, headers=headers)
        status = response.json()

        if status.get('done'):
            return True
        time.sleep(1)

    return False
```

### Issue 9: Web Search Not Working

**Checklist:**
- [ ] Web search enabled in settings
- [ ] Search engine API key configured (if required)
- [ ] `WEB_SEARCH_RESULT_COUNT` not set too high
- [ ] SearXNG accessible (if self-hosted)
- [ ] `WEB_SEARCH_TRUST_ENV=true` (if behind proxy)

### Issue 10: Permission Issues

**Knowledge Base Access:**
- Verify knowledge base visibility settings
- Check group assignments
- Confirm user has required permissions

---

## Environment Variables Reference

### Core RAG Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `RAG_EMBEDDING_MODEL` | Embedding model name | `sentence-transformers/all-MiniLM-L6-v2` |
| `RAG_TOP_K` | Number of chunks to retrieve | `5` |
| `RAG_SYSTEM_CONTEXT` | Inject context into system message | `true` |
| `CHUNK_SIZE` | Document chunk size (characters) | `1500` |
| `CHUNK_OVERLAP` | Chunk overlap (characters) | `150` |
| `VECTOR_DB` | Vector database backend | `chroma` |

### Vector Database Variables

| Variable | Description | Required For |
|----------|-------------|--------------|
| `VECTOR_DB` | `chroma`, `faiss`, or `pgvector` | All |
| `PGVECTOR_DB_URL` | PostgreSQL connection string | pgvector |
| `PGVECTOR_CREATE_EXTENSION` | Auto-create vector extension | pgvector |
| `PGVECTOR_INITIALIZE_MAX_VECTOR_LENGTH` | Vector dimension | pgvector |

### Upload Limits

| Variable | Description | Default |
|----------|-------------|---------|
| `RAG_FILE_MAX_SIZE` | Maximum file size | Unlimited |
| `RAG_FILE_MAX_COUNT` | Max files per chat upload | Unlimited |
| `FOLDER_MAX_FILE_COUNT` | Max files in folder upload | `100` |
| `RAG_ALLOWED_FILE_EXTENSIONS` | Allowed file types | All |

### Web Search Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `ENABLE_WEB_SEARCH` | Enable web search | `false` |
| `WEB_SEARCH_ENGINE` | Search engine | `searxng` |
| `SEARXNG_QUERY_URL` | SearXNG endpoint | - |
| `WEB_SEARCH_RESULT_COUNT` | Results per query | `3` |
| `WEB_SEARCH_CONCURRENT_REQUESTS` | Parallel requests | `10` |
| `WEB_SEARCH_TRUST_ENV` | Trust proxy env | `false` |

### Content Extraction

| Variable | Description |
|----------|-------------|
| `TIKA_SERVER_URL` | Apache Tika endpoint |
| `CONTENT_EXTRACTION_ENGINE` | `tika`, `docling`, or `default` |

### Advanced Variables

| Variable | Description |
|----------|-------------|
| `RAG_EMBEDDING_ENGINE` | `sentence-transformers`, `openai`, `external` |
| `RAG_RERANKING_ENGINE` | `default`, `external` |
| `RAG_RERANKING_MODEL` | Re-ranker model name |
| `CHUNK_MIN_SIZE_TARGET` | Minimum chunk size target |
| `ENABLE_RAG_HYBRID_SEARCH` | Enable hybrid search |

**Note:** Must include `ENABLE_` prefix
| `RAG_FULL_CONTEXT` | Enable full context mode |
| `RAG_BYPASS_EMBEDDING` | Bypass embedding entirely |

---

## Best Practices Summary

### Document Processing
1. ✅ Use Apache Tika for complex documents
2. ✅ Set appropriate chunk sizes (500-2000 tokens)
3. ✅ Configure chunk overlap (10-20%)
4. ✅ Set minimum chunk size to avoid fragments

### Embedding
1. ✅ Use quality embedding models (BGE, MPNet)
2. ✅ Reindex after changing embedding models
3. ✅ Consider GPU acceleration for large document sets

### Vector Database
1. ✅ Use Chroma for development
2. ✅ Use pgvector for production
3. ✅ Avoid FAISS for persistent storage

### Query Optimization
1. ✅ Enable hybrid search for better recall
2. ✅ Use re-ranking for precision
3. ✅ Set RAG_SYSTEM_CONTEXT=true for multi-turn chats
4. ✅ Use adequate context windows (8192+ tokens)

### Troubleshooting
1. ✅ Always check content extraction preview
2. ✅ Verify embedding model is loaded
3. ✅ Monitor token limits
4. ✅ Test with high-quality models (GPT-4) to isolate issues

---

## Quick Reference Card

```
┌────────────────────────────────────────────────────────────┐
│              OpenWebUI RAG Quick Reference                 │
├────────────────────────────────────────────────────────────┤
│                                                            │
│  SETUP                                                     │
│  ├── Vector DB: chroma (dev) / pgvector (prod)            │
│  ├── Embedding: all-MiniLM-L6-v2 (fast) / bge-large (best)│
│  └── Extraction: Tika for complex docs                     │
│                                                            │
│  CHUNKING                                                  │
│  ├── Size: 500-2000 tokens                                 │
│  ├── Overlap: 10-20%                                       │
│  └── Min Size: 50-60% of chunk size                        │
│                                                            │
│  RETRIEVAL                                                 │
│  ├── Top K: 5-10 chunks                                    │
│  ├── Enable: Hybrid Search + Re-ranking                   │
│  └── Context: RAG_SYSTEM_CONTEXT=true                     │
│                                                            │
│  CONTEXT                                                   │
│  ├── Local models: 4096-8192 tokens                       │
│  ├── Web search: 8192-16384 tokens                        │
│  └── Cloud models: 128k+ tokens                           │
│                                                            │
└────────────────────────────────────────────────────────────┘
```

---

**OpenWebUI Version:** v0.8.3
**Last Updated:** 2026-02-18

*Generated from OpenWebUI documentation as of February 2026*
*Sources: docs.openwebui.com, GitHub discussions, community tutorials*
