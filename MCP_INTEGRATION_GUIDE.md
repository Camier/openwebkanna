# OpenWebUI MCP Integration - Complete Guide

## Executive Summary

OpenWebUI v0.8.3 natively supports **MCP (Model Context Protocol)** via Streamable HTTP transport. This guide provides everything needed to integrate MCP servers with your thesis research setup.

**Key Points:**
- MCP support introduced in v0.6.31, fully stable in v0.8.3
- **MCPO proxy required** for stdio-based MCP servers
- **`WEBUI_SECRET_KEY` mandatory** for OAuth persistence
- Use **MCP (Streamable HTTP)** connection type, not OpenAPI

---

## 1. Quick Start (5 Minutes)

### Step 1: Verify Prerequisites

```bash
# Check your .env has these set
WEBUI_SECRET_KEY=your-32-byte-or-longer-secret-key  # REQUIRED
OPENAI_API_BASE_URL=http://cliproxyapi:8317/v1      # Your existing setup
```

### Step 2: Add MCPO to Docker Compose

```yaml
# Add to docker-compose.yml
services:
  openwebui:
    # ... your existing config ...
    depends_on:
      - mcpo  # Add this dependency

  mcpo:
    image: ghcr.io/open-webui/mcpo:latest
    container_name: mcpo
    ports:
      - "8000:8000"
    volumes:
      - ./mcp-config.json:/app/config.json
    environment:
      - MCPO_HOST=0.0.0.0
      - MCPO_PORT=8000
      - MCPO_LOG_LEVEL=INFO
    networks:
      - openwebui_network
```

### Step 3: Create MCP Configuration

```bash
# Create mcp-config.json
cat > mcp-config.json << 'EOF'
{
  "mcpServers": {
    "time": {
      "command": "uvx",
      "args": ["mcp-server-time", "--local-timezone=Europe/Paris"]
    }
  }
}
EOF
```

### Step 4: Restart and Configure

```bash
docker compose down
docker compose up -d

# Wait 30 seconds for services to start
sleep 30
```

### Step 5: Add to OpenWebUI

1. Open OpenWebUI → Admin Settings → External Tools
2. Click **+ Add Server**
3. **Type**: MCP (Streamable HTTP)
4. **URL**: `http://mcpo:8000`
5. **Auth**: None (for local MCPO)
6. Save and test

---

## 2. Official MCP Documentation

### Primary Documentation
- **MCP Overview**: https://docs.openwebui.com/features/extensibility/mcp/
- **MCP Servers**: https://docs.openwebui.com/features/extensibility/plugin/tools/openapi-servers/mcp
- **Tools Overview**: https://docs.openwebui.com/features/extensibility/plugin/tools/

### MCPO Proxy (Official)
- **GitHub**: https://github.com/open-webui/mcpo
- **Purpose**: Converts stdio/SSE MCP servers to HTTP endpoints
- **Install**: `uvx mcpo` or `docker pull ghcr.io/open-webui/mcpo`

---

## 3. Recommended MCP Servers for Academic Research

### 🔴 Critical: Zotero MCP Server

**Why**: Connect your citation library to OpenWebUI

**Installation**:
```bash
# Option 1: Using pip (recommended)
pip install mcp-zotero

# Option 2: Using uvx
uvx mcp-zotero

# Get your Zotero credentials:
# 1. User ID: https://www.zotero.org/settings/keys (look for "Your userID")
# 2. API Key: https://www.zotero.org/settings/keys/new
```

**Configuration** (`mcp-config.json`):
```json
{
  "mcpServers": {
    "zotero": {
      "command": "uvx",
      "args": ["mcp-zotero"],
      "env": {
        "ZOTERO_API_KEY": "your-api-key",
        "ZOTERO_USER_ID": "your-user-id"
      }
    }
  }
}
```

**Use Cases**:
- Search your entire Zotero library
- Extract citations for thesis chapters
- Find papers you forgot you had
- Add papers directly from OpenWebUI

### 🟡 High Priority: File System Server

**Why**: Access your research papers and data

**Installation**:
```bash
# Using npx
npx -y @modelcontextprotocol/server-filesystem /path/to/your/thesis/data
```

**Configuration**:
```json
{
  "mcpServers": {
    "thesis-files": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-filesystem",
        "/LAB/@thesis/openwebui/data"
      ]
    }
  }
}
```

**Security Note**: Only expose directories you want AI to access

### 🟢 Medium Priority: Web Search (Exa)

**Why**: Find latest research papers

**Installation**:
```bash
# Sign up at exa.ai for API key
pip install exa-mcp-server

# Configuration
{
  "mcpServers": {
    "exa-search": {
      "command": "uvx",
      "args": ["exa-mcp-server"],
      "env": {
        "EXA_API_KEY": "your-exa-api-key"
      }
    }
  }
}
```

### 🟢 Medium Priority: PostgreSQL MCP Server

**Why**: Query your research database directly

**Installation**:
```bash
# Docker image available
docker run -p 8080:8080 \
  -e POSTGRES_CONNECTION_STRING="postgresql://user:pass@postgres:5432/openwebui" \
  call518/mcp-postgresql-ops:latest
```

---

## 4. Complete Configuration Examples

### Example 1: Minimal Academic Setup

```json
{
  "mcpServers": {
    "zotero": {
      "command": "uvx",
      "args": ["mcp-zotero"],
      "env": {
        "ZOTERO_API_KEY": "${ZOTERO_API_KEY}",
        "ZOTERO_USER_ID": "${ZOTERO_USER_ID}"
      }
    },
    "thesis-files": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-filesystem",
        "/LAB/@thesis/openwebui/data/pdfs"
      ]
    }
  }
}
```

### Example 2: Full Research Stack

```json
{
  "mcpServers": {
    "zotero": {
      "command": "uvx",
      "args": ["mcp-zotero"],
      "env": {
        "ZOTERO_API_KEY": "${ZOTERO_API_KEY}",
        "ZOTERO_USER_ID": "${ZOTERO_USER_ID}"
      }
    },
    "filesystem": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-filesystem",
        "/LAB/@thesis/openwebui/data"
      ]
    },
    "memory": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-memory"]
    },
    "exa-search": {
      "command": "uvx",
      "args": ["exa-mcp-server"],
      "env": {
        "EXA_API_KEY": "${EXA_API_KEY}"
      }
    }
  }
}
```

### Example 3: With Environment Variables

Create `.env.mcp`:
```bash
ZOTERO_API_KEY=your-zotero-api-key
ZOTERO_USER_ID=your-zotero-user-id
EXA_API_KEY=your-exa-api-key
```

Load in docker-compose:
```yaml
services:
  mcpo:
    image: ghcr.io/open-webui/mcpo:latest
    env_file:
      - .env.mcp
    volumes:
      - ./mcp-config.json:/app/config.json
```

---

## 5. Docker Compose Integration

### Complete Setup for Thesis Work

```yaml
version: '3.8'

services:
  # Your existing OpenWebUI service
  openwebui:
    image: ghcr.io/open-webui/open-webui:v0.8.3
    container_name: openwebui
    ports:
      - "3000:8080"
    volumes:
      - openwebui_data:/app/backend/data
    environment:
      - WEBUI_SECRET_KEY=${WEBUI_SECRET_KEY}
      - OPENAI_API_BASE_URL=http://cliproxyapi:8317/v1
      - OPENAI_API_KEY=${OPENAI_API_KEY}
    depends_on:
      - postgres
      - mcpo
    networks:
      - openwebui_network

  # Your existing PostgreSQL
  postgres:
    image: pgvector/pgvector:pg16
    # ... your existing config ...
    networks:
      - openwebui_network

  # Your existing CLIProxyAPI
  cliproxyapi:
    # ... your existing config ...
    networks:
      - openwebui_network

  # NEW: MCPO for MCP servers
  mcpo:
    image: ghcr.io/open-webui/mcpo:latest
    container_name: mcpo
    ports:
      - "8000:8000"
    volumes:
      - ./mcp-config.json:/app/config.json:ro
      - ./data:/data:ro  # Read-only access to your data
    environment:
      - MCPO_HOST=0.0.0.0
      - MCPO_PORT=8000
      - MCPO_LOG_LEVEL=INFO
      # Load your MCP secrets
      - ZOTERO_API_KEY=${ZOTERO_API_KEY}
      - ZOTERO_USER_ID=${ZOTERO_USER_ID}
      - EXA_API_KEY=${EXA_API_KEY}
    networks:
      - openwebui_network
    restart: unless-stopped

networks:
  openwebui_network:
    driver: bridge

volumes:
  openwebui_data:
  postgres_data:
```

---

## 6. Troubleshooting Guide

### Issue: "Error decrypting tokens"

**Cause**: Missing or changed `WEBUI_SECRET_KEY`

**Solution**:
```bash
# Add to .env
WEBUI_SECRET_KEY=$(openssl rand -hex 32)

# Restart
docker compose restart openwebui
```

### Issue: "Failed to connect to MCP server"

**Solutions**:
1. Check MCPO is running: `docker compose logs mcpo`
2. Verify URL: Should be `http://mcpo:8000` (not localhost)
3. Check connection type: Must be "MCP (Streamable HTTP)"
4. Test: `curl http://localhost:8000/openapi.json`

### Issue: MCP tools not appearing in chat

**Solutions**:
1. Enable tools in model settings:
   - Admin Settings → Models → [Your Model]
   - Set "Function Calling" to "Native"
2. Enable MCP tools in chat:
   - Click `+` next to input box
   - Select your MCP tools
3. Check MCPO logs: `docker compose logs mcpo | tail -50`

### Issue: OAuth tokens lost on restart

**Cause**: `WEBUI_SECRET_KEY` not set or changed

**Solution**: Set persistent secret key in `.env` and never change it

### Issue: Filesystem MCP can't access files

**Cause**: Wrong path or Docker volume not mounted

**Solution**:
```json
{
  "mcpServers": {
    "filesystem": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-filesystem",
        "/data"  // Must match Docker volume mount
      ]
    }
  }
}
```

And mount in docker-compose:
```yaml
volumes:
  - ./data:/data:ro  # Read-only for safety
```

---

## 7. Security Best Practices

### Critical Rules

1. **Never commit secrets to git**
   ```bash
   # Add to .gitignore
   echo ".env.mcp" >> .gitignore
   echo "mcp-config.json" >> .gitignore
   ```

2. **Use read-only mounts for filesystem**
   ```yaml
   volumes:
     - ./data:/data:ro  # :ro = read-only
   ```

3. **Restrict MCP server access**
   - Only enable MCP servers you actively use
   - Disable unused servers in OpenWebUI settings

4. **Keep MCPO updated**
   ```bash
   docker compose pull mcpo
   docker compose up -d mcpo
   ```

5. **Monitor MCP logs**
   ```bash
   docker compose logs -f mcpo
   ```

---

## 8. Usage Examples

### Example 1: Search Zotero Library

```
"Search my Zotero library for papers about Sceletium tortuosum
 alkaloids published after 2020"
```

### Example 2: Analyze Research Data

```
"Read the file /data/extractions/sceletium-analysis.csv and
 calculate the average mesembrine concentration across all samples"
```

### Example 3: Find Recent Papers

```
"Search for recent papers on ethnopharmacology of South African
 medicinal plants published in the last 6 months"
```

### Example 4: Generate Citations

```
"From my Zotero library, generate a bibliography of all papers
 about Sceletium in Vancouver format"
```

---

## 9. Next Steps

### Immediate Actions (Today)

1. [ ] Add MCPO to docker-compose.yml
2. [ ] Create mcp-config.json with time server (test)
3. [ ] Restart and verify connection
4. [ ] Get Zotero API credentials
5. [ ] Add Zotero MCP server

### This Week

6. [ ] Test Zotero integration
7. [ ] Add filesystem server for data access
8. [ ] Configure web search (Exa)
9. [ ] Document your MCP setup

### Ongoing

10. [ ] Monitor MCP logs periodically
11. [ ] Update MCP servers monthly
12. [ ] Back up MCP configuration

---

## 10. Quick Reference

### Essential Commands

```bash
# View MCPO logs
docker compose logs mcpo

# Restart MCP only
docker compose restart mcpo

# Test MCP connection
curl http://localhost:8000/openapi.json

# Edit MCP config
nano mcp-config.json
docker compose restart mcpo

# Update MCPO
docker compose pull mcpo
docker compose up -d mcpo
```

### Key URLs

- **OpenWebUI**: http://localhost:3000
- **MCPO OpenAPI**: http://localhost:8000/openapi.json
- **Zotero Keys**: https://www.zotero.org/settings/keys
- **Exa API**: https://exa.ai

---

## Summary

**MCP brings external tools into OpenWebUI**:
- ✅ Zotero: Citation management
- ✅ Filesystem: Access your research data
- ✅ Web Search: Find latest papers
- ✅ Databases: Query your research DB

**Critical Requirements**:
1. `WEBUI_SECRET_KEY` must be set
2. MCPO proxy for stdio servers
3. Use "MCP (Streamable HTTP)" connection type
4. Never commit API keys to git

**Start with**: Zotero MCP + Filesystem MCP

**Ready to proceed?** Follow the Quick Start guide above.
