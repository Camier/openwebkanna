# Archived Scripts

This directory contains deprecated scripts that are no longer actively maintained but kept for reference.

## vLLM Scripts (Deprecated)

**Status:** Deprecated as of February 2026
**Replaced by:** CLIProxyAPI-first architecture
**Migration:** Use CLIProxyAPI scripts (`start-cliproxyapi.sh`, `stop-cliproxyapi.sh`, etc.)

### Archived Scripts
- `start-vllm.sh` - Start vLLM server (deprecated)
- `stop-vllm.sh` - Stop vLLM server (deprecated)
- `restart-vllm.sh` - Restart vLLM server (deprecated)
- `check-vllm.sh` - Check vLLM health (deprecated)

### Why Deprecated?

vLLM was the original local LLM serving solution for this RAG deployment. However, it has been replaced by **CLIProxyAPI** for the following reasons:

1. **Better OAuth Integration**: CLIProxyAPI provides seamless OAuth-backed model aliases (Z.ai, MiniMax, etc.)
2. **Production-Ready**: CLIProxyAPI is a managed, battle-tested solution
3. **Simplified Architecture**: No need to manage Python environments, CUDA drivers, and model downloads
4. **Multi-Provider Support**: Single interface for multiple upstream providers

### For Reference Only

These scripts remain here for:
- Historical reference
- Fallback if needed during migration
- Understanding previous architecture

### Current Architecture

```
OpenWebUI → CLIProxyAPI → OAuth Providers
              (port 8317)    - Z.ai glm-5
                             - MiniMax MiniMax-M2.5
                             - (optional local vLLM)
```

See [README.md](../README.md) for current deployment instructions.
