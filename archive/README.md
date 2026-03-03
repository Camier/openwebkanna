# Archived Scripts

This directory contains deprecated scripts that are no longer actively maintained but kept for reference.

## vLLM Scripts (Deprecated)

**Status:** Deprecated as of February 2026
**Replaced by:** LiteLLM-first architecture
**Migration:** Use standard stack operations (`deploy.sh`, `status.sh`, baseline tests). CLIProxyAPI remains optional legacy sidecar only.

### Archived Scripts
- `start-vllm.sh` - Start vLLM server (deprecated)
- `stop-vllm.sh` - Stop vLLM server (deprecated)
- `restart-vllm.sh` - Restart vLLM server (deprecated)
- `check-vllm.sh` - Check vLLM health (deprecated)

### Why Deprecated?

vLLM was the original local LLM serving solution for this RAG deployment. It is now non-default, with LiteLLM as the reference upstream and CLIProxyAPI kept only for legacy sidecar workflows.

1. **LiteLLM-first Operations**: the repository defaults are standardized on an OpenAI-compatible LiteLLM endpoint.
2. **Simplified Runtime Expectations**: fallback paths remain available but are no longer required for normal operation.
3. **Legacy OAuth Optionality**: CLIProxyAPI alias workflows remain possible without being the primary control plane.
4. **Reduced Drift**: docs/tests now validate the same default upstream path.

### For Reference Only

These scripts remain here for:
- Historical reference
- Fallback if needed during migration
- Understanding previous architecture

### Current Architecture

```
OpenWebUI → LiteLLM → providers
         \-> optional CLIProxyAPI sidecar (legacy)
         \-> optional local vLLM fallback
```

See [README.md](../README.md) for current deployment instructions.
