"""Native text query encoding runtime for the one-collection RAG.

This runtime preserves the active dense+sparse query contract used by the
`rag_evidence` text lane while removing the live import dependency on the
external `wow/thesis_graph` package.

Dense query encoding:
- Hugging Face Transformers `AutoModel.from_pretrained(..., trust_remote_code=True)`
- existing Nemotron query-model contract via `forward_queries(...)`

Sparse query encoding:
- `fastembed.SparseTextEmbedding`

The lane remains Qdrant-native; this module only owns query vector creation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from qdrant_client import models


@dataclass(frozen=True)
class NativeTextQueryEncoderConfig:
    model_name_or_path: str | None
    device: str = "cuda"
    attn_implementation: str = "sdpa"
    sparse_model_name: str = "Qdrant/bm25"
    trust_remote_code: bool = True


@dataclass
class QueryEncodingBundle:
    dense: list[float] | None = None
    sparse: models.SparseVector | None = None
    diagnostics: dict[str, Any] | None = None


class QueryEncodingRuntimeError(RuntimeError):
    pass


@dataclass
class NativeTextQueryEncoderRuntime:
    config: NativeTextQueryEncoderConfig
    _model: Any | None = field(default=None, init=False, repr=False)
    _sparse_model: Any | None = field(default=None, init=False, repr=False)
    _cache: dict[str, QueryEncodingBundle] = field(default_factory=dict, init=False, repr=False)

    def encode_query(self, text: str) -> QueryEncodingBundle:
        query = (text or "").strip()
        if not query:
            return QueryEncodingBundle(
                diagnostics={
                    "warnings": ["Empty query received by native text query runtime."],
                    "dense_enabled": False,
                    "sparse_enabled": False,
                }
            )

        cached = self._cache.get(query)
        if cached is not None:
            return cached

        diagnostics: dict[str, Any] = {
            "runtime": "native_text_query_runtime",
            "model_name_or_path": self.config.model_name_or_path,
            "device": self.config.device,
            "attn_implementation": self.config.attn_implementation,
            "sparse_model_name": self.config.sparse_model_name,
            "warnings": [],
        }

        dense = self._encode_dense(query, diagnostics)
        sparse = self._encode_sparse(query, diagnostics)
        diagnostics["dense_enabled"] = dense is not None
        diagnostics["sparse_enabled"] = sparse is not None

        bundle = QueryEncodingBundle(
            dense=dense,
            sparse=sparse,
            diagnostics=diagnostics,
        )
        self._cache[query] = bundle
        return bundle

    def _encode_dense(self, query: str, diagnostics: dict[str, Any]) -> list[float] | None:
        if not self.config.model_name_or_path:
            raise QueryEncodingRuntimeError(
                "Dense query model path is not configured. "
                "Set MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH or NEMOTRON_MODEL_PATH."
            )

        try:
            import torch
            from transformers import AutoModel
        except Exception as exc:
            raise QueryEncodingRuntimeError(
                "Dense query runtime requires torch and transformers in the active environment."
            ) from exc

        if self._model is None:
            dtype = torch.float32
            if self.config.device.startswith("cuda") and torch.cuda.is_available():
                dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16
            model_kwargs: dict[str, Any] = {
                "trust_remote_code": self.config.trust_remote_code,
                "low_cpu_mem_usage": True,
                "torch_dtype": dtype,
                "attn_implementation": self.config.attn_implementation,
            }
            if self.config.device == "auto":
                model_kwargs["device_map"] = "auto"
            try:
                model = AutoModel.from_pretrained(
                    self.config.model_name_or_path,
                    **model_kwargs,
                ).eval()
            except Exception as exc:
                raise QueryEncodingRuntimeError(
                    "Dense query model load failed for "
                    f"{self.config.model_name_or_path!r}: {exc}"
                ) from exc
            if self.config.device != "auto":
                try:
                    model = model.to(self.config.device)
                except Exception as exc:
                    raise QueryEncodingRuntimeError(
                        "Dense query model move failed for "
                        f"device={self.config.device!r}: {exc}"
                    ) from exc
            self._model = model

        if not hasattr(self._model, "forward_queries"):
            raise QueryEncodingRuntimeError(
                "Dense query model does not expose forward_queries(...). "
                "The configured model is not compatible with the indexed text vectors."
            )

        with torch.inference_mode():
            query_embeddings = self._model.forward_queries([query], batch_size=1)
        multi = self._to_tensor_2d(query_embeddings[0], torch).detach().cpu()
        multi = torch.nn.functional.normalize(multi, p=2, dim=-1)
        dense = torch.nn.functional.normalize(multi.mean(dim=0, keepdim=True), p=2, dim=-1).squeeze(0)
        diagnostics["dense_vector_dim"] = int(dense.shape[-1])
        return dense.tolist()

    def _encode_sparse(self, query: str, diagnostics: dict[str, Any]) -> models.SparseVector | None:
        try:
            from fastembed import SparseTextEmbedding
        except Exception as exc:
            diagnostics["warnings"].append(
                "Sparse query runtime unavailable because fastembed is not installed."
            )
            return None

        try:
            if self._sparse_model is None:
                self._sparse_model = SparseTextEmbedding(self.config.sparse_model_name)
            sparse = next(self._sparse_model.query_embed([query]))
        except Exception as exc:
            diagnostics["warnings"].append(f"Sparse query encoding failed: {exc}")
            return None

        indices = sparse.indices.tolist() if hasattr(sparse.indices, "tolist") else list(sparse.indices)
        values = sparse.values.tolist() if hasattr(sparse.values, "tolist") else list(sparse.values)
        diagnostics["sparse_terms"] = len(indices)
        return models.SparseVector(
            indices=[int(value) for value in indices],
            values=[float(value) for value in values],
        )

    @staticmethod
    def _to_tensor_2d(value: Any, torch_module: Any) -> Any:
        tensor = value if isinstance(value, torch_module.Tensor) else torch_module.tensor(value)
        if tensor.dim() == 1:
            tensor = tensor.unsqueeze(0)
        return tensor.float()
