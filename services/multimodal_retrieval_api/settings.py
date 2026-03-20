from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


def _env_path(name: str, default: Path | None = None) -> Path | None:
    raw = os.getenv(name)
    if raw:
        return Path(raw).expanduser().resolve()
    return default.resolve() if default is not None else None


def _env_value_from_file(path: Path | None, name: str) -> str | None:
    if path is None or not path.exists():
        return None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export ") :].strip()
        key, separator, value = line.partition("=")
        if not separator or key.strip() != name:
            continue
        normalized_value = value.strip()
        if (
            len(normalized_value) >= 2
            and normalized_value[0] == normalized_value[-1]
            and normalized_value[0] in {"'", '"'}
        ):
            normalized_value = normalized_value[1:-1]
        return normalized_value
    return None


@dataclass(frozen=True)
class ServiceSettings:
    app_name: str
    repo_root: Path
    extractions_root: Path
    qdrant_url: str
    qdrant_api_key: str | None
    compat_env_file: Path | None
    text_query_model_path: str | None
    text_query_device: str
    text_query_attn_implementation: str
    text_sparse_model_name: str
    rag_collection_name: str
    text_dense_vector_name: str
    text_sparse_vector_name: str
    vision_vector_name: str

    @classmethod
    def load(cls) -> "ServiceSettings":
        service_dir = Path(__file__).resolve().parent
        repo_root = _env_path(
            "MULTIMODAL_RETRIEVAL_API_REPO_ROOT",
            service_dir.parents[1],
        )
        if repo_root is None:
            raise RuntimeError("MULTIMODAL_RETRIEVAL_API_REPO_ROOT is required.")

        compat_env_file = _env_path("MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE")
        qdrant_url = os.getenv("MULTIMODAL_RETRIEVAL_API_QDRANT_URL")
        qdrant_api_key = os.getenv("MULTIMODAL_RETRIEVAL_API_QDRANT_API_KEY")
        if not qdrant_api_key:
            qdrant_api_key = os.getenv("QDRANT_API_KEY") or _env_value_from_file(compat_env_file, "QDRANT_API_KEY")
        text_query_model_path = os.getenv("MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH")
        if not text_query_model_path:
            text_query_model_path = os.getenv("NEMOTRON_MODEL_PATH") or _env_value_from_file(
                compat_env_file,
                "NEMOTRON_MODEL_PATH",
            )
        text_query_device = (
            os.getenv("MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_DEVICE")
            or os.getenv("NEMOTRON_DEVICE")
            or _env_value_from_file(compat_env_file, "NEMOTRON_DEVICE")
            or "cuda"
        )
        text_query_attn_implementation = (
            os.getenv("MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_ATTN_IMPLEMENTATION")
            or os.getenv("NEMOTRON_ATTN_IMPLEMENTATION")
            or _env_value_from_file(compat_env_file, "NEMOTRON_ATTN_IMPLEMENTATION")
            or "sdpa"
        )
        text_sparse_model_name = (
            os.getenv("MULTIMODAL_RETRIEVAL_API_TEXT_SPARSE_MODEL_NAME")
            or os.getenv("SPARSE_MODEL_NAME")
            or _env_value_from_file(compat_env_file, "SPARSE_MODEL_NAME")
            or "Qdrant/bm25"
        )

        extractions_root = _env_path(
            "MULTIMODAL_RETRIEVAL_API_EXTRACTIONS_ROOT",
            repo_root / "data" / "extractions",
        )
        if extractions_root is None:
            raise RuntimeError("MULTIMODAL_RETRIEVAL_API_EXTRACTIONS_ROOT is required.")

        return cls(
            app_name=os.getenv("MULTIMODAL_RETRIEVAL_API_NAME", "multimodal-retrieval-api"),
            repo_root=repo_root,
            extractions_root=extractions_root,
            qdrant_url=qdrant_url or "http://127.0.0.1:6335",
            qdrant_api_key=qdrant_api_key,
            compat_env_file=compat_env_file,
            text_query_model_path=text_query_model_path,
            text_query_device=text_query_device,
            text_query_attn_implementation=text_query_attn_implementation,
            text_sparse_model_name=text_sparse_model_name,
            rag_collection_name=os.getenv("MULTIMODAL_RETRIEVAL_API_RAG_COLLECTION", "rag_evidence"),
            text_dense_vector_name=os.getenv("MULTIMODAL_RETRIEVAL_API_TEXT_DENSE_VECTOR_NAME", "text_dense"),
            text_sparse_vector_name=os.getenv("MULTIMODAL_RETRIEVAL_API_TEXT_SPARSE_VECTOR_NAME", "text_sparse"),
            vision_vector_name=os.getenv("MULTIMODAL_RETRIEVAL_API_VISION_VECTOR_NAME", "vision_li"),
        )
