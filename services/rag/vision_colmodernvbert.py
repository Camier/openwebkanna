"""ColModernVBERT late-interaction runtime adapter.

Preferred runtime path:
- Hugging Face Transformers:
  - `ColModernVBertProcessor`
  - `ColModernVBertForRetrieval`

Fallback runtime path:
- Official ModernVBERT model card guidance, which currently uses
  `colpali` branch `vbert` with:
  - `colpali_engine.models.ColModernVBert`
  - `colpali_engine.models.ColModernVBertProcessor`

It produces multivector embeddings suitable for Qdrant's documented
late-interaction storage/query flow.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from .colmodernvbert_profile import ColModernVBERTImagePolicy, get_default_colmodernvbert_profile


@dataclass(frozen=True)
class ColModernVBERTConfig:
    """Runtime settings for ColModernVBERT inference."""

    model_name: str = get_default_colmodernvbert_profile().model_name
    batch_size: int = get_default_colmodernvbert_profile().embed_batch_size
    torch_dtype: str = "float32"
    attn_implementation: str | None = None
    device: str | None = None
    image_policy: ColModernVBERTImagePolicy = ColModernVBERTImagePolicy()


class ColModernVBERTVisionEmbedder:
    """Embed-and-upsert adapter implementing the vision embedder protocol."""

    def __init__(self, runtime: "ColModernVBERTLateInteractionRuntime" | None = None) -> None:
        self.runtime = runtime or ColModernVBERTLateInteractionRuntime()

    def embed_images(self, image_uris: list[str]) -> list[list[list[float]]]:
        """Return image multivectors for figure indexing."""

        return self.runtime.embed_images(image_uris)


class ColModernVBERTLateInteractionRuntime:
    """Late-interaction image/text encoder for ColModernVBERT."""

    def __init__(self, config: ColModernVBERTConfig | None = None) -> None:
        self.config = config or ColModernVBERTConfig()
        self._processor: Any | None = None
        self._model: Any | None = None
        self._torch: Any | None = None
        self._pil_image: Any | None = None
        self._runtime_kind: str | None = None

    def embed_images(self, image_uris: list[str]) -> list[list[list[float]]]:
        """Return one multivector per input image URI."""

        if not image_uris:
            return []

        processor, model, torch, pil_image, runtime_kind = self._load_runtime()
        outputs: list[list[list[float]]] = []
        for batch in _chunked(image_uris, self.config.batch_size):
            images = [self._prepare_image(self._open_image(uri, pil_image)) for uri in batch]
            inputs = (
                processor.process_images(images)
                if runtime_kind == "colpali_engine"
                else processor(images=images, return_tensors="pt")
            )
            inputs = self._move_to_device(inputs, model.device)
            with torch.inference_mode():
                raw_embeddings = model(**inputs)
            outputs.extend(self._to_multivectors(raw_embeddings))
        return outputs

    def embed_queries(self, texts: list[str]) -> list[list[list[float]]]:
        """Return one multivector query embedding per input text string."""

        if not texts:
            return []

        processor, model, torch, _, runtime_kind = self._load_runtime()
        outputs: list[list[list[float]]] = []
        for batch in _chunked(texts, self.config.batch_size):
            inputs = (
                processor.process_queries(batch)
                if runtime_kind == "colpali_engine"
                else processor(text=batch, return_tensors="pt")
            )
            inputs = self._move_to_device(inputs, model.device)
            with torch.inference_mode():
                raw_embeddings = model(**inputs)
            outputs.extend(self._to_multivectors(raw_embeddings))
        return outputs

    def _load_runtime(self) -> tuple[Any, Any, Any, Any, str]:
        if (
            self._model is not None
            and self._processor is not None
            and self._torch is not None
            and self._pil_image is not None
            and self._runtime_kind is not None
        ):
            return self._processor, self._model, self._torch, self._pil_image, self._runtime_kind

        try:
            import torch
            from PIL import Image
        except ImportError as exc:
            raise RuntimeError(
                "ColModernVBERT requires torch and Pillow in the active environment."
            ) from exc

        device = self.config.device or ("cuda" if torch.cuda.is_available() else "cpu")
        try:
            from transformers import ColModernVBertForRetrieval, ColModernVBertProcessor

            processor = ColModernVBertProcessor.from_pretrained(self.config.model_name)
            model_kwargs: dict[str, Any] = {
                "torch_dtype": _resolve_torch_dtype(torch, self.config.torch_dtype, device=device),
            }
            if self.config.attn_implementation:
                model_kwargs["attn_implementation"] = self.config.attn_implementation
            model = ColModernVBertForRetrieval.from_pretrained(
                self.config.model_name,
                **model_kwargs,
            )
            runtime_kind = "transformers"
        except Exception:
            try:
                from colpali_engine.models import ColModernVBert, ColModernVBertProcessor
            except ImportError as exc:
                raise RuntimeError(
                    "ColModernVBERT requires either Hugging Face Transformers with "
                    "ColModernVBert support or the official colpali_engine vbert branch."
                ) from exc

            processor = ColModernVBertProcessor.from_pretrained(self.config.model_name)
            model = ColModernVBert.from_pretrained(
                self.config.model_name,
                torch_dtype=_resolve_torch_dtype(torch, self.config.torch_dtype, device=device),
                trust_remote_code=True,
            )
            runtime_kind = "colpali_engine"

        model = model.to(device).eval()

        self._processor = processor
        self._model = model
        self._torch = torch
        self._pil_image = Image
        self._runtime_kind = runtime_kind
        return processor, model, torch, Image, runtime_kind

    def _open_image(self, image_uri: str, pil_image: Any) -> Any:
        parsed = urlparse(image_uri)
        if parsed.scheme in {"", "file"}:
            local_path = Path(parsed.path if parsed.scheme == "file" else image_uri)
            return pil_image.open(local_path)
        raise ValueError(f"Unsupported image URI for ColModernVBERT indexing: {image_uri}")

    def _prepare_image(self, image: Any) -> Any:
        image_policy = self.config.image_policy
        if image_policy.convert_to_rgb:
            image = image.convert("RGB")
        return _resize_image(
            image,
            max_long_edge=image_policy.max_long_edge,
            min_short_edge=image_policy.min_short_edge,
            upscale_small_images=image_policy.upscale_small_images,
        )

    def _move_to_device(self, batch_inputs: Any, device: str) -> Any:
        if hasattr(batch_inputs, "to"):
            return batch_inputs.to(device)
        moved: dict[str, Any] = {}
        for key, value in batch_inputs.items():
            moved[key] = value.to(device) if hasattr(value, "to") else value
        return moved

    def _to_multivectors(self, raw_embeddings: Any) -> list[list[list[float]]]:
        tensor = getattr(raw_embeddings, "embeddings", raw_embeddings)
        if hasattr(tensor, "detach"):
            tensor = tensor.detach().cpu()
        return [item.tolist() for item in tensor]


def probe_colmodernvbert_runtime(model_name: str) -> list[str]:
    """Return non-fatal readiness warnings for the ColModernVBERT runtime."""

    warnings: list[str] = []
    try:
        import torch  # noqa: F401
    except ImportError:
        warnings.append("ColModernVBERT visual lane requires torch in the active environment.")
        return warnings

    try:
        from PIL import Image  # noqa: F401
    except ImportError:
        warnings.append("ColModernVBERT visual lane requires Pillow in the active environment.")

    try:
        from transformers import ColModernVBertForRetrieval, ColModernVBertProcessor  # noqa: F401
    except ImportError:
        try:
            from colpali_engine.models import ColModernVBert, ColModernVBertProcessor  # noqa: F401
        except ImportError:
            warnings.append(
                "ColModernVBERT visual lane requires either Hugging Face Transformers with "
                "ColModernVBert support or the official colpali_engine vbert branch."
            )

    if not model_name:
        warnings.append("ColModernVBERT visual lane model name is empty.")
    return warnings


def _resolve_torch_dtype(torch: Any, dtype_name: str, *, device: str) -> Any:
    normalized = (dtype_name or "float32").lower()
    if device == "cpu":
        return torch.float32
    if normalized in {"bf16", "bfloat16"}:
        return torch.bfloat16
    if normalized in {"fp16", "float16", "half"}:
        return torch.float16
    return torch.float32


def _resize_image(image: Any, *, max_long_edge: int, min_short_edge: int, upscale_small_images: bool) -> Any:
    width, height = image.size
    long_edge = max(width, height)
    short_edge = min(width, height)

    if long_edge > max_long_edge:
        scale = max_long_edge / float(long_edge)
    elif upscale_small_images and short_edge < min_short_edge:
        scale = min_short_edge / float(short_edge)
    else:
        scale = 1.0

    if scale == 1.0:
        return image

    resized = (max(1, int(round(width * scale))), max(1, int(round(height * scale))))
    return image.resize(resized)


def _chunked(items: list[str], batch_size: int) -> list[list[str]]:
    return [items[i : i + batch_size] for i in range(0, len(items), batch_size)]
