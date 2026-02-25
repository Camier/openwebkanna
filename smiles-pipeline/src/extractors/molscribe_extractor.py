"""
MolScribe OCSR extractor.

Primary backend for molecular structure extraction from images.
Uses Swin Transformer architecture with character-level attention.

References:
- Qian et al. (2023): MolScribe JCIM, DOI: 10.1021/acs.jcim.2c01480
- GitHub: https://github.com/thomas0809/MolScribe
- License: MIT (commercial-safe)
"""

from typing import Dict, Any, Optional, List, Union
from pathlib import Path
import numpy as np
from PIL import Image


class MolScribeExtractor:
    """
    MolScribe-based optical chemical structure recognition.

    Provides GPU-accelerated molecular structure extraction with
    ~93% F1 accuracy on single molecules.
    """

    def __init__(
        self,
        model_name: str = "swin_base_char_aux_1m.pth",
        device: str = "cuda",
        confidence_threshold: float = 0.5,
        batch_size: int = 32,
    ):
        """
        Initialize MolScribe extractor.

        Args:
            model_name: Model checkpoint name
            device: "cuda" or "cpu"
            confidence_threshold: Minimum confidence to accept
            batch_size: Batch size for GPU inference
        """
        self.model_name = model_name
        self.device = self._select_device(device)
        self.confidence_threshold = confidence_threshold
        self.batch_size = batch_size
        self.predictor = None

    def _select_device(self, requested: str) -> str:
        """Select available device with fallback."""
        if requested == "cuda":
            try:
                import torch

                if torch.cuda.is_available():
                    return "cuda"
            except ImportError:
                pass
        return "cpu"

    def load_model(self):
        """Load MolScribe model weights."""
        from molscribe import MolScribe
        from pathlib import Path

        # Auto-download checkpoint from HuggingFace
        checkpoint_dir = Path.home() / ".cache" / "molscribe"
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        checkpoint_path = checkpoint_dir / self.model_name

        if not checkpoint_path.exists():
            from torch.hub import load_state_dict_from_url
            url = f"https://huggingface.co/yujieq/MolScribe/resolve/main/{self.model_name}"
            print(f"Downloading MolScribe checkpoint to {checkpoint_path}...")
            load_state_dict_from_url(url, model_dir=str(checkpoint_dir), map_location="cpu")

        self.predictor = MolScribe(model_path=str(checkpoint_path), device=self.device)

    def is_loaded(self) -> bool:
        """Check if model is loaded."""
        return self.predictor is not None

    def extract(
        self,
        image: Union[str, Path, Image.Image, np.ndarray],
        return_confidence: bool = True,
    ) -> Dict[str, Any]:
        """
        Extract SMILES from a chemical structure image.

        Args:
            image: Image path, PIL Image, or numpy array
            return_confidence: Include confidence scores

        Returns:
            Dictionary with extraction results:
            {
                'smiles': str or None,
                'confidence': float or None,
                'success': bool,
                'error': str or None,
                'backend': 'molscribe',
            }
        """
        result = {
            "smiles": None,
            "confidence": None,
            "success": False,
            "error": None,
            "backend": "molscribe",
        }

        # Lazy load model on first use
        if not self.is_loaded():
            try:
                self.load_model()
            except Exception as e:
                result["error"] = f"model_load_failed: {str(e)}"
                return result

        try:
            # Load image if path
            if isinstance(image, (str, Path)):
                image = Image.open(image).convert("RGB")
                image = np.array(image)  # Convert to numpy for MolScribe
            elif isinstance(image, Image.Image):
                image = np.array(image)  # Convert PIL to numpy
            # If already numpy, use as-is

            # Run prediction (MolScribe expects numpy arrays)
            prediction = self.predictor.predict_image(np.array(image) if not isinstance(image, np.ndarray) else image)

            smiles = prediction.get("smiles")
            confidence = prediction.get("confidence", 0.0)

            # Apply confidence threshold
            if smiles and confidence >= self.confidence_threshold:
                result["smiles"] = smiles
                result["confidence"] = confidence
                result["success"] = True
            elif smiles:
                result["smiles"] = smiles
                result["confidence"] = confidence
                result["success"] = True
                result["low_confidence"] = True
            else:
                result["error"] = "no_structure_detected"

        except Exception as e:
            result["error"] = f"extraction_failed: {str(e)}"

        return result

    def extract_batch(
        self,
        images: List[Union[str, Path, Image.Image]],
    ) -> List[Dict[str, Any]]:
        """
        Extract SMILES from multiple images in batch.

        Args:
            images: List of images (paths or PIL Images)

        Returns:
            List of extraction results
        """
        results = []

        # Process in batches
        for i in range(0, len(images), self.batch_size):
            batch = images[i : i + self.batch_size]
            batch_results = self._process_batch(batch)
            results.extend(batch_results)

        return results

    def _process_batch(
        self,
        images: List[Union[str, Path, Image.Image]],
    ) -> List[Dict[str, Any]]:
        """Process a single batch of images."""
        if not self.is_loaded():
            self.load_model()

        # Load all images
        loaded_images = []
        for img in images:
            if isinstance(img, (str, Path)):
                loaded_images.append(Image.open(img).convert("RGB"))
            elif isinstance(img, Image.Image):
                loaded_images.append(img.convert("RGB"))
            else:
                loaded_images.append(Image.fromarray(img))

        # Batch prediction
        try:
            predictions = self.predictor.predict_images(loaded_images)
        except Exception:
            # Fallback to individual prediction
            return [self.extract(img) for img in images]

        # Format results
        results = []
        for pred in predictions:
            results.append(
                {
                    "smiles": pred.get("smiles"),
                    "confidence": pred.get("confidence", 0.0),
                    "success": pred.get("smiles") is not None,
                    "error": None,
                    "backend": "molscribe",
                }
            )

        return results

    def get_model_info(self) -> Dict[str, Any]:
        """Get model metadata."""
        return {
            "backend": "molscribe",
            "model": self.model_name,
            "device": self.device,
            "confidence_threshold": self.confidence_threshold,
            "batch_size": self.batch_size,
            "loaded": self.is_loaded(),
            "license": "MIT",
            "accuracy": {
                "single_molecule_f1": 0.93,
                "stereochemistry": 0.85,
            },
        }


def extract_with_molscribe(
    image: Union[str, Path, Image.Image],
    **kwargs,
) -> Dict[str, Any]:
    """
    Convenience function for quick extraction.

    Args:
        image: Image to process
        **kwargs: Passed to MolScribeExtractor constructor

    Returns:
        Extraction result dictionary
    """
    extractor = MolScribeExtractor(**kwargs)
    return extractor.extract(image)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python molscribe_extractor.py <image_path>")
        sys.exit(1)

    image_path = sys.argv[1]
    result = extract_with_molscribe(image_path)

    print(f"Backend: {result['backend']}")
    print(f"Success: {result['success']}")
    print(f"SMILES: {result['smiles']}")
    print(f"Confidence: {result.get('confidence', 'N/A')}")
    if result.get("error"):
        print(f"Error: {result['error']}")
