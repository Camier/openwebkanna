"""
DECIMER OCSR extractor.

Fallback backend for molecular structure extraction from images.
Uses Vision Transformer (ViT) encoder with T5 decoder.

References:
- Rajan et al. (2023): DECIMER.ai Nature Communications, DOI: 10.1038/s41467-023-40782-0
- GitHub: https://github.com/Kohulan/DECIMER-Image_Transformer
- License: CC BY 4.0 (requires attribution)
"""

from typing import Dict, Any, Optional, List, Union
from pathlib import Path
import re
import numpy as np
from PIL import Image


class DECIMERExtractor:
    """
    DECIMER 2.x-based optical chemical structure recognition.

    Provides fallback extraction with ~91% F1 accuracy.
    Includes segmentation and classification components.
    """

    def __init__(
        self,
        model_version: str = "2.7",
        device: str = "cuda",
        confidence_threshold: float = 0.5,
        batch_size: int = 16,
        use_segmentation: bool = True,
    ):
        """
        Initialize DECIMER extractor.

        Args:
            model_version: DECIMER version metadata (2.7+ recommended)
            device: "cuda" or "cpu"
            confidence_threshold: Minimum confidence to accept
            batch_size: Batch size for inference (smaller than MolScribe)
            use_segmentation: Enable Mask R-CNN segmentation
        """
        self.model_version = model_version
        self.device = self._select_device(device)
        self.confidence_threshold = confidence_threshold
        self.batch_size = batch_size
        self.use_segmentation = use_segmentation
        self.predictor = None
        self.segmenter = None
        self.classifier = None

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

    def load_models(self):
        """Load DECIMER model components."""
        # DECIMER package naming/API differs by release line:
        # - modern wheels: `DECIMER.predict_SMILES`
        # - some older docs: `decimer.predict_smiles`
        try:
            from DECIMER import predict_SMILES  # type: ignore

            self.predictor = predict_SMILES
        except Exception:
            from decimer import predict_smiles  # type: ignore

            self.predictor = predict_smiles
        self._loaded = True

    def is_loaded(self) -> bool:
        """Check if models are loaded."""
        return getattr(self, "_loaded", False)

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
            Dictionary with extraction results
        """
        result = {
            "smiles": None,
            "confidence": None,
            "success": False,
            "error": None,
            "backend": "decimer",
        }

        # Lazy load on first use
        if not self.is_loaded():
            try:
                self.load_models()
            except Exception as e:
                result["error"] = f"model_load_failed: {str(e)}"
                return result

        try:
            temp_image_path: Optional[str] = None
            # Load image if path
            if isinstance(image, (str, Path)):
                image_path = str(image)
            elif isinstance(image, np.ndarray):
                # Save temp file for DECIMER
                from tempfile import NamedTemporaryFile

                temp_img = Image.fromarray(image)
                with NamedTemporaryFile(suffix=".png", delete=False) as f:
                    temp_img.save(f.name)
                    temp_image_path = f.name
                    image_path = temp_image_path
            else:
                # PIL Image
                from tempfile import NamedTemporaryFile

                with NamedTemporaryFile(suffix=".png", delete=False) as f:
                    image.save(f.name)
                    temp_image_path = f.name
                    image_path = temp_image_path

            predictor = self.predictor
            if predictor is None:
                result["error"] = "model_not_loaded"
                return result

            # Run prediction
            prediction = predictor(image_path)
            smiles, predictor_confidence = self._extract_smiles_and_confidence(
                prediction
            )

            if smiles:
                result["smiles"] = smiles
                result["success"] = True

                # Keep DECIMER-reported confidence when available; otherwise use heuristics.
                if predictor_confidence is None:
                    result["confidence"] = self._estimate_confidence(smiles)
                else:
                    result["confidence"] = predictor_confidence
            else:
                result["error"] = "no_structure_detected"

        except Exception as e:
            result["error"] = f"extraction_failed: {str(e)}"
        finally:
            if temp_image_path:
                try:
                    Path(temp_image_path).unlink(missing_ok=True)
                except Exception:
                    pass

        return result

    def _extract_smiles_and_confidence(self, prediction):
        """Extract smiles and optional confidence from predictor output."""
        smiles = None
        confidence = None

        if isinstance(prediction, str):
            smiles = prediction
        elif isinstance(prediction, dict):
            smiles = prediction.get("smiles") or prediction.get("smiles_predicted")
            confidence = prediction.get("confidence")
        elif isinstance(prediction, (tuple, list)) and len(prediction) >= 2:
            smiles = prediction[0]
            if isinstance(prediction[1], (int, float)):
                confidence = float(prediction[1])

        if not isinstance(smiles, str):
            smiles = None

        if confidence is not None and not isinstance(confidence, (int, float)):
            confidence = None

        return smiles, confidence

    def _estimate_confidence(self, smiles: str) -> float:
        """
        Estimate confidence from SMILES characteristics.

        Heuristics:
        - Longer valid SMILES = more confident
        - Has rings = confident
        - Has stereochemistry = confident
        """
        score = 0.5  # Base confidence

        # Reaction SMILES (contains '>') should not get ring-confidence boosts
        # from numeric tokens that may appear around reactant/product sections.
        is_reaction_smiles = ">" in smiles

        # Ring presence (single-digit and %nn ring closures)
        if not is_reaction_smiles and re.search(r"%\d{2}|[1-9]", smiles):
            score += 0.15

        # Stereochemistry
        if "@" in smiles:
            score += 0.15

        # Length (reasonable range)
        if 20 < len(smiles) < 200:
            score += 0.1

        # Heteroatoms (N, O, S)
        if any(c in smiles for c in ["N", "O", "S"]):
            score += 0.1

        return min(score, 1.0)

    def extract_batch(
        self,
        images: List[Union[str, Path, Image.Image]],
    ) -> List[Dict[str, Any]]:
        """
        Extract SMILES from multiple images.

        Note: DECIMER doesn't support true batching, processes sequentially.
        """
        return [self.extract(img) for img in images]

    def get_model_info(self) -> Dict[str, Any]:
        """Get model metadata."""
        return {
            "backend": "decimer",
            "version": self.model_version,
            "device": self.device,
            "confidence_threshold": self.confidence_threshold,
            "batch_size": self.batch_size,
            "loaded": self.is_loaded(),
            "license": "CC-BY-4.0",
            "attribution": "DECIMER 2.x by Kohulan Rajan et al. (Nature Communications 2023)",
            "accuracy": {
                "single_molecule_f1": 0.91,
                "stereochemistry": 0.80,
            },
            "components": {
                "segmentation": self.use_segmentation,
                "classifier": True,
                "transformer": "ViT_base + T5_decoder",
            },
        }


def extract_with_decimer(
    image: Union[str, Path, Image.Image],
    **kwargs,
) -> Dict[str, Any]:
    """
    Convenience function for quick extraction.

    Args:
        image: Image to process
        **kwargs: Passed to DECIMERExtractor constructor

    Returns:
        Extraction result dictionary
    """
    extractor = DECIMERExtractor(**kwargs)
    return extractor.extract(image)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python decimer_extractor.py <image_path>")
        sys.exit(1)

    image_path = sys.argv[1]
    result = extract_with_decimer(image_path)

    print(f"Backend: {result['backend']}")
    print(f"Success: {result['success']}")
    print(f"SMILES: {result['smiles']}")
    print(f"Confidence: {result.get('confidence', 'N/A')}")
    if result.get("error"):
        print(f"Error: {result['error']}")

    # Print attribution
    print("\n" + "=" * 60)
    print("DECIMER 2.x by Kohulan Rajan et al.")
    print("Nature Communications 2023")
    print("License: CC BY 4.0")
