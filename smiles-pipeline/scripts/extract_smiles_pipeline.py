#!/usr/bin/env python3
"""
SMILES extraction pipeline - New architecture.

Orchestrates OCSR extraction, 3-layer validation, and enrichment
using the modular validator/extractor/enricher components.

Usage:
    python extract_smiles_pipeline.py \\
        --input-dir data/extractions \\
        --output-dir smiles-pipeline/data/raw \\
        --backend-order molscribe,decimer \\
        --gpu

References:
- ARCHITECTURE_v1.md: Complete system design
- config/backends.yaml: OCSR backend configuration
- config/validation_rules.yaml: Validation thresholds
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from tqdm import tqdm

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from validators.indigo_validator import IndigoValidator
from validators.chemical_validator import ChemicalValidator
from validators.domain_validator import DomainValidator
from extractors.molscribe_extractor import MolScribeExtractor
from extractors.decimer_extractor import DECIMERExtractor
from enrichers.fingerprint_generator import FingerprintGenerator
from enrichers.property_calculator import PropertyCalculator


class ExtractionPipeline:
    """
    Complete SMILES extraction and validation pipeline.

    Orchestrates:
    1. OCSR extraction (MolScribe → DECIMER → MolVec)
    2. 3-layer validation (syntax → chemical → domain)
    3. Enrichment (properties + fingerprints)
    """

    def __init__(
        self,
        backend_order: List[str] = None,
        use_gpu: bool = True,
        confidence_threshold: float = 0.5,
        min_validation_score: int = 70,
    ):
        """
        Initialize pipeline.

        Args:
            backend_order: OCSR backend priority list
            use_gpu: Enable GPU acceleration
            confidence_threshold: Minimum OCSR confidence
            min_validation_score: Minimum validation score to accept
        """
        self.backend_order = backend_order or ["molscribe", "decimer"]
        self.use_gpu = use_gpu
        self.confidence_threshold = confidence_threshold
        self.min_validation_score = min_validation_score

        # Initialize components (lazy loaded)
        self._extractors = {}
        self._indigo_validator = None
        self._chemical_validator = None
        self._domain_validator = None
        self._property_calculator = None
        self._fingerprint_generator = None

        # Statistics
        self.stats = {
            "images_processed": 0,
            "extractions_successful": 0,
            "syntax_valid": 0,
            "chemically_valid": 0,
            "domain_valid": 0,
            "high_confidence": 0,
            "errors": [],
        }

    def _get_extractor(self, backend: str):
        """Get or create extractor for backend."""
        if backend not in self._extractors:
            device = "cuda" if self.use_gpu else "cpu"

            if backend == "molscribe":
                self._extractors[backend] = MolScribeExtractor(
                    device=device,
                    confidence_threshold=self.confidence_threshold,
                )
            elif backend == "decimer":
                self._extractors[backend] = DECIMERExtractor(
                    device=device,
                    confidence_threshold=self.confidence_threshold,
                )
            else:
                raise ValueError(f"Unknown backend: {backend}")

        return self._extractors[backend]

    @property
    def indigo(self):
        """Lazy load Indigo validator."""
        if self._indigo_validator is None:
            self._indigo_validator = IndigoValidator()
        return self._indigo_validator

    @property
    def chemical(self):
        """Lazy load chemical validator."""
        if self._chemical_validator is None:
            self._chemical_validator = ChemicalValidator()
        return self._chemical_validator

    @property
    def domain(self):
        """Lazy load domain validator."""
        if self._domain_validator is None:
            self._domain_validator = DomainValidator()
        return self._domain_validator

    @property
    def property_calc(self):
        """Lazy load property calculator."""
        if self._property_calculator is None:
            self._property_calculator = PropertyCalculator()
        return self._property_calculator

    @property
    def fingerprint_gen(self):
        """Lazy load fingerprint generator."""
        if self._fingerprint_generator is None:
            self._fingerprint_generator = FingerprintGenerator()
        return self._fingerprint_generator

    def extract_from_image(
        self,
        image_path: Path,
    ) -> Dict[str, Any]:
        """
        Extract SMILES from single image using backend chain.

        Args:
            image_path: Path to chemical structure image

        Returns:
            Extraction result with SMILES and metadata
        """
        result = {
            "smiles": None,
            "canonical_smiles": None,
            "backend_used": None,
            "confidence": 0.0,
            "success": False,
            "error": None,
            "image_path": str(image_path),
        }

        # Try backends in order
        for backend in self.backend_order:
            try:
                extractor = self._get_extractor(backend)
                extraction = extractor.extract(image_path)

                if extraction["success"] and extraction["smiles"]:
                    result["smiles"] = extraction["smiles"]
                    result["backend_used"] = backend
                    result["confidence"] = extraction.get("confidence", 0.0)
                    result["success"] = True
                    break

            except Exception as e:
                result["errors"] = result.get("errors", [])
                result["errors"].append(
                    {
                        "backend": backend,
                        "error": str(e),
                    }
                )

        if not result["success"]:
            result["error"] = "all_backends_failed"

        return result

    def validate_molecule(
        self,
        smiles: str,
    ) -> Dict[str, Any]:
        """
        Apply 3-layer validation to SMILES.

        Args:
            smiles: Raw SMILES string

        Returns:
            Validation results for all 3 layers
        """
        validation = {
            "level_1_syntax": None,
            "level_2_chemical": None,
            "level_3_domain": None,
            "is_valid": False,
            "confidence_score": 0,
        }

        # Level 1: Syntax (Indigo)
        syntax_result = self.indigo.validate_syntax(smiles)
        validation["level_1_syntax"] = syntax_result

        if not syntax_result["is_valid"]:
            return validation

        # Use canonical SMILES for subsequent validation
        canonical = syntax_result.get("canonical_smiles", smiles)

        # Level 2: Chemical Plausibility
        chemical_result = self.chemical.validate(smiles, canonical)
        validation["level_2_chemical"] = chemical_result

        if not chemical_result["is_valid"]:
            return validation

        # Level 3: Domain (Ethnopharmacology)
        domain_result = self.domain.validate(smiles, canonical)
        validation["level_3_domain"] = domain_result

        # Compute confidence score
        score = self._compute_confidence_score(
            syntax_result,
            chemical_result,
            domain_result,
        )
        validation["confidence_score"] = score
        validation["is_valid"] = score >= self.min_validation_score

        return validation

    def _compute_confidence_score(
        self,
        syntax: Dict[str, Any],
        chemical: Dict[str, Any],
        domain: Dict[str, Any],
    ) -> int:
        """
        Compute overall confidence score (0-100).

        Based on validation results and gold standard matches.
        """
        score = 0

        # Base scores by validation level
        if syntax["is_valid"]:
            score += 30
        if chemical["is_valid"]:
            score += 30
        if domain["is_valid"]:
            score += 20

        # Bonuses
        if domain.get("gold_standard_matches"):
            for match in domain["gold_standard_matches"]:
                if match["similarity"] >= 0.99:
                    score += 20  # Exact match
                    break
                elif match["similarity"] >= 0.7:
                    score += 10  # Similar
                    break

        if domain.get("matches_sceletium_scaffold"):
            score += 10

        if domain.get("methoxy_count", 0) >= 2:
            score += 5

        # Apply penalties
        if syntax.get("wildcard_count", 0) > 0:
            score -= 5

        return max(0, min(100, score))

    def enrich_molecule(
        self,
        smiles: str,
        canonical_smiles: str = None,
    ) -> Dict[str, Any]:
        """
        Compute properties and fingerprints.

        Args:
            smiles: SMILES string
            canonical_smiles: Pre-computed canonical SMILES

        Returns:
            Enrichment data (properties + fingerprints)
        """
        enrichment = {
            "properties": None,
            "fingerprints": None,
            "lipinski": None,
        }

        # Compute properties
        prop_result = self.property_calc.calculate(smiles, canonical_smiles)
        if prop_result["success"]:
            enrichment["properties"] = prop_result["properties"]
            enrichment["lipinski"] = prop_result["lipinski"]

        # Generate fingerprints
        fp_result = self.fingerprint_gen.generate(smiles, canonical_smiles)
        if fp_result["success"]:
            enrichment["fingerprints"] = fp_result["fingerprints"]

        return enrichment

    def process_image(
        self,
        image_path: Path,
        include_enrichment: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Process single image through entire pipeline.

        Args:
            image_path: Path to chemical structure image
            include_enrichment: Compute properties/fingerprints

        Returns:
            Complete molecule record or None if all validation failed
        """
        self.stats["images_processed"] += 1

        # Step 1: Extract SMILES
        extraction = self.extract_from_image(image_path)

        if not extraction["success"]:
            self.stats["errors"].append(
                {
                    "image": str(image_path),
                    "error": extraction["error"],
                }
            )
            return None

        self.stats["extractions_successful"] += 1

        # Step 2: Validate
        validation = self.validate_molecule(extraction["smiles"])

        if validation["level_1_syntax"]["is_valid"]:
            self.stats["syntax_valid"] += 1
        if validation.get("level_2_chemical", {}).get("is_valid"):
            self.stats["chemically_valid"] += 1
        if validation.get("level_3_domain", {}).get("is_valid"):
            self.stats["domain_valid"] += 1
        if validation["confidence_score"] >= self.min_validation_score:
            self.stats["high_confidence"] += 1

        if not validation["is_valid"]:
            return None

        # Step 3: Enrich (optional)
        canonical = validation["level_1_syntax"].get("canonical_smiles")
        enrichment = {}

        if include_enrichment:
            enrichment = self.enrich_molecule(extraction["smiles"], canonical)

        # Build complete record
        record = {
            "id": f"mol_{self.stats['images_processed']:06d}",
            "smiles": extraction["smiles"],
            "canonical_smiles": canonical,
            "backend_used": extraction["backend_used"],
            "extraction_confidence": extraction["confidence"],
            "validation": validation,
            "confidence_score": validation["confidence_score"],
            "compound_class": validation["level_3_domain"].get(
                "compound_class", "unknown"
            ),
            "source": {
                "image_path": str(image_path),
                "extracted_at": datetime.now().isoformat(),
            },
        }

        # Add enrichment
        if enrichment.get("properties"):
            record["properties"] = enrichment["properties"]
        if enrichment.get("fingerprints"):
            record["fingerprints"] = enrichment["fingerprints"]
        if enrichment.get("lipinski"):
            record["lipinski"] = enrichment["lipinski"]

        # Add gold standard matches
        if validation["level_3_domain"].get("gold_standard_matches"):
            record["gold_standard_matches"] = validation["level_3_domain"][
                "gold_standard_matches"
            ]

        return record

    def process_directory(
        self,
        input_dir: Path,
        output_dir: Path,
        limit: int = None,
    ) -> List[Dict[str, Any]]:
        """
        Process all images in directory.

        Args:
            input_dir: Directory with chemical structure images
            output_dir: Output directory for results
            limit: Maximum number of images to process

        Returns:
            List of validated molecule records
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Find all images
        image_extensions = {".jpg", ".jpeg", ".png", ".gif", ".bmp"}
        image_files = [
            f for f in input_dir.rglob("*") if f.suffix.lower() in image_extensions
        ]

        if limit:
            image_files = image_files[:limit]

        print(f"Found {len(image_files)} images to process")

        # Process images
        molecules = []
        for image_path in tqdm(image_files, desc="Processing images"):
            try:
                record = self.process_image(image_path)
                if record:
                    molecules.append(record)
            except Exception as e:
                self.stats["errors"].append(
                    {
                        "image": str(image_path),
                        "error": str(e),
                    }
                )

        # Save results
        output_file = output_dir / "molecules.jsonl"
        with open(output_file, "w") as f:
            for mol in molecules:
                f.write(json.dumps(mol) + "\n")

        # Save statistics
        stats_file = output_dir / "extraction_stats.json"
        with open(stats_file, "w") as f:
            json.dump(self.stats, f, indent=2)

        print(f"\nProcessed {len(molecules)} validated molecules")
        print(f"Results saved to: {output_file}")

        return molecules


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="SMILES extraction pipeline with 3-layer validation"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory with chemical structure images",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for results",
    )
    parser.add_argument(
        "--backend-order",
        type=str,
        default="molscribe,decimer",
        help="Comma-separated OCSR backend order",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Enable GPU acceleration",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Force CPU mode",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of images to process",
    )
    parser.add_argument(
        "--confidence-threshold",
        type=float,
        default=0.5,
        help="Minimum OCSR confidence",
    )
    parser.add_argument(
        "--min-validation-score",
        type=int,
        default=70,
        help="Minimum validation score (0-100)",
    )

    args = parser.parse_args()

    # Parse backend order
    backend_order = [b.strip() for b in args.backend_order.split(",")]

    # Determine GPU mode
    use_gpu = args.gpu or (not args.no_gpu)

    # Create pipeline
    pipeline = ExtractionPipeline(
        backend_order=backend_order,
        use_gpu=use_gpu,
        confidence_threshold=args.confidence_threshold,
        min_validation_score=args.min_validation_score,
    )

    # Run pipeline
    molecules = pipeline.process_directory(
        args.input_dir,
        args.output_dir,
        limit=args.limit,
    )

    # Print summary
    print("\n" + "=" * 60)
    print("Pipeline Summary")
    print("=" * 60)
    print(f"Images processed: {pipeline.stats['images_processed']}")
    print(f"Extractions successful: {pipeline.stats['extractions_successful']}")
    print(f"Syntax valid: {pipeline.stats['syntax_valid']}")
    print(f"Chemically valid: {pipeline.stats['chemically_valid']}")
    print(f"Domain valid: {pipeline.stats['domain_valid']}")
    print(
        f"High confidence (≥{args.min_validation_score}): {pipeline.stats['high_confidence']}"
    )
    print(f"Errors: {len(pipeline.stats['errors'])}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
