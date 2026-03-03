#!/usr/bin/env python3
"""
SMILES extraction pipeline - New architecture.

Orchestrates OCSR extraction, deterministic standardization,
3-layer validation, and enrichment
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
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from uuid import uuid4
from tqdm import tqdm
import yaml

try:
    from PIL import Image, ImageFilter, UnidentifiedImageError
except ImportError:  # pragma: no cover - optional dependency fallback
    Image = None
    ImageFilter = None
    UnidentifiedImageError = Exception

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

IndigoValidator = None
INDIGO_IMPORT_ERROR = None
StandardizationValidator = None
STANDARDIZATION_IMPORT_ERROR = None
ChemicalValidator = None
CHEMICAL_IMPORT_ERROR = None
DomainValidator = None
DOMAIN_IMPORT_ERROR = None
MolScribeExtractor = None
MOLSCRIBE_IMPORT_ERROR = None
DECIMERExtractor = None
DECIMER_IMPORT_ERROR = None

FingerprintGenerator = None
FINGERPRINT_IMPORT_ERROR = None
PropertyCalculator = None
PROPERTY_CALC_IMPORT_ERROR = None

SUPPORTED_BACKENDS = {"molscribe", "decimer"}


def load_symbol(module_names: List[str], symbol_name: str):
    """Load symbol from the first available module path."""
    last_error = None
    for module_name in module_names:
        try:
            module = import_module(module_name)
            return getattr(module, symbol_name), None
        except Exception as exc:
            last_error = exc
    return None, last_error


def load_yaml_config(path: Path) -> Tuple[Dict[str, Any], Optional[str]]:
    """Load a YAML config file and return (config, error_message)."""
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
            if not isinstance(data, dict):
                return {}, f"Config must be a YAML mapping: {path}"
            return data, None
    except FileNotFoundError:
        return {}, f"Missing config file: {path}"
    except yaml.YAMLError as exc:
        return {}, f"Invalid YAML in {path}: {exc}"


class ExtractionPipeline:
    """
    Complete SMILES extraction and validation pipeline.

    Orchestrates:
    1. OCSR extraction (MolScribe → DECIMER)
    2. Deterministic standardization (RDKit MolStandardize)
    3. 3-layer validation (syntax → chemical → domain)
    4. Enrichment (properties + fingerprints)
    """

    def __init__(
        self,
        backend_order: Optional[List[str]] = None,
        use_gpu: bool = True,
        confidence_threshold: Optional[float] = None,
        min_validation_score: Optional[int] = None,
        gold_standards_file: Optional[str] = None,
        backends_config: Optional[Dict[str, Any]] = None,
        validation_rules: Optional[Dict[str, Any]] = None,
        config_sources: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize pipeline.

        Args:
            backend_order: OCSR backend priority list (defaults to config priorities)
            use_gpu: Enable GPU acceleration
            confidence_threshold: Minimum OCSR confidence (defaults to config)
            min_validation_score: Minimum validation score to accept (defaults to config)
            gold_standards_file: Optional path override for gold standards file
            backends_config: Parsed contents of config/backends.yaml
            validation_rules: Parsed contents of config/validation_rules.yaml
            config_sources: Optional source metadata for loaded config paths
        """
        self.backends_config = backends_config or {}
        self.validation_rules = validation_rules or {}
        self.config_sources = config_sources or {}

        requested_backends = backend_order or self._resolve_backend_order()
        self.unsupported_backends = [
            b for b in requested_backends if b not in SUPPORTED_BACKENDS
        ]
        self.backend_order = [b for b in requested_backends if b in SUPPORTED_BACKENDS]
        if not self.backend_order:
            self.backend_order = ["molscribe", "decimer"]
        self.use_gpu = use_gpu
        self._confidence_threshold_overridden = confidence_threshold is not None
        self.confidence_threshold = (
            confidence_threshold
            if confidence_threshold is not None
            else self._resolve_confidence_threshold()
        )
        self.min_validation_score = (
            min_validation_score
            if min_validation_score is not None
            else self._resolve_min_validation_score()
        )
        self.min_extraction_confidence = self._resolve_min_extraction_confidence()
        (
            self.low_confidence_review_enabled,
            self.low_confidence_review_threshold,
        ) = self._resolve_low_confidence_review_policy()
        self.image_prefilter_config = self._resolve_image_prefilter_config()
        self.image_prefilter_enabled = bool(self.image_prefilter_config.get("enabled"))
        self.gold_standards_file = gold_standards_file
        self.run_id = self._generate_run_id()

        # Initialize components (lazy loaded)
        self._extractors = {}
        self._indigo_validator = None
        self._chemical_validator = None
        self._domain_validator = None
        self._standardizer = None
        self._property_calculator = None
        self._fingerprint_generator = None

        # Statistics
        self.stats = {
            "run_id": self.run_id,
            "started_at_utc": datetime.utcnow().isoformat() + "Z",
            "images_processed": 0,
            "extractions_successful": 0,
            "syntax_valid": 0,
            "chemically_valid": 0,
            "domain_valid": 0,
            "high_confidence": 0,
            "low_extraction_confidence_rejected": 0,
            "manual_review_candidates": 0,
            "images_discovered": 0,
            "images_prefilter_rejected": 0,
            "images_prefilter_soft_review": 0,
            "rejections_by_stage": {},
            "rejections_by_reason_code": {},
            "errors": [],
            "effective_config": self._build_effective_config_snapshot(),
            "runtime_manifest": self._build_runtime_manifest(),
        }

    def _resolve_backend_order(self) -> List[str]:
        """Resolve backend priority from config with unsupported backends filtered."""
        candidates = []
        for backend, settings in self.backends_config.items():
            if backend in {"routing", "fallback", "monitoring", "vlm"}:
                continue
            if not isinstance(settings, dict):
                continue
            if not settings.get("enabled", False):
                continue
            candidates.append((backend, settings.get("priority", 999)))

        if not candidates:
            return ["molscribe", "decimer"]

        candidates.sort(key=lambda x: x[1])
        ordered = [name for name, _ in candidates]
        supported = [name for name in ordered if name in SUPPORTED_BACKENDS]
        return supported or ["molscribe", "decimer"]

    def _resolve_confidence_threshold(self) -> float:
        """Resolve OCSR confidence threshold from backend/fallback config."""
        thresholds = []
        for backend in self.backend_order:
            backend_cfg = self.backends_config.get(backend, {})
            value = backend_cfg.get("confidence_threshold")
            if isinstance(value, (int, float)):
                thresholds.append(float(value))

        if thresholds:
            return min(thresholds)
        fallback = self.backends_config.get("fallback", {})
        if isinstance(fallback.get("min_confidence"), (int, float)):
            return float(fallback["min_confidence"])
        return 0.5

    def _resolve_min_validation_score(self) -> int:
        """Resolve high-confidence score threshold from validation rules."""
        scoring = self.validation_rules.get("confidence_scoring", {})
        thresholds = scoring.get("thresholds", {})
        high = thresholds.get("high_confidence")
        if isinstance(high, (int, float)):
            return int(high)
        return 70

    def _resolve_min_extraction_confidence(self) -> float:
        """Resolve minimum extraction confidence gate from backend fallback policy."""
        fallback = self.backends_config.get("fallback", {})
        min_conf = fallback.get("min_confidence")
        if isinstance(min_conf, (int, float)):
            return float(min_conf)
        return 0.4

    def _resolve_low_confidence_review_policy(self) -> tuple[bool, float]:
        """Resolve manual review policy for low-confidence extractions."""
        low_conf_review = self.backends_config.get("fallback", {}).get(
            "low_confidence_review", {}
        )
        enabled = bool(low_conf_review.get("enabled", True))
        threshold = low_conf_review.get("threshold")
        if isinstance(threshold, (int, float)):
            return enabled, float(threshold)
        return enabled, 0.6

    def _resolve_image_prefilter_config(self) -> Dict[str, Any]:
        """Resolve lightweight image prefilter policy from validation rules."""
        prefilter = self.validation_rules.get("image_prefilter", {})
        if not isinstance(prefilter, dict):
            prefilter = {}
        return {
            "enabled": bool(prefilter.get("enabled", False)),
            "min_width": int(prefilter.get("min_width", 200)),
            "min_height": int(prefilter.get("min_height", 150)),
            "max_aspect_ratio": float(prefilter.get("max_aspect_ratio", 5.0)),
            "dark_threshold": int(prefilter.get("dark_threshold", 48)),
            "min_dark_pixel_ratio": float(prefilter.get("min_dark_pixel_ratio", 0.003)),
            "max_dark_pixel_ratio": float(prefilter.get("max_dark_pixel_ratio", 0.35)),
            "edge_threshold": int(prefilter.get("edge_threshold", 32)),
            "min_edge_pixel_ratio": float(prefilter.get("min_edge_pixel_ratio", 0.005)),
            "max_edge_pixel_ratio": float(prefilter.get("max_edge_pixel_ratio", 0.45)),
            "min_entropy": float(prefilter.get("min_entropy", 0.2)),
            "max_entropy": float(prefilter.get("max_entropy", 7.8)),
        }

    def _build_indigo_config(self) -> Dict[str, Any]:
        """Map level_1_syntax config to IndigoValidator settings."""
        level_1 = self.validation_rules.get("level_1_syntax", {})
        runtime_options = level_1.get("runtime_options", {})
        set_options = level_1.get("set_options", {})
        return {
            "wildcard_threshold": level_1.get("wildcard_threshold", 2),
            "placeholder_values": level_1.get("placeholder_values", []),
            "runtime_options": {
                "timeout_ms": runtime_options.get("timeout_ms", 0),
                "reset_options_each_call": runtime_options.get(
                    "reset_options_each_call", False
                ),
            },
            "set_options": set_options if isinstance(set_options, dict) else {},
        }

    def _build_chemical_config(self) -> Dict[str, Any]:
        """Map level_2_chemical config to ChemicalValidator settings."""
        level_2 = self.validation_rules.get("level_2_chemical", {})
        mw = level_2.get("molecular_weight", {})
        atom = level_2.get("atom_count", {})
        fragments = level_2.get("fragments", {})
        rings = level_2.get("ring_systems", {})
        charge = level_2.get("charge", {})
        return {
            "molecular_weight": {
                "min": mw.get("min", 100),
                "max": mw.get("max", 1000),
                "warn_min": mw.get("warn_min", 150),
                "warn_max": mw.get("warn_max", 800),
            },
            "atom_count": {
                "min": atom.get("min", 5),
                "max": atom.get("max", 200),
            },
            "fragments": {
                "max": fragments.get("max", 5),
                "warn_max": fragments.get("warn_max", 3),
            },
            "ring_systems": {
                "max": rings.get("max_rings", 10),
            },
            "charge": {
                "max_absolute": charge.get("max_absolute_charge", 3),
            },
        }

    def _build_domain_config(self) -> Dict[str, Any]:
        """Expose level_3_domain policy to DomainValidator."""
        return self.validation_rules.get("level_3_domain", {})

    def _build_effective_config_snapshot(self) -> Dict[str, Any]:
        """Capture effective runtime settings for reproducibility."""
        return {
            "sources": self.config_sources,
            "backend_order": self.backend_order,
            "unsupported_backends": self.unsupported_backends,
            "confidence_threshold": self.confidence_threshold,
            "min_extraction_confidence": self.min_extraction_confidence,
            "min_validation_score": self.min_validation_score,
            "low_confidence_review_enabled": self.low_confidence_review_enabled,
            "low_confidence_review_threshold": self.low_confidence_review_threshold,
            "image_prefilter": self.image_prefilter_config,
            "gold_standards_file": self.gold_standards_file,
        }

    def _generate_run_id(self) -> str:
        """Generate unique run id for provenance across records and stats."""
        return f"smiles_run_{datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')}_{uuid4().hex[:8]}"

    def _safe_package_version(self, package_name: str) -> Optional[str]:
        """Return installed package version when available."""
        try:
            return version(package_name)
        except PackageNotFoundError:
            return None
        except Exception:
            return None

    def _first_available_version(self, candidates: List[str]) -> Optional[str]:
        """Return first discovered package version from candidate names."""
        for package_name in candidates:
            found = self._safe_package_version(package_name)
            if found:
                return found
        return None

    def _sha256_file(self, path_str: Optional[str]) -> Optional[str]:
        """Return SHA256 for a file path (None when missing/unreadable)."""
        if not path_str:
            return None
        file_path = Path(path_str)
        if not file_path.exists():
            return None
        hasher = hashlib.sha256()
        try:
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(8192), b""):
                    hasher.update(chunk)
            return hasher.hexdigest()
        except Exception:
            return None

    def _build_runtime_manifest(self) -> Dict[str, Any]:
        """Capture runtime/toolkit manifest for reproducibility audits."""
        dep_versions = {
            "python": platform.python_version(),
            "numpy": self._safe_package_version("numpy"),
            "rdkit": self._first_available_version(["rdkit", "rdkit-pypi"]),
            "epam-indigo": self._first_available_version(
                ["epam.indigo", "epam-indigo"]
            ),
            "molscribe": self._safe_package_version("molscribe"),
            "decimer": self._safe_package_version("decimer"),
            "pillow": self._safe_package_version("Pillow"),
            "pyyaml": self._safe_package_version("PyYAML"),
            "tqdm": self._safe_package_version("tqdm"),
        }
        config_hashes = {
            key: self._sha256_file(path)
            for key, path in (self.config_sources or {}).items()
        }
        return {
            "generated_at_utc": datetime.utcnow().isoformat() + "Z",
            "run_id": self.run_id,
            "platform": {
                "system": platform.system(),
                "release": platform.release(),
                "machine": platform.machine(),
            },
            "dependencies": dep_versions,
            "config_hashes": config_hashes,
            "indigo_runtime_policy": self._build_indigo_config().get("runtime_options"),
        }

    def _extractor_runtime_info(self, backend: str, extractor: Any) -> Dict[str, Any]:
        """Return backend runtime metadata without forcing extra model loads."""
        backend_cfg = self.backends_config.get(backend, {})
        info = {
            "backend": backend,
            "configured_version": backend_cfg.get("version"),
            "configured_model": backend_cfg.get("model"),
            "configured_confidence_threshold": backend_cfg.get("confidence_threshold"),
        }
        if extractor is not None and hasattr(extractor, "get_model_info"):
            try:
                model_info = extractor.get_model_info() or {}
                if isinstance(model_info, dict):
                    info["runtime_model_info"] = model_info
            except Exception:
                pass
        return info

    def _build_rejection_event(
        self,
        stage: str,
        reason_code: str,
        reason_detail: str,
        image_path: Path,
        **extra: Any,
    ) -> Dict[str, Any]:
        """Build normalized rejection payload for downstream QA/reporting."""
        event = {
            "image": str(image_path),
            "stage": stage,
            "reason_code": reason_code,
            "reason_detail": reason_detail,
            # Backward-compatible alias expected by older tooling/docs.
            "error": reason_code,
        }
        for key, value in extra.items():
            if value is not None:
                event[key] = value
        return event

    def _record_rejection(self, event: Dict[str, Any]) -> None:
        """Persist rejection event and keep aggregate counters in sync."""
        self.stats["errors"].append(event)
        stage = event.get("stage", "unknown")
        reason_code = event.get("reason_code", "unknown")
        self.stats["rejections_by_stage"][stage] = (
            self.stats["rejections_by_stage"].get(stage, 0) + 1
        )
        self.stats["rejections_by_reason_code"][reason_code] = (
            self.stats["rejections_by_reason_code"].get(reason_code, 0) + 1
        )

    def _derive_validation_rejection(
        self, image_path: Path, validation: Dict[str, Any], extraction: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Map failed validation output to a normalized rejection event."""
        syntax = validation.get("level_1_syntax") or {}
        if syntax and not syntax.get("is_valid", False):
            return self._build_rejection_event(
                stage="level_1_syntax",
                reason_code=str(syntax.get("error") or "syntax_validation_failed"),
                reason_detail="Failed Indigo syntax gate",
                image_path=image_path,
                backend_used=extraction.get("backend_used"),
                extraction_confidence=extraction.get("confidence"),
            )

        standardization = validation.get("standardization") or {}
        if standardization and not standardization.get("success", False):
            return self._build_rejection_event(
                stage="standardization",
                reason_code=str(
                    standardization.get("error") or "standardization_failed"
                ),
                reason_detail="Failed deterministic RDKit standardization",
                image_path=image_path,
                backend_used=extraction.get("backend_used"),
                extraction_confidence=extraction.get("confidence"),
            )

        chemical = validation.get("level_2_chemical") or {}
        if chemical and not chemical.get("is_valid", False):
            return self._build_rejection_event(
                stage="level_2_chemical",
                reason_code=str(
                    chemical.get("rejection_reason")
                    or chemical.get("error")
                    or "chemical_validation_failed"
                ),
                reason_detail="Failed RDKit chemical plausibility gate",
                image_path=image_path,
                backend_used=extraction.get("backend_used"),
                extraction_confidence=extraction.get("confidence"),
            )

        domain = validation.get("level_3_domain") or {}
        if domain and not domain.get("is_valid", False):
            return self._build_rejection_event(
                stage="level_3_domain",
                reason_code=str(
                    domain.get("rejection_reason")
                    or domain.get("error")
                    or "domain_validation_failed"
                ),
                reason_detail="Failed domain-specific gate",
                image_path=image_path,
                backend_used=extraction.get("backend_used"),
                extraction_confidence=extraction.get("confidence"),
            )

        return self._build_rejection_event(
            stage="confidence_routing",
            reason_code="below_min_validation_score",
            reason_detail=(
                f"Score {validation.get('confidence_score', 0)} below "
                f"min {self.min_validation_score}"
            ),
            image_path=image_path,
            backend_used=extraction.get("backend_used"),
            extraction_confidence=extraction.get("confidence"),
            confidence_score=validation.get("confidence_score"),
            min_validation_score=self.min_validation_score,
            triage_label=validation.get("triage_label"),
        )

    def _get_extractor(self, backend: str):
        """Get or create extractor for backend."""
        global MolScribeExtractor, MOLSCRIBE_IMPORT_ERROR
        global DECIMERExtractor, DECIMER_IMPORT_ERROR
        if backend not in self._extractors:
            backend_cfg = self.backends_config.get(backend, {})
            device = "cuda" if self.use_gpu else "cpu"
            batch_size = backend_cfg.get("batch_size")
            backend_conf = backend_cfg.get(
                "confidence_threshold", self.confidence_threshold
            )
            if self._confidence_threshold_overridden:
                backend_conf = self.confidence_threshold

            if backend == "molscribe":
                if MolScribeExtractor is None and MOLSCRIBE_IMPORT_ERROR is None:
                    MolScribeExtractor, MOLSCRIBE_IMPORT_ERROR = load_symbol(
                        [
                            "extractors.molscribe_extractor",
                            "src.extractors.molscribe_extractor",
                        ],
                        "MolScribeExtractor",
                    )
                if MolScribeExtractor is None:
                    raise RuntimeError(
                        "MolScribe extractor dependency is unavailable. "
                        f"Original import error: {MOLSCRIBE_IMPORT_ERROR}"
                    ) from MOLSCRIBE_IMPORT_ERROR
                self._extractors[backend] = MolScribeExtractor(
                    model_name=backend_cfg.get("model", "swin_base_char_aux_1m.pth"),
                    device=device,
                    confidence_threshold=backend_conf,
                    batch_size=batch_size if isinstance(batch_size, int) else 32,
                )
            elif backend == "decimer":
                if DECIMERExtractor is None and DECIMER_IMPORT_ERROR is None:
                    DECIMERExtractor, DECIMER_IMPORT_ERROR = load_symbol(
                        [
                            "extractors.decimer_extractor",
                            "src.extractors.decimer_extractor",
                        ],
                        "DECIMERExtractor",
                    )
                if DECIMERExtractor is None:
                    raise RuntimeError(
                        "DECIMER extractor dependency is unavailable. "
                        f"Original import error: {DECIMER_IMPORT_ERROR}"
                    ) from DECIMER_IMPORT_ERROR
                segmentation_cfg = backend_cfg.get("segmentation", {})
                self._extractors[backend] = DECIMERExtractor(
                    model_version=str(backend_cfg.get("version", "2.7")),
                    device=device,
                    confidence_threshold=backend_conf,
                    batch_size=batch_size if isinstance(batch_size, int) else 16,
                    use_segmentation=segmentation_cfg.get("enabled", True),
                )
            else:
                raise ValueError(f"Unknown backend: {backend}")

        return self._extractors[backend]

    def _evaluate_image_prefilter(self, image_path: Path) -> Dict[str, Any]:
        """
        Lightweight image gate to skip obvious non-structure figures.

        The gate is intentionally conservative: images are rejected only when
        they violate hard geometric/content bounds often seen in charts/photos.
        """
        if not self.image_prefilter_enabled:
            return {"is_candidate": True, "reason_code": "prefilter_disabled"}

        if Image is None or ImageFilter is None:
            return {"is_candidate": True, "reason_code": "prefilter_dependency_missing"}

        cfg = self.image_prefilter_config
        try:
            with Image.open(image_path) as img:
                gray = img.convert("L")
                width, height = gray.size
                if width <= 0 or height <= 0:
                    return {
                        "is_candidate": False,
                        "reason_code": "invalid_image_dimensions",
                        "reason_detail": f"invalid dimensions {width}x{height}",
                    }

                total_pixels = width * height
                aspect_ratio = max(width / height, height / width)
                hist = gray.histogram()
                dark_threshold = max(0, min(int(cfg["dark_threshold"]), 255))
                dark_pixels = sum(hist[: dark_threshold + 1])
                dark_ratio = dark_pixels / total_pixels
                entropy = float(gray.entropy())

                edge_img = gray.filter(ImageFilter.FIND_EDGES)
                edge_hist = edge_img.histogram()
                edge_threshold = max(0, min(int(cfg["edge_threshold"]), 255))
                edge_pixels = sum(edge_hist[edge_threshold:])
                edge_ratio = edge_pixels / total_pixels
        except Exception as exc:
            if isinstance(exc, (UnidentifiedImageError, OSError)):
                reason_code = "image_open_failed"
            else:
                reason_code = "prefilter_runtime_error"
            return {
                "is_candidate": False,
                "reason_code": reason_code,
                "reason_detail": str(exc),
            }

        metrics = {
            "width": width,
            "height": height,
            "aspect_ratio": round(aspect_ratio, 4),
            "dark_pixel_ratio": round(dark_ratio, 6),
            "edge_pixel_ratio": round(edge_ratio, 6),
            "entropy": round(entropy, 6),
        }

        if width < cfg["min_width"] or height < cfg["min_height"]:
            return {
                "is_candidate": False,
                "reason_code": "image_too_small",
                "reason_detail": (
                    f"size={width}x{height} below "
                    f"{cfg['min_width']}x{cfg['min_height']}"
                ),
                "metrics": metrics,
            }
        if aspect_ratio > cfg["max_aspect_ratio"]:
            return {
                "is_candidate": False,
                "reason_code": "aspect_ratio_out_of_bounds",
                "reason_detail": (
                    f"aspect_ratio={aspect_ratio:.3f} > {cfg['max_aspect_ratio']}"
                ),
                "metrics": metrics,
            }
        if dark_ratio < cfg["min_dark_pixel_ratio"]:
            return {
                "is_candidate": False,
                "reason_code": "too_sparse_for_structure",
                "reason_detail": (
                    f"dark_pixel_ratio={dark_ratio:.6f} < {cfg['min_dark_pixel_ratio']}"
                ),
                "metrics": metrics,
            }
        if dark_ratio > cfg["max_dark_pixel_ratio"]:
            return {
                "is_candidate": False,
                "reason_code": "too_dense_for_structure",
                "reason_detail": (
                    f"dark_pixel_ratio={dark_ratio:.6f} > {cfg['max_dark_pixel_ratio']}"
                ),
                "metrics": metrics,
            }
        if edge_ratio < cfg["min_edge_pixel_ratio"]:
            return {
                "is_candidate": False,
                "reason_code": "insufficient_edges",
                "reason_detail": (
                    f"edge_pixel_ratio={edge_ratio:.6f} < {cfg['min_edge_pixel_ratio']}"
                ),
                "metrics": metrics,
            }
        if edge_ratio > cfg["max_edge_pixel_ratio"]:
            return {
                "is_candidate": False,
                "reason_code": "edge_noise_too_high",
                "reason_detail": (
                    f"edge_pixel_ratio={edge_ratio:.6f} > {cfg['max_edge_pixel_ratio']}"
                ),
                "metrics": metrics,
            }
        if entropy < cfg["min_entropy"] or entropy > cfg["max_entropy"]:
            return {
                "is_candidate": False,
                "reason_code": "entropy_out_of_bounds",
                "reason_detail": (
                    f"entropy={entropy:.3f} outside "
                    f"[{cfg['min_entropy']}, {cfg['max_entropy']}]"
                ),
                "metrics": metrics,
            }

        return {
            "is_candidate": True,
            "reason_code": "prefilter_pass",
            "metrics": metrics,
        }

    @property
    def indigo(self):
        """Lazy load Indigo validator."""
        global IndigoValidator, INDIGO_IMPORT_ERROR
        if IndigoValidator is None and INDIGO_IMPORT_ERROR is None:
            IndigoValidator, INDIGO_IMPORT_ERROR = load_symbol(
                ["validators.indigo_validator", "src.validators.indigo_validator"],
                "IndigoValidator",
            )
        if IndigoValidator is None:
            missing_name = getattr(INDIGO_IMPORT_ERROR, "name", "indigo")
            raise RuntimeError(
                "Indigo validator dependency is unavailable "
                f"(missing module: {missing_name}). "
                "Install the Python package `epam.indigo` in this environment."
            ) from INDIGO_IMPORT_ERROR
        if self._indigo_validator is None:
            self._indigo_validator = IndigoValidator(config=self._build_indigo_config())
        return self._indigo_validator

    @property
    def standardizer(self):
        """Lazy load deterministic standardizer."""
        global StandardizationValidator, STANDARDIZATION_IMPORT_ERROR
        if StandardizationValidator is None and STANDARDIZATION_IMPORT_ERROR is None:
            StandardizationValidator, STANDARDIZATION_IMPORT_ERROR = load_symbol(
                [
                    "validators.standardization_validator",
                    "src.validators.standardization_validator",
                ],
                "StandardizationValidator",
            )
        if StandardizationValidator is None:
            raise RuntimeError(
                "Standardization validator dependency is unavailable. "
                f"Original import error: {STANDARDIZATION_IMPORT_ERROR}"
            ) from STANDARDIZATION_IMPORT_ERROR
        if self._standardizer is None:
            self._standardizer = StandardizationValidator()
        return self._standardizer

    @property
    def chemical(self):
        """Lazy load chemical validator."""
        global ChemicalValidator, CHEMICAL_IMPORT_ERROR
        if ChemicalValidator is None and CHEMICAL_IMPORT_ERROR is None:
            ChemicalValidator, CHEMICAL_IMPORT_ERROR = load_symbol(
                ["validators.chemical_validator", "src.validators.chemical_validator"],
                "ChemicalValidator",
            )
        if ChemicalValidator is None:
            raise RuntimeError(
                "Chemical validator dependency is unavailable. "
                f"Original import error: {CHEMICAL_IMPORT_ERROR}"
            ) from CHEMICAL_IMPORT_ERROR
        if self._chemical_validator is None:
            self._chemical_validator = ChemicalValidator(
                config=self._build_chemical_config()
            )
        return self._chemical_validator

    @property
    def domain(self):
        """Lazy load domain validator."""
        global DomainValidator, DOMAIN_IMPORT_ERROR
        if DomainValidator is None and DOMAIN_IMPORT_ERROR is None:
            DomainValidator, DOMAIN_IMPORT_ERROR = load_symbol(
                ["validators.domain_validator", "src.validators.domain_validator"],
                "DomainValidator",
            )
        if DomainValidator is None:
            raise RuntimeError(
                "Domain validator dependency is unavailable. "
                f"Original import error: {DOMAIN_IMPORT_ERROR}"
            ) from DOMAIN_IMPORT_ERROR
        if self._domain_validator is None:
            self._domain_validator = DomainValidator(
                gold_standards_file=self.gold_standards_file,
                config=self._build_domain_config(),
            )
        return self._domain_validator

    @property
    def property_calc(self):
        """Lazy load property calculator."""
        global PropertyCalculator, PROPERTY_CALC_IMPORT_ERROR
        if PropertyCalculator is None and PROPERTY_CALC_IMPORT_ERROR is None:
            PropertyCalculator, PROPERTY_CALC_IMPORT_ERROR = load_symbol(
                ["enrichers.property_calculator", "src.enrichers.property_calculator"],
                "PropertyCalculator",
            )
        if PropertyCalculator is None:
            raise RuntimeError(
                "Property calculator dependency is unavailable. "
                f"Original import error: {PROPERTY_CALC_IMPORT_ERROR}"
            ) from PROPERTY_CALC_IMPORT_ERROR
        if self._property_calculator is None:
            self._property_calculator = PropertyCalculator()
        return self._property_calculator

    @property
    def fingerprint_gen(self):
        """Lazy load fingerprint generator."""
        global FingerprintGenerator, FINGERPRINT_IMPORT_ERROR
        if FingerprintGenerator is None and FINGERPRINT_IMPORT_ERROR is None:
            FingerprintGenerator, FINGERPRINT_IMPORT_ERROR = load_symbol(
                [
                    "enrichers.fingerprint_generator",
                    "src.enrichers.fingerprint_generator",
                ],
                "FingerprintGenerator",
            )
        if FingerprintGenerator is None:
            raise RuntimeError(
                "Fingerprint generator dependency is unavailable. "
                f"Original import error: {FINGERPRINT_IMPORT_ERROR}"
            ) from FINGERPRINT_IMPORT_ERROR
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
            "backend_runtime": None,
            "cross_backend_agreement": None,
            "agreement_checked": False,
            "agreement_reference_backend": None,
            "agreement_reference_confidence": None,
            "agreement_reference_runtime": None,
            "success": False,
            "error": None,
            "image_path": str(image_path),
        }

        # Try backends in order. Prefer the first extraction that meets the
        # configured minimum confidence, but keep the best low-confidence
        # candidate as a fallback for downstream rejection analytics.
        best_candidate: Optional[Dict[str, Any]] = None
        for backend in self.backend_order:
            try:
                extractor = self._get_extractor(backend)
                extraction = extractor.extract(image_path)

                if extraction["success"] and extraction["smiles"]:
                    current = {
                        "smiles": extraction["smiles"],
                        "backend_used": backend,
                        "confidence": float(extraction.get("confidence", 0.0) or 0.0),
                        "backend_runtime": self._extractor_runtime_info(
                            backend, extractor
                        ),
                        "success": True,
                    }

                    if (
                        best_candidate is None
                        or current["confidence"] > best_candidate["confidence"]
                    ):
                        best_candidate = current

                    if current["confidence"] >= float(self.min_extraction_confidence):
                        result["smiles"] = current["smiles"]
                        result["backend_used"] = current["backend_used"]
                        result["confidence"] = current["confidence"]
                        result["backend_runtime"] = current["backend_runtime"]
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

        if not result["success"] and best_candidate is not None:
            result["smiles"] = best_candidate["smiles"]
            result["backend_used"] = best_candidate["backend_used"]
            result["confidence"] = best_candidate["confidence"]
            result["backend_runtime"] = best_candidate["backend_runtime"]
            result["success"] = True

        if not result["success"]:
            result["error"] = "all_backends_failed"
            return result

        # Optional consensus pass: compare primary prediction with one additional backend.
        consensus_cfg = self.backends_config.get("routing", {}).get(
            "consensus_check", {}
        )
        if consensus_cfg.get("enabled", False):
            check_below = float(consensus_cfg.get("check_when_confidence_below", 0.7))
            check_always = bool(consensus_cfg.get("check_always", False))
            if check_always or float(result["confidence"]) < check_below:
                primary_backend = result["backend_used"]
                primary_smiles = result["smiles"]
                secondary_smiles = None
                secondary_backend = None
                secondary_conf = None
                for backend in self.backend_order:
                    if backend == primary_backend:
                        continue
                    try:
                        extractor = self._get_extractor(backend)
                        extraction = extractor.extract(image_path)
                        if extraction["success"] and extraction["smiles"]:
                            secondary_backend = backend
                            secondary_smiles = extraction["smiles"]
                            secondary_conf = extraction.get("confidence", 0.0)
                            result["agreement_reference_runtime"] = (
                                self._extractor_runtime_info(backend, extractor)
                            )
                            break
                    except Exception:
                        continue

                if secondary_smiles:
                    primary_canonical = self.indigo.canonicalize(primary_smiles)
                    secondary_canonical = self.indigo.canonicalize(secondary_smiles)
                    agreed = False
                    if primary_canonical and secondary_canonical:
                        agreed = primary_canonical == secondary_canonical
                    else:
                        agreed = primary_smiles == secondary_smiles
                    result["cross_backend_agreement"] = agreed
                    result["agreement_checked"] = True
                    result["agreement_reference_backend"] = secondary_backend
                    result["agreement_reference_confidence"] = secondary_conf

        return result

    def validate_molecule(
        self,
        smiles: str,
        extraction_confidence: Optional[float] = None,
        cross_backend_agreement: Optional[bool] = None,
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
            "standardization": None,
            "level_2_chemical": None,
            "level_3_domain": None,
            "is_valid": False,
            "confidence_score": 0,
            "confidence_breakdown": {},
            "triage_label": "REJECT_OR_DEEP_REVIEW",
        }

        # Level 1: Syntax (Indigo)
        syntax_result = self.indigo.validate_syntax(smiles)
        validation["level_1_syntax"] = syntax_result

        if not syntax_result["is_valid"]:
            return validation

        # Use canonical SMILES for subsequent validation
        canonical = syntax_result.get("canonical_smiles", smiles)

        # Deterministic standardization (between Level 1 and Level 2)
        standardization_result = self.standardizer.standardize(canonical)
        validation["standardization"] = standardization_result

        if not standardization_result["success"]:
            return validation

        standardized_smiles = standardization_result["standardized_smiles"]

        # Level 2: Chemical Plausibility (on standardized structure)
        chemical_result = self.chemical.validate(
            standardized_smiles, standardized_smiles
        )
        validation["level_2_chemical"] = chemical_result

        if not chemical_result["is_valid"]:
            return validation

        # Level 3: Domain (Ethnopharmacology) on standardized structure
        domain_result = self.domain.validate(standardized_smiles, standardized_smiles)
        validation["level_3_domain"] = domain_result

        # Compute confidence score
        score, breakdown = self._compute_confidence_score(
            syntax_result,
            chemical_result,
            domain_result,
            extraction_confidence=extraction_confidence,
            cross_backend_agreement=cross_backend_agreement,
        )
        validation["confidence_score"] = score
        validation["confidence_breakdown"] = breakdown
        validation["triage_label"] = self._label_confidence(score)
        validation["is_valid"] = score >= self.min_validation_score

        return validation

    def _compute_confidence_score(
        self,
        syntax: Dict[str, Any],
        chemical: Dict[str, Any],
        domain: Dict[str, Any],
        extraction_confidence: Optional[float] = None,
        cross_backend_agreement: Optional[bool] = None,
    ) -> tuple[int, Dict[str, Any]]:
        """
        Compute overall confidence score (0-100).

        Uses configurable base/bonus/penalty policy from validation_rules.yaml.
        """
        scoring_cfg = self.validation_rules.get("confidence_scoring", {})
        base_cfg = scoring_cfg.get("base_scores", {})
        bonus_cfg = scoring_cfg.get("bonuses", {})
        penalty_cfg = scoring_cfg.get("penalties", {})

        breakdown = {
            "base": {},
            "bonuses": {},
            "penalties": {},
            "requirements": {},
            "inputs": {
                "extraction_confidence": extraction_confidence,
                "cross_backend_agreement": cross_backend_agreement,
                "has_gold_match": bool(domain.get("gold_standard_matches")),
                "matches_sceletium_scaffold": bool(
                    domain.get("matches_sceletium_scaffold")
                ),
            },
        }
        score = 0

        # Base scores by validation level
        if syntax["is_valid"]:
            points = int(base_cfg.get("level_1_syntax_pass", 30))
            score += points
            breakdown["base"]["level_1_syntax_pass"] = points
        if chemical["is_valid"]:
            points = int(base_cfg.get("level_2_chemical_pass", 30))
            score += points
            breakdown["base"]["level_2_chemical_pass"] = points
        if domain["is_valid"]:
            points = int(base_cfg.get("level_3_domain_pass", 20))
            score += points
            breakdown["base"]["level_3_domain_pass"] = points

        # Bonuses
        has_exact_match = False
        has_similar_match = False
        gold_similarity_threshold = (
            self.validation_rules.get("level_3_domain", {})
            .get("gold_standard_comparison", {})
            .get("similarity_threshold", 0.7)
        )
        if domain.get("gold_standard_matches"):
            for match in domain["gold_standard_matches"]:
                if match["similarity"] >= 0.99:
                    has_exact_match = True
                    break
                elif match["similarity"] >= float(gold_similarity_threshold):
                    has_similar_match = True
                    break

        positive_indicators = {}
        if has_exact_match:
            points = int(bonus_cfg.get("exact_gold_standard_match", 20))
            score += points
            breakdown["bonuses"]["exact_gold_standard_match"] = points
            positive_indicators["exact_gold_standard_match"] = True
        elif has_similar_match:
            # Similar-match bonus is kept for backwards compatibility.
            points = int(bonus_cfg.get("similar_gold_standard_match", 10))
            score += points
            breakdown["bonuses"]["similar_gold_standard_match"] = points

        if domain.get("matches_sceletium_scaffold"):
            points = int(bonus_cfg.get("scaffold_match", 10))
            score += points
            breakdown["bonuses"]["scaffold_match"] = points
            positive_indicators["scaffold_match"] = True

        if domain.get("methoxy_count", 0) >= 1:
            points = int(bonus_cfg.get("methoxy_groups_present", 5))
            score += points
            breakdown["bonuses"]["methoxy_groups_present"] = points
            positive_indicators["methoxy_groups_present"] = True

        if domain.get("mw_in_range"):
            points = int(bonus_cfg.get("mw_in_expected_range", 5))
            score += points
            breakdown["bonuses"]["mw_in_expected_range"] = points

        if cross_backend_agreement is True:
            points = int(bonus_cfg.get("cross_backend_agreement", 5))
            score += points
            breakdown["bonuses"]["cross_backend_agreement"] = points

        if "@" in (syntax.get("canonical_smiles") or ""):
            points = int(bonus_cfg.get("stereochemistry_present", 5))
            score += points
            breakdown["bonuses"]["stereochemistry_present"] = points

        # Apply penalties
        low_conf_threshold = (
            self.backends_config.get("fallback", {})
            .get("low_confidence_review", {})
            .get("threshold", 0.6)
        )
        if extraction_confidence is not None and extraction_confidence < float(
            low_conf_threshold
        ):
            points = int(penalty_cfg.get("low_confidence_backend", -15))
            score += points
            breakdown["penalties"]["low_confidence_backend"] = points

        if chemical.get("properties", {}).get("num_fragments", 1) > 1:
            points = int(penalty_cfg.get("multiple_fragments", -10))
            score += points
            breakdown["penalties"]["multiple_fragments"] = points

        if syntax.get("wildcard_count", 0) > 0:
            points = int(penalty_cfg.get("wildcards_present", -5))
            score += points
            breakdown["penalties"]["wildcards_present"] = points

        if domain.get("matches_sceletium_scaffold") and "@" not in (
            syntax.get("canonical_smiles") or ""
        ):
            points = int(penalty_cfg.get("no_stereochemistry", -5))
            score += points
            breakdown["penalties"]["no_stereochemistry"] = points

        mw = float(chemical.get("properties", {}).get("molecular_weight", 0.0) or 0.0)
        domain_mw_cfg = (
            self.validation_rules.get("level_3_domain", {})
            .get("sceletium_alkaloids", {})
            .get("expected_mw_range", {})
        )
        mw_min = float(domain_mw_cfg.get("min", 250))
        mw_max = float(domain_mw_cfg.get("max", 400))
        borderline_margin = float(scoring_cfg.get("borderline_mw_margin", 10.0))
        lower_edge = mw_min - borderline_margin
        upper_edge = mw_max + borderline_margin
        if mw and (lower_edge <= mw <= mw_min or mw_max <= mw <= upper_edge):
            points = int(penalty_cfg.get("borderline_mw", -5))
            score += points
            breakdown["penalties"]["borderline_mw"] = points

        requirements_cfg = scoring_cfg.get("requirements", {})
        positive_list = requirements_cfg.get("positive_indicators", [])
        has_positive_indicator = any(
            positive_indicators.get(name, False) for name in positive_list
        )
        require_positive = bool(
            requirements_cfg.get("require_positive_indicator", False)
        )
        high_threshold = int(
            scoring_cfg.get("thresholds", {}).get("high_confidence", 70)
        )
        if require_positive and not has_positive_indicator and score >= high_threshold:
            adjusted_score = high_threshold - 1
            breakdown["requirements"]["positive_indicator_required"] = True
            breakdown["requirements"]["score_capped_from"] = score
            breakdown["requirements"]["score_capped_to"] = adjusted_score
            score = adjusted_score

        return max(0, min(100, score)), breakdown

    def _label_confidence(self, score: int) -> str:
        """Assign configured triage label from confidence thresholds."""
        scoring_cfg = self.validation_rules.get("confidence_scoring", {})
        thresholds = scoring_cfg.get("thresholds", {})
        labels = scoring_cfg.get("labels", {})

        high_threshold = int(thresholds.get("high_confidence", 70))
        medium_threshold = int(thresholds.get("medium_confidence", 50))

        if score >= high_threshold:
            return labels.get("high", "READY")
        if score >= medium_threshold:
            return labels.get("medium", "MANUAL_REVIEW")
        return labels.get("low", "REJECT_OR_DEEP_REVIEW")

    def enrich_molecule(
        self,
        smiles: str,
        canonical_smiles: Optional[str] = None,
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
        prefilter_result: Optional[Dict[str, Any]] = None,
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
            self._record_rejection(
                self._build_rejection_event(
                    stage="extraction",
                    reason_code=str(extraction.get("error") or "extraction_failed"),
                    reason_detail="No backend returned a valid SMILES output",
                    image_path=image_path,
                    backend_order=self.backend_order,
                    backend_errors=extraction.get("errors"),
                )
            )
            return None

        self.stats["extractions_successful"] += 1

        extraction_confidence = float(extraction.get("confidence", 0.0))
        if (
            extraction_confidence < self.min_extraction_confidence
            and extraction.get("cross_backend_agreement") is not True
        ):
            self.stats["low_extraction_confidence_rejected"] += 1
            self._record_rejection(
                self._build_rejection_event(
                    stage="extraction",
                    reason_code="low_extraction_confidence",
                    reason_detail=(
                        f"confidence={extraction_confidence:.3f} < "
                        f"min_extraction_confidence={self.min_extraction_confidence:.3f}"
                    ),
                    image_path=image_path,
                    confidence=extraction_confidence,
                    min_extraction_confidence=self.min_extraction_confidence,
                    backend_used=extraction.get("backend_used"),
                    cross_backend_agreement=extraction.get("cross_backend_agreement"),
                )
            )
            return None

        requires_manual_review = False
        review_reasons: List[str] = []
        if (prefilter_result or {}).get("soft_failed"):
            requires_manual_review = True
            review_reasons.append("image_too_small")
            self.stats["manual_review_candidates"] += 1

        if (
            self.low_confidence_review_enabled
            and extraction_confidence < self.low_confidence_review_threshold
        ):
            requires_manual_review = True
            review_reasons.append("low_extraction_confidence")
            self.stats["manual_review_candidates"] += 1

        # Step 2: Validate
        validation = self.validate_molecule(
            extraction["smiles"],
            extraction_confidence=extraction_confidence,
            cross_backend_agreement=extraction.get("cross_backend_agreement"),
        )

        if validation["level_1_syntax"]["is_valid"]:
            self.stats["syntax_valid"] += 1
        if (validation.get("level_2_chemical") or {}).get("is_valid"):
            self.stats["chemically_valid"] += 1
        if (validation.get("level_3_domain") or {}).get("is_valid"):
            self.stats["domain_valid"] += 1
        if validation["confidence_score"] >= self.min_validation_score:
            self.stats["high_confidence"] += 1

        if not validation["is_valid"]:
            if validation.get("triage_label") == "MANUAL_REVIEW":
                requires_manual_review = True
                review_reasons.append("below_min_validation_score")
                self.stats["manual_review_candidates"] += 1
            else:
                self._record_rejection(
                    self._derive_validation_rejection(
                        image_path, validation, extraction
                    )
                )
                return None

        # Step 3: Enrich (optional)
        canonical = validation["level_1_syntax"].get("canonical_smiles")
        standardization = validation.get("standardization") or {}
        standardized_smiles = standardization.get("standardized_smiles") or canonical
        enrichment = {}

        if include_enrichment:
            enrichment = self.enrich_molecule(standardized_smiles, standardized_smiles)

        # Build complete record
        record = {
            "id": f"mol_{self.stats['images_processed']:06d}",
            "smiles": extraction["smiles"],
            "canonical_smiles": canonical,
            "standardized_smiles": standardized_smiles,
            "standard_inchi": standardization.get("standard_inchi"),
            "standard_inchikey": standardization.get("standard_inchikey"),
            "standardization_policy_version": standardization.get("policy_version"),
            "backend_used": extraction["backend_used"],
            "extraction_confidence": extraction_confidence,
            "backend_runtime": extraction.get("backend_runtime"),
            "cross_backend_agreement": extraction.get("cross_backend_agreement"),
            "agreement_checked": extraction.get("agreement_checked"),
            "agreement_reference_backend": extraction.get(
                "agreement_reference_backend"
            ),
            "agreement_reference_confidence": extraction.get(
                "agreement_reference_confidence"
            ),
            "agreement_reference_runtime": extraction.get(
                "agreement_reference_runtime"
            ),
            "validation": validation,
            "confidence_score": validation["confidence_score"],
            "triage_label": validation.get("triage_label"),
            "ready_for_kb": (
                validation.get("triage_label") == "READY" and not requires_manual_review
            ),
            "compound_class": validation["level_3_domain"].get(
                "compound_class", "unknown"
            ),
            "requires_manual_review": requires_manual_review,
            "review_reasons": review_reasons,
            "source": {
                "image_path": str(image_path),
                "extracted_at": datetime.now().isoformat(),
                "run_id": self.run_id,
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

        prefilter_context = prefilter_result or {}
        if prefilter_context.get("soft_failed"):
            record["prefilter"] = {
                "reason_code": prefilter_context.get("reason_code"),
                "reason_detail": prefilter_context.get("reason_detail"),
                "metrics": prefilter_context.get("metrics"),
            }

        return record

    def process_directory(
        self,
        input_dir: Path,
        output_dir: Path,
        limit: Optional[int] = None,
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

        self.stats["images_discovered"] = len(image_files)
        print(f"Found {len(image_files)} images to process")

        # Process images
        molecules = []
        for image_path in tqdm(image_files, desc="Processing images"):
            prefilter_result = self._evaluate_image_prefilter(image_path)
            if not prefilter_result.get("is_candidate", True):
                reason_code = str(
                    prefilter_result.get("reason_code") or "prefilter_rejected"
                )
                if reason_code == "image_too_small":
                    prefilter_result = dict(prefilter_result)
                    prefilter_result["is_candidate"] = True
                    prefilter_result["soft_failed"] = True
                    self.stats["images_prefilter_soft_review"] += 1
                else:
                    self.stats["images_prefilter_rejected"] += 1
                    self._record_rejection(
                        self._build_rejection_event(
                            stage="prefilter",
                            reason_code=reason_code,
                            reason_detail=str(
                                prefilter_result.get("reason_detail")
                                or "Rejected by image prefilter"
                            ),
                            image_path=image_path,
                            prefilter_metrics=prefilter_result.get("metrics"),
                        )
                    )
                    continue
            try:
                record = self.process_image(
                    image_path, prefilter_result=prefilter_result
                )
                if record:
                    molecules.append(record)
            except Exception as e:
                self._record_rejection(
                    self._build_rejection_event(
                        stage="runtime",
                        reason_code="unexpected_exception",
                        reason_detail=str(e),
                        image_path=image_path,
                    )
                )

        # Save tiered outputs for downstream import workflows
        output_file = output_dir / "molecules.jsonl"
        high_confidence_file = output_dir / "molecules_high_confidence.jsonl"
        manual_review_file = output_dir / "molecules_manual_review.jsonl"

        high_confidence_molecules = [
            mol for mol in molecules if bool(mol.get("ready_for_kb"))
        ]
        manual_review_molecules = [
            mol
            for mol in molecules
            if bool(mol.get("requires_manual_review"))
            or mol.get("triage_label") == "MANUAL_REVIEW"
        ]

        with open(output_file, "w") as f:
            for mol in molecules:
                f.write(json.dumps(mol) + "\n")

        with open(high_confidence_file, "w") as f:
            for mol in high_confidence_molecules:
                f.write(json.dumps(mol) + "\n")

        with open(manual_review_file, "w") as f:
            for mol in manual_review_molecules:
                f.write(json.dumps(mol) + "\n")

        # Save statistics
        stats_file = output_dir / "extraction_stats.json"
        self.stats["completed_at_utc"] = datetime.utcnow().isoformat() + "Z"
        self.stats["images_accepted"] = len(molecules)
        self.stats["high_confidence"] = len(high_confidence_molecules)
        self.stats["manual_review_candidates"] = len(manual_review_molecules)
        self.stats["acceptance_rate"] = (
            round(len(molecules) / self.stats["images_processed"], 6)
            if self.stats["images_processed"]
            else 0.0
        )
        with open(stats_file, "w") as f:
            json.dump(self.stats, f, indent=2)

        print(f"\nProcessed {len(molecules)} validated molecules")
        print(f"Results saved to: {output_file}")
        print(
            "Tiered outputs: "
            f"high_confidence={len(high_confidence_molecules)} -> {high_confidence_file.name}, "
            f"manual_review={len(manual_review_molecules)} -> {manual_review_file.name}"
        )

        return molecules


def run_preflight(pipeline: ExtractionPipeline) -> Dict[str, Any]:
    """
    Execute runtime dependency checks before full extraction.

    Returns a structured report with pass/warn/fail outcomes.
    """
    checks: List[Dict[str, Any]] = []

    def _add_check(name: str, status: str, detail: str, **extra: Any) -> None:
        entry = {"name": name, "status": status, "detail": detail}
        for key, value in extra.items():
            if value is not None:
                entry[key] = value
        checks.append(entry)

    runtime_manifest = pipeline._build_runtime_manifest()
    dependencies = runtime_manifest.get("dependencies", {})

    def _parse_version_triplet(raw: Optional[str]) -> Optional[tuple[int, int, int]]:
        """Parse major.minor.patch triplet from a version string."""
        if not raw:
            return None
        cleaned = raw.strip()
        core = cleaned.split("+", 1)[0]
        for marker in ("rc", "a", "b", "dev", "post"):
            if marker in core:
                core = core.split(marker, 1)[0]
        parts = core.split(".")
        if len(parts) < 3:
            return None
        try:
            return int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            return None

    for cfg_name, cfg_path in (pipeline.config_sources or {}).items():
        cfg_file = Path(cfg_path)
        if cfg_file.exists():
            _add_check(
                name=f"config:{cfg_name}",
                status="pass",
                detail=f"Loaded config file: {cfg_file}",
            )
        else:
            _add_check(
                name=f"config:{cfg_name}",
                status="fail",
                detail=f"Missing config file: {cfg_file}",
            )

    numpy_version = dependencies.get("numpy")
    if numpy_version:
        numpy_major = str(numpy_version).split(".", 1)[0]
        if numpy_major.isdigit() and int(numpy_major) >= 2:
            _add_check(
                name="dependency:numpy",
                status="warn",
                detail=(
                    f"NumPy {numpy_version} detected; ensure RDKit build is ABI-compatible "
                    "(many wheels still target NumPy <2)."
                ),
                version=numpy_version,
            )
        else:
            _add_check(
                name="dependency:numpy",
                status="pass",
                detail=f"NumPy version looks compatible: {numpy_version}",
                version=numpy_version,
            )
    else:
        _add_check(
            name="dependency:numpy",
            status="warn",
            detail="NumPy version not detected via package metadata.",
        )

    # Indigo version contract
    indigo_version = dependencies.get("epam-indigo")
    indigo_min_supported = (1, 40, 0)
    indigo_max_validated_exclusive = (1, 41, 0)
    if not indigo_version:
        _add_check(
            name="dependency:indigo_version",
            status="warn",
            detail=(
                "Indigo package version not detected via metadata; "
                "validated range is >=1.40.0,<1.41.0."
            ),
        )
    else:
        parsed_indigo_version = _parse_version_triplet(str(indigo_version))
        if parsed_indigo_version is None:
            _add_check(
                name="dependency:indigo_version",
                status="warn",
                detail=(
                    f"Could not parse Indigo version '{indigo_version}'; "
                    "validated range is >=1.40.0,<1.41.0."
                ),
                version=indigo_version,
            )
        elif parsed_indigo_version < indigo_min_supported:
            _add_check(
                name="dependency:indigo_version",
                status="fail",
                detail=(
                    f"Indigo {indigo_version} is below minimum supported 1.40.0; "
                    "upgrade required."
                ),
                version=indigo_version,
            )
        elif parsed_indigo_version >= indigo_max_validated_exclusive:
            _add_check(
                name="dependency:indigo_version",
                status="warn",
                detail=(
                    f"Indigo {indigo_version} is above validated range (<1.41.0); "
                    "run regression tests before production use."
                ),
                version=indigo_version,
            )
        else:
            _add_check(
                name="dependency:indigo_version",
                status="pass",
                detail=f"Indigo version within validated contract: {indigo_version}",
                version=indigo_version,
            )

    # Indigo readiness
    try:
        indigo = pipeline.indigo
        can = indigo.canonicalize("C[C@H](O)N")
        if can:
            _add_check(
                name="dependency:indigo",
                status="pass",
                detail="Indigo loaded and canonicalization succeeded.",
                canonical_smiles=can,
            )
        else:
            _add_check(
                name="dependency:indigo",
                status="fail",
                detail="Indigo loaded but canonicalization returned empty result.",
            )
        if hasattr(indigo, "get_runtime_policy"):
            policy = indigo.get_runtime_policy()
            if policy.get("option_errors"):
                _add_check(
                    name="indigo:runtime_options",
                    status="warn",
                    detail="Some Indigo runtime options failed to apply.",
                    option_errors=policy.get("option_errors"),
                )
            else:
                _add_check(
                    name="indigo:runtime_options",
                    status="pass",
                    detail="Indigo runtime options applied cleanly.",
                    policy=policy,
                )
    except Exception as exc:
        _add_check(
            name="dependency:indigo",
            status="fail",
            detail=str(exc),
        )

    def _truncate(text: str, limit: int = 220) -> str:
        if not text:
            return ""
        cleaned = " ".join(text.strip().split())
        if len(cleaned) <= limit:
            return cleaned
        return cleaned[: limit - 3] + "..."

    def _extract_json(stdout: str) -> Dict[str, Any]:
        for line in reversed((stdout or "").splitlines()):
            candidate = line.strip()
            if not candidate.startswith("{") or not candidate.endswith("}"):
                continue
            try:
                return json.loads(candidate)
            except json.JSONDecodeError:
                continue
        return {}

    def _run_probe(code: str, timeout_sec: int = 20) -> Dict[str, Any]:
        try:
            proc = subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True,
                text=True,
                timeout=timeout_sec,
            )
            return {
                "returncode": proc.returncode,
                "stdout": proc.stdout or "",
                "stderr": proc.stderr or "",
                "timed_out": False,
            }
        except subprocess.TimeoutExpired as exc:
            return {
                "returncode": -9,
                "stdout": exc.stdout or "",
                "stderr": exc.stderr or "",
                "timed_out": True,
            }

    # RDKit-backed readiness checks run in isolated subprocesses to avoid
    # interpreter crashes from binary ABI mismatches.
    probe_smiles = "OCCN(CC)CC"
    src_dir = str(Path(__file__).resolve().parent.parent / "src")

    rdkit_core_probe = _run_probe(
        "import json,sys\n"
        f"probe = {probe_smiles!r}\n"
        "try:\n"
        "    from rdkit import Chem\n"
        "    mol = Chem.MolFromSmiles(probe)\n"
        "    if mol is None:\n"
        "        print(json.dumps({'ok': False, 'detail': 'MolFromSmiles returned None'}))\n"
        "        sys.exit(1)\n"
        "    print(json.dumps({'ok': True, 'detail': 'RDKit core import probe succeeded'}))\n"
        "except Exception as exc:\n"
        "    print(json.dumps({'ok': False, 'detail': str(exc)}))\n"
        "    sys.exit(1)\n"
    )
    rdkit_core_payload = _extract_json(rdkit_core_probe["stdout"])
    rdkit_core_ok = rdkit_core_probe["returncode"] == 0 and bool(
        rdkit_core_payload.get("ok")
    )
    if rdkit_core_ok:
        _add_check(
            name="dependency:rdkit_core",
            status="pass",
            detail=rdkit_core_payload.get(
                "detail", "RDKit subprocess probe completed."
            ),
        )
    else:
        core_detail_parts = []
        if rdkit_core_payload.get("detail"):
            core_detail_parts.append(str(rdkit_core_payload.get("detail")))
        if rdkit_core_probe.get("timed_out"):
            core_detail_parts.append("probe timed out")
        core_detail_parts.append(f"exit_code={rdkit_core_probe['returncode']}")
        stderr_excerpt = _truncate(rdkit_core_probe.get("stderr", ""))
        if stderr_excerpt:
            core_detail_parts.append(f"stderr={stderr_excerpt}")
        _add_check(
            name="dependency:rdkit_core",
            status="fail",
            detail="; ".join(core_detail_parts),
        )

    if rdkit_core_ok:
        probe_specs = [
            (
                "dependency:rdkit_standardization",
                "from validators.standardization_validator import StandardizationValidator\n"
                f"res = StandardizationValidator().standardize({probe_smiles!r})\n"
                "ok = bool(res.get('success'))\n"
                "detail = 'standardization_ok' if ok else str(res.get('error') or 'standardization_failed')\n",
            ),
            (
                "dependency:rdkit_chemical_validation",
                "from validators.chemical_validator import ChemicalValidator\n"
                f"res = ChemicalValidator().validate({probe_smiles!r}, {probe_smiles!r})\n"
                "ok = isinstance(res, dict)\n"
                "detail = 'chemical_validator_executed' if ok else 'chemical_validator_failed'\n",
            ),
            (
                "dependency:domain_validation",
                "from validators.domain_validator import DomainValidator\n"
                f"res = DomainValidator().validate({probe_smiles!r}, {probe_smiles!r})\n"
                "ok = isinstance(res, dict) and ('compound_class' in res)\n"
                "detail = 'domain_validator_executed' if ok else 'domain_validator_failed'\n",
            ),
            (
                "dependency:property_calculator",
                "from enrichers.property_calculator import PropertyCalculator\n"
                f"res = PropertyCalculator().calculate({probe_smiles!r}, {probe_smiles!r})\n"
                "ok = bool(res.get('success'))\n"
                "detail = 'property_calculator_ok' if ok else str(res.get('error') or 'property_calculator_failed')\n",
            ),
            (
                "dependency:fingerprint_generator",
                "from enrichers.fingerprint_generator import FingerprintGenerator\n"
                f"res = FingerprintGenerator().generate({probe_smiles!r}, {probe_smiles!r})\n"
                "ok = bool(res.get('success'))\n"
                "detail = 'fingerprint_generator_ok' if ok else str(res.get('error') or 'fingerprint_generation_failed')\n",
            ),
        ]

        for check_name, probe_body in probe_specs:
            probe_code = (
                "import json,sys\n"
                f"sys.path.insert(0, {src_dir!r})\n"
                "try:\n"
                + "\n".join(f"    {line}" for line in probe_body.splitlines())
                + "\n"
                "    print(json.dumps({'ok': bool(ok), 'detail': str(detail)}))\n"
                "    sys.exit(0 if ok else 1)\n"
                "except Exception as exc:\n"
                "    print(json.dumps({'ok': False, 'detail': str(exc)}))\n"
                "    sys.exit(1)\n"
            )
            probe_result = _run_probe(probe_code)
            payload = _extract_json(probe_result["stdout"])
            ok = probe_result["returncode"] == 0 and bool(payload.get("ok"))
            detail_parts = []
            if payload.get("detail"):
                detail_parts.append(str(payload.get("detail")))
            if probe_result.get("timed_out"):
                detail_parts.append("probe timed out")
            if not ok:
                detail_parts.append(f"exit_code={probe_result['returncode']}")
                stderr_excerpt = _truncate(probe_result.get("stderr", ""))
                if stderr_excerpt:
                    detail_parts.append(f"stderr={stderr_excerpt}")
            _add_check(
                name=check_name,
                status="pass" if ok else "fail",
                detail="; ".join(detail_parts) if detail_parts else "probe completed",
            )
    else:
        blocked_detail = (
            "Skipped: rdkit_core subprocess probe failed in this environment."
        )
        _add_check(
            name="dependency:rdkit_standardization",
            status="fail",
            detail=blocked_detail,
        )
        _add_check(
            name="dependency:rdkit_chemical_validation",
            status="fail",
            detail=blocked_detail,
        )
        _add_check(
            name="dependency:domain_validation",
            status="fail",
            detail=blocked_detail,
        )
        _add_check(
            name="dependency:property_calculator",
            status="fail",
            detail=blocked_detail,
        )
        _add_check(
            name="dependency:fingerprint_generator",
            status="fail",
            detail=blocked_detail,
        )

    has_failures = any(check["status"] == "fail" for check in checks)
    has_warnings = any(check["status"] == "warn" for check in checks)
    overall_status = "fail" if has_failures else ("warn" if has_warnings else "pass")

    return {
        "status": overall_status,
        "generated_at_utc": datetime.utcnow().isoformat() + "Z",
        "run_id": pipeline.run_id,
        "checks": checks,
        "runtime_manifest": runtime_manifest,
    }


def main():
    """Main entry point."""
    default_config_dir = Path(__file__).resolve().parent.parent / "config"
    parser = argparse.ArgumentParser(
        description="SMILES extraction pipeline with 3-layer validation"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=False,
        help="Directory with chemical structure images",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=False,
        help="Output directory for results",
    )
    parser.add_argument(
        "--backend-order",
        type=str,
        default=None,
        help="Comma-separated OCSR backend order (defaults to enabled+priority in backends.yaml)",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=default_config_dir,
        help="Directory containing backends.yaml and validation_rules.yaml",
    )
    parser.add_argument(
        "--backends-config-file",
        type=Path,
        default=None,
        help="Optional override path to backends.yaml",
    )
    parser.add_argument(
        "--validation-rules-file",
        type=Path,
        default=None,
        help="Optional override path to validation_rules.yaml",
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
        default=None,
        help="Minimum OCSR confidence (defaults to config)",
    )
    parser.add_argument(
        "--min-validation-score",
        type=int,
        default=None,
        help="Minimum validation score (0-100, defaults to config)",
    )
    parser.add_argument(
        "--gold-standards-file",
        type=str,
        default=None,
        help=(
            "Path to gold_standards.json "
            "(overrides SMILES_GOLD_STANDARDS_FILE and default repo path)"
        ),
    )
    parser.add_argument(
        "--preflight",
        action="store_true",
        help="Run runtime dependency/config checks and exit",
    )

    args = parser.parse_args()
    if not args.preflight and (args.input_dir is None or args.output_dir is None):
        parser.error(
            "--input-dir and --output-dir are required unless --preflight is used"
        )

    backends_config_file = args.backends_config_file or (
        args.config_dir / "backends.yaml"
    )
    validation_rules_file = args.validation_rules_file or (
        args.config_dir / "validation_rules.yaml"
    )
    backends_config, backends_config_error = load_yaml_config(backends_config_file)
    validation_rules, validation_rules_error = load_yaml_config(validation_rules_file)
    config_errors = [
        err for err in (backends_config_error, validation_rules_error) if err is not None
    ]
    for err in config_errors:
        print(err, file=sys.stderr)
    if config_errors and not args.preflight:
        return 1

    # Parse backend order
    backend_order = (
        [b.strip() for b in args.backend_order.split(",") if b.strip()]
        if args.backend_order
        else None
    )

    # Determine GPU mode
    use_gpu = args.gpu or (not args.no_gpu)

    # Create pipeline
    pipeline = ExtractionPipeline(
        backend_order=backend_order,
        use_gpu=use_gpu,
        confidence_threshold=args.confidence_threshold,
        min_validation_score=args.min_validation_score,
        gold_standards_file=args.gold_standards_file,
        backends_config=backends_config,
        validation_rules=validation_rules,
        config_sources={
            "backends_config_file": str(backends_config_file),
            "validation_rules_file": str(validation_rules_file),
        },
    )

    print(f"Using backends config: {backends_config_file}")
    print(f"Using validation rules: {validation_rules_file}")
    print(f"Effective backend order: {','.join(pipeline.backend_order)}")
    if pipeline.unsupported_backends:
        print(
            "Ignored unsupported backends: " + ",".join(pipeline.unsupported_backends)
        )
    print(f"Effective OCSR confidence threshold: {pipeline.confidence_threshold}")
    print(f"Effective min validation score: {pipeline.min_validation_score}")

    if args.preflight:
        report = run_preflight(pipeline)
        print("\n" + "=" * 60)
        print("Preflight Report")
        print("=" * 60)
        print(json.dumps(report, indent=2))
        return 0 if report["status"] == "pass" else 1

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
        f"High confidence (≥{pipeline.min_validation_score}): {pipeline.stats['high_confidence']}"
    )
    print(f"Errors: {len(pipeline.stats['errors'])}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
