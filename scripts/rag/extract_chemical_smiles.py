#!/usr/bin/env python3

from __future__ import annotations

import argparse
import base64
import html
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import cv2
import numpy as np

try:
    from rdkit import Chem

    RDKIT_AVAILABLE = True
except Exception:
    Chem = None
    RDKIT_AVAILABLE = False


REPO_ROOT = Path(__file__).resolve().parents[2]
EXTRACTIONS_DIR = REPO_ROOT / "data" / "extractions"
DEFAULT_OUTPUT_NAME = "smiles_extracted.json"
DEFAULT_TIMEOUT_SECONDS = 1800
DEFAULT_MIN_CONFIDENCE = 0.5
DEFAULT_BACKEND = "molgrapher"
DEFAULT_DETECTOR = "moldetv2-general"
DEFAULT_DETECTOR_CONFIDENCE = 0.5
DEFAULT_DETECTOR_IOU = 0.5
DEFAULT_MOLDETV2_REPO = "UniParser/MolDetv2"
DEFAULT_MOLDETV2_GENERAL_FILE = "moldet_v2_yolo11n_640_general.onnx"
DEFAULT_MOLDETV2_DOC_FILE = "moldet_v2_yolo11n_960_doc.onnx"
DEFAULT_DETECTOR_PADDING = 0.03
DEFAULT_MOLSCRIBE_CKPT_REPO = "yujieq/MolScribe"
DEFAULT_MOLSCRIBE_CKPT_FILE = "swin_base_char_aux_1m680k.pth"
MOLGRAPHER_RUNNER_CODE = r"""
import json
import os
import sys
from pathlib import Path

import molgrapher
from molgrapher.models.molgrapher_model import MolgrapherModel

manifest_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
user_workdir = Path(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else None

if user_workdir:
    os.chdir(user_workdir)
else:
    os.chdir(Path(molgrapher.__file__).resolve().parents[1])

image_paths = json.loads(manifest_path.read_text(encoding="utf-8"))
model = MolgrapherModel()
predictions = model.predict_batch(image_paths)
output_path.write_text(json.dumps(predictions, ensure_ascii=True), encoding="utf-8")
"""
MOLSCRIBE_RUNNER_CODE = r"""
import json
import sys
from pathlib import Path

import torch
from huggingface_hub import hf_hub_download
from molscribe import MolScribe

manifest_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
checkpoint_path = sys.argv[3] if len(sys.argv) > 3 else ""
checkpoint_repo = sys.argv[4] if len(sys.argv) > 4 else "yujieq/MolScribe"
checkpoint_file = sys.argv[5] if len(sys.argv) > 5 else "swin_base_char_aux_1m680k.pth"

image_paths = json.loads(manifest_path.read_text(encoding="utf-8"))
resolved_checkpoint = checkpoint_path or hf_hub_download(checkpoint_repo, checkpoint_file)
model = MolScribe(resolved_checkpoint, device=torch.device("cpu"))
predictions = model.predict_image_files(image_paths, return_atoms_bonds=True, return_confidence=True)
output_path.write_text(json.dumps(predictions, ensure_ascii=True), encoding="utf-8")
"""


@dataclass
class ChemicalImage:
    block_id: str
    page: int | None
    bbox: list[float] | None
    image_name: str
    description: str
    image_value: str


@dataclass
class ChemicalCrop:
    block_id: str
    page: int | None
    bbox: list[float] | None
    image_name: str
    description: str
    image_bytes: bytes
    source_image_name: str
    detector: str
    detector_confidence: float | None
    detector_bbox: list[float] | None
    crop_index: int


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def strip_tags(value: str) -> str:
    text = re.sub(r"<[^>]+>", " ", value)
    return " ".join(html.unescape(text).split())


def extract_img_alt(value: str) -> str:
    match = re.search(r'<img[^>]*\salt="([^"]+)"', value, flags=re.IGNORECASE)
    return html.unescape(match.group(1)).strip() if match else ""


def normalize_text(value: str) -> str:
    lowered = value.casefold()
    lowered = re.sub(r"[^0-9a-z]+", " ", lowered)
    return " ".join(lowered.split())


def build_description(block: dict[str, Any]) -> str:
    html_value = str(block.get("html") or "")
    parts = [extract_img_alt(html_value), strip_tags(html_value)]
    deduped: list[str] = []
    seen: set[str] = set()
    for part in parts:
        normalized = normalize_text(part)
        if part and normalized and normalized not in seen:
            deduped.append(part)
            seen.add(normalized)
    return " ".join(deduped).strip()


def decode_image_bytes(value: str) -> bytes:
    cleaned = value.strip()
    if "," in cleaned and cleaned.lower().startswith("data:image/"):
        cleaned = cleaned.split(",", 1)[1]
    return base64.b64decode(cleaned)


def looks_like_image_block(block: dict[str, Any]) -> bool:
    block_id = str(block.get("id") or "")
    return "ChemicalBlock" in block_id and isinstance(block.get("images"), dict) and bool(block.get("images"))


def collect_chemical_images(paper_dir: Path) -> tuple[list[ChemicalImage], int]:
    normalized_path = paper_dir / "normalized.json"
    payload = json.loads(normalized_path.read_text(encoding="utf-8"))
    blocks = ((payload.get("raw") or {}).get("chunks") or {}).get("blocks") or []
    entries: list[ChemicalImage] = []
    total_images = 0
    for block in blocks:
        if not looks_like_image_block(block):
            continue
        description = build_description(block)
        for image_name, image_value in (block.get("images") or {}).items():
            if not isinstance(image_name, str) or not isinstance(image_value, str):
                continue
            total_images += 1
            entries.append(
                ChemicalImage(
                    block_id=str(block.get("id") or image_name),
                    page=block.get("page") if isinstance(block.get("page"), int) else None,
                    bbox=block.get("bbox") if isinstance(block.get("bbox"), list) else None,
                    image_name=image_name,
                    description=description,
                    image_value=image_value,
                )
            )
    return entries, total_images


def ensure_hf_hub_available() -> Any:
    try:
        from huggingface_hub import hf_hub_download
    except Exception as exc:
        raise RuntimeError("huggingface_hub is required for MolDetv2 model resolution") from exc
    return hf_hub_download


def resolve_detector_model_path(*, detector: str, detector_model_path: str | None, detector_model_repo: str) -> str | None:
    if detector == "none":
        return None
    if detector_model_path:
        return detector_model_path
    hf_hub_download = ensure_hf_hub_available()
    filename = {
        "moldetv2-general": DEFAULT_MOLDETV2_GENERAL_FILE,
        "moldetv2-doc": DEFAULT_MOLDETV2_DOC_FILE,
    }.get(detector)
    if not filename:
        raise ValueError(f"Unsupported detector: {detector}")
    return hf_hub_download(repo_id=detector_model_repo, filename=filename)


def decode_cv2_image(image_bytes: bytes) -> np.ndarray | None:
    array = np.frombuffer(image_bytes, dtype=np.uint8)
    if array.size == 0:
        return None
    return cv2.imdecode(array, cv2.IMREAD_COLOR)


def clamp_bbox(x1: float, y1: float, x2: float, y2: float, *, width: int, height: int) -> list[int] | None:
    left = max(0, min(width - 1, int(np.floor(x1))))
    top = max(0, min(height - 1, int(np.floor(y1))))
    right = max(0, min(width, int(np.ceil(x2))))
    bottom = max(0, min(height, int(np.ceil(y2))))
    if right - left < 4 or bottom - top < 4:
        return None
    return [left, top, right, bottom]


def detect_moldetv2_boxes(
    image: np.ndarray,
    model: Any,
    *,
    input_size: int,
    conf_threshold: float,
    iou_threshold: float,
) -> list[dict[str, Any]]:
    height, width = image.shape[:2]
    detections: list[dict[str, Any]] = []
    results = model.predict(image, imgsz=input_size, conf=conf_threshold, iou=iou_threshold, verbose=False)
    for result in results:
        for box in result.boxes:
            confidence = float(box.conf)
            coords = box.xyxy[0].tolist()
            bbox = clamp_bbox(float(coords[0]), float(coords[1]), float(coords[2]), float(coords[3]), width=width, height=height)
            if bbox is None:
                continue
            detections.append({"bbox": bbox, "confidence": confidence})
    detections.sort(key=lambda item: (-float(item["confidence"]), item["bbox"][1], item["bbox"][0]))
    return detections


def encode_png_bytes(image: np.ndarray) -> bytes:
    ok, encoded = cv2.imencode(".png", image)
    if not ok:
        raise RuntimeError("Failed to encode detector crop as PNG")
    return encoded.tobytes()


def crop_with_padding(image: np.ndarray, bbox: list[int], *, padding_ratio: float) -> tuple[bytes, list[int]]:
    height, width = image.shape[:2]
    left, top, right, bottom = bbox
    pad_x = int(round((right - left) * padding_ratio))
    pad_y = int(round((bottom - top) * padding_ratio))
    padded_left = max(0, left - pad_x)
    padded_top = max(0, top - pad_y)
    padded_right = min(width, right + pad_x)
    padded_bottom = min(height, bottom + pad_y)
    crop = image[padded_top:padded_bottom, padded_left:padded_right]
    if crop.size == 0:
        raise RuntimeError("Detector produced an empty crop")
    return encode_png_bytes(crop), [padded_left, padded_top, padded_right, padded_bottom]


def build_detector_crops(
    entries: list[ChemicalImage],
    *,
    detector: str,
    detector_model_path: str | None,
    detector_confidence: float,
    detector_iou: float,
    detector_padding: float,
) -> tuple[list[ChemicalCrop], int]:
    if detector == "none":
        return (
            [
                ChemicalCrop(
                    block_id=entry.block_id,
                    page=entry.page,
                    bbox=entry.bbox,
                    image_name=entry.image_name,
                    description=entry.description,
                    image_bytes=decode_image_bytes(entry.image_value),
                    source_image_name=entry.image_name,
                    detector="none",
                    detector_confidence=None,
                    detector_bbox=None,
                    crop_index=0,
                )
                for entry in entries
            ],
            0,
        )

    try:
        from ultralytics import YOLO
    except Exception as exc:
        raise RuntimeError("ultralytics is required for MolDetv2 detector inference") from exc

    if not detector_model_path:
        raise RuntimeError("Detector model path was not resolved")
    input_size = 640 if detector == "moldetv2-general" else 960
    model = YOLO(detector_model_path, task="detect")
    crops: list[ChemicalCrop] = []
    skipped_no_detections = 0
    for entry in entries:
        image = decode_cv2_image(decode_image_bytes(entry.image_value))
        if image is None:
            skipped_no_detections += 1
            continue
        detections = detect_moldetv2_boxes(
            image,
            model,
            input_size=input_size,
            conf_threshold=detector_confidence,
            iou_threshold=detector_iou,
        )
        if not detections:
            skipped_no_detections += 1
            continue
        for crop_index, detection in enumerate(detections):
            crop_bytes, padded_bbox = crop_with_padding(
                image,
                detection["bbox"],
                padding_ratio=detector_padding,
            )
            crop_name = f"{Path(entry.image_name).stem}__det{crop_index:02d}.png"
            crops.append(
                ChemicalCrop(
                    block_id=entry.block_id,
                    page=entry.page,
                    bbox=entry.bbox,
                    image_name=crop_name,
                    description=entry.description,
                    image_bytes=crop_bytes,
                    source_image_name=entry.image_name,
                    detector=detector,
                    detector_confidence=float(detection["confidence"]),
                    detector_bbox=[float(value) for value in padded_bbox],
                    crop_index=crop_index,
                )
            )
    return crops, skipped_no_detections


def canonicalize_smiles(value: str) -> tuple[str | None, str]:
    if not value:
        return None, "backend_empty"
    if not RDKIT_AVAILABLE or Chem is None:
        return value, "backend_unvalidated"
    mol = Chem.MolFromSmiles(value)
    if mol is None:
        return None, "invalid_smiles"
    return Chem.MolToSmiles(mol), "validated"


def parse_confidence(value: Any) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        cleaned = value.strip()
        if not cleaned:
            return None
        try:
            return float(cleaned)
        except ValueError:
            return None
    return None


def stage_images(entries: list[ChemicalCrop], temp_dir: Path) -> tuple[list[str], dict[str, ChemicalCrop]]:
    staged_paths: list[str] = []
    mapping: dict[str, ChemicalCrop] = {}
    for index, entry in enumerate(entries):
        suffix = Path(entry.image_name).suffix or ".png"
        staged_name = f"{index:04d}_{Path(entry.image_name).stem}{suffix}"
        staged_path = temp_dir / staged_name
        staged_path.write_bytes(entry.image_bytes)
        staged_paths.append(str(staged_path))
        mapping[str(staged_path)] = entry
    return staged_paths, mapping


def resolve_python_bin(value: str | None, *, backend: str) -> str:
    if value:
        return value
    env_var = {
        "molgrapher": "MOLGRAPHER_PYTHON_BIN",
        "molscribe": "MOLSCRIBE_PYTHON_BIN",
    }[backend]
    env_value = os.environ.get(env_var, "").strip()
    if env_value:
        return env_value
    return sys.executable


def run_molgrapher(
    *,
    image_paths: list[str],
    python_bin: str,
    timeout_seconds: int,
    workdir: str | None,
) -> list[dict[str, Any]]:
    with tempfile.TemporaryDirectory(prefix="molgrapher-run-") as temp_root:
        temp_root_path = Path(temp_root)
        manifest_path = temp_root_path / "images.json"
        output_path = temp_root_path / "predictions.json"
        manifest_path.write_text(json.dumps(image_paths, ensure_ascii=True), encoding="utf-8")
        result = subprocess.run(
            [
                python_bin,
                "-c",
                MOLGRAPHER_RUNNER_CODE,
                str(manifest_path),
                str(output_path),
                workdir or "",
            ],
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            cwd=str(REPO_ROOT),
        )
        if result.returncode != 0:
            stderr = result.stderr.strip() or "MolGrapher exited with a non-zero status"
            raise RuntimeError(stderr)
        if not output_path.exists():
            raise RuntimeError("MolGrapher completed without writing predictions.json")
        return json.loads(output_path.read_text(encoding="utf-8"))


def run_molscribe(
    *,
    image_paths: list[str],
    python_bin: str,
    timeout_seconds: int,
    checkpoint_path: str | None,
    checkpoint_repo: str,
    checkpoint_file: str,
) -> list[dict[str, Any]]:
    with tempfile.TemporaryDirectory(prefix="molscribe-run-") as temp_root:
        temp_root_path = Path(temp_root)
        manifest_path = temp_root_path / "images.json"
        output_path = temp_root_path / "predictions.json"
        manifest_path.write_text(json.dumps(image_paths, ensure_ascii=True), encoding="utf-8")
        result = subprocess.run(
            [
                python_bin,
                "-c",
                MOLSCRIBE_RUNNER_CODE,
                str(manifest_path),
                str(output_path),
                checkpoint_path or "",
                checkpoint_repo,
                checkpoint_file,
            ],
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            cwd=str(REPO_ROOT),
        )
        if result.returncode != 0:
            stderr = result.stderr.strip() or "MolScribe exited with a non-zero status"
            raise RuntimeError(stderr)
        if not output_path.exists():
            raise RuntimeError("MolScribe completed without writing predictions.json")
        return json.loads(output_path.read_text(encoding="utf-8"))


def run_backend(
    *,
    backend: str,
    image_paths: list[str],
    python_bin: str,
    timeout_seconds: int,
    workdir: str | None,
    molscribe_checkpoint_path: str | None,
    molscribe_checkpoint_repo: str,
    molscribe_checkpoint_file: str,
) -> list[dict[str, Any]]:
    if backend == "molgrapher":
        return run_molgrapher(
            image_paths=image_paths,
            python_bin=python_bin,
            timeout_seconds=timeout_seconds,
            workdir=workdir,
        )
    if backend == "molscribe":
        return run_molscribe(
            image_paths=image_paths,
            python_bin=python_bin,
            timeout_seconds=timeout_seconds,
            checkpoint_path=molscribe_checkpoint_path,
            checkpoint_repo=molscribe_checkpoint_repo,
            checkpoint_file=molscribe_checkpoint_file,
        )
    raise ValueError(f"Unsupported backend: {backend}")


def extract_backend_fields(backend: str, prediction: dict[str, Any]) -> tuple[str, float | None, dict[str, Any]]:
    if backend == "molgrapher":
        return (
            str(prediction.get("smi") or "").strip(),
            parse_confidence(prediction.get("conf")),
            {
                "abbreviations": prediction.get("abbreviations", []),
                "abbreviations_ocr": prediction.get("abbreviations_ocr", []),
                "annotator": prediction.get("annotator", {}),
                "file_info": prediction.get("file-info", {}),
            },
        )
    if backend == "molscribe":
        return (
            str(prediction.get("smiles") or "").strip(),
            parse_confidence(prediction.get("confidence")),
            {
                "molfile": prediction.get("molfile"),
                "atoms": prediction.get("atoms", []),
                "bonds": prediction.get("bonds", []),
            },
        )
    raise ValueError(f"Unsupported backend: {backend}")


def build_picture_record(
    entry: ChemicalCrop,
    prediction: dict[str, Any],
    *,
    backend: str,
    min_confidence: float | None,
) -> dict[str, Any]:
    raw_smiles, confidence, backend_fields = extract_backend_fields(backend, prediction)
    molecule_smiles, status = canonicalize_smiles(raw_smiles)
    if molecule_smiles and min_confidence is not None and confidence is not None and confidence < min_confidence:
        molecule_smiles = None
        status = "filtered_low_confidence"
    candidate = {
        "backend": backend,
        "smiles": raw_smiles or None,
        "canonical_smiles": molecule_smiles,
        "status": status,
        "confidence": confidence,
        **backend_fields,
    }
    return {
        "block_id": entry.block_id,
        "page": entry.page,
        "bbox": entry.bbox,
        "image": entry.image_name,
        "source_image": entry.source_image_name,
        "crop_index": entry.crop_index,
        "description": entry.description,
        "backend": backend,
        "detector": entry.detector,
        "detector_confidence": entry.detector_confidence,
        "detector_bbox": entry.detector_bbox,
        "backends_tried": [backend],
        "attempts": 1,
        "success": bool(molecule_smiles),
        "error": None if molecule_smiles else status,
        "status": status,
        "molecule_smiles": molecule_smiles,
        "raw_smiles": raw_smiles or None,
        "confidence": confidence,
        "raw_candidates": [candidate],
        **backend_fields,
    }


def build_empty_payload(
    paper_dir: Path,
    warning: str,
    *,
    backend: str,
    detector: str,
    detector_model_path: str | None,
    python_bin: str,
    min_confidence: float | None,
    total_images: int,
    skipped_no_detections: int,
) -> dict[str, Any]:
    relative_paper_dir = str(paper_dir.relative_to(REPO_ROOT))
    return {
        "paper": paper_dir.name,
        "paper_dir": relative_paper_dir,
        "generated_at": utc_now_iso(),
        "pictures_detected_total": total_images,
        "pictures_found": 0,
        "pictures_processed": 0,
        "pictures_with_smiles": 0,
        "pictures_skipped_no_detections": skipped_no_detections,
        "rdkit_available": RDKIT_AVAILABLE,
        "rdkit_status": "available" if RDKIT_AVAILABLE else "unavailable",
        "backend": backend,
        "backends_order": [backend],
        "backend_python_bin": python_bin,
        "detector": detector,
        "detector_model_path": detector_model_path,
        "min_confidence": min_confidence,
        "results": [],
        "picture_data": [],
        "smiles_extracted": [],
        "valid_smiles": [],
        "invalid_smiles": [],
        "images_found_warning": warning,
    }


def process_paper(
    paper_dir: Path,
    *,
    backend: str,
    detector: str,
    python_bin: str,
    output_name: str,
    timeout_seconds: int,
    min_confidence: float | None,
    overwrite: bool,
    dry_run: bool,
    backend_workdir: str | None,
    detector_model_path: str | None,
    detector_confidence: float,
    detector_iou: float,
    detector_padding: float,
    molscribe_checkpoint_path: str | None,
    molscribe_checkpoint_repo: str,
    molscribe_checkpoint_file: str,
) -> dict[str, Any]:
    output_path = paper_dir / output_name
    if output_path.exists() and not overwrite and not dry_run:
        return {"status": "skipped", "paper": paper_dir.name, "reason": f"{output_name} already exists"}

    source_entries, total_images = collect_chemical_images(paper_dir)
    detector_entries, skipped_no_detections = build_detector_crops(
        source_entries,
        detector=detector,
        detector_model_path=detector_model_path,
        detector_confidence=detector_confidence,
        detector_iou=detector_iou,
        detector_padding=detector_padding,
    )
    if not detector_entries:
        payload = build_empty_payload(
            paper_dir,
            "No detector-backed molecule crops found in ChemicalBlock images",
            backend=backend,
            detector=detector,
            detector_model_path=detector_model_path,
            python_bin=python_bin,
            min_confidence=min_confidence,
            total_images=total_images,
            skipped_no_detections=skipped_no_detections,
        )
        if not dry_run:
            output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
        return {
            "status": "success",
            "paper": paper_dir.name,
            "pictures_found": 0,
            "pictures_with_smiles": 0,
            "output_path": str(output_path.relative_to(REPO_ROOT)),
        }

    with tempfile.TemporaryDirectory(prefix=f"{paper_dir.name[:24]}-{backend}-") as temp_root:
        temp_root_path = Path(temp_root)
        staged_paths, mapping = stage_images(detector_entries, temp_root_path)
        predictions = run_backend(
            backend=backend,
            image_paths=staged_paths,
            python_bin=python_bin,
            timeout_seconds=timeout_seconds,
            workdir=backend_workdir,
            molscribe_checkpoint_path=molscribe_checkpoint_path,
            molscribe_checkpoint_repo=molscribe_checkpoint_repo,
            molscribe_checkpoint_file=molscribe_checkpoint_file,
        )

    picture_data: list[dict[str, Any]] = []
    valid_smiles: list[str] = []
    invalid_smiles: list[str] = []
    seen_valid: set[str] = set()
    seen_invalid: set[str] = set()

    for staged_path, prediction in zip(staged_paths, predictions, strict=False):
        entry = mapping.get(staged_path)
        if entry is None:
            continue
        record = build_picture_record(
            entry,
            prediction if isinstance(prediction, dict) else {},
            backend=backend,
            min_confidence=min_confidence,
        )
        picture_data.append(record)
        if record["molecule_smiles"]:
            value = str(record["molecule_smiles"])
            if value not in seen_valid:
                valid_smiles.append(value)
                seen_valid.add(value)
        elif record["raw_smiles"]:
            value = str(record["raw_smiles"])
            if value not in seen_invalid:
                invalid_smiles.append(value)
                seen_invalid.add(value)

    relative_paper_dir = str(paper_dir.relative_to(REPO_ROOT))
    payload = {
        "paper": paper_dir.name,
        "paper_dir": relative_paper_dir,
        "generated_at": utc_now_iso(),
        "pictures_detected_total": total_images,
        "pictures_found": len(detector_entries),
        "pictures_processed": len(picture_data),
        "pictures_with_smiles": sum(1 for item in picture_data if item.get("molecule_smiles")),
        "pictures_skipped_no_detections": skipped_no_detections,
        "rdkit_available": RDKIT_AVAILABLE,
        "rdkit_status": "available" if RDKIT_AVAILABLE else "unavailable",
        "backend": backend,
        "backends_order": [backend],
        "backend_python_bin": python_bin,
        "detector": detector,
        "detector_model_path": detector_model_path,
        "min_confidence": min_confidence,
        "results": picture_data,
        "picture_data": picture_data,
        "smiles_extracted": valid_smiles,
        "valid_smiles": valid_smiles,
        "invalid_smiles": invalid_smiles,
    }

    if not dry_run:
        output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")

    return {
        "status": "success",
        "paper": paper_dir.name,
        "pictures_found": len(detector_entries),
        "pictures_with_smiles": payload["pictures_with_smiles"],
        "output_path": str(output_path.relative_to(REPO_ROOT)),
    }


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract SMILES from ChemicalBlock images and write smiles_extracted.json per paper."
    )
    parser.add_argument("--paper", help="Substring filter applied to extraction directory names.")
    parser.add_argument("--limit", type=int, help="Limit the number of extraction directories processed.")
    parser.add_argument(
        "--backend",
        choices=("molgrapher", "molscribe"),
        default=DEFAULT_BACKEND,
        help=f"OCSR backend to run. Default: {DEFAULT_BACKEND}",
    )
    parser.add_argument(
        "--detector",
        choices=("none", "moldetv2-general", "moldetv2-doc"),
        default=DEFAULT_DETECTOR,
        help=f"Detection stage applied before OCSR. Default: {DEFAULT_DETECTOR}",
    )
    parser.add_argument(
        "--python-bin",
        default=None,
        help="Python interpreter for the selected backend. Defaults to a backend-specific env var or the current interpreter.",
    )
    parser.add_argument(
        "--detector-model-path",
        default=os.environ.get("MOLDETV2_MODEL_PATH", "").strip() or None,
        help="Optional local MolDetv2 ONNX path. If omitted, the detector downloads from Hugging Face Hub.",
    )
    parser.add_argument(
        "--detector-model-repo",
        default=os.environ.get("MOLDETV2_MODEL_REPO", DEFAULT_MOLDETV2_REPO),
        help=f"Hugging Face repo for MolDetv2 weights. Default: {DEFAULT_MOLDETV2_REPO}",
    )
    parser.add_argument(
        "--detector-confidence",
        type=float,
        default=DEFAULT_DETECTOR_CONFIDENCE,
        help=f"Minimum MolDetv2 confidence retained before OCSR. Default: {DEFAULT_DETECTOR_CONFIDENCE}",
    )
    parser.add_argument(
        "--detector-iou",
        type=float,
        default=DEFAULT_DETECTOR_IOU,
        help=f"NMS IoU threshold for MolDetv2 boxes. Default: {DEFAULT_DETECTOR_IOU}",
    )
    parser.add_argument(
        "--detector-padding",
        type=float,
        default=DEFAULT_DETECTOR_PADDING,
        help=f"Relative padding applied around MolDetv2 crops before OCSR. Default: {DEFAULT_DETECTOR_PADDING}",
    )
    parser.add_argument(
        "--backend-workdir",
        default=os.environ.get("MOLGRAPHER_WORKDIR", "").strip() or None,
        help="Optional working directory for MolGrapher. Ignored for MolScribe.",
    )
    parser.add_argument(
        "--molscribe-checkpoint-path",
        default=os.environ.get("MOLSCRIBE_CKPT_PATH", "").strip() or None,
        help="Optional local MolScribe checkpoint path. If omitted, the backend downloads from Hugging Face Hub.",
    )
    parser.add_argument(
        "--molscribe-checkpoint-repo",
        default=os.environ.get("MOLSCRIBE_CKPT_REPO", DEFAULT_MOLSCRIBE_CKPT_REPO),
        help=f"Hugging Face repo for MolScribe checkpoints. Default: {DEFAULT_MOLSCRIBE_CKPT_REPO}",
    )
    parser.add_argument(
        "--molscribe-checkpoint-file",
        default=os.environ.get("MOLSCRIBE_CKPT_FILE", DEFAULT_MOLSCRIBE_CKPT_FILE),
        help=f"Checkpoint filename inside the Hugging Face repo. Default: {DEFAULT_MOLSCRIBE_CKPT_FILE}",
    )
    parser.add_argument("--output-name", default=DEFAULT_OUTPUT_NAME, help=f"Output filename. Default: {DEFAULT_OUTPUT_NAME}")
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TIMEOUT_SECONDS, help="Backend timeout per paper.")
    parser.add_argument(
        "--min-confidence",
        type=float,
        default=DEFAULT_MIN_CONFIDENCE,
        help=f"Drop otherwise valid backend predictions below this confidence. Default: {DEFAULT_MIN_CONFIDENCE}",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing smiles_extracted.json files.")
    parser.add_argument("--dry-run", action="store_true", help="Run extraction without writing output files.")
    return parser.parse_args(argv)


def iter_paper_dirs(*, paper_filter: str | None, limit: int | None) -> list[Path]:
    paper_dirs = sorted(path for path in EXTRACTIONS_DIR.iterdir() if path.is_dir() and (path / "normalized.json").exists())
    if paper_filter:
        lowered = paper_filter.casefold()
        paper_dirs = [path for path in paper_dirs if lowered in path.name.casefold()]
    if limit:
        paper_dirs = paper_dirs[:limit]
    return paper_dirs


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    python_bin = resolve_python_bin(args.python_bin, backend=args.backend)
    detector_model_path = resolve_detector_model_path(
        detector=args.detector,
        detector_model_path=args.detector_model_path,
        detector_model_repo=args.detector_model_repo,
    )
    paper_dirs = iter_paper_dirs(paper_filter=args.paper, limit=args.limit)

    if not paper_dirs:
        print("No extraction directories matched the current filters.", file=sys.stderr)
        return 1

    print(f"[extract] repo_root={REPO_ROOT}")
    print(f"[extract] backend={args.backend}")
    print(f"[extract] detector={args.detector}")
    print(f"[extract] python_bin={python_bin}")
    if detector_model_path:
        print(f"[extract] detector_model_path={detector_model_path}")
    print(f"[extract] papers={len(paper_dirs)}")

    failures = 0
    for index, paper_dir in enumerate(paper_dirs, start=1):
        print(f"[{index}/{len(paper_dirs)}] {paper_dir.name[:72]}...", end=" ")
        try:
            result = process_paper(
                paper_dir,
                backend=args.backend,
                detector=args.detector,
                python_bin=python_bin,
                output_name=args.output_name,
                timeout_seconds=args.timeout_seconds,
                min_confidence=args.min_confidence,
                overwrite=args.overwrite,
                dry_run=args.dry_run,
                backend_workdir=args.backend_workdir,
                detector_model_path=detector_model_path,
                detector_confidence=args.detector_confidence,
                detector_iou=args.detector_iou,
                detector_padding=args.detector_padding,
                molscribe_checkpoint_path=args.molscribe_checkpoint_path,
                molscribe_checkpoint_repo=args.molscribe_checkpoint_repo,
                molscribe_checkpoint_file=args.molscribe_checkpoint_file,
            )
        except Exception as exc:
            failures += 1
            print(f"ERROR {exc}")
            continue

        status = result.get("status")
        if status == "success":
            print(
                f"ok pictures={result.get('pictures_found', 0)} "
                f"with_smiles={result.get('pictures_with_smiles', 0)}"
            )
        else:
            print(f"skip {result.get('reason', 'unknown')}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
