#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


PROD_MAX_DIR = Path("/LAB/@thesis/openwebui/prod_max")
DEFAULT_MAX_NEW_TOKENS = 512
DEFAULT_VLM_MODEL = "HuggingFaceTB/SmolVLM-256M-Instruct"
DEFAULT_SMILES_BACKEND_ORDER = ["molscribe", "decimer", "imago", "vlm"]
DEFAULT_IMAGO_TIMEOUT_SECONDS = 30

MOLECULE_PROMPT = """You are a chemistry-aware image analyst. This image may contain a chemical structure.

Return ONLY a JSON object with this exact schema:
{
  \"description\": \"short description of what is shown\",
  \"smiles_candidates\": [\"<SMILES string>\", ...],
  \"notes\": \"optional notes\"
}

Rules:
- If no clear molecular structure is present, return [] for smiles_candidates.
- Keep each SMILES candidate as a plain string, no markdown.
- Prefer valid, complete SMILES and avoid extra explanations.
"""

SMILES_FALLBACK_PROMPT = """Extract only SMILES strings from this image.
Return strict JSON: {"smiles_candidates": ["..."]}
If none are visible, return: {"smiles_candidates": []}
"""

SMILES_PATTERNS = [
    re.compile(r"(?im)^\s*SMILES\s*:\s*([^\n\r]+)$"),
    re.compile(
        r"(?i)SMILES(?:\s*string)?\s*[:=]\s*([A-Za-z0-9@+\\-\\[\\]\\(\\)=#%\\./:]+)"
    ),
    re.compile(r"(?i)\b([A-Z][A-Za-z0-9@+\\-\\[\\]\\(\\)=#%\\./:]{3,})\b"),
]

PLACEHOLDER_VALUES = {
    "",
    "n/a",
    "na",
    "none",
    "not available",
    "unknown",
    "no smiles",
    "not found",
    "failed",
    "...",
    "…",
}


def _normalize_smiles(candidate: str) -> str:
    if not isinstance(candidate, str):
        return ""
    value = candidate.strip().strip("`\"'")
    value = re.sub(r"\s+", "", value)
    return value


def _dedupe(values: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for value in values:
        normalized = _normalize_smiles(value)
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        out.append(normalized)
    return out


def _is_placeholder(candidate: str) -> bool:
    return candidate.strip().lower() in PLACEHOLDER_VALUES


def _as_list_of_strings(value) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    return [str(value).strip()]


def _parse_json_payload(raw: str) -> Dict[str, object]:
    block = _extract_json_block(raw)
    if not block:
        return {}
    try:
        parsed = json.loads(block)
        if isinstance(parsed, dict):
            return parsed
    except Exception:
        return {}
    return {}


def _extract_json_block(raw: str) -> Optional[str]:
    fenced = re.search(r"```\s*json\s*(\{.*?\})\s*```", raw, flags=re.S | re.I)
    if fenced:
        return fenced.group(1).strip()

    start = raw.find("{")
    while start != -1:
        depth = 0
        in_string = False
        escape = False
        for index in range(start, len(raw)):
            ch = raw[index]
            if in_string:
                if escape:
                    escape = False
                elif ch == "\\":
                    escape = True
                elif ch == '"':
                    in_string = False
                continue

            if ch == '"':
                in_string = True
                continue
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    return raw[start : index + 1].strip()
        start = raw.find("{", start + 1)

    return None


def _parse_candidate_from_json(payload: Dict[str, object]) -> Tuple[str, List[str]]:
    description = ""
    description_raw = payload.get("description")
    if isinstance(description_raw, str):
        description = description_raw.strip()

    smiles_values: List[str] = []
    for key in ("smiles_candidates", "smiles", "candidate_smiles", "results"):
        if key in payload:
            smiles_values.extend(_as_list_of_strings(payload.get(key)))

    if not smiles_values and isinstance(payload.get("result"), dict):
        nested = payload.get("result")
        if isinstance(nested, dict):
            for key in ("smiles_candidates", "smiles"):
                smiles_values.extend(_as_list_of_strings(nested.get(key)))

    return description, smiles_values


def _extract_payload(raw: str) -> Tuple[str, List[str], Dict[str, object]]:
    payload = _parse_json_payload(raw)
    if payload:
        description, smiles = _parse_candidate_from_json(payload)
        smiles = _dedupe(
            [
                _normalize_smiles(s)
                for s in smiles
                if not _is_placeholder(_normalize_smiles(s))
            ]
        )
        if smiles:
            return description, smiles, payload

    smiles = _iter_smiles(raw)
    return "", smiles, {}


def _iter_smiles(raw: str) -> List[str]:
    candidates: List[str] = []
    for pattern in SMILES_PATTERNS:
        candidates.extend(pattern.findall(raw))
    return _dedupe(candidates)


def _parse_env_bool(value: Optional[str], default: bool = False) -> bool:
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on", "y"}


def _parse_env_int(value: Optional[str], default: int) -> int:
    if value is None:
        return default
    try:
        parsed = int(value)
        return parsed if parsed >= 0 else default
    except ValueError:
        return default


def _parse_backend_order(raw: str) -> List[str]:
    entries = [entry.strip().lower() for entry in raw.split(",") if entry.strip()]
    normalized: List[str] = []
    for entry in entries:
        if (
            entry in {"molscribe", "decimer", "imago", "vlm"}
            and entry not in normalized
        ):
            normalized.append(entry)
    return normalized


class SmilesValidator:
    def __init__(self) -> None:
        self._chem: Any = None
        self.available = False
        self.reason = ""

        disable_rdkit = os.getenv("SMILES_DISABLE_RDKIT", "").strip().lower()
        if disable_rdkit in {"1", "true", "yes", "on"}:
            self.reason = "disabled_by_env"
            return

        try:
            import numpy as np

            major = int(str(np.__version__).split(".", maxsplit=1)[0])
            if major >= 2:
                self.reason = "numpy_major_gte_2"
                return
        except Exception:
            self.reason = "numpy_version_unknown"

        try:
            from rdkit import Chem
            from rdkit import RDLogger

            if _parse_env_bool(os.getenv("SMILES_SILENCE_RDKIT_ERRORS"), True):
                RDLogger.DisableLog("rdApp.error")

            self._chem = Chem
            self.available = True
            self.reason = "available"
        except Exception:
            self._chem = None
            self.available = False
            self.reason = "rdkit_import_failed"

    def normalize_and_validate(self, smiles: str) -> Optional[str]:
        candidate = _normalize_smiles(smiles)
        if not candidate or _is_placeholder(candidate):
            return None
        if all(ch == "." for ch in candidate):
            return None

        if self._chem is None:
            return candidate

        try:
            mol_from_smiles = getattr(self._chem, "MolFromSmiles", None)
            mol_to_smiles = getattr(self._chem, "MolToSmiles", None)
            if not callable(mol_from_smiles) or not callable(mol_to_smiles):
                return None

            molecule = mol_from_smiles(candidate)
            if molecule is None:
                return None
            canonical = mol_to_smiles(molecule)
            if not isinstance(canonical, str):
                return None
            return canonical
        except Exception:
            return None


def find_images_for_paper(paper_dir: Path) -> List[Path]:
    images: List[Path] = []
    for assets_dir in paper_dir.glob("_assets/*/images"):
        images.extend(assets_dir.glob("*.jpg"))
        images.extend(assets_dir.glob("*.jpeg"))
        images.extend(assets_dir.glob("*.png"))
        images.extend(assets_dir.glob("*.JPG"))
        images.extend(assets_dir.glob("*.JPEG"))
        images.extend(assets_dir.glob("*.PNG"))
    return sorted(images, key=lambda p: p.name)


def _build_image_metadata_map(paper_dir: Path) -> Dict[str, Dict[str, object]]:
    normalized_file = paper_dir / "normalized.json"
    if not normalized_file.exists():
        return {}

    mapping: Dict[str, Dict[str, object]] = {}
    try:
        with normalized_file.open() as handle:
            payload = json.load(handle)
        blocks = payload.get("raw", {}).get("chunks", {}).get("blocks", [])
        for block in blocks:
            images = block.get("images", {}) or {}
            if not isinstance(images, dict) or not images:
                continue

            page = block.get("page", 0)
            bbox = block.get("bbox")
            for image_name in images.keys():
                mapping[image_name] = {
                    "page": page,
                    "bbox": bbox,
                }
    except Exception:
        return {}

    return mapping


def _classify_smiles(
    candidate_lists: Sequence[str],
    validator: SmilesValidator,
) -> Tuple[List[str], List[str], List[str]]:
    cleaned: List[str] = []
    valid_smiles: List[str] = []
    invalid_smiles: List[str] = []
    for value in candidate_lists:
        normalized = _normalize_smiles(value)
        if not normalized or _is_placeholder(normalized):
            continue
        cleaned.append(normalized)

    for candidate in _dedupe(cleaned):
        canonical = validator.normalize_and_validate(candidate)
        if canonical:
            valid_smiles.append(canonical)
        else:
            invalid_smiles.append(candidate)

    return cleaned, _dedupe(valid_smiles), _dedupe(invalid_smiles)


def _prediction_to_candidates(
    raw_prediction: Any,
) -> Tuple[List[str], Dict[str, object]]:
    metadata: Dict[str, object] = {}

    if raw_prediction is None:
        return [], metadata

    if isinstance(raw_prediction, str):
        text = raw_prediction.strip()
        if not text:
            return [], metadata
        payload = _parse_json_payload(text)
        if payload:
            metadata["parsed_payload"] = payload
            _desc, values = _parse_candidate_from_json(payload)
            if _desc:
                metadata["description"] = _desc
            return _dedupe([_normalize_smiles(v) for v in values]), metadata

        return _iter_smiles(text), {"raw_output": text}

    if isinstance(raw_prediction, dict):
        metadata["parsed_payload"] = raw_prediction
        _, smiles = _parse_candidate_from_json(raw_prediction)
        return _dedupe([_normalize_smiles(v) for v in smiles]), metadata

    if isinstance(raw_prediction, (list, tuple)):
        candidates: List[str] = []
        if (
            len(raw_prediction) == 2
            and isinstance(raw_prediction[0], str)
            and not isinstance(raw_prediction[1], (str, bytes))
        ):
            candidates.extend(_as_list_of_strings(raw_prediction[0]))
            metadata["return_confidence"] = raw_prediction[1]
            return _dedupe([_normalize_smiles(v) for v in candidates]), metadata
        for item in raw_prediction:
            if isinstance(item, str):
                candidates.append(item)
            elif isinstance(item, dict):
                _, smiles = _parse_candidate_from_json(item)
                candidates.extend(smiles)
        return _dedupe([_normalize_smiles(v) for v in candidates]), metadata

    return [], metadata


def _resolve_device(raw_device: str) -> str:
    if not raw_device or raw_device == "auto":
        try:
            import torch

            return "cuda" if torch.cuda.is_available() else "cpu"
        except Exception:
            return "cpu"
    return raw_device


def _initialize_molscribe_state(
    args: argparse.Namespace, state: Dict[str, Any]
) -> None:
    if state.get("initialized"):
        return
    state["initialized"] = True
    state["available"] = False
    state["backend"] = "molscribe"
    state["attempted"] = False

    checkpoint = args.molscribe_checkpoint or os.getenv("SMILES_MOLSCRIBE_CHECKPOINT")
    if not checkpoint:
        state["reason"] = "missing_checkpoint"
        return

    model_path = Path(checkpoint).expanduser()
    if not model_path.exists():
        state["reason"] = "checkpoint_not_found"
        return

    repo = args.molscribe_repo or os.getenv("SMILES_MOLSCRIBE_REPO")
    if repo:
        repo_path = Path(repo).expanduser()
        if not repo_path.exists():
            state["reason"] = "repo_not_found"
            return
        if str(repo_path) not in sys.path:
            sys.path.insert(0, str(repo_path))

    try:
        import importlib
        import torch

        molscribe_module = importlib.import_module("molscribe")
        MolScribe = getattr(molscribe_module, "MolScribe")
    except Exception as error:
        state["reason"] = f"import_failed: {error}"
        return

    device = _resolve_device(args.smiles_device)
    num_workers = max(1, int(getattr(args, "molscribe_num_workers", 1)))
    try:
        start = time.perf_counter()
        model = MolScribe(
            str(model_path),
            device=torch.device(device),
            num_workers=num_workers,
        )
        state["load_time_s"] = round(time.perf_counter() - start, 6)
        state["available"] = True
        state["reason"] = "ready"
        state["model"] = model
        state["num_workers"] = num_workers
        state["device"] = device
    except Exception as error:
        state["reason"] = f"init_failed: {error}"


def _initialize_decimer_state(state: Dict[str, Any]) -> None:
    if state.get("initialized"):
        return
    state["initialized"] = True
    state["available"] = False
    state["backend"] = "decimer"

    try:
        import DECIMER
    except Exception as error:
        state["reason"] = f"import_failed: {error}"
        return

    state["available"] = True
    state["reason"] = "ready"
    state["module"] = DECIMER


def _initialize_imago_state(args: argparse.Namespace, state: Dict[str, Any]) -> None:
    if state.get("initialized"):
        return
    state["initialized"] = True
    state["available"] = False
    state["backend"] = "imago"

    command = (
        args.imago_command or os.getenv("SMILES_IMAGO_COMMAND") or "imago"
    ).strip()
    if not command:
        state["reason"] = "missing_command"
    else:
        binary = shlex.split(command)[0]
        if shutil.which(binary) is not None:
            state["available"] = True
            state["reason"] = "ready"
            state["command"] = command
        else:
            state["reason"] = f"command_missing: {binary}"

    if not state["available"]:
        try:
            import importlib

            importlib.import_module("imago")
        except Exception as error:
            state.setdefault("module_import_reason", f"import_failed: {error}")
        else:
            state["available"] = True
            state["reason"] = "module_available"
            state["module_path"] = "imago"


def _run_molscribe(
    image_path: Path, args: argparse.Namespace, state: Dict[str, Any]
) -> Tuple[List[str], Dict[str, object]]:
    if args.molscribe_use_subprocess:
        command = shlex.split(args.molscribe_command)
        cmd = command + ["--image", str(image_path)]
        if args.molscribe_checkpoint:
            cmd += ["--checkpoint", str(args.molscribe_checkpoint)]
        if args.molscribe_repo:
            cmd += ["--repo", str(args.molscribe_repo)]

        try:
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
                timeout=args.molscribe_timeout,
            )
        except subprocess.TimeoutExpired:
            return [], {
                "status": "timeout",
                "command": args.molscribe_command,
                "error": f"timeout_{args.molscribe_timeout}s",
            }
        except Exception as error:
            return [], {
                "status": "command_error",
                "command": args.molscribe_command,
                "error": str(error),
            }

        combined_output = "\n".join(
            [
                part.strip()
                for part in ((process.stdout or ""), (process.stderr or ""))
                if part and part.strip()
            ]
        )
        payload = _parse_json_payload(combined_output)
        if process.returncode != 0:
            return [], {
                "status": "command_failed",
                "command": args.molscribe_command,
                "return_code": process.returncode,
                "stderr": (process.stderr or "").strip(),
                "parsed_payload": payload,
            }

        smiles_values: List[str] = []
        if payload:
            smiles_values.extend(_as_list_of_strings(payload.get("smiles_candidates")))
            if not smiles_values:
                smiles_values.extend(_as_list_of_strings(payload.get("smiles")))

        metadata = {
            "command": args.molscribe_command,
            "return_code": process.returncode,
            "parsed_payload": payload,
        }
        if "description" in payload and isinstance(payload["description"], str):
            metadata["description"] = payload["description"]

        return _dedupe([_normalize_smiles(value) for value in smiles_values]), metadata

    _initialize_molscribe_state(args, state)
    if not state.get("available"):
        return [], {"status": "backend_unavailable", "reason": state.get("reason")}

    model = state.get("model")
    if model is None:
        return [], {"status": "backend_unavailable", "reason": state.get("reason")}

    try:
        prediction = model.predict_image_file(
            str(image_path),
            return_confidence=args.molscribe_return_confidence,
        )
        candidates, metadata = _prediction_to_candidates(prediction)
        if isinstance(prediction, dict):
            metadata["confidence"] = prediction.get("confidence")
        return candidates, metadata
    except Exception as error:
        return [], {"status": "error", "error": str(error)}


def _run_decimer(
    image_path: Path, args: argparse.Namespace, state: Dict[str, Any]
) -> Tuple[List[str], Dict[str, object]]:
    if not state.get("initialized"):
        _initialize_decimer_state(state)
    if not state.get("available"):
        return [], {"status": "backend_unavailable", "reason": state.get("reason")}

    module = state.get("module")
    if module is None:
        return [], {"status": "backend_unavailable", "reason": state.get("reason")}

    try:
        prediction = module.predict_SMILES(
            str(image_path),
            confidence=args.decimer_return_confidence,
            hand_drawn=args.decimer_hand_drawn,
        )
        candidates, metadata = _prediction_to_candidates(prediction)
        metadata["hand_drawn"] = args.decimer_hand_drawn
        return candidates, metadata
    except Exception as error:
        return [], {"status": "error", "error": str(error)}


def _run_imago(
    image_path: Path, args: argparse.Namespace, state: Dict[str, Any]
) -> Tuple[List[str], Dict[str, object]]:
    if not state.get("initialized"):
        _initialize_imago_state(args, state)
    if not state.get("available"):
        return [], {"status": "backend_unavailable", "reason": state.get("reason")}

    command = state.get("command")
    if command:
        args_list = shlex.split(command) + [str(image_path)]
        try:
            process = subprocess.run(
                args_list,
                capture_output=True,
                text=True,
                check=False,
                timeout=args.imago_timeout,
            )
            if process.returncode != 0 and not process.stdout:
                return [], {
                    "status": "command_failed",
                    "return_code": process.returncode,
                    "stderr": (process.stderr or "").strip(),
                }

            output = "\n".join(
                [
                    part.strip()
                    for part in ((process.stdout or ""), (process.stderr or ""))
                    if part and part.strip()
                ]
            )
            candidates, metadata = _prediction_to_candidates(output)
            metadata.update(
                {
                    "command": command,
                    "return_code": process.returncode,
                }
            )
            return candidates, metadata
        except subprocess.TimeoutExpired:
            return [], {
                "status": "timeout",
                "command": command,
                "error": f"timeout_{args.imago_timeout}s",
            }
        except Exception as error:
            return [], {
                "status": "command_error",
                "command": command,
                "error": str(error),
            }

    return [], {
        "status": "module_fallback_not_implemented",
        "reason": state.get("reason"),
    }


def _initialize_vlm_state(args: argparse.Namespace, state: Dict[str, Any]) -> None:
    if state.get("initialized"):
        return

    state["initialized"] = True
    state["available"] = False
    state["backend"] = "vlm"
    model_name = args.vlm_model or DEFAULT_VLM_MODEL

    try:
        import torch
        from transformers import AutoModelForVision2Seq, AutoProcessor
    except Exception as error:
        state["reason"] = f"import_failed: {error}"
        return

    try:
        device = _resolve_device(args.smiles_device)
        dtype = torch.float16 if str(device).startswith("cuda") else torch.float32
        processor = AutoProcessor.from_pretrained(model_name)
        model = AutoModelForVision2Seq.from_pretrained(model_name, torch_dtype=dtype)
        model = model.to(device)
        model.eval()
        state["available"] = True
        state["reason"] = "ready"
        state["device"] = device
        state["model_name"] = model_name
        state["processor"] = processor
        state["model"] = model
    except Exception as error:
        state["reason"] = f"init_failed: {error}"


def _generate_vlm_output(
    image_path: Path,
    prompt: str,
    max_new_tokens: int,
    state: Dict[str, Any],
) -> str:
    import torch
    from PIL import Image

    processor = state["processor"]
    model = state["model"]
    device = state["device"]

    image = Image.open(image_path).convert("RGB")
    messages = [
        {
            "role": "user",
            "content": [
                {"type": "image"},
                {"type": "text", "text": prompt},
            ],
        }
    ]
    rendered = processor.apply_chat_template(
        messages,
        tokenize=False,
        add_generation_prompt=True,
    )
    inputs = processor(
        text=[rendered], images=[image], return_tensors="pt", padding=True
    )
    inputs = {name: tensor.to(device) for name, tensor in inputs.items()}

    with torch.no_grad():
        generated_ids = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=False,
        )

    decoded = processor.batch_decode(
        [out[len(inp) :] for inp, out in zip(inputs["input_ids"], generated_ids)],
        skip_special_tokens=True,
        clean_up_tokenization_spaces=True,
    )[0]
    return decoded


def _run_vlm(
    image_path: Path, args: argparse.Namespace, state: Dict[str, Any]
) -> Tuple[List[str], Dict[str, object]]:
    _initialize_vlm_state(args, state)
    if not state.get("available"):
        return [], {"status": "backend_unavailable", "reason": state.get("reason")}

    try:
        primary = _generate_vlm_output(
            image_path=image_path,
            prompt=MOLECULE_PROMPT,
            max_new_tokens=args.max_new_tokens,
            state=state,
        )
        description, candidates, payload = _extract_payload(primary)
        metadata: Dict[str, object] = {
            "model": state.get("model_name"),
            "device": state.get("device"),
            "description": description,
            "raw_output": primary,
        }
        if payload:
            metadata["parsed_payload"] = payload

        if candidates:
            return candidates, metadata

        fallback = _generate_vlm_output(
            image_path=image_path,
            prompt=SMILES_FALLBACK_PROMPT,
            max_new_tokens=args.max_new_tokens,
            state=state,
        )
        description_fb, candidates_fb, payload_fb = _extract_payload(fallback)
        metadata["fallback_output"] = fallback
        if description_fb:
            metadata["description"] = description_fb
        if payload_fb:
            metadata["fallback_payload"] = payload_fb
        return candidates_fb, metadata
    except Exception as error:
        return [], {"status": "error", "error": str(error)}


def _attempt_backends(
    image_path: Path,
    args: argparse.Namespace,
    validator: SmilesValidator,
    backend_order: List[str],
    backend_states: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    attempts: List[Dict[str, object]] = []
    final_status = "failed"
    selected_smiles: List[str] = []
    selected_valid: List[str] = []
    selected_invalid: List[str] = []
    backend_used: Optional[str] = None
    backends_tried: List[str] = []

    for backend in backend_order:
        start = time.perf_counter()
        backends_tried.append(backend)

        if backend == "molscribe":
            raw_candidates, attempt_metadata = _run_molscribe(
                image_path,
                args,
                backend_states.setdefault("molscribe", {}),
            )
        elif backend == "decimer":
            raw_candidates, attempt_metadata = _run_decimer(
                image_path,
                args,
                backend_states.setdefault("decimer", {}),
            )
        elif backend == "imago":
            raw_candidates, attempt_metadata = _run_imago(
                image_path,
                args,
                backend_states.setdefault("imago", {}),
            )
        elif backend == "vlm":
            raw_candidates, attempt_metadata = _run_vlm(
                image_path,
                args,
                backend_states.setdefault("vlm", {}),
            )
        else:
            raw_candidates, attempt_metadata = [], {"status": "unknown_backend"}

        duration_s = round(time.perf_counter() - start, 6)

        candidates, valid_smiles, invalid_smiles = _classify_smiles(
            raw_candidates, validator
        )

        if candidates:
            if valid_smiles:
                status = "valid"
                final_status = "valid"
                backend_used = backend
                selected_smiles = candidates
                selected_valid = valid_smiles
                selected_invalid = invalid_smiles
            else:
                status = "invalid_candidates"
                if final_status == "failed":
                    final_status = "invalid_candidates"
                    selected_smiles = candidates
                    selected_valid = []
                    selected_invalid = invalid_smiles
        else:
            if attempt_metadata.get("status") in {
                "command_failed",
                "error",
                "backend_unavailable",
                "timeout",
                "command_error",
            }:
                status = "error"
                if final_status == "failed":
                    final_status = "error"
            else:
                status = "no_candidates"
                if final_status == "failed":
                    final_status = "no_candidates"

        attempt = {
            "backend": backend,
            "status": status,
            "duration_s": duration_s,
            "attempted_smiles": candidates,
            "valid_smiles": valid_smiles,
            "invalid_smiles": invalid_smiles,
            "metadata": attempt_metadata,
        }
        attempts.append(attempt)

        if final_status == "valid":
            break

    payload = {
        "attempts": attempts,
        "backends_tried": backends_tried,
        "backend": backend_used,
        "status": final_status,
        "description": "",
        "smiles_extracted": selected_smiles,
        "valid_smiles": selected_valid,
        "invalid_smiles": selected_invalid,
    }

    if payload["smiles_extracted"] and final_status != "valid":
        # If we got invalid candidates but no valid smiles, keep the first non-empty list.
        payload["smiles_extracted"] = _dedupe(payload["smiles_extracted"])

    return payload


def process_image(
    image_path: Path,
    picture_meta: Dict[str, object],
    args: argparse.Namespace,
    validator: SmilesValidator,
    backend_states: Dict[str, Dict[str, Any]],
) -> Dict:
    from PIL import Image  # kept local to avoid import overhead when not needed

    result: Dict[str, object] = {
        "image": str(image_path),
        "description": "",
        "smiles_extracted": [],
        "valid_smiles": [],
        "invalid_smiles": [],
        "status": "failed",
        "page": picture_meta.get("page", 0),
        "bbox": picture_meta.get("bbox", ""),
        "response": "",
        "success": False,
        "backend": None,
        "backends_tried": [],
        "attempts": [],
    }

    # Validate the file is loadable before backend execution.
    try:
        with Image.open(image_path) as handle:
            handle.verify()
    except Exception as error:
        result.update(
            {
                "status": "error",
                "response": "",
                "success": False,
                "error": f"invalid_image: {error}",
            }
        )
        return result

    attempt_payload = _attempt_backends(
        image_path,
        args,
        validator,
        _parse_backend_order(
            ",".join(args.backend_order)
            if isinstance(args.backend_order, list)
            else args.backend_order
        ),
        backend_states,
    )

    result.update(
        {
            "description": attempt_payload.get("description", ""),
            "smiles_extracted": attempt_payload.get("smiles_extracted", []),
            "valid_smiles": attempt_payload.get("valid_smiles", []),
            "invalid_smiles": attempt_payload.get("invalid_smiles", []),
            "status": attempt_payload.get("status", "failed"),
            "response": "",
            "success": bool(attempt_payload.get("status") == "valid"),
            "backend": attempt_payload.get("backend"),
            "backends_tried": attempt_payload.get("backends_tried", []),
            "attempts": attempt_payload.get("attempts", []),
        }
    )

    return result


def process_paper(
    paper_dir: Path,
    args: argparse.Namespace,
    validator: SmilesValidator,
    dry_run: bool,
    write_molecules_jsonl: bool,
    backend_states: Dict[str, Dict[str, Any]],
) -> Dict:
    paper_name = paper_dir.name
    images = find_images_for_paper(paper_dir)
    if args.max_images_per_paper and args.max_images_per_paper > 0:
        images = images[: args.max_images_per_paper]
    print(f"  Found {len(images)} images")

    image_meta = _build_image_metadata_map(paper_dir)
    results: List[Dict] = []

    for index, image in enumerate(images, start=1):
        print(f"  Processing image {index}/{len(images)}: {image.name[:30]}...")
        results.append(
            process_image(
                image,
                image_meta.get(image.name, {}),
                args,
                validator,
                backend_states,
            )
        )

    all_smiles: List[str] = []
    all_valid_smiles: List[str] = []
    all_invalid_smiles: List[str] = []
    picture_data: List[Dict] = []

    backend_usage: List[str] = []

    for result in results:
        all_smiles.extend(result.get("smiles_extracted", []))
        all_valid_smiles.extend(result.get("valid_smiles", []))
        all_invalid_smiles.extend(result.get("invalid_smiles", []))

        backend = result.get("backend")
        if backend and backend not in backend_usage:
            backend_usage.append(backend)

        picture_data.append(
            {
                "image": result.get("image", ""),
                "page": result.get("page", 0),
                "bbox": result.get("bbox", ""),
                "description": result.get("description", ""),
                "status": result.get("status", ""),
                "molecule_smiles": (
                    result.get("valid_smiles", [None])[0]
                    if result.get("valid_smiles")
                    else None
                ),
                "raw_candidates": result.get("smiles_extracted", []),
                "success": result.get("success", False),
                "error": result.get("error", ""),
                "backend": result.get("backend"),
                "backends_tried": result.get("backends_tried", []),
                "attempts": result.get("attempts", []),
            }
        )

    backend_used: Optional[str] = None
    if len(backend_usage) == 1:
        backend_used = backend_usage[0]
    elif len(backend_usage) > 1:
        backend_used = "mixed"

    output = {
        "paper": paper_name,
        "paper_dir": str(paper_dir),
        "pictures_found": len(images),
        "pictures_processed": len(images),
        "pictures_with_smiles": len([r for r in results if r.get("valid_smiles")]),
        "rdkit_available": validator.available,
        "rdkit_status": validator.reason,
        "backend": backend_used,
        "backends_order": args.backend_order,
        "results": results,
        "picture_data": picture_data,
        "smiles_extracted": _dedupe(all_smiles),
        "valid_smiles": _dedupe(all_valid_smiles),
        "invalid_smiles": _dedupe(all_invalid_smiles),
    }

    if not dry_run:
        output_file = paper_dir / "smiles_extracted.json"
        with output_file.open("w") as handle:
            json.dump(output, handle, indent=2)
        print(f"  Saved to: {output_file}")

        if write_molecules_jsonl:
            molecules_file = paper_dir / "molecules.jsonl"
            with molecules_file.open("w") as handle:
                for item in picture_data:
                    valid_list = []
                    for entry in item.get("raw_candidates", []):
                        canonical = validator.normalize_and_validate(str(entry))
                        if canonical:
                            valid_list.append(canonical)

                    for smiles in _dedupe(valid_list):
                        record = {
                            "paper": paper_name,
                            "paper_dir": str(paper_dir),
                            "image": item.get("image", ""),
                            "page": item.get("page", 0),
                            "bbox": item.get("bbox", ""),
                            "status": item.get("status", ""),
                            "description": item.get("description", ""),
                            "smiles": smiles,
                        }
                        handle.write(json.dumps(record, ensure_ascii=False) + "\n")
            print(f"  Saved to: {molecules_file}")

    return output


def _resolve_papers(args: argparse.Namespace, prod_root: Path) -> List[Path]:
    if not prod_root.exists() or not prod_root.is_dir():
        return []

    all_papers = [item for item in sorted(prod_root.iterdir()) if item.is_dir()]
    if args.all:
        return all_papers

    if args.paper:
        return [item for item in all_papers if args.paper.lower() in item.name.lower()]

    return []


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract SMILES from molecule images")

    parser.add_argument(
        "--paper", type=str, help="Paper directory name (partial match)"
    )
    parser.add_argument("--all", action="store_true", help="Process all papers")

    parser.add_argument(
        "--backend-order",
        type=str,
        default=os.getenv(
            "SMILES_BACKEND_ORDER",
            ",".join(DEFAULT_SMILES_BACKEND_ORDER),
        ),
        help="Comma-separated backend order: molscribe,decimer,imago,vlm",
    )

    parser.add_argument(
        "--smiles-device", type=str, default=os.getenv("SMILES_DEVICE", "auto")
    )

    parser.add_argument(
        "--molscribe-repo", type=str, default=os.getenv("SMILES_MOLSCRIBE_REPO")
    )
    parser.add_argument(
        "--molscribe-checkpoint",
        type=str,
        default=os.getenv("SMILES_MOLSCRIBE_CHECKPOINT"),
        help="MolScribe checkpoint path (required when using molscribe)",
    )
    parser.add_argument(
        "--molscribe-num-workers",
        type=int,
        default=_parse_env_int(os.getenv("SMILES_MOLSCRIBE_NUM_WORKERS"), 1),
        help="Number of workers passed to MolScribe",
    )
    parser.add_argument(
        "--molscribe-return-confidence",
        action="store_true",
        default=_parse_env_bool(os.getenv("MOLSCRIBE_RETURN_CONFIDENCE"), False),
        help="Request confidence field from MolScribe",
    )
    parser.add_argument(
        "--molscribe-use-subprocess",
        action=argparse.BooleanOptionalAction,
        default=_parse_env_bool(os.getenv("SMILES_MOLSCRIBE_USE_SUBPROCESS"), True),
        help="Run MolScribe using external command",
    )
    parser.add_argument(
        "--molscribe-command",
        type=str,
        default=os.getenv(
            "SMILES_MOLSCRIBE_COMMAND",
            "micromamba run -n molscribe-legacy python scripts/molscribe_predict.py",
        ),
        help="Command used when MolScribe subprocess mode is enabled",
    )
    parser.add_argument(
        "--molscribe-timeout",
        type=int,
        default=_parse_env_int(os.getenv("SMILES_MOLSCRIBE_TIMEOUT"), 240),
        help="Timeout in seconds for MolScribe subprocess invocation",
    )

    parser.add_argument(
        "--decimer-hand-drawn",
        action="store_true",
        default=_parse_env_bool(os.getenv("DECIMER_HAND_DRAWN"), False),
        help="Use DECIMER hand-drawn variant",
    )
    parser.add_argument(
        "--decimer-return-confidence",
        action="store_true",
        default=_parse_env_bool(os.getenv("DECIMER_RETURN_CONFIDENCE"), False),
        help="Request confidence output from DECIMER",
    )

    parser.add_argument(
        "--imago-command",
        type=str,
        default=os.getenv("SMILES_IMAGO_COMMAND", "imago"),
        help="CLI command used to call imago backend",
    )
    parser.add_argument(
        "--imago-timeout",
        type=int,
        default=_parse_env_int(
            os.getenv("SMILES_IMAGO_TIMEOUT"), DEFAULT_IMAGO_TIMEOUT_SECONDS
        ),
        help="Timeout in seconds for imago CLI invocation",
    )

    parser.add_argument(
        "--vlm-model",
        type=str,
        default=os.getenv("SMILES_VLM_MODEL", DEFAULT_VLM_MODEL),
    )
    parser.add_argument("--max-new-tokens", type=int, default=DEFAULT_MAX_NEW_TOKENS)

    parser.add_argument(
        "--dry-run", action="store_true", help="Do not write output files"
    )
    parser.add_argument("--limit", type=int, help="Limit number of papers")
    parser.add_argument(
        "--max-images-per-paper",
        type=int,
        default=_parse_env_int(os.getenv("SMILES_MAX_IMAGES_PER_PAPER"), 0),
        help="Process at most N images per paper (0 = all)",
    )
    parser.add_argument(
        "--data-root",
        type=str,
        default=str(PROD_MAX_DIR),
        help="Root folder for extracted papers",
    )
    parser.add_argument(
        "--write-molecules-jsonl",
        action="store_true",
        help="Write validated molecule rows to molecules.jsonl",
    )
    args = parser.parse_args()

    args.backend_order = _parse_backend_order(args.backend_order)
    if not args.backend_order:
        args.backend_order = DEFAULT_SMILES_BACKEND_ORDER.copy()

    prod_root = Path(args.data_root)

    if not args.paper and not args.all:
        parser.error("Either --paper or --all is required")

    paper_dirs = _resolve_papers(args, prod_root)
    if args.limit:
        paper_dirs = paper_dirs[: args.limit]

    if not paper_dirs:
        print("No matching paper directories found")
        return

    print(f"Found {len(paper_dirs)} paper(s)")

    validator = SmilesValidator()
    results: List[Dict] = []
    backend_states: Dict[str, Dict[str, Any]] = {}

    for index, paper_dir in enumerate(paper_dirs, start=1):
        images = find_images_for_paper(paper_dir)
        print(f"\n[{index}/{len(paper_dirs)}] {paper_dir.name[:50]}")

        if not images:
            output = {
                "paper": paper_dir.name,
                "paper_dir": str(paper_dir),
                "pictures_found": 0,
                "pictures_processed": 0,
                "pictures_with_smiles": 0,
                "rdkit_available": validator.available,
                "rdkit_status": validator.reason,
                "backend": None,
                "backends_order": args.backend_order,
                "results": [],
                "picture_data": [],
                "smiles_extracted": [],
                "valid_smiles": [],
                "invalid_smiles": [],
                "images_found_warning": "No images found in _assets/*/images",
            }
            if not args.dry_run:
                output_file = paper_dir / "smiles_extracted.json"
                with output_file.open("w") as handle:
                    json.dump(output, handle, indent=2)
                print(f"  No images found, saved summary to: {output_file}")
                if args.write_molecules_jsonl:
                    molecules_file = paper_dir / "molecules.jsonl"
                    molecules_file.write_text("")
                    print(
                        f"  No images found, saved empty molecules file to: {molecules_file}"
                    )
            results.append(output)
            continue

        results.append(
            process_paper(
                paper_dir,
                args,
                validator,
                dry_run=args.dry_run,
                write_molecules_jsonl=args.write_molecules_jsonl,
                backend_states=backend_states,
            )
        )

    total_images = sum(result.get("pictures_processed", 0) for result in results)
    total_candidates = sum(
        len(result.get("smiles_extracted", [])) for result in results
    )
    total_valid = sum(len(result.get("valid_smiles", [])) for result in results)
    total_invalid = sum(len(result.get("invalid_smiles", [])) for result in results)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print(f"  Papers processed: {len(results)}")
    print(f"  Images processed: {total_images}")
    print(f"  SMILES candidates: {total_candidates}")
    print(f"  Valid SMILES: {total_valid}")
    print(f"  Invalid SMILES: {total_invalid}")
    print(f"  RDKit available: {validator.available}")


if __name__ == "__main__":
    main()
