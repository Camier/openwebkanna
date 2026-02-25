#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


CHECKPOINT_CANDIDATES = [
    "swin_base_char_aux_1m680k.pth",
    "swin_base_char_aux_1m.pth",
]


def _normalize(candidate: str) -> str:
    return "".join(candidate.strip().strip("`\"'").split())


def _to_smiles_list(prediction) -> list[str]:
    if prediction is None:
        return []
    if isinstance(prediction, str):
        value = _normalize(prediction)
        return [value] if value else []
    if isinstance(prediction, dict):
        values = []
        if prediction.get("smiles"):
            values.append(str(prediction.get("smiles")))
        if prediction.get("smiles_candidates"):
            raw = prediction.get("smiles_candidates")
            if isinstance(raw, list):
                values.extend([str(x) for x in raw])
            else:
                values.append(str(raw))
        out = []
        seen = set()
        for value in values:
            n = _normalize(value)
            if n and n not in seen:
                seen.add(n)
                out.append(n)
        return out
    if isinstance(prediction, (list, tuple)):
        out = []
        seen = set()
        for item in prediction:
            for value in _to_smiles_list(item):
                if value not in seen:
                    seen.add(value)
                    out.append(value)
        return out
    return []


def _resolve_checkpoint(checkpoint: str | None, repo: str) -> str:
    if checkpoint:
        path = Path(checkpoint).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"checkpoint_not_found: {path}")
        return str(path)

    from huggingface_hub import hf_hub_download

    for name in CHECKPOINT_CANDIDATES:
        try:
            return hf_hub_download(repo_id=repo, filename=name)
        except Exception:
            continue
    raise RuntimeError("checkpoint_download_failed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", type=str, required=True)
    parser.add_argument("--checkpoint", type=str)
    parser.add_argument("--repo", type=str, default="yujieq/MolScribe")
    parser.add_argument("--device", type=str, default="cpu")
    parser.add_argument("--return-confidence", action="store_true")
    args = parser.parse_args()

    import importlib
    import torch

    molscribe_module = importlib.import_module("molscribe")
    MolScribe = getattr(molscribe_module, "MolScribe")

    image_path = Path(args.image)
    if not image_path.exists():
        print(json.dumps({"status": "error", "error": f"missing_image: {image_path}"}))
        raise SystemExit(2)

    checkpoint = _resolve_checkpoint(args.checkpoint, args.repo)
    model = MolScribe(checkpoint, device=torch.device(args.device))
    prediction = model.predict_image_file(
        str(image_path),
        return_confidence=args.return_confidence,
    )
    smiles = _to_smiles_list(prediction)

    payload = {
        "status": "ok",
        "backend": "molscribe",
        "model": "MolScribe",
        "checkpoint": checkpoint,
        "smiles_candidates": smiles,
    }
    print(json.dumps(payload, ensure_ascii=False))


if __name__ == "__main__":
    main()
