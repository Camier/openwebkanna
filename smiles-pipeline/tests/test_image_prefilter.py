from pathlib import Path
import json
import subprocess
import sys

import pytest
from PIL import Image, ImageDraw

from extract_smiles_pipeline import ExtractionPipeline


def _build_pipeline_for_prefilter() -> ExtractionPipeline:
    return ExtractionPipeline(
        backend_order=["molscribe"],
        use_gpu=False,
        backends_config={},
        validation_rules={
            "image_prefilter": {
                "enabled": True,
                "min_width": 200,
                "min_height": 150,
                "max_aspect_ratio": 5.0,
                "dark_threshold": 48,
                "min_dark_pixel_ratio": 0.003,
                "max_dark_pixel_ratio": 0.35,
                "edge_threshold": 32,
                "min_edge_pixel_ratio": 0.005,
                "max_edge_pixel_ratio": 0.45,
                "min_entropy": 0.2,
                "max_entropy": 7.8,
            }
        },
    )


def test_prefilter_rejects_too_small_image(tmp_path: Path):
    pipeline = _build_pipeline_for_prefilter()
    image_path = tmp_path / "tiny.png"
    Image.new("RGB", (120, 80), color="white").save(image_path)

    result = pipeline._evaluate_image_prefilter(image_path)
    assert result["is_candidate"] is False
    assert result["reason_code"] == "image_too_small"


def test_prefilter_accepts_structure_like_image(tmp_path: Path):
    pipeline = _build_pipeline_for_prefilter()
    image_path = tmp_path / "structure_like.png"
    image = Image.new("L", (300, 220), color=255)
    draw = ImageDraw.Draw(image)
    for y in range(20, 200, 10):
        draw.line((30, y, 260, y), fill=0, width=1)
    for x in range(30, 260, 12):
        draw.line((x, 20, x, 200), fill=0, width=1)
    image.convert("RGB").save(image_path)

    result = pipeline._evaluate_image_prefilter(image_path)
    assert result["is_candidate"] is True
    assert result["reason_code"] == "prefilter_pass"
    assert "metrics" in result


@pytest.mark.integration
def test_pipeline_preflight_command_runs():
    repo_root = Path(__file__).resolve().parents[2]
    script_path = repo_root / "smiles-pipeline" / "scripts" / "extract_smiles_pipeline.py"
    config_dir = repo_root / "smiles-pipeline" / "config"
    completed = subprocess.run(
        [
            sys.executable,
            str(script_path),
            "--preflight",
            "--config-dir",
            str(config_dir),
            "--no-gpu",
        ],
        capture_output=True,
        text=True,
    )

    # Preflight exits non-zero on warn/fail by design, so parse report from stdout.
    assert completed.stdout
    lines = completed.stdout.splitlines()
    start_idx = None
    for idx, line in enumerate(lines):
        if line.strip() == "{":
            start_idx = idx
            break
    assert start_idx is not None
    report = json.loads("\n".join(lines[start_idx:]))
    assert "status" in report
    assert "checks" in report
    assert any(c["name"] == "dependency:indigo" for c in report["checks"])
