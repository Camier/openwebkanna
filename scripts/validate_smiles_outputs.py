#!/usr/bin/env python3
"""Compatibility wrapper for legacy root path.

Forwards execution to smiles-pipeline/archive/legacy-v1/validate_smiles_outputs.py.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def main() -> None:
    target = (
        Path(__file__).resolve().parents[1]
        / "smiles-pipeline"
        / "archive"
        / "legacy-v1"
        / "validate_smiles_outputs.py"
    )
    if not target.exists():
        raise SystemExit(f"Missing target script: {target}")

    os.execv(sys.executable, [sys.executable, str(target), *sys.argv[1:]])


if __name__ == "__main__":
    main()
