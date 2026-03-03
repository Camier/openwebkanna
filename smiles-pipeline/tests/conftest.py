from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SMILES_PIPELINE_ROOT = REPO_ROOT / "smiles-pipeline"
SRC_DIR = SMILES_PIPELINE_ROOT / "src"
SCRIPTS_DIR = SMILES_PIPELINE_ROOT / "scripts"

for candidate in (str(SRC_DIR), str(SCRIPTS_DIR)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)
