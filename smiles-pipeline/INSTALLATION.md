# SMILES Pipeline - Installation Guide

Installation and setup guide for the SMILES extraction pipeline with MolScribe and DECIMER OCSR backends.

## Prerequisites

- **Linux** (tested on Ubuntu/Debian-based systems)
- **micromamba** or **conda** for environment management
- **NVIDIA GPU** (optional but recommended for ~10x speedup)
- **15GB+ free disk space** (model checkpoints are ~1GB each)

## Quick Start

### 1. Create Environment

```bash
cd /LAB/@thesis/openwebui/smiles-pipeline/envs

# Using micromamba (recommended)
micromamba env create -f environment-smiles-extraction.yml --yes

# Using conda
conda env create -f environment-smiles-extraction.yml

# Takes ~10-15 minutes depending on internet speed
```

### 2. Verify Installation

```bash
# Activate environment
micromamba activate smiles-extraction

# Test all imports
python -c "
import torch
import tensorflow
import molscribe
import DECIMER
import rdkit
import indigo
print('✓ All packages installed successfully')
print(f'  PyTorch: {torch.__version__}')
print(f'  TensorFlow: {tensorflow.__version__}')
print(f'  RDKit: {rdkit.__version__}')
print(f'  DECIMER: {DECIMER.__version__}')
"
```

### 3. Test OCSR Extraction

```bash
# Run on 5 sample images (CPU mode)
cd /LAB/@thesis/openwebui
./smiles-pipeline/run_extraction.sh --limit 5 --no-gpu

# With GPU acceleration
./smiles-pipeline/run_extraction.sh --limit 5 --gpu
```

## Configuration

### GPU vs CPU

**Default behavior**: Uses GPU if available (via `--gpu` flag)

- **GPU**: ~150-300ms per image
- **CPU**: ~1-2 seconds per image

```bash
# Force CPU mode (useful for debugging)
./smiles-pipeline/run_extraction.sh --no-gpu

# Use GPU (default)
./smiles-pipeline/run_extraction.sh --gpu
```

### Backend Selection

Edit `smiles-pipeline/config/backends.yaml` to change backend priority:

```yaml
backend_order:
  - molscribe  # Primary (93% F1)
  - decimer    # Fallback (91% F1)
  # - molvec   # CPU-only fallback (not implemented)
```

## Model Checkpoints

### MolScribe

- **Size**: ~1.06 GB
- **Location**: `~/.cache/molscribe/swin_base_char_aux_1m.pth`
- **Download**: Automatic on first use from HuggingFace
- **Manual download** (if needed):
  ```bash
  wget -P ~/.cache/molscribe \
    https://huggingface.co/yujieq/MolScribe/resolve/main/swin_base_char_aux_1m.pth
  ```

### DECIMER

- **Size**: Downloaded automatically via pystow
- **Location**: `~/.pystow/DECIMER-V2/`
- **Automatic download** on first import

## Troubleshooting

### Common Issues

#### 1. "torch version mismatch" error
```bash
# Reinstall PyTorch for your CUDA version
micromamba run -n smiles-extraction pip install --upgrade torch torchaudio
```

#### 2. "Import 'molscribe' could not be resolved"
This is expected until the environment is created. LSP errors are normal during development.

**Solution**: Activate the environment:
```bash
micromamba activate smiles-extraction
```

#### 3. Out of memory (OOM) on GPU
Reduce batch size in `config/backends.yaml`:
```yaml
molscribe:
  batch_size: 8  # Default: 32
```

#### 4. Slow extraction speed
- **Check GPU usage**: `nvidia-smi`
- **Ensure GPU mode**: Use `--gpu` flag
- **Reduce batch size** if OOM

#### 5. RDKit API errors
Current environment uses RDKit 2024.03.5. If you see API errors:
```bash
# The pipeline uses CalcNumStereoCenters (new API)
# Not CalcNumChiralCenters (deprecated)
```

### Environment Diagnostics

```bash
# Check installed packages
micromamba run -n smiles-extraction pip list | grep -E 'torch|tensorflow|molscribe|decimer|rdkit'

# Test GPU availability
micromamba run -n smiles-extraction python -c "import torch; print('CUDA:', torch.cuda.is_available())"

# Check model cache
ls -lh ~/.cache/molscribe/
ls -lh ~/.pystow/DECIMER-V2/ 2>/dev/null || echo "DECIMER models will download on first use"
```

## Updating

### Update Dependencies

```bash
cd /LAB/@thesis/openwebui/smiles-pipeline/envs

# Update environment file
# Edit environment-smiles-extraction.yml

# Re-create environment
micromamba env remove -n smiles-extraction --yes
micromamba env create -f environment-smiles-extraction.yml --yes
```

### Check for New MolScribe Version

```bash
# Check GitHub for updates
curl -s https://api.github.com/repos/thomas0809/MolScribe/commits?sha=main | jq '.[0].sha, .[0].commit.message'

# MolScribe is installed from GitHub main branch (auto-updates on reinstall)
```

### Check for New DECIMER Version

```bash
pip index versions decimer
```

## Performance Benchmarks

### GPU (RTX 3060/3080/4090)

| Component | Time per Image |
|-----------|---------------|
| MolScribe extraction | 150-300ms |
| DECIMER extraction | 200-500ms |
| 3-layer validation | <10ms |
| Property enrichment | <20ms |
| **Total** | **~200-350ms** |

### CPU (8-core)

| Component | Time per Image |
|-----------|---------------|
| MolScribe extraction | 1-2s |
| DECIMER extraction | 2-5s |
| 3-layer validation | <10ms |
| Property enrichment | <20ms |
| **Total** | **~1-3s** |

## Next Steps

After installation:

1. **Run extraction**:
   ```bash
   ./smiles-pipeline/run_extraction.sh --gpu
   ```

2. **Review results**:
   ```bash
   head -n 5 smiles-pipeline/data/validated/molecules.jsonl | jq
   ```

3. **Import to OpenWebUI**:
   ```bash
   ./smiles-pipeline/import-smiles-to-kb.sh --input-dir smiles-pipeline/data/validated
   ```

4. **Generate QA report**:
   ```bash
   micromamba run -n smiles-extraction python smiles-pipeline/scripts/generate_qa_report.py \
     --input-dir smiles-pipeline/data/validated \
     --output-file smiles-pipeline/data/summary/qa_report.json
   ```

## Support

- **OCSR Accuracy Issues**: Check image quality (blurry, low resolution)
- **Installation Problems**: See troubleshooting section
- **Pipeline Errors**: Check `smiles-pipeline/data/validated/extraction_stats.json`

## License

- **Pipeline code**: MIT
- **MolScribe**: MIT
- **DECIMER**: CC BY 4.0 (requires attribution in publications)

---

*Last updated: 2026-02-25*
*Status: Production-ready ✅*
