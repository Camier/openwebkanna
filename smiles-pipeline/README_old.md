# SMILES Extraction Pipeline

**Optical Chemical Structure Recognition (OCSR) pipeline** for extracting molecular structures from ethnopharmacological research papers.

## Overview

This pipeline extracts SMILES (Simplified Molecular Input Line Entry System) strings from images of chemical structures in academic PDFs, validates them, and imports them into OpenWebUI Knowledge Bases for RAG retrieval.

## Quick Start

### Prerequisites

```bash
# Install Conda environment
conda env create -f smiles-pipeline/envs/environment-smiles-extraction.yml

# Or use mise (if available)
mise use conda@latest
```

### Full Pipeline (Recommended)

```bash
# Process all papers in data/pdfs
mise run smiles:pipeline

# Or run directly
./smiles-pipeline/scripts/extract_smiles_from_images.py \
  --input-dir data/extractions \
  --output-dir smiles-pipeline/data \
  --backend-order molscribe,decimer,imago,vlm
```

### Publish to Knowledge Base

```bash
# Import high-confidence molecules
./smiles-pipeline/import-smiles-to-kb.sh

# Convert to markdown format (better readability)
mise run smiles:publish-kb-markdown
```

## Pipeline Architecture

```
PDF Images → OCSR Backends → Validation → Filtering → Knowledge Base
     ↓              ↓              ↓           ↓            ↓
  Figure    MolScribe      RDKit     Confidence   OpenWebUI
  blocks    DECIMER        Canonical   Score      KB Upload
            Imago          SMILES
            VLM
```

## Components

### Scripts (`smiles-pipeline/scripts/`)

| Script | Purpose |
|--------|---------|
| `extract_smiles_from_images.py` | Main extraction pipeline (1,345 lines) |
| `filter_high_confidence_smiles.py` | Filter by confidence score |
| `smiles_qa_summary.py` | Quality assurance metrics |
| `molscribe_predict.py` | MolScribe subprocess wrapper |
| `validate_smiles_outputs.py` | RDKit validation |
| `jsonl_to_markdown_smiles.py` | Convert to markdown for KB |

### Shell Scripts

| Script | Purpose |
|--------|---------|
| `smiles-pipeline/import-smiles-to-kb.sh` | Upload to OpenWebUI KB |
| `smiles-pipeline/dedupe-smiles-kb.sh` | Remove duplicate KBs |

### Environments (`smiles-pipeline/envs/`)

| File | Purpose |
|------|---------|
| `environment-smiles-extraction.yml` | Conda environment |
| `requirements-smiles-extraction.txt` | Pip requirements |
| `environment-molscribe-legacy.yml` | Legacy MolScribe env |

### Documentation (`smiles-pipeline/docs/`)

| File | Purpose |
|------|---------|
| `SMILES_PIPELINE.md` | Complete pipeline guide |
| `SMILES_OCSR_REFERENCE_MATRIX.md` | Backend comparison matrix |

## Backend Configuration

### OCSR Backend Order

```bash
# Default order (best to fastest)
SMILES_BACKEND_ORDER=molscribe,decimer,imago,vlm

# Custom order
SMILES_BACKEND_ORDER=decimer,molscribe

# Single backend
SMILES_BACKEND_ORDER=molscribe
```

### Backend Characteristics

| Backend | Accuracy | Speed | Model Size | Hardware |
|---------|----------|-------|------------|----------|
| **MolScribe** | High | Medium | ~100MB | CPU/GPU |
| **DECIMER** | High | Slow | ~500MB | CPU/GPU |
| **Imago** | Medium | Fast | ~50MB | CPU only |
| **VLM** | Variable | Slow | API-based | Cloud |

## Data Flow

### Input Structure

```
data/extractions/{paper_id}/
└── _assets/{doc_id}/images/
    ├── {hash1}.jpg
    ├── {hash2}.jpg
    └── ...
```

### Output Structure

```
smiles-pipeline/data/
├── {paper_id}/
│   ├── smiles_extracted.json           # Raw extraction results
│   ├── molecules.jsonl                 # Validated molecules
│   └── molecules_high_confidence.jsonl # Filtered for import
└── summary/
    └── smiles_qa_summary.json          # Aggregate metrics
```

### Molecule JSONL Format

```json
{
  "doc_id": "sha256_hash",
  "image_name": "aspirin_structure.jpg",
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "confidence": 0.95,
  "backend": "molscribe",
  "validation_status": "valid",
  "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "page": 3,
  "block_id": "/page/2/Figure/1"
}
```

## Quality Control

### Confidence Filtering

```bash
# Default high-confidence threshold
CONFIDENCE_THRESHOLD=0.8

# Custom threshold
CONFIDENCE_THRESHOLD=0.9 ./smiles-pipeline/scripts/filter_high_confidence_smiles.py
```

### Validation Checks

1. **RDKit canonicalization** - Invalid SMILES rejected
2. **Placeholder rejection** - n/a, unknown, failed, etc.
3. **Confidence score** - Must meet threshold
4. **Structure matching** - Cross-backend consistency

### QA Metrics

```bash
# Generate QA summary
./smiles-pipeline/scripts/smiles_qa_summary.py \
  --input-dir smiles-pipeline/data \
  --output-file smiles-pipeline/data/summary/qa_report.json

# Metrics include:
# - Extraction success rate
# - Validation pass rate
# - Confidence distribution
# - Backend comparison
```

## OpenWebUI Integration

### Import to Knowledge Base

```bash
# Import all high-confidence molecules
./smiles-pipeline/import-smiles-to-kb.sh \
  --input-dir smiles-pipeline/data \
  --kb-name "Molecular Structures" \
  --format jsonl

# Or markdown format (better readability)
./smiles-pipeline/import-smiles-to-kb.sh --format markdown
```

### Knowledge Base Structure

```
OpenWebUI KB: "Molecular Structures"
├── [Paper 1] Aspirin synthesis
│   ├── CC(=O)Oc1ccccc1C(=O)O  # Aspirin
│   └── C1=CC=C(C=C1)C(=O)O    # Salicylic acid
├── [Paper 2] Sceletium alkaloids
│   ├── CN1CCC23C=CC(=CC2CC1)C3  # Mesembrine
│   └── ...
└── ...
```

### Deduplication

```bash
# Remove duplicate KBs (keep one)
./smiles-pipeline/dedupe-smiles-kb.sh --dry-run
./smiles-pipeline/dedupe-smiles-kb.sh
```

## Advanced Usage

### Batch Processing

```bash
# Process specific papers only
./smiles-pipeline/scripts/extract_smiles_from_images.py \
  --paper-ids "2024_Kang_Glucose,2023_Gericke_Sceletium"

# Parallel processing (experimental)
SMILES_PARALLEL=4 ./smiles-pipeline/scripts/extract_smiles_from_images.py
```

### Custom Validation

```python
from validate_smiles_outputs import SmilesValidator

validator = SmilesValidator()
result = validator.validate("CC(=O)Oc1ccccc1C(=O)O")

if result.is_valid:
    print(f"Valid SMILES: {result.canonical_smiles}")
    print(f"Molecular weight: {result.molecular_weight}")
    print(f"logP: {result.logp}")
```

### A/B Testing Backends

```bash
# Compare backend performance
for backend in molscribe decimer imago; do
    SMILES_BACKEND_ORDER=$backend \
    ./smiles-pipeline/scripts/extract_smiles_from_images.py \
      --paper-ids "test_set_10_papers"
done

# Analyze results
./smiles-pipeline/scripts/smiles_qa_summary.py --compare-backends
```

## Troubleshooting

### Common Issues

**Issue: MolScribe not found**
```bash
# Check Conda env
conda activate smiles-extraction
which molscribe
```

**Issue: Low extraction rate**
```bash
# Try different backend order
SMILES_BACKEND_ORDER=decimer,molscribe,imago
```

**Issue: Invalid SMILES in output**
```bash
# Re-validate with stricter thresholds
./smiles-pipeline/scripts/validate_smiles_outputs.py \
  --strict-mode \
  --min-confidence 0.9
```

### Performance Tuning

**GPU Acceleration:**
```bash
# MolScribe with GPU
CUDA_VISIBLE_DEVICES=0 ./smiles-pipeline/scripts/extract_smiles_from_images.py
```

**Memory Optimization:**
```bash
# Process one paper at a time
SMILES_BATCH_SIZE=1 ./smiles-pipeline/scripts/extract_smiles_from_images.py
```

## Performance Benchmarks

| Metric | Value |
|--------|-------|
| MolScribe throughput | ~2 images/sec (GPU) |
| DECIMER throughput | ~0.5 images/sec |
| Validation rate | ~95% of extractions |
| High-confidence rate | ~70-80% of valid |
| Typical paper | 10-30 structures |

## Research Context

This pipeline supports ethnopharmacological research on:
- **Sceletium tortuosum** (Kanna) alkaloids
- Structure-activity relationships
- Traditional medicinal plant chemistry
- Molecular mechanism elucidation

## Maintenance

### Update Dependencies

```bash
# Update Conda env
conda env update -f smiles-pipeline/envs/environment-smiles-extraction.yml

# Update pip packages
pip install -r smiles-pipeline/envs/requirements-smiles-extraction.txt --upgrade
```

### Backup Data

```bash
# Backup extracted molecules
tar -czf smiles-backup-$(date +%Y%m%d).tar.gz \
  smiles-pipeline/data/molecules*.jsonl
```

## References

- **MolScribe:** https://github.com/HongjianNi/MolScribe
- **DECIMER:** https://github.com/Kohulan/DECIMER-Image_Transformer
- **Imago:** https://github.com/epam/imago
- **RDKit:** https://www.rdkit.org/

---

**For integration with RAG system:** See [../../OPERATIONS.md](../../OPERATIONS.md)
