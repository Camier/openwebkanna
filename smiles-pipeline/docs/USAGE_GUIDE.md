# SMILES Pipeline - Usage Guide

Operational procedures for running the SMILES extraction pipeline.

---

## Quick Reference

### Common Commands

| Task | Command |
|------|---------|
| Install environment | `micromamba env create -f envs/environment-smiles-extraction.yml` |
| Extract from images (GPU) | `python scripts/extract_smiles_pipeline.py --input-dir data/extractions --gpu` |
| Extract from images (CPU) | `python scripts/extract_smiles_pipeline.py --input-dir data/extractions --no-gpu` |
| Test on 5 images | `python scripts/extract_smiles_pipeline.py --limit 5 --no-gpu` |
| Import to KB | `./import-smiles-to-kb.sh --input-dir data/validated` |
| Generate QA report | `python scripts/generate_qa_report.py --input-dir data/validated` |
| View extraction stats | `cat data/validated/extraction_stats.json \| jq` |

### Typical Workflow

```bash
# 1. Process images
cd /LAB/@thesis/openwebui
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --output-dir smiles-pipeline/data/validated \
  --gpu

# 2. Review results
cat smiles-pipeline/data/validated/extraction_stats.json | jq

# 3. Convert to Markdown
python smiles-pipeline/scripts/jsonl_to_markdown.py \
  --input smiles-pipeline/data/validated/molecules.jsonl \
  --output smiles-pipeline/data/validated/markdown

# 4. Import to OpenWebUI
./smiles-pipeline/import-smiles-to-kb.sh \
  --input-dir smiles-pipeline/data/validated

# 5. Generate QA report
python smiles-pipeline/scripts/generate_qa_report.py \
  --input-dir smiles-pipeline/data/validated \
  --output-file smiles-pipeline/data/summary/qa_report.json
```

**Total time** (1000 images, GPU): ~6 minutes

---

## Step-by-Step Procedures

### Procedure 1: Extract SMILES from Images

**Input**: Chemical structure images (PNG/JPG/GIF/BMP)
**Output**: `molecules.jsonl` with extracted and validated molecules

```bash
# Basic usage
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir <path_to_images> \
  --output-dir <path_to_output>

# Examples

# Extract from all images in data/extractions (GPU)
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --output-dir smiles-pipeline/data/validated \
  --gpu

# Extract from subset (first 100 images)
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --limit 100 \
  --gpu

# CPU-only mode (for debugging)
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --no-gpu
```

**Options**:

| Flag | Description | Default |
|------|-------------|---------|
| `--input-dir` | Directory containing images | Required |
| `--output-dir` | Output directory | `smiles-pipeline/data/validated` |
| `--gpu` | Use GPU acceleration | True |
| `--no-gpu` | Force CPU mode | False |
| `--limit` | Process only N images | All |
| `--backend-order` | Backend priority | `molscribe,decimer` |
| `--confidence-threshold` | Min confidence | 0.3 |
| `--min-validation-score` | Min validation score | 70 |

### Procedure 2: Convert to Markdown for OpenWebUI

**Input**: `molecules.jsonl`
**Output**: Markdown files in `markdown/` directory

```bash
python smiles-pipeline/scripts/jsonl_to_markdown.py \
  --input smiles-pipeline/data/validated/molecules.jsonl \
  --output smiles-pipeline/data/validated/markdown

# Options
--split-by-paper    # Create one MD file per paper (default)
--single-file       # Combine all molecules into one MD file
--include-invalid   # Include molecules that failed validation
```

### Procedure 3: Import to OpenWebUI Knowledge Base

**Input**: Markdown files
**Output**: OpenWebUI KB with searchable molecules

```bash
./smiles-pipeline/import-smiles-to-kb.sh \
  --input-dir smiles-pipeline/data/validated \
  --kb-name "Sceletium Alkaloids"

# Auto-detects OpenWebUI API key from .env
# Progress: [1/100] Importing mesembrine.md → Done
```

**Prerequisites**:
- OpenWebUI running (port 3000)
- API key in `.env` (`OPENWEBUI_API_KEY`)
- Admin access to OpenWebUI

### Procedure 4: Generate QA Report

**Input**: `molecules.jsonl` and `extraction_stats.json`
**Output**: `qa_report.json` with analytics

```bash
python smiles-pipeline/scripts/generate_qa_report.py \
  --input-dir smiles-pipeline/data/validated \
  --output-file smiles-pipeline/data/summary/qa_report.json
```

**Report includes**:
- Total molecules extracted
- Validation pass rates
- Backend performance breakdown
- Property distributions (MW, LogP, TPSA)
- Gold standard matches
- Confidence score distribution

### Procedure 5: Troubleshoot Extraction

**When extraction fails or produces low-quality results**:

```bash
# 1. Check environment
micromamba run -n smiles-extraction pip list | grep -E 'torch|tensorflow|molscribe|decimer'

# 2. Test OCSR backends individually
python /tmp/test_molscribe.py
python /tmp/test_decimer.py

# 3. View extraction errors
cat smiles-pipeline/data/validated/extraction_stats.json | jq '.errors[]'

# 4. Check GPU availability
micromamba run -n smiles-extraction python -c "import torch; print('CUDA:', torch.cuda.is_available())"

# 5. Run with verbose logging
python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --verbose \
  --log-file extraction.log
```

---

## Best Practices

### Image Preparation

**Ideal image characteristics**:
- Resolution: ≥300 DPI
- Format: PNG (lossless) or high-quality JPG
- Content: Single chemical structure per image
- Background: Clean, high contrast (black on white)

**Preprocessing** (if needed):
```bash
# Convert PDFs to images (300 DPI)
pdftoppm -png -r 300 paper.pdf output

# Extract figures from PDFs
python scripts/extract_figures_from_pdf.py paper.pdf --output images/
```

### Batch Processing

**For large datasets (100+ images)**:

```bash
# Process in batches of 100
for batch in {0..9}; do
  micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
    --input-dir data/extractions \
    --output-dir smiles-pipeline/data/validated/batch_${batch} \
    --limit 100 \
    --skip $((batch * 100))
done

# Merge results
cat smiles-pipeline/data/validated/batch_*/molecules.jsonl > smiles-pipeline/data/validated/merged_molecules.jsonl
```

### GPU Optimization

**Maximize GPU utilization**:

```yaml
# config/backends.yaml
batch_sizes:
  molscribe: 64    # Increase for GPU with more VRAM
  decimer: 32

gpu_enabled: true
```

**Monitor GPU usage**:
```bash
watch -n 1 nvidia-smi
```

### Validation Tuning

**Adjust validation strictness**:

```yaml
# config/validation_rules.yaml

# Relax MW limits for broader compound classes
level_2_chemical:
  mw_min: 50       # Default: 100
  mw_max: 1500     # Default: 1000

# Allow more fragments for complex molecules
  max_fragments: 10  # Default: 5

# Reduce domain specificity
level_3_domain:
  require_nitrogen: false   # Accept non-alkaloids
  scaffold_similarity_min: 0.5  # Default: 0.7
```

---

## Example Sessions

### Session 1: Full Pipeline Run

**Goal**: Extract all molecules from 129 papers (~1000 images)

```bash
# Step 1: Extract (15 minutes with GPU)
time micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --output-dir smiles-pipeline/data/validated \
  --gpu
# Processed 1000 images in 14:32 minutes (average 0.87s/image)

# Step 2: Review stats
cat smiles-pipeline/data/validated/extraction_stats.json | jq
# {
#   "images_processed": 1000,
#   "extractions_successful": 987,
#   "syntax_valid": 975,
#   "chemically_valid": 892,
#   "domain_valid": 456,
#   "high_confidence": 423
# }

# Step 3: Generate report
python smiles-pipeline/scripts/generate_qa_report.py \
  --input-dir smiles-pipeline/data/validated

# Step 4: Import to OpenWebUI
./smiles-pipeline/import-smiles-to-kb.sh \
  --input-dir smiles-pipeline/data/validated \
  --kb-name "Sceletium Compounds v2.0"
```

**Results**:
- 987 successful extractions (98.7% success rate)
- 456 domain-validated alkaloids
- 23 exact gold standard matches
- 156 similar matches (Tanimoto ≥0.7)

### Session 2: High-Confidence Extraction

**Goal**: Extract only highest-confidence molecules (≥90 score)

```bash
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --output-dir smiles-pipeline/data/validated/hq \
  --confidence-threshold 0.8 \
  --min-validation-score 90 \
  --gpu
```

**Results**:
- Fewer molecules extracted (stricter filtering)
- Higher precision (>95% correct)
- Fewer false positives

### Session 3: CPU-Only Processing

**Goal**: Process on CPU-only machine (no GPU)

```bash
# Slower but functional
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --no-gpu \
  --limit 50  # Process subset first to test
```

**Expected speed**: ~5 seconds/image → ~4 minutes for 50 images

---

## Monitoring and Logging

### Progress Monitoring

**Built-in progress bar**:
```
Processing images:  60%|██████    | 300/500 [02:34<01:42,  1.96it/s]
```

**Log files**:
```
smiles-pipeline/logs/extraction_{timestamp}.log
```

### Metrics to Watch

| Metric | Target | Alert If |
|--------|--------|----------|
| Extraction success rate | ≥95% | <90% |
| Syntax validation rate | ≥95% | <90% |
| Chemical validation rate | ≥85% | <75% |
| Domain validation rate | ≥50% | <30% |
| Average confidence | ≥0.7 | <0.5 |
| Processing speed (GPU) | ≥1 image/s | <0.5 image/s |

---

## Data Management

### File Organization

```
smiles-pipeline/data/validated/
├── molecules.jsonl              # Main output
├── extraction_stats.json        # Metrics
├── markdown/                    # OpenWebUI-ready
│   ├── paper_001.md
│   ├── paper_002.md
│   └── ...
├── low_confidence.jsonl         # Below threshold (manual review)
└── errors.log                   # Failed extractions
```

### Backup and archival

```bash
# Backup validated molecules
tar -czvf smiles-validated-$(date +%Y%m%d).tar.gz \
  smiles-pipeline/data/validated/

# Archive to cold storage
mv smiles-validated-*.tar.gz /data/archive/smiles/
```

---

## Troubleshooting Guide

### Problem: No molecules in output

**Symptoms**:
```json
{
  "images_processed": 100,
  "extractions_successful": 0
}
```

**Diagnosis**:
```bash
# Check if images exist
ls -la data/extractions/*/_assets/*/images/*.jpg | wc -l

# Test extraction on single image
python /tmp/test_molscribe.py

# Check logs
cat smiles-pipeline/logs/extraction_*.log | tail -50
```

**Solutions**:
1. Verify image paths are correct
2. Check GPU availability (if using `--gpu`)
3. Try `--no-gpu` mode
4. Increase batch size if OOM

### Problem: All validations failing

**Symptoms**:
```json
{
  "syntax_valid": 0,
  "chemically_valid": 0,
  "domain_valid": 0
}
```

**Diagnosis**:
```bash
# Check RDKit installation
micromamba run -n smiles-extraction python -c "from rdkit import Chem; print('OK')"

# Test validator
python -c "
from validators.indigo_validator import IndigoValidator
v = IndigoValidator()
print(v.validate('CN1CC2(C1)C3=CC=CC=C3'))
"
```

**Solutions**:
1. Reinstall RDKit: `pip install rdkit==2024.03.5`
2. Check for API changes in validators
3. Relax validation rules in config

### Problem: Slow extraction speed

**Symptoms**: >10 seconds per image

**Diagnosis**:
```bash
# Check GPU usage
nvidia-smi  # Should show Python process using GPU

# Monitor during extraction
watch -n 1 'nvidia-smi | grep "Default Process"'
```

**Solutions**:
1. Ensure GPU mode: `--gpu` flag
2. Reduce batch size if OOM
3. Close other GPU processes
4. Check thermal throttling

### Problem: Model download fails

**Symptoms**:
```
Error downloading checkpoint: Connection timeout
```

**Solutions**:
```bash
# Manual download (MolScribe)
wget -P ~/.cache/molscribe \
  https://huggingface.co/yujieq/MolScribe/resolve/main/swin_base_char_aux_1m.pth

# Manual download (DECIMER)
python -c "
from pystow import ensure
ensure('DECIMER-V2', 'https://zenodo.org/record/.../models.zip')
"
```

---

## Advanced Usage

### Custom Backend Chain

**Use only DECIMER for specific image types**:
```bash
python scripts/extract_smiles_pipeline.py \
  --input-dir data/hand_drawn \
  --backend-order decimer,molscribe \
  --gpu
```

### Confidence-Based Routing

**Route high-confidence to KB, low-confidence to manual review**:
```python
import json

with open('molecules.jsonl') as f:
    molecules = [json.loads(line) for line in f]

high_conf = [m for m in molecules if m['confidence_score'] >= 70]
low_conf = [m for m in molecules if m['confidence_score'] < 70]

# Save separately
with open('molecules_high_confidence.jsonl', 'w') as f:
    for m in high_conf:
        f.write(json.dumps(m) + '\n')

with open('molecules_manual_review.jsonl', 'w') as f:
    for m in low_conf:
        f.write(json.dumps(m) + '\n')
```

### Custom Validation Rules

**Create domain-specific validation**:
```python
from validators.domain_validator import DomainValidator

class CustomValidator(DomainValidator):
    def validate_custom(self, smiles):
        # Add your custom validation logic
        result = self.validate(smiles)
        result['custom_check'] = self.my_custom_function(smiles)
        return result

validator = CustomValidator()
```

---

## FAQ

**Q: Can I extract from multi-molecule images?**
A: Currently, the pipeline expects single molecules per image. Multi-molecule images may have reduced accuracy. For best results, split images first.

**Q: How do I add a new OCSR backend?**
A: Create a new extractor class in `src/extractors/`, add to `config/backends.yaml`, and update the backend chain in `scripts/extract_smiles_pipeline.py`.

**Q: Can I run this on Windows?**
A: The pipeline is Linux-tested. Windows may work with WSL2, but GPU support is limited.

**Q: What if my compounds aren't alkaloids?**
A: Adjust `config/validation_rules.yaml` to relax domain validation or disable Level 3 entirely.

**Q: How do I integrate with my own database?**
A: The `molecules.jsonl` output is standard JSON Lines format. Parse and import to any database that supports JSON.

---

**Last updated**: 2026-02-25
**Version**: 2.0.0
