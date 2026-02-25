# SMILES Pipeline v2.0 - Architecture Documentation

## System Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         SMILES Pipeline v2.0                            │
├─────────────────────────────────────────────────────────────────────────┤
│  INPUT LAYER                                                            │
│  ───────────                                                            │
│  Images (PNG/JPG/GIF/BMP) → Preprocessing → Normalization              │
│                                                                       │
│  EXTRACTION LAYER                                                       │
│  ────────────────                                                       │
│  ┌─────────────┐  ┌─────────────┐                                      │
│  │ MolScribe   │  │ DECIMER 2.7 │  ← Parallel OCSR backends           │
│  │ (PyTorch)   │  │ (TensorFlow)│                                      │
│  │ 93% F1      │  │ 91% F1      │  ← Fallback chain                   │
│  └─────────────┘  └─────────────┘                                      │
│         │                │                                              │
│         └────────┬───────┘                                             │
│                  │                                                      │
│  VALIDATION LAYER                                                       │
│  ─────────────────                                                      │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐                    │
│  │ Level 1     │→ │ Level 2     │→ │ Level 3     │                    │
│  │ Syntax      │  │ Chemical    │  │ Domain      │                    │
│  │ (Indigo)    │  │ (RDKit)     │  │ (Alkaloid)  │                    │
│  └─────────────┘  └─────────────┘  └─────────────┘                    │
│                                                                       │
│  ENRICHMENT LAYER                                                       │
│  ─────────────────                                                      │
│  Properties (MW, LogP, TPSA) + Fingerprints (ECFP4, MACCS)           │
│                                                                       │
│  OUTPUT LAYER                                                           │
│  ────────────                                                           │
│  JSONL → Markdown → OpenWebUI KB                                      │
└─────────────────────────────────────────────────────────────────────────┘
```

## Component Details

### 1. OCSR Extraction Layer

#### MolScribe Backend

**Architecture**: Swin Transformer with character-level attention

| Attribute | Value |
|-----------|-------|
| **Framework** | PyTorch 2.10+ |
| **Model** | `swin_base_char_aux_1m.pth` |
| **Size** | 1.06 GB |
| **Accuracy** | 93% F1 (single molecule) |
| **Speed (GPU)** | 150-300ms/image |
| **Speed (CPU)** | 1-2s/image |
| **License** | MIT |

**Installation**:
```bash
# Auto-installed from GitHub main branch
pip install git+https://github.com/thomas0809/MolScribe.git

# Checkpoint auto-downloads to ~/.cache/molscribe/
```

**Usage**:
```python
from extractors.molscribe_extractor import MolScribeExtractor

extractor = MolScribeExtractor(
    confidence_threshold=0.5,
    device="cuda"  # or "cpu"
)
result = extractor.extract(image_path)
# {'smiles': '...', 'confidence': 0.92, 'success': True}
```

#### DECIMER Backend

**Architecture**: Vision Transformer (ViT) + T5 decoder

| Attribute | Value |
|-----------|-------|
| **Framework** | TensorFlow 2.15 |
| **Version** | DECIMER 2.7.2 |
| **Accuracy** | 91% F1 (single molecule) |
| **Speed (GPU)** | 200-500ms/image |
| **Speed (CPU)** | 2-5s/image |
| **License** | CC BY 4.0 |

**Installation**:
```bash
pip install decimer>=2.7.0

# Models auto-download via pystow to ~/.pystow/DECIMER-V2/
```

**Usage**:
```python
from extractors.decimer_extractor import DECIMERExtractor

extractor = DECIMERExtractor()
result = extractor.extract(image_path)
```

### 2. Validation Layer

#### Level 1: Syntax Validation (Indigo)

**Purpose**: Ensure SMILES is syntactically valid and canonical

| Check | Description |
|-------|-------------|
| **Parseability** | RDKit/Indigo can parse SMILES |
| **Valence rules** | No valence violations |
| **Ring closures** | All rings properly closed |
| **Wildcards** | Max 2 wildcards allowed |
| **Placeholders** | Reject "n/a", "unknown", etc. |
| **Canonicalization** | Standardize with stereochemistry |

**Implementation**: `src/validators/indigo_validator.py`

```python
from validators.indigo_validator import IndigoValidator

validator = IndigoValidator()
result = validator.validate("CN1CC[C@]2(C1)C3=CC=CC=C3")
# {'is_valid': True, 'canonical_smiles': '...', 'wildcard_count': 0}
```

#### Level 2: Chemical Validation (RDKit)

**Purpose**: Ensure chemical plausibility

| Property | Valid Range | Warning Range |
|----------|-------------|---------------|
| **Molecular Weight** | 100-1000 Da | <150 or >800 |
| **Heavy Atoms** | 5-200 | - |
| **Fragments** | ≤5 | >3 |
| **Ring Systems** | ≤10 | - |
| **Formal Charge** | \|charge\| ≤3 | - |
| **Garbage Detection** | No C.C.C patterns | - |

**Implementation**: `src/validators/chemical_validator.py`

```python
from validators.chemical_validator import ChemicalValidator

validator = ChemicalValidator()
result = validator.validate("CN1CC[C@]2(C1)C3=CC=CC=C3")
# {
#   'is_valid': True,
#   'properties': {
#     'molecular_weight': 147.22,
#     'logp': 2.15,
#     'tpsa': 3.24,
#     ...
#   }
# }
```

#### Level 3: Domain Validation (Alkaloid Filter)

**Purpose**: Ensure relevance to ethnopharmacological research

| Check | Description |
|-------|-------------|
| **Nitrogen presence** | Required for alkaloid classification |
| **Sceletium features** | Methoxy groups, MW 250-400 Da |
| **Scaffold matching** | Mesembrine core, aryl piperidine |
| **Gold standard comparison** | Tanimoto similarity to 7 authenticated compounds |

**Gold Standards**:

| Compound | SMILES | Similarity Threshold |
|----------|--------|---------------------|
| Mesembrine | `CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC` | ≥0.7 |
| Mesembrenone | `CN1CC[C@]2([C@@H]1CC(=O)C=C2)C3=CC(=C(C=C3)OC)OC` | ≥0.7 |
| Tortuosamine (L/D) | `CNCC[C@@]1(CCC2=C(C1)C=CC=N2)C3=...` | ≥0.7 |

**Implementation**: `src/validators/domain_validator.py`

```python
from validators.domain_validator import DomainValidator

validator = DomainValidator()
result = validator.validate("CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC")
# {
#   'is_valid': True,
#   'compound_class': 'alkaloid',
#   'gold_standard_matches': [
#     {'id': 'mesembrine', 'similarity': 1.0, 'match_type': 'exact'}
#   ]
# }
```

### 3. Enrichment Layer

#### Property Calculator

**Properties computed**:
- Molecular weight (MW)
- LogP (octanol-water partition coefficient)
- TPSA (Topological Polar Surface Area)
- H-bond acceptors (HBA)
- H-bond donors (HBD)
- Rotatable bonds
- Stereocenters
- Lipinski Rule of 5 compliance

**Implementation**: `src/enrichers/property_calculator.py`

#### Fingerprint Generator

**Fingerprints generated**:
- **ECFP4** (Morgan, radius=2, 2048 bits) - Used for similarity search
- **MACCS keys** (166 bits) - Used for substructure matching

**Implementation**: `src/enrichers/fingerprint_generator.py`

### 4. Integration Layer

#### JSONL → Markdown Converter

Converts validated molecules to OpenWebUI-compatible Markdown format.

**Output format**:
```markdown
# Mesembrine

**SMILES**: `CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC`

## Properties

| Property | Value |
|----------|-------|
| MW | 289.41 Da |
| LogP | 2.0 |
| TPSA | 32.2 Å² |
| HBA | 4 |
| HBD | 0 |

## Validation

- ✓ Syntax valid (Indigo)
- ✓ Chemically valid (RDKit)
- ✓ Domain valid (Alkaloid)
- Gold standard match: Mesembrine (100% similarity)
```

**Implementation**: `src/utils/kb_converter.py`

## Data Flow

### Complete Pipeline

```python
# 1. Load image
image = load_image("paper_123_figure_4.png")

# 2. Extract (try backends in order)
for backend in config.backend_order:
    result = extractor.extract(image)
    if result['success']:
        break

# 3. Validate (3 layers)
if not indigo.validate(result['smiles'])['is_valid']:
    reject("syntax_invalid")

if not chemical.validate(result['smiles'])['is_valid']:
    reject("chemically_invalid")

if not domain.validate(result['smiles'])['is_valid']:
    reject("domain_invalid")

# 4. Enrich
properties = calculator.calculate(result['smiles'])
fingerprints = generator.compute(result['smiles'])

# 5. Output
output = {
    'smiles': result['smiles'],
    'validation': {...},
    'properties': properties,
    'fingerprints': fingerprints,
    'confidence_score': compute_confidence(...)
}
```

## Configuration Files

### backends.yaml

```yaml
backend_order:
  - molscribe  # Primary (93% F1)
  - decimer    # Fallback (91% F1)

confidence_thresholds:
  molscribe: 0.5
  decimer: 0.3

batch_sizes:
  molscribe: 32
  decimer: 16

gpu_enabled: true
```

### validation_rules.yaml

```yaml
level_1_syntax:
  max_wildcards: 2
  reject_placeholders: true
  placeholders:
    - "n/a"
    - "unknown"
    - "not available"

level_2_chemical:
  mw_min: 100
  mw_max: 1000
  mw_warn_min: 150
  mw_warn_max: 800
  max_heavy_atoms: 200
  max_fragments: 5
  max_rings: 10
  max_formal_charge: 3

level_3_domain:
  require_nitrogen: true
  min_methoxy_groups: 1
  mw_sceletium_min: 250
  mw_sceletium_max: 400
  scaffold_similarity_min: 0.7
```

## Performance Characteristics

### Memory Usage

| Component | Peak Memory |
|-----------|-------------|
| MolScribe (GPU) | 4-6 GB VRAM |
| MolScribe (CPU) | 2-3 GB RAM |
| DECIMER (GPU) | 3-5 GB VRAM |
| DECIMER (CPU) | 2-4 GB RAM |
| Validators | <100 MB |
| Enrichers | <50 MB |

### Computational Complexity

| Operation | Time Complexity | Space Complexity |
|-----------|----------------|------------------|
| MolScribe extraction | O(n²) attention | O(n) activations |
| DECIMER extraction | O(n²) attention | O(n) activations |
| Indigo validation | O(n) | O(1) |
| RDKit validation | O(n²) | O(n) |
| Domain validation | O(n²) × k | O(n) |
| Fingerprint gen | O(n) | O(1) |

Where n = molecule size, k = gold standard count.

## Error Handling

### Failure Modes

| Failure Type | Recovery Strategy |
|--------------|-------------------|
| OCSR backend fails | Try next backend in chain |
| Syntax invalid | Reject, log reason |
| Chemical invalid | Reject, but save to manual review queue |
| Domain invalid | Accept but mark as "non-alkaloid" |
| GPU OOM | Fallback to CPU |
| Model not found | Auto-download from HuggingFace/pystow |

### Retry Logic

```python
def extract_with_retry(image, backends, max_retries=2):
    for attempt in range(max_retries):
        for backend_name in backends:
            try:
                result = backend.extract(image)
                if result['success']:
                    return result
            except ModelLoadError:
                # Auto-download and retry
                backend.download_checkpoint()
            except GPUError:
                # Fallback to CPU
                backend.device = "cpu"
            except Exception as e:
                log_error(e)
    return None
```

## Testing Strategy

### Unit Tests

```bash
# Validators
pytest tests/test_indigo_validator.py
pytest tests/test_chemical_validator.py
pytest tests/test_domain_validator.py

# Extractors
pytest tests/test_molscribe_extractor.py
pytest tests/test_decimer_extractor.py

# Enrichers
pytest tests/test_property_calculator.py
pytest tests/test_fingerprint_generator.py
```

### Integration Tests

```bash
# Full pipeline (10 images)
python scripts/extract_smiles_pipeline.py --limit 10 --test-mode

# End-to-end (100 images)
python scripts/extract_smiles_pipeline.py --limit 100 --gpu
```

### Smoke Tests

```bash
# Verify imports
python -c "from extractors.molscribe import MolScribeExtractor"

# Test extraction
python /tmp/test_molscribe.py
python /tmp/test_decimer.py
```

## Security Considerations

| Concern | Mitigation |
|---------|------------|
| Malicious SMILES | Validation rejects invalid structures |
| Large molecules | MW and atom count limits |
| Model integrity | Checkpoints from trusted sources (HuggingFace, official) |
| API exposure | No external API by default |
| Data leakage | All processing local, no cloud upload |

## Future Enhancements

### Planned (v2.1)

- [ ] MolVec backend (CPU-only fallback)
- [ ] Automated test suite (GitHub Actions)
- [ ] Batch processing optimizations
- [ ] Progress API for OpenWebUI integration
- [ ] Error rate analytics dashboard

### Under Consideration

- [ ] Multi-molecule detection in single image
- [ ] Reaction extraction
- [ ] 3D structure prediction
- [ ] Literature citation extraction

---

**Last updated**: 2026-02-25
**Maintainer**: Thesis Research Team
