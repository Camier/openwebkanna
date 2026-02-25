# SMILES Pipeline Documentation

This directory contains comprehensive documentation for the SMILES extraction pipeline.

---

## Document Index

| Document | Audience | Purpose | Status |
|----------|----------|---------|--------|
| [ARCHITECTURE_v2.md](ARCHITECTURE_v2.md) | Developers, Architects | System design, component details, technical specifications | ✅ Complete |
| [USAGE_GUIDE.md](USAGE_GUIDE.md) | Operators, Researchers | Step-by-step procedures, examples, troubleshooting | ✅ Complete |
| [INSTALLATION.md](../INSTALLATION.md) | New users | Setup, environment creation, verification | ✅ Complete |
| [README.md](../README.md) | All users | Quick start, overview, quick reference | ✅ Complete |
| [MIGRATION.md](../MIGRATION.md) | v1 users | Migration from legacy v1 to v2 | ✅ Available |

---

## By Audience

### For New Users

**Start here**:
1. [README.md](../README.md) - Overview and quick start
2. [INSTALLATION.md](../INSTALLATION.md) - Setup environment
3. [USAGE_GUIDE.md](USAGE_GUIDE.md) § "Quick Reference" - Run first extraction

**Time to first result**: ~30 minutes

### For Operators (Daily Use)

**Key documents**:
- [USAGE_GUIDE.md](USAGE_GUIDE.md) - Procedures and best practices
- [USAGE_GUIDE.md](USAGE_GUIDE.md) § "Common Commands" - Quick reference
- [USAGE_GUIDE.md](USAGE_GUIDE.md) § "Troubleshooting Guide" - Problem resolution

### For Developers

**Technical documentation**:
1. [ARCHITECTURE_v2.md](ARCHITECTURE_v2.md) - System architecture and design
2. [ARCHITECTURE_v2.md](ARCHITECTURE_v2.md) § "Component Details" - API reference
3. [ARCHITECTURE_v2.md](ARCHITECTURE_v2.md) § "Testing Strategy" - Test guidelines

### For Researchers

**Focus on**:
- [README.md](../README.md) § "Gold Standard Compounds" - Validation targets
- [USAGE_GUIDE.md](USAGE_GUIDE.md) § "Example Sessions" - Typical workflows
- [USAGE_GUIDE.md](USAGE_GUIDE.md) § "Data Management" - Output organization

---

## Document Purposes

### ARCHITECTURE_v2.md

**What it covers**:
- System diagram and data flow
- Component specifications (MolScribe, DECIMER, validators)
- API signatures and usage
- Performance characteristics
- Configuration file formats
- Error handling strategies

**Use when**:
- Understanding how components interact
- Adding new backends or validators
- Debugging complex issues
- Optimizing performance

### USAGE_GUIDE.md

**What it covers**:
- Quick reference commands
- Step-by-step procedures
- Best practices
- Example sessions
- Troubleshooting guide
- FAQ

**Use when**:
- Running extraction pipeline
- Troubleshooting errors
- Learning typical workflows
- Tuning parameters

### README.md

**What it covers**:
- Project overview
- Quick start (3 commands)
- Architecture summary
- Installation status
- Quick links to other docs

**Use when**:
- First time using pipeline
- Need quick reminder
- Looking for documentation links

### INSTALLATION.md

**What it covers**:
- Prerequisites
- Environment creation
- Verification steps
- Troubleshooting installation issues
- Update procedures

**Use when**:
- Setting up on new machine
- Reinstalling environment
- Diagnosing import errors

---

## Quick Reference

### Installation

```bash
cd /LAB/@thesis/openwebui/smiles-pipeline/envs
micromamba env create -f environment-smiles-extraction.yml --yes
```

### Basic Extraction

```bash
micromamba run -n smiles-extraction python smiles-pipeline/scripts/extract_smiles_pipeline.py \
  --input-dir data/extractions \
  --output-dir smiles-pipeline/data/validated \
  --gpu
```

### View Results

```bash
cat smiles-pipeline/data/validated/extraction_stats.json | jq
```

### Import to OpenWebUI

```bash
./smiles-pipeline/import-smiles-to-kb.sh --input-dir smiles-pipeline/data/validated
```

---

## File Organization

```
smiles-pipeline/
├── README.md                    # Main documentation (this document)
├── INSTALLATION.md             # Setup guide
├── MIGRATION.md                # v1 → v2 migration
├── docs/                       # Detailed documentation
│   ├── README_DOCS.md          # This file (document index)
│   ├── ARCHITECTURE_v2.md      # Technical architecture
│   └── USAGE_GUIDE.md          # Operational procedures
├── config/                     # Configuration files
│   ├── backends.yaml
│   └── validation_rules.yaml
├── src/                        # Source code
├── scripts/                    # Pipeline scripts
└── data/                       # Data directory
```

---

## Document Maintenance

### When to Update

| Trigger | Action |
|---------|--------|
| New feature added | Update ARCHITECTURE_v2.md, USAGE_GUIDE.md |
| Bug fix with user impact | Update USAGE_GUIDE.md § Troubleshooting |
| New dependency | Update INSTALLATION.md |
| Performance optimization | Update ARCHITECTURE_v2.md § Performance |
| API change | Update ARCHITECTURE_v2.md § API Reference |

### Version Control

- All docs versioned with pipeline (v2.0.0)
- Changelog in README.md
- Breaking changes migrate to MIGRATION.md

### Style Guide

- Use markdown formatting
- Include code examples for all commands
- Use tables for comparisons
- Include troubleshooting sections
- Add "Last updated" date

---

## Finding Information

**Lookup table for common questions**:

| Question | Document | Section |
|----------|----------|---------|
| How do I install? | INSTALLATION.md | Quick Start |
| How do I run extraction? | USAGE_GUIDE.md | Procedure 1 |
| What's the architecture? | ARCHITECTURE_v2.md | System Overview |
| Why is it slow? | USAGE_GUIDE.md | Slow extraction speed |
| How do I add a backend? | ARCHITECTURE_v2.md | Adding New Backends |
| What are gold standards? | README.md | Gold Standard Compounds |
| How do I import to OpenWebUI? | USAGE_GUIDE.md | Procedure 3 |
| Where are outputs? | USAGE_GUIDE.md | Data Management |

---

## Support

### Documentation Issues

If documentation is unclear, incomplete, or incorrect:

1. Check if another doc addresses it (use index table)
2. Verify against code in `src/`
3. Create issue with specific section reference

### Requesting Additions

For missing documentation:

- **Procedures** → Add to USAGE_GUIDE.md
- **Technical details** → Add to ARCHITECTURE_v2.md
- **Setup steps** → Add to INSTALLATION.md
- **Overview** → Add to README.md

---

## Document Status

| Document | Completeness | Accuracy | Last Review |
|----------|--------------|----------|-------------|
| README.md | ✅ Complete | ✅ Verified | 2026-02-25 |
| INSTALLATION.md | ✅ Complete | ✅ Verified | 2026-02-25 |
| ARCHITECTURE_v2.md | ✅ Complete | ✅ Verified | 2026-02-25 |
| USAGE_GUIDE.md | ✅ Complete | ✅ Verified | 2026-02-25 |
| MIGRATION.md | ⚠️ Partial | Needs update | 2026-01-31 |

**Status legend**:
- ✅ Complete and tested
- ⚠️ Partial or needs update
- ❌ Incomplete or outdated

---

**Last updated**: 2026-02-25
**Maintainer**: Thesis Research Team
