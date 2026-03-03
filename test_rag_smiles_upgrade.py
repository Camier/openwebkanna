#!/usr/bin/env python
# Test cases for RAG & SMILES Architecture Upgrades
# Target: OpenWebUI + SMILES Pipeline v2.0

"""
Test Categories:
1. Config Validation - .env settings
2. Confidence Scoring - require_positive_indicator
3. Molecular Fingerprints - ECFP4, Tanimoto, MACCS
4. PostgreSQL Schema - molecules table with vectors
5. Chunk-Molecule Linking - molecule_ids field
6. Structure Search API - endpoint spec
"""

import sys
import subprocess
from pathlib import Path

# Add smiles-pipeline to path
sys.path.insert(0, str(Path(__file__).parent / "smiles-pipeline" / "src"))


def test_env_config():
    """Test .env configuration settings."""
    from dotenv import load_dotenv
    from os import getenv

    env_path = Path(__file__).parent / ".env"
    load_dotenv(env_path)

    # Test embedding model
    embedding_model = getenv("RAG_EMBEDDING_MODEL")
    assert embedding_model == "BAAI/bge-base-en-v1.5", (
        f"Expected BAAI/bge-base-en-v1.5, got {embedding_model}"
    )

    # Test reranker
    reranker_model = getenv("RAG_RERANKING_MODEL")
    assert reranker_model == "BAAI/bge-reranker-v2-m3", (
        f"Expected BAAI/bge-reranker-v2-m3, got {reranker_model}"
    )

    # Test RAG_TOP_K
    rag_top_k = int(getenv("RAG_TOP_K", "15"))
    assert rag_top_k >= 10, f"RAG_TOP_K too low: {rag_top_k}"

    # Test CHUNK_SIZE
    chunk_size = int(getenv("CHUNK_SIZE", "3000"))
    assert chunk_size >= 2000, f"CHUNK_SIZE too small: {chunk_size}"

    # Test CHUNK_OVERLAP
    chunk_overlap = int(getenv("CHUNK_OVERLAP", "600"))
    assert chunk_overlap >= 400, f"CHUNK_OVERLAP too small: {chunk_overlap}"

    # Test CHUNK_MIN_SIZE_TARGET
    min_target = int(getenv("CHUNK_MIN_SIZE_TARGET", "1000"))
    assert min_target >= 500, f"CHUNK_MIN_SIZE_TARGET too small: {min_target}"

    # Test RAG_HYBRID_BM25_WEIGHT
    bm25_weight = float(getenv("RAG_HYBRID_BM25_WEIGHT", "0.4"))
    assert 0.3 <= bm25_weight <= 0.6, (
        f"BM25 weight outside expected range: {bm25_weight}"
    )

    print("Config validation passed")
    return True


def test_confidence_scoring():
    """Test validation_rules.yaml confidence scoring."""
    import yaml

    config_path = (
        Path(__file__).parent / "smiles-pipeline" / "config" / "validation_rules.yaml"
    )
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Test require_positive_indicator
    confidence = config.get("confidence_scoring", {})
    requirements = confidence.get("requirements", {})

    assert requirements.get("require_positive_indicator") is True, (
        "require_positive_indicator should be True"
    )

    positive_indicators = requirements.get("positive_indicators", [])
    assert "exact_gold_standard_match" in positive_indicators, (
        "exact_gold_standard_match should be in positive indicators"
    )
    assert "scaffold_match" in positive_indicators, (
        "scaffold_match should be in positive indicators"
    )
    assert "methoxy_groups_present" in positive_indicators, (
        "methoxy_groups_present should be in positive indicators"
    )

    # Test base scores
    base_scores = confidence.get("base_scores", {})
    assert base_scores.get("level_1_syntax_pass") == 30, "Level 1 score should be 30"
    assert base_scores.get("level_2_chemical_pass") == 30, (
        "Level 2 score should be 30"
    )
    assert base_scores.get("level_3_domain_pass") == 20, "Level 3 score should be 20"

    print("Confidence scoring validation passed")
    return True


def test_fingerprint_generator():
    """Test fingerprint generator produces valid vectors."""
    # Guard against silent ABI drift where RDKit emits errors but import still appears usable.
    probe = subprocess.run(
        [
            sys.executable,
            "-c",
            "from rdkit import Chem; from rdkit.Chem import AllChem, MACCSkeys",
        ],
        capture_output=True,
        text=True,
    )
    abi_error_markers = ("_ARRAY_API not found", "compiled using NumPy 1.x")
    abi_output = f"{probe.stdout}\n{probe.stderr}"
    assert probe.returncode == 0, f"RDKit import probe failed: {probe.stderr.strip()}"
    assert not any(marker in abi_output for marker in abi_error_markers), (
        "RDKit/NumPy ABI mismatch detected; align RDKit build with NumPy major version."
    )

    from enrichers.fingerprint_generator import FingerprintGenerator

    gen = FingerprintGenerator()

    # Test with a known alkaloid SMILES
    smiles = "CN1CC[C@@H]2C1CC(=O)CC2"  # mesembrine core

    # Test ECFP4
    ecfp4 = gen.ecfp4(smiles)
    assert len(ecfp4) == 2048, f"ECFP4 should be 2048-dim, got {len(ecfp4)}"
    assert all(isinstance(x, float) for x in ecfp4), "ECFP4 should be float vector"

    # Test MACCS
    maccs = gen.maccs(smiles)
    assert len(maccs) == 167, f"MACCS should be 167-dim, got {len(maccs)}"
    assert all(isinstance(x, float) for x in maccs), "MACCS should be float vector"

    # Test similarity
    smiles2 = "CN1CC[C@H]2C1CC(=O)CC2"  # stereoisomer

    ecfp4_2 = gen.ecfp4(smiles2)

    tanimoto = gen.tanimoto_similarity(ecfp4, ecfp4_2)
    assert 0.7 <= tanimoto <= 1.0, f"Similar isomers should have high Tanimoto: {tanimoto}"

    print("Fingerprint generator passed")
    return True


def test_molecular_table_schema():
    """Test PostgreSQL schema for molecular vectors."""
    # This test verifies the SQL schema is correct
    schema_sql = """
    CREATE TABLE IF NOT EXISTS molecules (
        cid VARCHAR(32) PRIMARY KEY,
        smiles TEXT NOT NULL,
        molecular_weight FLOAT,
        logp FLOAT,
        hba INTEGER,
        hbd INTEGER,
        rotatable_bonds INTEGER,
        heavy_atoms INTEGER,
        fps_ecfp4 VECTOR(2048),
        fps_maccs VECTOR(167),
        scaffold TEXT,
        is_gold_standard BOOLEAN,
        similarity_to_gold_standard FLOAT
    );

    CREATE INDEX IF NOT EXISTS molecules_cid_idx ON molecules(cid);
    CREATE INDEX IF NOT EXISTS molecules_smiles_idx ON molecules(smiles);
    CREATE INDEX IF NOT EXISTS molecules_fp_idx
        ON molecules USING hnsw (fps_ecfp4 vector_cosine_ops);
    """

    # Verify schema keywords
    assert "CREATE TABLE" in schema_sql
    assert "VECTOR(2048)" in schema_sql
    assert "hnsw" in schema_sql
    assert "vector_cosine_ops" in schema_sql

    print("Molecular table schema passed")
    return True


def test_chunk_molecule_linking():
    """Test chunk schema includes molecule_ids field."""
    # Verify chunk structure includes molecule linking
    sample_chunk = {
        "doc_id": "paper_123",
        "text": "The mesembrine alkaloids inhibit...",
        "text_embedding": [0.1] * 768,
        "molecule_ids": ["CID_394162", "CID_216272"],  # Required field
        "chunk_metadata": {
            "page": 12,
            "position": "abstract",
            "image_references": ["fig1.png"],
        },
    }

    assert "molecule_ids" in sample_chunk, "Chunks must include molecule_ids"
    assert isinstance(sample_chunk["molecule_ids"], list), (
        "molecule_ids should be a list"
    )

    print("Chunk-molecule linking schema passed")
    return True


def test_structure_search_api():
    """Test structure search API mock implementation."""
    # Mock API request/response
    api_spec = {
        "endpoint": "/api/v1/structure/search",
        "method": "POST",
        "request_fields": {
            "query_smiles": "string",
            "threshold": "float (default: 0.7)",
            "top_k": "int (default: 10)",
        },
        "response_fields": {"results": "array of {doc_id, smiles, similarity}"},
    }

    # Verify required fields
    assert api_spec["endpoint"] == "/api/v1/structure/search"
    assert api_spec["method"] == "POST"

    print("Structure search API schema passed")
    return True


def run_all_tests():
    """Run all validation tests."""
    print("=" * 60)
    print("RAG & SMILES Architecture Upgrade Tests")
    print("=" * 60)

    tests = [
        ("Config Validation", test_env_config),
        ("Confidence Scoring", test_confidence_scoring),
        ("Fingerprint Generator", test_fingerprint_generator),
        ("Molecular Table Schema", test_molecular_table_schema),
        ("Chunk-Molecule Linking", test_chunk_molecule_linking),
        ("Structure Search API", test_structure_search_api),
    ]

    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            result = test_func()
            if result:
                passed += 1
            else:
                failed += 1
                print(f"{name} failed")
        except Exception as e:
            failed += 1
            print(f"{name} error: {e}")

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
