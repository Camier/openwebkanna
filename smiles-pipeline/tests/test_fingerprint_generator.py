from enrichers.fingerprint_generator import FingerprintGenerator


def test_fingerprint_generator_includes_maccs_key_and_metadata():
    smiles = "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC"
    generator = FingerprintGenerator()

    result = generator.generate(smiles)

    assert result["success"] is True
    assert "ecfp4" in result["fingerprints"]
    assert "maccs" in result["fingerprints"]
    assert "macces" not in result["fingerprints"]
    assert result["fingerprints"]["maccs"]["nbits"] == 167
    assert result["fingerprints"]["maccs"]["hex"]


def test_fingerprint_generator_accepts_legacy_macces_alias():
    smiles = "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC"
    generator = FingerprintGenerator(fingerprint_types=["macces"])

    result = generator.generate(smiles)

    assert result["success"] is True
    assert "maccs" in result["fingerprints"]
    assert result["fingerprints"]["maccs"]["nbits"] == 167
