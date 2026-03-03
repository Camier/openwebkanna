from validators.chemical_validator import ChemicalValidator


def test_chemical_validator_rejects_low_mw_fragment():
    validator = ChemicalValidator()
    result = validator.validate("CCO")
    assert result["is_valid"] is False
    assert "mw_too_low" in (result["rejection_reason"] or "")


def test_chemical_validator_accepts_reasonable_alkaloid_like_structure():
    validator = ChemicalValidator()
    mesembrine_like = "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC"
    result = validator.validate(mesembrine_like)
    assert result["is_valid"] is True
    assert result["rejection_reason"] is None
