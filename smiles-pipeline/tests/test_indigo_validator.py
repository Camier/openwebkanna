import pytest

pytest.importorskip("indigo")

from validators.indigo_validator import IndigoValidator


def test_indigo_valid_smiles_passes():
    validator = IndigoValidator()
    result = validator.validate_syntax("CCO")
    assert result["is_valid"] is True
    assert result["canonical_smiles"]
    assert result["error"] is None


def test_indigo_placeholder_rejected():
    validator = IndigoValidator()
    result = validator.validate_syntax("unknown")
    assert result["is_valid"] is False
    assert result["error"] == "placeholder_value"


def test_indigo_wildcard_threshold_enforced():
    validator = IndigoValidator(config={"wildcard_threshold": 0})
    result = validator.validate_syntax("CC*")
    assert result["is_valid"] is False
    assert result["error"] == "too_many_wildcards"
