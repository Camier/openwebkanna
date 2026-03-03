from validators.domain_validator import DomainValidator


def test_domain_validator_non_alkaloid_classification():
    validator = DomainValidator()
    result = validator.validate("CCO")
    assert result["compound_class"] == "non_alkaloid"
    assert result["is_valid"] is False
    assert result["rejection_reason"] == "no_nitrogen_not_alkaloid"


def test_domain_validator_alkaloid_classification():
    validator = DomainValidator()
    mesembrine_like = "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC"
    result = validator.validate(mesembrine_like)
    assert result["compound_class"] == "alkaloid"
