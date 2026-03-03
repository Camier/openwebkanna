from integrations.rag_fingerprint_storage import RAGFingerprintStorage


def test_resolve_channel_accepts_ecfp4_and_maccs():
    storage = RAGFingerprintStorage("postgresql://unused")

    channel_ecfp4, col_ecfp4 = storage._resolve_channel("ecfp4")
    channel_maccs, col_maccs = storage._resolve_channel("maccs")

    assert channel_ecfp4 == "ecfp4"
    assert col_ecfp4 == "molecule_fingerprint"
    assert channel_maccs == "maccs"
    assert col_maccs == "molecule_fingerprint_maccs"


def test_resolve_channel_rejects_unknown_type():
    storage = RAGFingerprintStorage("postgresql://unused")

    try:
        storage._resolve_channel("unknown")
        assert False, "Expected ValueError for unknown fingerprint type"
    except ValueError as exc:
        assert "Unsupported fingerprint_type" in str(exc)
