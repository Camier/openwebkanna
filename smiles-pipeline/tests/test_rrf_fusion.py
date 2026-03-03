from retrieval.rrf_fusion import compute_recall_at_k, fuse_rankings_rrf


def test_fuse_rankings_rrf_prioritizes_multi_channel_hits():
    fused = fuse_rankings_rrf(
        {
            "text_dense": [
                {"doc_id": "doc-a", "score": 0.9},
                {"doc_id": "doc-b", "score": 0.8},
            ],
            "bm25": [
                {"doc_id": "doc-b", "score": 16.0},
                {"doc_id": "doc-c", "score": 10.0},
            ],
            "smiles_fingerprint": [
                {"doc_id": "doc-b", "score": 0.85},
                {"doc_id": "doc-d", "score": 0.7},
            ],
        },
        top_k=3,
        rrf_k=60,
    )

    assert fused[0]["doc_id"] == "doc-b"
    assert fused[0]["rrf_score"] > fused[1]["rrf_score"]
    assert "text_dense" in fused[0]["channels"]
    assert "bm25" in fused[0]["channels"]
    assert "smiles_fingerprint" in fused[0]["channels"]


def test_compute_recall_at_k():
    ranked = ["doc-1", "doc-2", "doc-3", "doc-4"]
    relevant = ["doc-3", "doc-9"]
    recall_at_3 = compute_recall_at_k(ranked, relevant, 3)
    assert recall_at_3 == 0.5
