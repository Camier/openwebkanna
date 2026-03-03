from evaluate_retrieval_channels import evaluate_queries


def test_evaluate_queries_reports_modality_specific_pre_and_post_recall():
    queries = [
        {
            "query_id": "t1",
            "modality": "text",
            "relevant_doc_ids": ["doc-a", "doc-b"],
            "channels": {
                "text_dense": [{"doc_id": "doc-a"}],
                "bm25": [{"doc_id": "doc-b"}],
                "smiles_fingerprint": [],
            },
        },
        {
            "query_id": "s1",
            "modality": "smiles",
            "relevant_doc_ids": ["doc-z"],
            "channels": {
                "text_dense": [{"doc_id": "doc-x"}],
                "bm25": [{"doc_id": "doc-y"}],
                "smiles_fingerprint": [{"doc_id": "doc-z"}],
            },
        },
    ]

    report = evaluate_queries(
        queries=queries,
        k_values=[1, 5],
        rrf_k=60,
        weights={"smiles_fingerprint": 2.0},
        enable_smiles_embedding=False,
    )

    assert report["summary"]["text"]["Recall@1"]["pre_fusion"] == 0.5
    assert report["summary"]["text"]["Recall@1"]["post_fusion"] == 0.5
    assert report["summary"]["smiles"]["Recall@1"]["pre_fusion"] == 1.0
    assert report["summary"]["smiles"]["Recall@1"]["post_fusion"] == 1.0
