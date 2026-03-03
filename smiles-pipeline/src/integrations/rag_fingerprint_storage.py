"""
RAG Fingerprint Storage Integration

Integrates molecular fingerprint generation with the RAG pipeline.
Stores ECFP4 + MACCS fingerprints in PostgreSQL pgvector for structure similarity
search.

Phase 2, Task 2.2
"""

import logging
from datetime import datetime
from typing import Any, Dict, List, Literal, Optional, Tuple

import psycopg2
from psycopg2.extras import Json

from enrichers.fingerprint_generator import FingerprintGenerator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RAGFingerprintStorage:
    """
    Store and retrieve molecular fingerprints in PostgreSQL pgvector.

    Usage:
        storage = RAGFingerprintStorage(db_url)
        storage.store_fingerprint(doc_id, smiles, properties)
        results = storage.search_similar(query_fp, threshold=0.7)
    """

    def __init__(self, db_url: str):
        """
        Initialize fingerprint storage.

        Args:
            db_url: PostgreSQL connection URL with pgvector extension
        """
        self.db_url = db_url
        self.fp_generator = FingerprintGenerator()
        self._channel_columns = {
            "ecfp4": "molecule_fingerprint",
            "maccs": "molecule_fingerprint_maccs",
        }

    def _generate_vectors(
        self, smiles: str
    ) -> Tuple[List[float], List[float], str, str]:
        """Generate ECFP4 and MACCS vectors and convert to pgvector literals."""
        ecfp4 = self.fp_generator.ecfp4(smiles)
        maccs = self.fp_generator.maccs(smiles)

        ecfp4_vector = f"[{','.join(map(str, ecfp4))}]"
        maccs_vector = f"[{','.join(map(str, maccs))}]"
        return ecfp4, maccs, ecfp4_vector, maccs_vector

    def _resolve_channel(
        self, fingerprint_type: str
    ) -> Tuple[Literal["ecfp4", "maccs"], str]:
        """Resolve channel label to a concrete DB column."""
        channel = (fingerprint_type or "ecfp4").strip().lower()
        if channel not in self._channel_columns:
            supported = ", ".join(sorted(self._channel_columns.keys()))
            raise ValueError(
                f"Unsupported fingerprint_type '{fingerprint_type}'. "
                f"Supported: {supported}"
            )
        return channel, self._channel_columns[channel]

    def store_fingerprint(
        self,
        doc_id: str,
        smiles: str,
        properties: Optional[Dict[str, Any]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Generate ECFP4 + MACCS fingerprints and store in document_chunk table.

        Args:
            doc_id: Document ID in PostgreSQL
            smiles: SMILES string for fingerprint generation
            properties: Optional molecule properties (MW, LogP, etc.)
            metadata: Optional additional metadata

        Returns:
            Result dict with success status and fingerprint info

        Raises:
            ValueError: If SMILES is invalid
            psycopg2.Error: If database operation fails
        """
        try:
            ecfp4, maccs, ecfp4_vector, maccs_vector = self._generate_vectors(smiles)

            # Build metadata
            full_metadata = {
                "source": "smiles_pipeline",
                "generated_at": datetime.now().isoformat(),
                "fingerprint_types": ["ecfp4", "maccs"],
                "fingerprint_dims": {"ecfp4": len(ecfp4), "maccs": len(maccs)},
            }

            if properties:
                full_metadata["properties"] = properties
            if metadata:
                full_metadata.update(metadata)

            # Update database
            self._update_document(
                doc_id,
                smiles,
                ecfp4_vector,
                maccs_vector,
                full_metadata,
            )

            logger.info(
                "Stored fingerprints for doc %s (ecfp4=%s dims, maccs=%s dims)",
                doc_id,
                len(ecfp4),
                len(maccs),
            )

            return {
                "success": True,
                "doc_id": doc_id,
                "smiles": smiles,
                "fingerprint_dims": {"ecfp4": len(ecfp4), "maccs": len(maccs)},
                "metadata": full_metadata,
            }

        except ValueError as e:
            logger.error(f"Invalid SMILES for doc {doc_id}: {e}")
            raise
        except Exception as e:
            logger.error(f"Fingerprint storage failed: {e}")
            raise

    def store_fingerprints_bulk(
        self, documents: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Store fingerprints for multiple document_chunk in a single transaction.

        Args:
            document_chunk: List of dicts with keys:
                - doc_id (str): Document ID
                - smiles (str): SMILES string
                - properties (dict, optional): Molecule properties
                - metadata (dict, optional): Additional metadata

        Returns:
            Summary with success/failure counts

        Example:
            document_chunk = [
                {"doc_id": "doc_001", "smiles": "CN1CC...", "properties": {...}},
                {"doc_id": "doc_002", "smiles": "CCO...", "properties": {...}},
            ]
            result = storage.store_fingerprints_bulk(document_chunk)
        """
        success_count = 0
        failure_count = 0
        failures = []

        try:
            with psycopg2.connect(self.db_url) as conn:
                with conn.cursor() as cur:
                    for doc in documents:
                        doc_id = doc["doc_id"]
                        smiles = doc["smiles"]
                        properties = doc.get("properties", {})
                        metadata = doc.get("metadata", {})

                        try:
                            ecfp4, maccs, ecfp4_vector, maccs_vector = (
                                self._generate_vectors(smiles)
                            )

                            # Build metadata
                            full_metadata = {
                                "source": "smiles_pipeline",
                                "generated_at": datetime.now().isoformat(),
                                "fingerprint_types": ["ecfp4", "maccs"],
                                "fingerprint_dims": {
                                    "ecfp4": len(ecfp4),
                                    "maccs": len(maccs),
                                },
                                "properties": properties,
                            }
                            if metadata:
                                full_metadata.update(metadata)

                            # Update
                            cur.execute(
                                """
                                UPDATE document_chunk
                                SET molecule_fingerprint = %s::halfvec,
                                    molecule_fingerprint_maccs = %s::halfvec,
                                    molecule_smiles = %s,
                                    molecule_metadata = %s::jsonb
                                WHERE id = %s
                            """,
                                (
                                    ecfp4_vector,
                                    maccs_vector,
                                    smiles,
                                    Json(full_metadata),
                                    doc_id,
                                ),
                            )

                            success_count += 1

                        except Exception as e:
                            failure_count += 1
                            failures.append(
                                {"doc_id": doc_id, "smiles": smiles, "error": str(e)}
                            )
                            logger.warning(f"Failed for doc {doc_id}: {e}")

            logger.info(
                f"Bulk storage complete: {success_count} succeeded, "
                f"{failure_count} failed"
            )

            return {
                "success": True,
                "total": len(documents),
                "succeeded": success_count,
                "failed": failure_count,
                "failures": failures,
            }

        except Exception as e:
            logger.error(f"Bulk storage failed: {e}")
            return {
                "success": False,
                "error": str(e),
                "partial_results": {
                    "succeeded": success_count,
                    "failed": failure_count,
                },
            }

    def _update_document(
        self,
        doc_id: str,
        smiles: str,
        ecfp4_vector: str,
        maccs_vector: str,
        metadata: Dict[str, Any],
    ):
        """Update single document with molecular data."""
        with psycopg2.connect(self.db_url) as conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    UPDATE document_chunk
                    SET molecule_fingerprint = %s::halfvec,
                        molecule_fingerprint_maccs = %s::halfvec,
                        molecule_smiles = %s,
                        molecule_metadata = %s::jsonb
                    WHERE id = %s
                """,
                    (ecfp4_vector, maccs_vector, smiles, Json(metadata), doc_id),
                )

    def search_similar(
        self,
        query_smiles: str,
        threshold: float = 0.7,
        top_k: int = 10,
        fingerprint_type: Literal["ecfp4", "maccs"] = "ecfp4",
    ) -> List[Dict[str, Any]]:
        """
        Search for molecules similar to query SMILES.

        Args:
            query_smiles: Query SMILES string
            threshold: Minimum similarity threshold (0.0-1.0)
            top_k: Maximum results to return
            fingerprint_type: Channel to query ('ecfp4' or 'maccs')

        Returns:
            List of similar molecules with similarity scores
        """
        channel, fp_column = self._resolve_channel(fingerprint_type)

        # Generate query fingerprint
        if channel == "maccs":
            query_fp = self.fp_generator.maccs(query_smiles)
        else:
            query_fp = self.fp_generator.ecfp4(query_smiles)
        fp_vector = f"[{','.join(map(str, query_fp))}]"

        with psycopg2.connect(self.db_url) as conn:
            with conn.cursor() as cur:
                cur.execute(
                    f"""
                    SELECT
                        id,
                        molecule_smiles,
                        1 - ({fp_column} <=> %s::halfvec) as similarity,
                        metadata->>'title' as title,
                        molecule_metadata
                    FROM document_chunk
                    WHERE {fp_column} IS NOT NULL
                      AND 1 - ({fp_column} <=> %s::halfvec) >= %s
                    ORDER BY similarity DESC
                    LIMIT %s
                """,
                    (fp_vector, fp_vector, threshold, top_k),
                )

                results = []
                for row in cur.fetchall():
                    results.append(
                        {
                            "doc_id": row[0],
                            "smiles": row[1],
                            "similarity": float(row[2]),
                            "fingerprint_type": channel,
                            "document_title": row[3],
                            "molecule_metadata": row[4],
                        }
                    )

                return results

    def get_fingerprint(
        self, doc_id: str, fingerprint_type: Literal["ecfp4", "maccs"] = "ecfp4"
    ) -> Optional[List[float]]:
        """
        Retrieve stored fingerprint for a document.

        Args:
            doc_id: Document ID
            fingerprint_type: Channel to retrieve ('ecfp4' or 'maccs')

        Returns:
            Fingerprint vector as float list, or None if not found
        """
        _, fp_column = self._resolve_channel(fingerprint_type)
        with psycopg2.connect(self.db_url) as conn:
            with conn.cursor() as cur:
                cur.execute(
                    f"""
                    SELECT {fp_column}
                    FROM document_chunk
                    WHERE id = %s AND {fp_column} IS NOT NULL
                """,
                    (doc_id,),
                )

                result = cur.fetchone()
                if result:
                    if isinstance(result[0], list):
                        return [float(x) for x in result[0]]
                    # Parse pgvector format: "[0.0,1.0,0.0,...]"
                    fp_str = result[0].strip("[]")
                    return [float(x) for x in fp_str.split(",")]
                return None

    def get_statistics(self) -> Dict[str, Any]:
        """
        Get fingerprint storage statistics.

        Returns:
            Dict with counts and metadata
        """
        with psycopg2.connect(self.db_url) as conn:
            with conn.cursor() as cur:
                # Count document_chunk with fingerprints
                cur.execute("""
                    SELECT
                        COUNT(*) as total_docs,
                        COUNT(*) FILTER (
                            WHERE molecule_fingerprint IS NOT NULL
                               OR molecule_fingerprint_maccs IS NOT NULL
                        ) as docs_with_any_fps,
                        COUNT(molecule_fingerprint) as docs_with_ecfp4,
                        COUNT(molecule_fingerprint_maccs) as docs_with_maccs,
                        COUNT(DISTINCT molecule_smiles) as unique_smiles
                    FROM document_chunk
                """)

                row = cur.fetchone()

                if row is None:
                    return {
                        "total_documents": 0,
                        "documents_with_fingerprints": 0,
                        "documents_with_ecfp4": 0,
                        "documents_with_maccs": 0,
                        "unique_smiles": 0,
                        "coverage": 0,
                        "coverage_ecfp4": 0,
                        "coverage_maccs": 0,
                    }

                return {
                    "total_documents": row[0],
                    "documents_with_fingerprints": row[1],
                    "documents_with_ecfp4": row[2],
                    "documents_with_maccs": row[3],
                    "unique_smiles": row[4],
                    "coverage": row[1] / row[0] if row[0] > 0 else 0,
                    "coverage_ecfp4": row[2] / row[0] if row[0] > 0 else 0,
                    "coverage_maccs": row[3] / row[0] if row[0] > 0 else 0,
                }


def process_extraction_results(
    extraction_results: List[Dict[str, Any]], db_url: str, batch_size: int = 100
) -> Dict[str, Any]:
    """
    Process SMILES extraction results and store fingerprints.

    Convenience function for batch processing extraction pipeline output.

    Args:
        extraction_results: List from smiles-pipeline extraction
            Each item should have:
            - paper_id (str)
            - molecules (list): List of molecule dicts with smiles
        db_url: PostgreSQL URL
        batch_size: Batch size for bulk storage

    Returns:
        Processing summary
    """
    storage = RAGFingerprintStorage(db_url)

    # Flatten to document list
    document_chunk = []
    for result in extraction_results:
        paper_id = result.get("paper_id")
        molecules = result.get("molecules", [])

        for mol in molecules:
            smiles = mol.get("smiles")
            if not smiles:
                continue

            properties = {
                "confidence": mol.get("confidence_score"),
                "compound_class": mol.get("compound_class"),
                "gold_standard_matches": mol.get("gold_standard_matches", []),
            }

            metadata = {
                "paper_id": paper_id,
                "backend_used": mol.get("backend_used"),
                "validation": mol.get("validation"),
            }

            document_chunk.append(
                {
                    "doc_id": paper_id,
                    "smiles": smiles,
                    "properties": properties,
                    "metadata": metadata,
                }
            )

    # Store in batches
    all_results = []
    for i in range(0, len(document_chunk), batch_size):
        batch = document_chunk[i : i + batch_size]
        result = storage.store_fingerprints_bulk(batch)
        all_results.append(result)

    # Aggregate results
    total_success = sum(r.get("succeeded", 0) for r in all_results)
    total_failed = sum(r.get("failed", 0) for r in all_results)

    return {
        "total_processed": len(document_chunk),
        "succeeded": total_success,
        "failed": total_failed,
        "batches": len(all_results),
    }


if __name__ == "__main__":
    # Example usage
    import os

    db_url = os.getenv(
        "DATABASE_URL",
        "postgresql://openwebui:openwebui_pgvector_pass_2026@localhost:5432/openwebui",
    )

    storage = RAGFingerprintStorage(db_url)

    # Test single storage
    test_smiles = "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC"  # Mesembrine

    try:
        result = storage.store_fingerprint(
            doc_id="test_001",
            smiles=test_smiles,
            properties={"mw": 289.41, "logp": 2.0},
        )
        print(f"Storage result: {result}")
    except Exception as e:
        print(f"Storage failed: {e}")

    # Test similarity search
    results = storage.search_similar(query_smiles=test_smiles, threshold=0.5, top_k=5)
    print(f"Similar molecules: {results}")

    # Get statistics
    stats = storage.get_statistics()
    print(f"Database statistics: {stats}")
