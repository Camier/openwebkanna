"""
Structure Search API

FastAPI endpoints for molecular structure similarity search.
Uses ECFP4 fingerprints and pgvector for sub-second similarity queries.

Phase 2, Task 2.3
"""

import os
import logging
from typing import List, Optional
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel, Field
import psycopg2

from enrichers.fingerprint_generator import FingerprintGenerator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Database URL from environment
DATABASE_URL = os.getenv(
    "DATABASE_URL",
    "postgresql://openwebui:openwebui_pgvector_pass_2026@localhost:5432/openwebui",
)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan manager."""
    # Startup: verify database connection
    try:
        conn = psycopg2.connect(DATABASE_URL)
        conn.close()
        logger.info("Database connection verified")
    except Exception as e:
        logger.warning(f"Database connection check failed: {e}")

    yield

    # Shutdown: cleanup if needed


app = FastAPI(
    title="Structure Search API",
    description="Molecular structure similarity search using ECFP4 fingerprints",
    version="1.0.0",
    lifespan=lifespan,
)


# =============================================================================
# Request/Response Models
# =============================================================================


class StructureSearchRequest(BaseModel):
    """Request model for structure similarity search."""

    query_smiles: str = Field(..., description="Query SMILES string")
    threshold: float = Field(
        default=0.7,
        ge=0.0,
        le=1.0,
        description="Minimum similarity threshold (0.0-1.0)",
    )
    top_k: int = Field(
        default=10, ge=1, le=100, description="Maximum number of results to return"
    )


class StructureSearchResult(BaseModel):
    """Single search result."""

    doc_id: str = Field(..., description="Document ID")
    smiles: str = Field(..., description="SMILES string")
    similarity: float = Field(..., description="Similarity score (0.0-1.0)")
    document_title: Optional[str] = Field(None, description="Document title")
    molecule_metadata: Optional[dict] = Field(None, description="Molecule metadata")


class StructureSearchResponse(BaseModel):
    """Response model for structure search."""

    query_smiles: str
    results: List[StructureSearchResult]
    total_results: int


class MoleculeInfo(BaseModel):
    """Molecule information response."""

    doc_id: str
    smiles: str
    similarity: Optional[float] = None
    document_title: Optional[str] = None
    molecule_metadata: Optional[dict] = None


# =============================================================================
# API Endpoints
# =============================================================================


@app.post(
    "/v1/structure/search",
    response_model=StructureSearchResponse,
    tags=["Structure Search"],
)
async def search_similar_structures(request: StructureSearchRequest):
    """
    Search for molecules similar to query structure.

    Uses ECFP4 fingerprints with cosine similarity search via pgvector HNSW index.
    Returns molecules with similarity >= threshold, sorted by similarity descending.

    **Example:**
    ```json
    {
        "query_smiles": "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC",
        "threshold": 0.7,
        "top_k": 10
    }
    ```
    """
    try:
        # Step 1: Generate query fingerprint
        fp_generator = FingerprintGenerator()
        query_fp = fp_generator.ecfp4(request.query_smiles)

        # Convert to pgvector format: "[0.0,1.0,0.0,...]"
        fp_vector = f"[{','.join(map(str, query_fp))}]"

        # Step 2: Execute similarity search
        conn = psycopg2.connect(DATABASE_URL)
        cur = conn.cursor()

        cur.execute(
            """
            SELECT
                id,
                molecule_smiles,
                1 - (molecule_fingerprint <=> %s::halfvec) as similarity,
                metadata->>'title' as title,
                molecule_metadata
            FROM document_chunk
            WHERE molecule_fingerprint IS NOT NULL
              AND 1 - (molecule_fingerprint <=> %s::halfvec) >= %s
            ORDER BY similarity DESC
            LIMIT %s
        """,
            (fp_vector, fp_vector, request.threshold, request.top_k),
        )

        results = []
        for row in cur.fetchall():
            results.append(
                StructureSearchResult(
                    doc_id=row[0],
                    smiles=row[1],
                    similarity=float(row[2]),
                    document_title=row[3],
                    molecule_metadata=row[4],
                )
            )

        cur.close()
        conn.close()

        return StructureSearchResponse(
            query_smiles=request.query_smiles,
            results=results,
            total_results=len(results),
        )

    except ValueError as e:
        # Invalid SMILES
        raise HTTPException(status_code=400, detail=f"Invalid SMILES: {str(e)}")
    except psycopg2.Error as e:
        logger.error(f"Database error: {e}")
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        logger.error(f"Search failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get(
    "/v1/structure/{doc_id}", response_model=MoleculeInfo, tags=["Molecule Lookup"]
)
async def get_molecule(doc_id: str):
    """
    Retrieve molecule information by document ID.

    Returns stored SMILES, fingerprint metadata, and document information.
    """
    try:
        conn = psycopg2.connect(DATABASE_URL)
        cur = conn.cursor()

        cur.execute(
            """
            SELECT
                id,
                molecule_smiles,
                metadata->>'title' as title,
                molecule_metadata
            FROM document_chunk
            WHERE id = %s AND molecule_smiles IS NOT NULL
        """,
            (doc_id,),
        )

        row = cur.fetchone()
        cur.close()
        conn.close()

        if not row:
            raise HTTPException(
                status_code=404, detail=f"Molecule not found for document {doc_id}"
            )

        return MoleculeInfo(
            doc_id=row[0],
            smiles=row[1],
            document_title=row[2],
            molecule_metadata=row[3],
        )

    except HTTPException:
        raise
    except psycopg2.Error as e:
        logger.error(f"Database error: {e}")
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        logger.error(f"Lookup failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/v1/structure/stats", tags=["Statistics"])
async def get_statistics():
    """
    Get fingerprint database statistics.

    Returns counts of molecules, document_chunk, and coverage metrics.
    """
    try:
        conn = psycopg2.connect(DATABASE_URL)
        cur = conn.cursor()

        cur.execute(
            """
            SELECT
                COUNT(*) as total_docs,
                COUNT(molecule_fingerprint) as docs_with_fps,
                COUNT(DISTINCT molecule_smiles) as unique_smiles,
                COUNT(DISTINCT id) as unique_documents
            FROM document_chunk
        """
        )

        row = cur.fetchone()
        cur.close()
        conn.close()

        if row is None:
            return {
                "total_documents": 0,
                "documents_with_fingerprints": 0,
                "unique_smiles": 0,
                "unique_documents": 0,
                "coverage": 0.0,
            }

        return {
            "total_documents": row[0],
            "documents_with_fingerprints": row[1],
            "unique_smiles": row[2],
            "unique_documents": row[3],
            "coverage": row[1] / row[0] if row[0] > 0 else 0.0,
        }

    except Exception as e:
        logger.error(f"Statistics query failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/v1/structure/batch-search", tags=["Batch Search"])
async def batch_search(
    query_smiles_list: List[str],
    threshold: float = Query(default=0.7, ge=0.0, le=1.0),
    top_k: int = Query(default=5, ge=1, le=50),
):
    """
    Search for multiple query structures in parallel.

    Returns results for each query SMILES.
    Useful for comparing multiple compounds or scaffold hopping.
    """
    try:
        fp_generator = FingerprintGenerator()
        all_results = {}

        conn = psycopg2.connect(DATABASE_URL)
        cur = conn.cursor()

        for query_smiles in query_smiles_list:
            try:
                # Generate fingerprint
                query_fp = fp_generator.ecfp4(query_smiles)
                fp_vector = f"[{','.join(map(str, query_fp))}]"

                # Search
                cur.execute(
                    """
                    SELECT
                        id,
                        molecule_smiles,
                        1 - (molecule_fingerprint <=> %s::halfvec) as similarity,
                        metadata->>'title' as title
                    FROM document_chunk
                    WHERE molecule_fingerprint IS NOT NULL
                      AND 1 - (molecule_fingerprint <=> %s::halfvec) >= %s
                    ORDER BY similarity DESC
                    LIMIT %s
                """,
                    (fp_vector, fp_vector, threshold, top_k),
                )

                results = [
                    {
                        "doc_id": row[0],
                        "smiles": row[1],
                        "similarity": float(row[2]),
                        "document_title": row[3],
                    }
                    for row in cur.fetchall()
                ]

                all_results[query_smiles] = {
                    "query_smiles": query_smiles,
                    "results": results,
                    "total_results": len(results),
                }

            except Exception as e:
                all_results[query_smiles] = {
                    "query_smiles": query_smiles,
                    "error": str(e),
                    "results": [],
                }

        cur.close()
        conn.close()

        return {"queries": all_results}

    except Exception as e:
        logger.error(f"Batch search failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# =============================================================================
# Health Check
# =============================================================================


@app.get("/health", tags=["Health"])
async def health_check():
    """API health check."""
    try:
        conn = psycopg2.connect(DATABASE_URL)
        cur = conn.cursor()
        cur.execute("SELECT 1")
        cur.close()
        conn.close()

        return {"status": "healthy", "database": "connected"}
    except Exception as e:
        return {"status": "unhealthy", "database": "disconnected", "error": str(e)}


# =============================================================================
# Main
# =============================================================================


if __name__ == "__main__":
    import uvicorn

    # Run with: python -m api.structure_search
    uvicorn.run(
        app,
        host="0.0.0.0",
        port=8000,
        log_level="info",
    )
