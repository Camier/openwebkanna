"""
Structure Search API

FastAPI endpoints for molecular structure similarity search.
Uses fingerprint vectors and pgvector for sub-second similarity queries.

Phase 2, Task 2.3
"""

import os
import logging
from typing import Any, Dict, List, Optional
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel, Field
import psycopg2

from enrichers.fingerprint_generator import FingerprintGenerator
from retrieval.rrf_fusion import fuse_rankings_rrf

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Database URL from environment
DATABASE_URL = os.getenv(
    "DATABASE_URL",
    "postgresql://openwebui:openwebui_pgvector_pass_2026@localhost:5432/openwebui",
)


def _parse_enabled_fingerprint_types(raw: Optional[str]) -> tuple[str, ...]:
    """Parse SMILES fingerprint channels from env and keep supported values only."""
    allowed = {"ecfp4", "maccs"}
    if not raw:
        return ("ecfp4",)

    parsed: list[str] = []
    for item in raw.split(","):
        channel = item.strip().lower()
        if channel in allowed and channel not in parsed:
            parsed.append(channel)

    return tuple(parsed) if parsed else ("ecfp4",)


ENABLED_FINGERPRINT_TYPES = _parse_enabled_fingerprint_types(
    os.getenv("SMILES_FINGERPRINT_TYPES", "ecfp4")
)
ENABLED_FINGERPRINT_TYPES_LABEL = ", ".join(ENABLED_FINGERPRINT_TYPES)


def _validate_fingerprint_channel(channel: str) -> str:
    """Validate and normalize fingerprint channel against runtime-enabled channels."""
    normalized = (channel or "").strip().lower()
    if normalized not in ENABLED_FINGERPRINT_TYPES:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Invalid fingerprint_type '{channel}'. "
                f"Enabled values: {ENABLED_FINGERPRINT_TYPES_LABEL}"
            ),
        )
    return normalized


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
    description=(
        "Molecular structure similarity search using fingerprint vectors "
        f"(enabled: {ENABLED_FINGERPRINT_TYPES_LABEL})"
    ),
    version="1.0.0",
    lifespan=lifespan,
)


# =============================================================================
# Schema Helpers
# =============================================================================


def _get_document_chunk_columns(cur) -> set[str]:
    """Return available columns on `document_chunk` for runtime-safe SQL."""
    cur.execute(
        """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = 'public'
          AND table_name = 'document_chunk'
    """
    )
    return {row[0] for row in cur.fetchall()}


def _title_expr(columns: set[str]) -> str:
    """Return SQL expression for document title across schema variants."""
    if "vmetadata" in columns:
        return "vmetadata->>'title'"
    if "metadata" in columns:
        return "metadata->>'title'"
    return "NULL"


def _resolve_fingerprint_column(channel: str, columns: set[str]) -> str:
    """Resolve requested fingerprint channel to an existing DB column."""
    channel = _validate_fingerprint_channel(channel)
    if channel == "ecfp4":
        if "molecule_fingerprint" not in columns:
            raise HTTPException(
                status_code=503,
                detail=(
                    "ECFP4 fingerprint column 'molecule_fingerprint' is missing in "
                    "document_chunk. Apply SMILES DB migrations first."
                ),
            )
        return "molecule_fingerprint"

    if channel == "maccs":
        if "molecule_fingerprint_maccs" not in columns:
            raise HTTPException(
                status_code=400,
                detail=(
                    "MACCS channel unavailable: column "
                    "'molecule_fingerprint_maccs' is missing in document_chunk."
                ),
            )
        return "molecule_fingerprint_maccs"

    # Guard rail: should never happen because _validate_fingerprint_channel
    # enforces runtime-enabled values.
    raise HTTPException(status_code=500, detail="Unhandled fingerprint_type mapping")


# =============================================================================
# Request/Response Models
# =============================================================================


class StructureSearchRequest(BaseModel):
    """Request model for structure similarity search."""

    query_smiles: str = Field(..., description="Query SMILES string")
    fingerprint_type: str = Field(
        default="ecfp4",
        description=(
            "Fingerprint channel enabled at runtime via "
            f"SMILES_FINGERPRINT_TYPES (current: {ENABLED_FINGERPRINT_TYPES_LABEL})"
        ),
    )
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
    fingerprint_type: str = Field(..., description="Similarity channel used")
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


class RankedDoc(BaseModel):
    """Ranked retrieval hit from a single channel."""

    doc_id: str = Field(..., description="Document ID")
    score: Optional[float] = Field(None, description="Channel-native score")
    metadata: Optional[Dict[str, Any]] = Field(
        None, description="Optional hit metadata"
    )


class RetrievalFuseRequest(BaseModel):
    """Request model for multi-channel RRF fusion."""

    text_dense_results: List[RankedDoc] = Field(
        default_factory=list, description="Dense text retrieval ranking"
    )
    bm25_results: List[RankedDoc] = Field(
        default_factory=list, description="Sparse BM25 retrieval ranking"
    )
    smiles_fingerprint_results: List[RankedDoc] = Field(
        default_factory=list, description="SMILES fingerprint retrieval ranking"
    )
    smiles_embedding_results: List[RankedDoc] = Field(
        default_factory=list,
        description="Optional molecular embedding retrieval ranking",
    )
    include_smiles_embedding: bool = Field(
        default=False,
        description="Include smiles_embedding_results in fusion",
    )
    top_k: int = Field(default=20, ge=1, le=200)
    rrf_k: int = Field(default=60, ge=1, le=1000)
    channel_weights: Optional[Dict[str, float]] = Field(
        default=None,
        description="Optional channel weights (e.g. {'text_dense': 1.0})",
    )


class RetrievalFuseResponse(BaseModel):
    """RRF fusion response payload."""

    total_results: int
    included_channels: List[str]
    results: List[Dict[str, Any]]


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

    Uses enabled fingerprint channel (`SMILES_FINGERPRINT_TYPES`) with cosine
    similarity search via
    pgvector IVFFlat index. Returns molecules with similarity >= threshold, sorted
    by similarity descending.

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
        channel = _validate_fingerprint_channel(request.fingerprint_type)

        # Step 1: Generate query fingerprint
        fp_generator = FingerprintGenerator()
        if channel == "maccs":
            query_fp = fp_generator.maccs(request.query_smiles)
        else:
            query_fp = fp_generator.ecfp4(request.query_smiles)

        # Convert to pgvector format: "[0.0,1.0,0.0,...]"
        fp_vector = f"[{','.join(map(str, query_fp))}]"

        # Step 2: Execute similarity search
        conn = psycopg2.connect(DATABASE_URL)
        cur = conn.cursor()
        columns = _get_document_chunk_columns(cur)
        fp_column = _resolve_fingerprint_column(channel, columns)
        title_column_expr = _title_expr(columns)

        cur.execute(
            f"""
            SELECT
                id,
                molecule_smiles,
                1 - ({fp_column} <=> %s::halfvec) as similarity,
                {title_column_expr} as title,
                molecule_metadata
            FROM document_chunk
            WHERE {fp_column} IS NOT NULL
              AND 1 - ({fp_column} <=> %s::halfvec) >= %s
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
                    fingerprint_type=channel,
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

    except HTTPException:
        raise
    except ValueError as e:
        # Invalid SMILES
        raise HTTPException(status_code=400, detail=f"Invalid SMILES: {str(e)}")
    except psycopg2.Error as e:
        logger.error(f"Database error: {e}")
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        logger.error(f"Search failed: {e}")
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
        columns = _get_document_chunk_columns(cur)
        has_ecfp4 = "molecule_fingerprint" in columns
        has_maccs = "molecule_fingerprint_maccs" in columns

        any_fps_predicates = []
        if has_ecfp4:
            any_fps_predicates.append("molecule_fingerprint IS NOT NULL")
        if has_maccs:
            any_fps_predicates.append("molecule_fingerprint_maccs IS NOT NULL")

        any_fps_filter = " OR ".join(any_fps_predicates) if any_fps_predicates else "FALSE"
        ecfp4_count_expr = "COUNT(molecule_fingerprint)" if has_ecfp4 else "0"
        maccs_count_expr = "COUNT(molecule_fingerprint_maccs)" if has_maccs else "0"

        cur.execute(
            f"""
            SELECT
                COUNT(*) as total_docs,
                COUNT(*) FILTER (
                    WHERE {any_fps_filter}
                ) as docs_with_any_fps,
                {ecfp4_count_expr} as docs_with_ecfp4,
                {maccs_count_expr} as docs_with_maccs,
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
                "documents_with_ecfp4": 0,
                "documents_with_maccs": 0,
                "unique_smiles": 0,
                "unique_documents": 0,
                "coverage": 0.0,
            }

        return {
            "total_documents": row[0],
            "documents_with_fingerprints": row[1],
            "documents_with_ecfp4": row[2],
            "documents_with_maccs": row[3],
            "unique_smiles": row[4],
            "unique_documents": row[5],
            "coverage": row[1] / row[0] if row[0] > 0 else 0.0,
            "coverage_ecfp4": row[2] / row[0] if row[0] > 0 else 0.0,
            "coverage_maccs": row[3] / row[0] if row[0] > 0 else 0.0,
            "enabled_fingerprint_types": list(ENABLED_FINGERPRINT_TYPES),
        }

    except Exception as e:
        logger.error(f"Statistics query failed: {e}")
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
        columns = _get_document_chunk_columns(cur)
        title_column_expr = _title_expr(columns)

        cur.execute(
            f"""
            SELECT
                id,
                molecule_smiles,
                {title_column_expr} as title,
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


@app.post("/v1/structure/batch-search", tags=["Batch Search"])
async def batch_search(
    query_smiles_list: List[str],
    fingerprint_type: str = Query(default="ecfp4"),
    threshold: float = Query(default=0.7, ge=0.0, le=1.0),
    top_k: int = Query(default=5, ge=1, le=50),
):
    """
    Search for multiple query structures in parallel.

    Returns results for each query SMILES.
    Useful for comparing multiple compounds or scaffold hopping.
    """
    try:
        channel = _validate_fingerprint_channel(fingerprint_type)
        fp_generator = FingerprintGenerator()
        all_results = {}

        conn = psycopg2.connect(DATABASE_URL)
        cur = conn.cursor()
        columns = _get_document_chunk_columns(cur)
        title_column_expr = _title_expr(columns)
        fp_column = _resolve_fingerprint_column(channel, columns)

        for query_smiles in query_smiles_list:
            try:
                # Generate fingerprint
                if channel == "maccs":
                    query_fp = fp_generator.maccs(query_smiles)
                else:
                    query_fp = fp_generator.ecfp4(query_smiles)
                fp_vector = f"[{','.join(map(str, query_fp))}]"

                # Search
                cur.execute(
                    f"""
                    SELECT
                        id,
                        molecule_smiles,
                        1 - ({fp_column} <=> %s::halfvec) as similarity,
                        {title_column_expr} as title
                    FROM document_chunk
                    WHERE {fp_column} IS NOT NULL
                      AND 1 - ({fp_column} <=> %s::halfvec) >= %s
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
                        "fingerprint_type": channel,
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

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Batch search failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post(
    "/v1/retrieval/fuse",
    response_model=RetrievalFuseResponse,
    tags=["Retrieval Fusion"],
)
async def fuse_retrieval_rankings(request: RetrievalFuseRequest):
    """
    Fuse heterogeneous retrieval rankings with Reciprocal Rank Fusion (RRF).

    Intended channels:
    - text_dense
    - bm25
    - smiles_fingerprint
    - smiles_embedding (optional)
    """
    def _dump(item: RankedDoc) -> Dict[str, Any]:
        if hasattr(item, "model_dump"):
            return item.model_dump()
        return item.dict()

    channel_rankings: Dict[str, List[Dict[str, Any]]] = {
        "text_dense": [_dump(item) for item in request.text_dense_results],
        "bm25": [_dump(item) for item in request.bm25_results],
        "smiles_fingerprint": [_dump(item) for item in request.smiles_fingerprint_results],
    }
    if request.include_smiles_embedding:
        channel_rankings["smiles_embedding"] = [
            _dump(item) for item in request.smiles_embedding_results
        ]

    fused = fuse_rankings_rrf(
        channel_rankings,
        top_k=request.top_k,
        rrf_k=request.rrf_k,
        channel_weights=request.channel_weights,
    )

    return RetrievalFuseResponse(
        total_results=len(fused),
        included_channels=list(channel_rankings.keys()),
        results=fused,
    )


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
