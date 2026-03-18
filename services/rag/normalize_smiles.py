"""RDKit-backed SMILES normalization for exact chemistry retrieval.

This module is the chemistry identity gate for the one-collection RAG.
Exact molecule retrieval should only use normalized outputs from here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class NormalizedSmiles:
    """Normalized chemistry identity derived from a SMILES string."""

    raw_smiles: str
    valid: bool
    canonical_smiles: str | None
    inchikey: str | None
    review_status: str
    error: str | None = None

    @property
    def is_exact_matchable(self) -> bool:
        """Return True when this molecule can participate in exact retrieval."""
        return self.valid and bool(self.canonical_smiles) and bool(self.inchikey)


def normalize_smiles(raw_smiles: str) -> NormalizedSmiles:
    """Parse and canonicalize a SMILES string with RDKit.

    Rules:
    - returns a structured invalid result instead of raising on parse failure
    - computes canonical isomeric SMILES for exact identity
    - computes an InChIKey for stable equality filtering
    """

    raw_smiles = (raw_smiles or "").strip()
    if not raw_smiles:
        return NormalizedSmiles(
            raw_smiles=raw_smiles,
            valid=False,
            canonical_smiles=None,
            inchikey=None,
            review_status="invalid",
            error="empty_smiles",
        )

    try:
        from rdkit import Chem
        from rdkit.Chem import inchi
    except Exception as exc:  # pragma: no cover - environment dependent
        return NormalizedSmiles(
            raw_smiles=raw_smiles,
            valid=False,
            canonical_smiles=None,
            inchikey=None,
            review_status="unavailable",
            error=f"rdkit_unavailable:{exc}",
        )

    try:
        mol = Chem.MolFromSmiles(raw_smiles)
    except Exception as exc:  # pragma: no cover - rdkit edge case
        return NormalizedSmiles(
            raw_smiles=raw_smiles,
            valid=False,
            canonical_smiles=None,
            inchikey=None,
            review_status="invalid",
            error=f"parse_error:{exc}",
        )

    if mol is None:
        return NormalizedSmiles(
            raw_smiles=raw_smiles,
            valid=False,
            canonical_smiles=None,
            inchikey=None,
            review_status="invalid",
            error="parse_failed",
        )

    try:
        canonical_smiles = Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=True,
        )
        inchikey = inchi.MolToInchiKey(mol)
    except Exception as exc:  # pragma: no cover - rdkit edge case
        return NormalizedSmiles(
            raw_smiles=raw_smiles,
            valid=False,
            canonical_smiles=None,
            inchikey=None,
            review_status="invalid",
            error=f"normalization_error:{exc}",
        )

    return NormalizedSmiles(
        raw_smiles=raw_smiles,
        valid=True,
        canonical_smiles=canonical_smiles,
        inchikey=inchikey,
        review_status="parsed",
        error=None,
    )


def to_payload_dict(normalized: NormalizedSmiles) -> dict[str, Any]:
    """Convert a normalized record into a payload-ready dictionary."""
    return {
        "raw_smiles": normalized.raw_smiles,
        "canonical_smiles": normalized.canonical_smiles,
        "inchikey": normalized.inchikey,
        "review_status": normalized.review_status,
    }
