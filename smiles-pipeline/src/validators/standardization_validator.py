"""
Deterministic molecular standardization (post-syntax stage).

Applies RDKit MolStandardize cleanup in a fixed sequence so downstream
validation, deduplication, and persistence are reproducible.
"""

from typing import Any, Dict

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


class StandardizationValidator:
    """
    Standardize SMILES with a deterministic RDKit policy.

    Policy sequence:
    1. Cleanup
    2. Fragment parent extraction
    3. Uncharge
    4. Canonical tautomer selection
    """

    POLICY_VERSION = "rdkit_molstandardize_v1"

    def __init__(self):
        self._uncharger = rdMolStandardize.Uncharger()
        self._tautomer_enumerator = rdMolStandardize.TautomerEnumerator()

    def standardize(self, smiles: str) -> Dict[str, Any]:
        """
        Standardize a SMILES string and derive InChI identifiers.

        Args:
            smiles: Parsed/canonical SMILES from Level 1

        Returns:
            {
                "success": bool,
                "standardized_smiles": str | None,
                "standard_inchi": str | None,
                "standard_inchikey": str | None,
                "policy_version": str,
                "error": str | None,
                "warnings": list[str],
            }
        """
        result = {
            "success": False,
            "standardized_smiles": None,
            "standard_inchi": None,
            "standard_inchikey": None,
            "policy_version": self.POLICY_VERSION,
            "error": None,
            "warnings": [],
        }

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "rdkit_parse_failure"
            return result

        try:
            cleaned = rdMolStandardize.Cleanup(mol)
            parent = rdMolStandardize.FragmentParent(cleaned)
            uncharged = self._uncharger.uncharge(parent)
            standardized = self._tautomer_enumerator.Canonicalize(uncharged)
        except Exception as exc:
            result["error"] = f"standardization_failed: {exc}"
            return result

        try:
            result["standardized_smiles"] = Chem.MolToSmiles(
                standardized,
                canonical=True,
                isomericSmiles=True,
            )
            result["success"] = True
        except Exception as exc:
            result["error"] = f"smiles_serialization_failed: {exc}"
            return result

        # InChI generation is non-fatal for standardization success.
        try:
            result["standard_inchi"] = Chem.MolToInchi(standardized)
            result["standard_inchikey"] = Chem.MolToInchiKey(standardized)
        except Exception as exc:
            result["warnings"].append(f"inchi_generation_failed: {exc}")

        return result
