"""
Domain-specific validation for ethnopharmacology (Level 3).

Validates molecules against Sceletium alkaloid characteristics and
compares to authenticated gold standard references.

References:
- PubChem Sceletium alkaloids
- Gericke (2008): Sceletium botanical review
"""

from typing import Dict, Any, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
import json
from pathlib import Path


class DomainValidator:
    """
    Level 3 domain validator for Sceletium alkaloids.

    Validates against:
    1. Alkaloid requirements (nitrogen content)
    2. Sceletium-specific features (methoxy groups, MW range)
    3. Gold standard authenticated alkaloids
    """

    def __init__(self, gold_standards_file: str = None):
        """
        Initialize domain validator.

        Args:
            gold_standards_file: Path to gold_standards.json
        """
        self.gold_standards = self._load_gold_standards(gold_standards_file)

    def _load_gold_standards(self, filepath: str) -> List[Dict[str, Any]]:
        """Load gold standard reference molecules."""
        if filepath is None:
            filepath = Path(__file__).parent.parent / "config" / "gold_standards.json"

        try:
            with open(filepath) as f:
                data = json.load(f)
                return data.get("gold_standards", [])
        except FileNotFoundError:
            return []

    def validate(self, smiles: str, canonical_smiles: str = None) -> Dict[str, Any]:
        """
        Perform Level 3 domain validation.

        Args:
            smiles: SMILES string
            canonical_smiles: Pre-computed canonical SMILES

        Returns:
            Validation results with compound classification and gold standard matches
        """
        result = {
            "is_valid": True,
            "compound_class": "unknown",
            "warnings": [],
            "gold_standard_matches": [],
            "matches_sceletium_scaffold": False,
        }

        mol = Chem.MolFromSmiles(canonical_smiles or smiles)
        if mol is None:
            result["is_valid"] = False
            return result

        # Check alkaloid requirements
        alkaloid_result = self._check_alkaloid(mol)
        result["compound_class"] = alkaloid_result["compound_class"]

        # Check Sceletium-specific features
        sceletium_result = self._check_sceletium_features(mol)
        result["warnings"].extend(sceletium_result.get("warnings", []))
        result["matches_sceletium_scaffold"] = sceletium_result[
            "matches_sceletium_scaffold"
        ]
        result["methoxy_count"] = sceletium_result["methoxy_count"]

        # Compare to gold standards
        matches = self._compare_to_gold_standards(mol)
        result["gold_standard_matches"] = matches

        return result

    def _check_alkaloid(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Check if molecule qualifies as alkaloid.

        Alkaloids must contain nitrogen.
        """
        result = {"compound_class": "unknown"}

        nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "N")

        if nitrogen_count == 0:
            result["compound_class"] = "non_alkaloid"
            result["is_alkaloid"] = False
            result["nitrogen_count"] = 0
        else:
            result["compound_class"] = "alkaloid"
            result["is_alkaloid"] = True
            result["nitrogen_count"] = nitrogen_count

        return result

    def _check_sceletium_features(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Check for Sceletium alkaloid hallmarks.

        Features:
        - 3,4-dimethoxy phenyl pattern
        - MW in 250-400 Da range
        """
        result = {
            "matches_sceletium_scaffold": False,
            "methoxy_count": 0,
            "warnings": [],
        }

        # Count methoxy groups (COc pattern)
        methoxy_smarts = "COc"
        methoxy_pattern = Chem.MolFromSmarts(methoxy_smarts)

        if methoxy_pattern:
            matches = mol.GetSubstructMatches(methoxy_pattern)
            result["methoxy_count"] = len(matches)

            if len(matches) == 0:
                result["warnings"].append("no_methoxy_groups")
            elif len(matches) >= 2:
                result["matches_sceletium_scaffold"] = True

        # Check MW range for mesembrine-type
        mw = Descriptors.MolWt(mol)
        if 250 <= mw <= 400:
            result["mw_in_range"] = True
        else:
            result["mw_in_range"] = False
            result["warnings"].append(f"mw_outside_sceletium_range: {mw:.1f}")

        return result

    def _compare_to_gold_standards(self, mol: Chem.Mol) -> List[Dict[str, Any]]:
        """
        Compare molecule to authenticated gold standard alkaloids.

        Uses Tanimoto similarity on ECFP4 fingerprints.
        """
        if not self.gold_standards:
            return []

        results = []
        mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

        for ref in self.gold_standards:
            ref_mol = Chem.MolFromSmiles(ref["smiles_isomeric"])
            if ref_mol is None:
                continue

            ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)

            similarity = DataStructs.TanimotoSimilarity(mol_fp, ref_fp)

            if similarity >= 0.7:  # Threshold for "similar"
                results.append(
                    {
                        "id": ref["id"],
                        "similarity": similarity / 1.0,
                        "match_type": ("exact" if similarity >= 0.99 else "similar"),
                        "cid": ref.get("pubchem_cid"),
                    }
                )

        # Sort by similarity descending
        results.sort(key=lambda x: x["similarity"], reverse=True)

        return results

    def get_compound_summary(
        self, smiles: str, validation_results: Dict[str, Any]
    ) -> str:
        """
        Generate human-readable compound summary.

        Args:
            smiles: Original SMILES
            validation_results: Results from validate()

        Returns:
            Summary string
        """
        parts = []

        # Compound class
        compound_class = validation_results.get("compound_class", "unknown")
        parts.append(f"**Class**: {compound_class}")

        # Gold standard matches
        matches = validation_results.get("gold_standard_matches", [])
        if matches:
            parts.append("**Matches**:")
            for match in matches[:3]:  # Top 3
                parts.append(
                    f"  - {match['id']}: {match['match_type']} "
                    f"(Tanimoto: {match['similarity']:.2f})"
                )

        # Sceletium features
        if validation_results.get("matches_sceletium_scaffold"):
            parts.append("**Sceletium scaffold**: Yes")

        return "\n".join(parts)


def validate_ethnopharmacology_domain(
    smiles: str, canonical_smiles: str = None
) -> Dict[str, Any]:
    """
    Convenience function for quick domain validation.

    Args:
        smiles: SMILES string
        canonical_smiles: Pre-computed canonical SMILES

    Returns:
        Validation result dictionary
    """
    validator = DomainValidator()
    return validator.validate(smiles, canonical_smiles)


if __name__ == "__main__":
    # Test examples
    test_cases = [
        (
            "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC",
            True,
        ),  # Mesembrine
        (
            "CN1CCC2(C1)C=CC(=O)C2c3ccc(OC)c(OC)c3",
            True,
        ),  # Mesembrenone
        ("CCO", False),  # Ethanol (not alkaloid)
    ]

    validator = DomainValidator()
    print("Domain Validator Test Results")
    print("=" * 60)

    for smiles, expected_alkaloid in test_cases:
        result = validator.validate(smiles)
        is_alkaloid = result["compound_class"] == "alkaloid"
        status = "✓" if is_alkaloid == expected_alkaloid else "✗"
        print(f"{status} {smiles[:40]:<40} -> {result['compound_class']}")

        if result["gold_standard_matches"]:
            for match in result["gold_standard_matches"]:
                print(f"   Match: {match['id']} ({match['similarity']:.2f})")
