"""
Chemical plausibility validation (Level 2).

Uses RDKit to validate molecular properties and filter extraction artifacts.
Implements standard cheminformatics sanity checks used in OCSR pipelines.

References:
- RDKit Documentation: https://www.rdkit.org/docs/
- MolVS: https://molvs.readthedocs.io/
- Krasnov et al. (2024): Digital Discovery, DOI: 10.1039/D3DD00228D
"""

from typing import Dict, Any, List
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import re


class ChemicalValidator:
    """
    Level 2 chemical plausibility validator.

    Filters extraction artifacts using molecular property thresholds
    based on natural product alkaloid characteristics.
    """

    def __init__(self, config: Dict[str, Any] = None):
        """
        Initialize validator with configuration.

        Args:
            config: Validation thresholds (uses defaults if not provided)
        """
        self.config = config or self._default_config()

    def _default_config(self) -> Dict[str, Any]:
        """Default validation thresholds."""
        return {
            "molecular_weight": {
                "min": 100,
                "max": 1000,
                "warn_min": 150,
                "warn_max": 800,
            },
            "atom_count": {
                "min": 5,
                "max": 200,
            },
            "fragments": {
                "max": 5,
                "warn_max": 3,
            },
            "ring_systems": {
                "max": 10,
            },
            "charge": {
                "max_absolute": 3,
            },
        }

    def validate(self, smiles: str, canonical_smiles: str = None) -> Dict[str, Any]:
        """
        Perform Level 2 chemical plausibility validation.

        Args:
            smiles: SMILES string (canonical preferred)
            canonical_smiles: Pre-computed canonical SMILES (optional)

        Returns:
            Dictionary with validation results and computed properties
        """
        result = {
            "is_valid": True,
            "warnings": [],
            "rejection_reason": None,
            "properties": {},
        }

        # Parse SMILES with RDKit
        mol = Chem.MolFromSmiles(canonical_smiles or smiles)

        if mol is None:
            result["is_valid"] = False
            result["rejection_reason"] = "rdkit_parse_failure"
            return result

        # Compute properties
        result["properties"] = self._compute_properties(mol, smiles)
        props = result["properties"]

        # Apply validation rules
        check_results = [
            self._check_molecular_weight(props),
            self._check_atom_count(mol, props),
            self._check_fragments(mol, props),
            self._check_garbage_patterns(smiles),
            self._check_rings(mol, props),
            self._check_charge(mol, props),
        ]

        for check_result in check_results:
            if not check_result["is_valid"]:
                result["is_valid"] = False
                result["rejection_reason"] = check_result["reason"]
                return result

            result["warnings"].extend(check_result.get("warnings", []))

        return result

    def _compute_properties(self, mol: Chem.Mol, smiles: str) -> Dict[str, Any]:
        """Compute molecular properties."""
        return {
            "molecular_weight": Descriptors.MolWt(mol),
            "exact_mass": Descriptors.ExactMolWt(mol),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "num_fragments": len(Chem.GetMolFrags(mol, asMols=True)),
            "num_rings": mol.GetRingInfo().NumRings(),
            "num_aromatic_rings": sum(
                1
                for ring in mol.GetRingInfo().AtomRings()
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            ),
            "formal_charge": Chem.GetFormalCharge(mol),
            "num_hba": rdMolDescriptors.CalcNumHBA(mol),
            "num_hbd": rdMolDescriptors.CalcNumHBD(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
        }

    def _check_molecular_weight(self, props: Dict[str, Any]) -> Dict[str, Any]:
        """Check molecular weight thresholds."""
        mw = props["molecular_weight"]
        mw_config = self.config["molecular_weight"]

        if mw < mw_config["min"]:
            return {
                "is_valid": False,
                "reason": f"mw_too_low: {mw:.1f} Da",
            }

        if mw > mw_config["max"]:
            return {
                "is_valid": False,
                "reason": f"mw_too_high: {mw:.1f} Da",
            }

        warnings = []
        if mw < mw_config["warn_min"]:
            warnings.append(f"low_mw: {mw:.1f} Da")
        if mw > mw_config["warn_max"]:
            warnings.append(f"high_mw: {mw:.1f} Da")

        return {"is_valid": True, "warnings": warnings}

    def _check_atom_count(self, mol: Chem.Mol, props: Dict[str, Any]) -> Dict[str, Any]:
        """Check heavy atom count thresholds."""
        num_atoms = props["num_heavy_atoms"]
        atom_config = self.config["atom_count"]

        if num_atoms < atom_config["min"]:
            return {
                "is_valid": False,
                "reason": f"too_few_atoms: {num_atoms}",
            }

        if num_atoms > atom_config["max"]:
            return {
                "is_valid": False,
                "reason": f"too_many_atoms: {num_atoms}",
            }

        return {"is_valid": True}

    def _check_fragments(self, mol: Chem.Mol, props: Dict[str, Any]) -> Dict[str, Any]:
        """Check fragment count (disconnected components)."""
        num_frags = props["num_fragments"]
        frag_config = self.config["fragments"]

        if num_frags > frag_config["max"]:
            return {
                "is_valid": False,
                "reason": f"too_many_fragments: {num_frags}",
            }

        warnings = []
        if num_frags > frag_config["warn_max"]:
            warnings.append(f"multiple_fragments: {num_frags}")

        return {"is_valid": True, "warnings": warnings}

    def _check_garbage_patterns(self, smiles: str) -> Dict[str, Any]:
        """Detect common OCSR garbage patterns."""
        smiles_clean = smiles.replace(".", "").upper()

        # Carbon-only garbage (C.C.C patterns)
        if smiles_clean and all(c in "C0123456789" for c in smiles_clean):
            return {
                "is_valid": False,
                "reason": "carbon_only_garbage",
            }

        # Single element only (O.O.O, N.N.N)
        for element in ["O", "N", "S", "P", "F", "CL", "BR", "I"]:
            pattern = f"^{re.escape(element)}[0-9]*(\\.{re.escape(element)}[0-9]*)*$"
            if re.match(pattern, smiles.upper()):
                return {
                    "is_valid": False,
                    "reason": f"single_element_only: {element}",
                }

        # Repetitive carbonyl garbage C(=O)(=O)(=O)
        if "(=O)(=O)" in smiles:
            return {
                "is_valid": False,
                "reason": "impossible_bonding",
            }

        return {"is_valid": True}

    def _check_rings(self, mol: Chem.Mol, props: Dict[str, Any]) -> Dict[str, Any]:
        """Check ring system complexity."""
        num_rings = props["num_rings"]
        ring_config = self.config["ring_systems"]

        if num_rings > ring_config["max"]:
            return {
                "is_valid": False,
                "reason": f"too_many_rings: {num_rings}",
            }

        return {"is_valid": True}

    def _check_charge(self, mol: Chem.Mol, props: Dict[str, Any]) -> Dict[str, Any]:
        """Check formal charge."""
        charge = props["formal_charge"]
        charge_config = self.config["charge"]

        if abs(charge) > charge_config["max_absolute"]:
            return {
                "is_valid": False,
                "reason": f"excessive_charge: {charge}",
            }

        return {"is_valid": True}


def validate_chemical_plausibility(
    smiles: str, canonical_smiles: str = None
) -> Dict[str, Any]:
    """
    Convenience function for quick chemical validation.

    Args:
        smiles: SMILES string to validate
        canonical_smiles: Pre-computed canonical SMILES (optional)

    Returns:
        Validation result dictionary
    """
    validator = ChemicalValidator()
    return validator.validate(smiles, canonical_smiles)


if __name__ == "__main__":
    # Test examples
    test_cases = [
        # Valid natural products
        (
            "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC",
            True,
        ),  # Mesembrine (MW 289)
        (
            "CN1CCC2(C1)C=CC(=O)C2c3ccc(OC)c(OC)c3",
            True,
        ),  # Mesembrenone
        # Invalid
        ("C.C.C.C.C.C.C.C.C.C", False),  # Carbon garbage
        ("C1CO1", False),  # Too small (MW ~58)
    ]

    validator = ChemicalValidator()
    print("Chemical Validator Test Results")
    print("=" * 60)

    for smiles, expected_valid in test_cases:
        result = validator.validate(smiles)
        status = "✓" if result["is_valid"] == expected_valid else "✗"
        print(f"{status} {smiles[:40]:<40} -> {result['is_valid']}")
        if result["properties"]:
            props = result["properties"]
            print(
                f"   MW: {props['molecular_weight']:.1f}, "
                f"Atoms: {props['num_heavy_atoms']}, "
                f"Frags: {props['num_fragments']}"
            )
        if result.get("rejection_reason"):
            print(f"   Rejected: {result['rejection_reason']}")
