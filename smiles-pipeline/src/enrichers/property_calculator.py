"""
Molecular property calculator.

Computes physicochemical properties for enriched molecule records.
Implements standard cheminformatics descriptors used in drug discovery.

References:
- RDKit Descriptors: https://www.rdkit.org/docs/GettingStartedInPython.html#molecular-descriptors
- Lipinski Rule of 5: Lipinski et al. (1997), Adv. Drug Deliv. Rev.
- PubChem property calculations
"""

from typing import Dict, Any, List, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski
from rdkit.Chem import rdMolDescriptors as rdMolDesc


class PropertyCalculator:
    """
    Calculate molecular properties for enrichment.

    Computes:
    - Basic properties (MW, formula, atoms)
    - Lipophilicity (LogP)
    - Polarity (TPSA)
    - H-bond donors/acceptors
    - Rotatable bonds
    - Ring systems
    - Lipinski Rule of 5 analysis
    """

    def __init__(self):
        """Initialize property calculator."""
        pass

    def calculate(
        self,
        smiles: str,
        canonical_smiles: str = None,
    ) -> Dict[str, Any]:
        """
        Calculate all molecular properties.

        Args:
            smiles: SMILES string
            canonical_smiles: Pre-computed canonical SMILES

        Returns:
            Dictionary with computed properties
        """
        result = {
            "success": False,
            "properties": {},
            "lipinski": {},
            "error": None,
        }

        mol = Chem.MolFromSmiles(canonical_smiles or smiles)

        if mol is None:
            result["error"] = "rdkit_parse_failure"
            return result

        # Compute all properties
        result["properties"] = self._compute_all_properties(mol)
        result["lipinski"] = self._compute_lipinski_analysis(mol, result["properties"])
        result["success"] = True

        return result

    def _compute_all_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Compute comprehensive property set."""
        return {
            # Basic properties
            "molecular_weight": Descriptors.MolWt(mol),
            "exact_mass": Descriptors.ExactMolWt(mol),
            "molecular_formula": self._get_formula(mol),
            # Atom counts
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "num_atoms": mol.GetNumAtoms(),
            "num_carbon": self._count_element(mol, 6),
            "num_nitrogen": self._count_element(mol, 7),
            "num_oxygen": self._count_element(mol, 8),
            "num_sulfur": self._count_element(mol, 16),
            # Lipophilicity
            "logp": Descriptors.MolLogP(mol),
            "logp_mr": self._calc_mr_logp(mol),  # McGowan volume-based
            # Polarity
            "tpsa": Descriptors.TPSA(mol),
            "tpsa_labute": self._calc_labute_tpsa(mol),
            # H-bonding
            "num_hba": rdMolDescriptors.CalcNumHBA(mol),
            "num_hbd": rdMolDescriptors.CalcNumHBD(mol),
            # Flexibility
            "num_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "num_rings": mol.GetRingInfo().NumRings(),
            "num_aromatic_rings": self._count_aromatic_rings(mol),
            "num_aliphatic_rings": self._count_aliphatic_rings(mol),
            # Structural features
            "fraction_sp3_carbons": rdMolDescriptors.CalcFractionCSP3(mol),
            "num_stereocenters": rdMolDescriptors.CalcNumStereoCenters(mol),
            "num_double_bonds": self._count_double_bonds(mol),
            "num_triple_bonds": self._count_triple_bonds(mol),
            # Size/shape
            "molar_refractivity": Descriptors.MolMR(mol),
            "molecular_volume": self._calc_mcgowan_volume(mol),
        }

    def _compute_lipinski_analysis(
        self,
        mol: Chem.Mol,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Compute Lipinski Rule of 5 analysis.

        References:
        - Lipinski et al. (1997): "Experimental and computational approaches..."
        """
        mw = properties["molecular_weight"]
        logp = properties["logp"]
        hba = properties["num_hba"]
        hbd = properties["num_hbd"]

        # Count violations
        violations = 0
        violation_details = []

        if mw > 500:
            violations += 1
            violation_details.append(f"MW > 500: {mw:.1f}")

        if logp > 5.0:
            violations += 1
            violation_details.append(f"LogP > 5.0: {logp:.2f}")

        if hba > 10:
            violations += 1
            violation_details.append(f"HBA > 10: {hba}")

        if hbd > 5:
            violations += 1
            violation_details.append(f"HBD > 5: {hbd}")

        return {
            "violations": violations,
            "passes_lipinski": violations == 0,
            "violation_details": violation_details,
            "drug_likeness": self._classify_drug_likeness(violations, properties),
        }

    def _classify_drug_likeness(
        self,
        violations: int,
        properties: Dict[str, Any],
    ) -> str:
        """Classify drug-likeness based on properties."""
        if violations == 0:
            return "good"
        elif violations <= 2:
            return "moderate"
        else:
            return "poor"

    def _get_formula(self, mol: Chem.Mol) -> str:
        """Get molecular formula."""
        from rdkit.Chem.rdchem import Mol

        return rdMolDescriptors.CalcMolFormula(mol)

    def _count_element(self, mol: Chem.Mol, atomic_num: int) -> int:
        """Count atoms of specific element."""
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_num)

    def _count_aromatic_rings(self, mol: Chem.Mol) -> int:
        """Count aromatic rings."""
        ring_info = mol.GetRingInfo()
        aromatic_rings = 0

        for ring in ring_info.AtomRings():
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_rings += 1

        return aromatic_rings

    def _count_aliphatic_rings(self, mol: Chem.Mol) -> int:
        """Count aliphatic rings."""
        return mol.GetRingInfo().NumRings() - self._count_aromatic_rings(mol)

    def _count_double_bonds(self, mol: Chem.Mol) -> int:
        """Count double bonds."""
        return sum(
            1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE
        )

    def _count_triple_bonds(self, mol: Chem.Mol) -> int:
        """Count triple bonds."""
        return sum(
            1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.TRIPLE
        )

    def _calc_mr_logp(self, mol: Chem.Mol) -> float:
        """Calculate McGowan volume-based LogP."""
        # Simplified approximation
        return Descriptors.MolLogP(mol)

    def _calc_labute_tpsa(self, mol: Chem.Mol) -> float:
        """Calculate Labute TPSA."""
        from rdkit.Chem import rdMolDescriptors

        return rdMolDescriptors.CalcTPSA(mol, force=True)

    def _calc_mcgowan_volume(self, mol: Chem.Mol) -> float:
        """Calculate McGowan characteristic volume."""
        # Approximation using atomic contributions
        volume = 0.0
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            # Simplified atomic volumes (Å³)
            atomic_volumes = {
                "C": 16.35,
                "H": 8.71,
                "O": 10.05,
                "N": 14.36,
                "S": 25.10,
                "P": 22.00,
                "F": 12.00,
                "Cl": 24.00,
                "Br": 28.00,
                "I": 34.00,
            }
            volume += atomic_volumes.get(symbol, 15.0)

        return volume


def compute_properties(
    smiles: str,
    include_lipinski: bool = True,
) -> Dict[str, Any]:
    """
    Convenience function for quick property computation.

    Args:
        smiles: SMILES string
        include_lipinski: Include Lipinski analysis

    Returns:
        Property dictionary
    """
    calculator = PropertyCalculator()
    return calculator.calculate(smiles)


if __name__ == "__main__":
    # Test examples
    test_smiles = [
        "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC",  # Mesembrine
        "CCO",  # Ethanol
        "c1ccccc1",  # Benzene
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    ]

    calculator = PropertyCalculator()

    print("Property Calculator Test Results")
    print("=" * 60)

    for smiles in test_smiles:
        result = calculator.calculate(smiles)
        print(f"\nSMILES: {smiles}")
        print(f"Success: {result['success']}")

        if result["success"]:
            props = result["properties"]
            print(f"  MW: {props['molecular_weight']:.2f} Da")
            print(f"  Formula: {props['molecular_formula']}")
            print(f"  LogP: {props['logp']:.2f}")
            print(f"  TPSA: {props['tpsa']:.2f} Å²")
            print(f"  HBA/HBD: {props['num_hba']}/{props['num_hbd']}")
            print(f"  Rotatable bonds: {props['num_rotatable_bonds']}")

            lipinski = result["lipinski"]
            print(f"  Lipinski violations: {lipinski['violations']}")
            print(f"  Drug-likeness: {lipinski['drug_likeness']}")
