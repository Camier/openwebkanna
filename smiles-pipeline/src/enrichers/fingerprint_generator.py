"""
Molecular fingerprint generator.

Generates molecular fingerprints for similarity search and clustering.
Supports ECFP4 (Morgan) and MACCS keys fingerprints.

References:
- RDKit fingerprints: https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprints
- ECFP4: Rogers & Hahn (2010), J. Chem. Inf. Model.
- MACCS keys: PubChem fingerprint system
"""

from typing import Dict, Any, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit import DataStructs
import numpy as np


class FingerprintGenerator:
    """
    Generate molecular fingerprints for similarity search.

    Supports:
    - ECFP4 (Morgan fingerprints, radius=2)
    - MACCS keys (166-bit structural keys)
    """

    def __init__(
        self,
        fingerprint_types: List[str] = None,
        ecfp_radius: int = 2,
        ecfp_nbits: int = 2048,
    ):
        """
        Initialize fingerprint generator.

        Args:
            fingerprint_types: List of fingerprint types to generate
            ecfp_radius: ECFP radius (2 = ECFP4, 3 = ECFP6)
            ecfp_nbits: Number of bits for ECFP
        """
        self.fingerprint_types = fingerprint_types or ["ecfp4", "macces"]
        self.ecfp_radius = ecfp_radius
        self.ecfp_nbits = ecfp_nbits

    def generate(
        self,
        smiles: str,
        canonical_smiles: str = None,
    ) -> Dict[str, Any]:
        """
        Generate fingerprints for a molecule.

        Args:
            smiles: SMILES string
            canonical_smiles: Pre-computed canonical SMILES

        Returns:
            Dictionary with fingerprints in hex format
        """
        result = {
            "success": False,
            "fingerprints": {},
            "error": None,
        }

        mol = Chem.MolFromSmiles(canonical_smiles or smiles)

        if mol is None:
            result["error"] = "rdkit_parse_failure"
            return result

        # Generate requested fingerprints
        for fp_type in self.fingerprint_types:
            if fp_type.lower() == "ecfp4":
                fp = self._generate_ecfp(mol)
                result["fingerprints"]["ecfp4"] = {
                    "hex": self._fp_to_hex(fp),
                    "nbits": self.ecfp_nbits,
                    "radius": self.ecfp_radius,
                }

            elif fp_type.lower() == "macces":
                fp = self._generate_maccs(mol)
                result["fingerprints"]["macces"] = {
                    "hex": self._fp_to_hex(fp),
                    "nbits": 167,  # MACCS is fixed at 167 bits (1-based)
                }

        result["success"] = True
        return result

    def _generate_ecfp(self, mol: Chem.Mol):
        """Generate ECFP4 (Morgan) fingerprint."""
        return AllChem.GetMorganFingerprintAsBitVect(
            mol,
            self.ecfp_radius,
            nBits=self.ecfp_nbits,
        )

    def _generate_maccs(self, mol: Chem.Mol):
        """Generate MACCS keys fingerprint."""
        return MACCSkeys.GenMACCSKeys(mol)

    def _fp_to_hex(self, fp) -> str:
        """Convert RDKit fingerprint to hex string for storage."""
        arr = np.zeros((1,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr.tobytes().hex()

    def _hex_to_fp(
        self,
        hex_str: str,
        nbits: int,
    ) -> "DataStructs.ExplicitBitVect":
        """Convert hex string back to RDKit fingerprint."""
        byte_data = bytes.fromhex(hex_str)
        fp = DataStructs.ExplicitBitVect(nbits)
        DataStructs.SetBitVectFromBinary(fp, byte_data)
        return fp

    def calculate_similarity(
        self,
        fp1_hex: str,
        fp2_hex: str,
        nbits: int,
        method: str = "tanimoto",
    ) -> float:
        """
        Calculate similarity between two fingerprints.

        Args:
            fp1_hex: First fingerprint (hex string)
            fp2_hex: Second fingerprint (hex string)
            nbits: Number of bits in fingerprint
            method: Similarity metric ("tanimoto", "dice", "cosine")

        Returns:
            Similarity score (0.0 to 1.0)
        """
        fp1 = self._hex_to_fp(fp1_hex, nbits)
        fp2 = self._hex_to_fp(fp2_hex, nbits)

        if method.lower() == "tanimoto":
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        elif method.lower() == "dice":
            return DataStructs.DiceSimilarity(fp1, fp2)
        elif method.lower() == "cosine":
            return DataStructs.CosineSimilarity(fp1, fp2)
        else:
            raise ValueError(f"Unknown similarity method: {method}")

    def find_similar_molecules(
        self,
        query_fp_hex: str,
        candidate_fps: List[Dict[str, Any]],
        nbits: int,
        threshold: float = 0.7,
        top_k: int = 10,
    ) -> List[Dict[str, Any]]:
        """
        Find molecules similar to query.

        Args:
            query_fp_hex: Query fingerprint (hex)
            candidate_fps: List of candidate molecules with fingerprints
            nbits: Number of bits in fingerprint
            threshold: Minimum similarity threshold
            top_k: Maximum number of results to return

        Returns:
            List of similar molecules with similarity scores
        """
        results = []

        for candidate in candidate_fps:
            candidate_fp = candidate.get("fingerprints", {}).get("ecfp4", {}).get("hex")

            if not candidate_fp:
                continue

            similarity = self.calculate_similarity(
                query_fp_hex,
                candidate_fp,
                nbits,
            )

            if similarity >= threshold:
                results.append(
                    {
                        "smiles": candidate.get("smiles"),
                        "similarity": similarity,
                        "id": candidate.get("id"),
                    }
                )

        # Sort by similarity descending
        results.sort(key=lambda x: x["similarity"], reverse=True)

        # Return top_k
        return results[:top_k]


def generate_fingerprints(
    smiles: str,
    fingerprint_types: List[str] = None,
) -> Dict[str, Any]:
    """
    Convenience function for quick fingerprint generation.

    Args:
        smiles: SMILES string
        fingerprint_types: Types to generate

    Returns:
        Fingerprint dictionary
    """
    generator = FingerprintGenerator(fingerprint_types)
    return generator.generate(smiles)


if __name__ == "__main__":
    # Test examples
    test_smiles = [
        "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC",  # Mesembrine
        "CCO",  # Ethanol
        "c1ccccc1",  # Benzene
    ]

    generator = FingerprintGenerator()

    print("Fingerprint Generator Test Results")
    print("=" * 60)

    for smiles in test_smiles:
        result = generator.generate(smiles)
        print(f"\nSMILES: {smiles}")
        print(f"Success: {result['success']}")

        if result["success"]:
            for fp_type, fp_data in result["fingerprints"].items():
                print(f"  {fp_type}: {fp_data['nbits']} bits")
                print(f"    Hex: {fp_data['hex'][:40]}...")

        # Test similarity (mesembrine vs mesembrine = 1.0)
        if len(test_smiles) >= 2:
            fp1 = generator.generate(test_smiles[0])
            fp2 = generator.generate(test_smiles[0])  # Same molecule

            if fp1["success"] and fp2["success"]:
                sim = generator.calculate_similarity(
                    fp1["fingerprints"]["ecfp4"]["hex"],
                    fp2["fingerprints"]["ecfp4"]["hex"],
                    2048,
                )
                print(f"\nSelf-similarity (should be 1.0): {sim:.3f}")
