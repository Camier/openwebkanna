"""
Indigo-based SMILES syntax validation (Level 1).

This module provides syntax validation using the Indigo toolkit from EPAM.
Indigo is chosen over RDKit for canonicalization because it better preserves
stereochemistry information - critical for natural product alkaloids.

References:
- Indigo Toolkit: https://lifescience.opensource.epam.com/indigo/
- Indigo API: https://lifescience.opensource.epam.com/indigo/manual/

Validation checks:
1. SMILES parseability (valence rules, ring closures)
2. No excessive wildcards (* atoms)
3. No placeholder text ("n/a", "unknown", etc.)
4. Canonicalization with stereochemistry preservation
"""

from typing import Dict, Any, Optional
from indigo import Indigo, IndigoException


class IndigoValidator:
    """
    Level 1 syntax validator using Indigo toolkit.

    Provides SMILES validation, canonicalization, and syntax error detection.
    Thread-safe: each instance has its own Indigo context.
    """

    # Placeholder values that indicate extraction failures
    DEFAULT_PLACEHOLDER_VALUES = frozenset(
        [
            "n/a",
            "na",
            "unknown",
            "failed",
            "none",
            "null",
            "",
            "not available",
            "not determined",
            "mixture",
        ]
    )

    # Maximum allowed wildcards for non-R-group papers
    DEFAULT_MAX_WILDCARDS = 2

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize Indigo context and optional syntax policy."""
        config = config or {}
        self.indigo = Indigo()
        self.runtime_options = config.get("runtime_options", {})
        self.extra_set_options = config.get("set_options", {})
        self.reset_options_each_call = bool(
            self.runtime_options.get("reset_options_each_call", False)
        )
        self.max_wildcards = config.get(
            "wildcard_threshold",
            self.DEFAULT_MAX_WILDCARDS,
        )
        placeholders = config.get("placeholder_values", list(self.DEFAULT_PLACEHOLDER_VALUES))
        self.placeholder_values = frozenset(str(v).lower().strip() for v in placeholders)
        self._runtime_option_errors = []
        self._apply_runtime_options()

    def _set_indigo_option(self, option_name: str, option_value: Any) -> Optional[str]:
        """Safely apply a single Indigo option and return optional error string."""
        try:
            if isinstance(option_value, bool):
                self.indigo.setOption(option_name, option_value)
            elif isinstance(option_value, int):
                self.indigo.setOption(option_name, option_value)
            elif isinstance(option_value, float):
                self.indigo.setOption(option_name, option_value)
            elif isinstance(option_value, str):
                self.indigo.setOption(option_name, option_value)
            else:
                self.indigo.setOption(option_name, str(option_value))
            return None
        except Exception as exc:
            return f"{option_name}: {exc}"

    def _apply_runtime_options(self) -> None:
        """Apply configured Indigo runtime options (best-effort, non-fatal)."""
        self._runtime_option_errors = []

        timeout_ms = self.runtime_options.get("timeout_ms")
        if timeout_ms is not None:
            err = self._set_indigo_option("timeout", int(timeout_ms))
            if err:
                self._runtime_option_errors.append(err)

        # Optional raw Indigo setOption map for advanced runtime tuning.
        if isinstance(self.extra_set_options, dict):
            for option_name, option_value in self.extra_set_options.items():
                if not option_name:
                    continue
                err = self._set_indigo_option(str(option_name), option_value)
                if err:
                    self._runtime_option_errors.append(err)

    def validate_syntax(self, smiles: str) -> Dict[str, Any]:
        """
        Perform Level 1 syntax validation on a SMILES string.

        Args:
            smiles: Raw SMILES string from OCSR backend

        Returns:
            Dictionary with validation results:
            {
                'is_valid': bool,
                'canonical_smiles': str or None,
                'error': str or None,
                'wildcard_count': int or None,
                'molecular_formula': str or None
            }
        """
        result = {
            "is_valid": False,
            "canonical_smiles": None,
            "error": None,
            "wildcard_count": None,
            "molecular_formula": None,
        }

        # Check for placeholder values first (fast path)
        if self._is_placeholder(smiles):
            result["error"] = "placeholder_value"
            return result

        # Attempt to parse and canonicalize
        try:
            if self.reset_options_each_call:
                self._apply_runtime_options()
            mol = self.indigo.loadMolecule(smiles)

            # Count wildcards
            wildcard_count = self._count_wildcards(mol)
            result["wildcard_count"] = wildcard_count

            # Reject excessive wildcards (unless R-group paper)
            if wildcard_count > self.max_wildcards:
                result["error"] = "too_many_wildcards"
                return result

            # Get canonical SMILES (preserves stereochemistry)
            canonical = mol.canonicalSmiles()
            result["canonical_smiles"] = canonical

            # Get molecular formula for downstream validation
            result["molecular_formula"] = self._get_formula(mol)
            if self._runtime_option_errors:
                result["runtime_option_warnings"] = list(self._runtime_option_errors)

            # Additional valence check
            if not self._check_valence(mol):
                result["error"] = "valence_error"
                return result

            result["is_valid"] = True

        except IndigoException as e:
            result["error"] = f"parse_error: {str(e)}"
        except Exception as e:
            result["error"] = f"unexpected_error: {str(e)}"

        return result

    def _is_placeholder(self, smiles: str) -> bool:
        """Check if SMILES is a placeholder value."""
        return smiles.lower().strip() in self.placeholder_values

    def _count_wildcards(self, mol) -> int:
        """Count wildcard atoms (*) in molecule."""
        count = 0
        for atom in mol.iterateAtoms():
            if atom.symbol() == "*":
                count += 1
        return count

    def _get_formula(self, mol) -> str:
        """Get molecular formula from molecule."""
        return mol.molecularFormula()

    def _check_valence(self, mol) -> bool:
        """
        Check if all atoms satisfy valence rules.

        Indigo's loadSMILES already does basic valence checking,
        but we add an additional check for aromaticity issues.
        """
        try:
            # Attempt aromaticity perception
            mol.aromatize()
            return True
        except IndigoException:
            return False

    def canonicalize(self, smiles: str) -> Optional[str]:
        """
        Canonicalize a SMILES string without full validation.

        Use when you only need canonicalization, not validation.

        Args:
            smiles: Input SMILES string

        Returns:
            Canonical SMILES or None if parsing fails
        """
        try:
            if self.reset_options_each_call:
                self._apply_runtime_options()
            mol = self.indigo.loadMolecule(smiles)
            return mol.canonicalSmiles()
        except Exception:
            return None

    def smiles_to_cml(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to CML (Chemical Markup Language) format.

        Useful for visualization and interoperability.

        Args:
            smiles: Input SMILES string

        Returns:
            CML string or None if conversion fails
        """
        try:
            if self.reset_options_each_call:
                self._apply_runtime_options()
            mol = self.indigo.loadMolecule(smiles)
            return mol.cml()
        except Exception:
            return None

    def smiles_to_inchi(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to InChI format.

        Args:
            smiles: Input SMILES string

        Returns:
            InChI string or None if conversion fails
        """
        try:
            if self.reset_options_each_call:
                self._apply_runtime_options()
            mol = self.indigo.loadMolecule(smiles)
            return mol.inchi()
        except Exception:
            return None

    def smiles_to_inchikey(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to InChIKey.

        Args:
            smiles: Input SMILES string

        Returns:
            InChIKey string or None if conversion fails
        """
        try:
            if self.reset_options_each_call:
                self._apply_runtime_options()
            mol = self.indigo.loadMolecule(smiles)
            return mol.inchiKey()
        except Exception:
            return None

    def get_stereochemistry_info(self, smiles: str) -> Dict[str, Any]:
        """
        Extract stereochemistry information from SMILES.

        Args:
            smiles: Input SMILES string

        Returns:
            Dictionary with stereochemistry details
        """
        result = {
            "has_stereochemistry": False,
            "chiral_centers": 0,
            "double_bond_stereo": 0,
            "is_absolute": False,
        }

        try:
            if self.reset_options_each_call:
                self._apply_runtime_options()
            mol = self.indigo.loadMolecule(smiles)

            # Count chiral centers
            chiral_count = 0
            cdb_count = 0  # Cis-double bond count

            for atom in mol.iterateAtoms():
                if atom.piElectrons() == 0 and atom.bondCount() == 4:
                    # Potential chiral center
                    if mol.cipChirality(atom.id()) != "":
                        chiral_count += 1

            for bond in mol.iterateBonds():
                if bond.bondOrder() == 2:
                    # Double bond
                    stereo = bond.stereo()
                    if stereo > 0:  # CIS or TRANS
                        cdb_count += 1

            result["chiral_centers"] = chiral_count
            result["double_bond_stereo"] = cdb_count
            result["has_stereochemistry"] = (chiral_count > 0) or (cdb_count > 0)
            result["is_absolute"] = chiral_count > 0  # Has defined chiral centers

        except Exception:
            pass

        return result

    def get_runtime_policy(self) -> Dict[str, Any]:
        """Return active runtime policy and option application warnings."""
        return {
            "timeout_ms": self.runtime_options.get("timeout_ms"),
            "reset_options_each_call": self.reset_options_each_call,
            "set_options": dict(self.extra_set_options)
            if isinstance(self.extra_set_options, dict)
            else {},
            "option_errors": list(self._runtime_option_errors),
        }


def validate_smiles_syntax(smiles: str) -> Dict[str, Any]:
    """
    Convenience function for quick syntax validation.

    Args:
        smiles: SMILES string to validate

    Returns:
        Validation result dictionary
    """
    validator = IndigoValidator()
    return validator.validate_syntax(smiles)


if __name__ == "__main__":
    # Test examples
    test_cases = [
        # Valid molecules
        ("CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC", True),  # Mesembrine
        ("CCO", True),  # Ethanol
        ("C.C.C", True),  # Three separate molecules (valid but suspicious)
        # Invalid
        ("n/a", False),  # Placeholder
        ("unknown", False),  # Placeholder
        ("C(*)*", False),  # Too many wildcards (2 is ok, but let's test)
        ("invalid_smiles", False),  # Syntax error
        ("", False),  # Empty
    ]

    validator = IndigoValidator()
    print("Indigo Validator Test Results")
    print("=" * 60)

    for smiles, expected_valid in test_cases:
        result = validator.validate_syntax(smiles)
        status = "✓" if result["is_valid"] == expected_valid else "✗"
        print(
            f"{status} {smiles[:40]:<40} -> {result['is_valid']} ({result.get('error', 'ok')})"
        )
