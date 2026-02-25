"""
OpenWebUI Knowledge Base format converter.

Converts validated molecule JSONL to OpenWebUI KB formats:
- Markdown (human-readable documents)
- JSONL (machine-parseable metadata)

References:
- OpenWebUI KB API: https://docs.openwebui.com/
- OpenAPI KB endpoints: POST /api/v1/knowledge/
"""

import json
import markdown
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime


class KBConverter:
    """
    Convert molecules to OpenWebUI KB format.

    Supports:
    - Markdown documents with formatted properties
    - JSONL for batch import
    - Document metadata for retrieval
    """

    def __init__(
        self,
        output_format: str = "markdown",
        include_properties: bool = True,
        include_fingerprints: bool = False,
    ):
        """
        Initialize KB converter.

        Args:
            output_format: "markdown" or "jsonl"
            include_properties: Include property tables
            include_fingerprints: Include fingerprint data (increases size)
        """
        self.output_format = output_format
        self.include_properties = include_properties
        self.include_fingerprints = include_fingerprints

    def convert(
        self,
        molecule: Dict[str, Any],
        output_path: Path = None,
    ) -> str:
        """
        Convert single molecule to KB format.

        Args:
            molecule: Validated molecule dictionary
            output_path: Optional file path to write

        Returns:
            Formatted document string
        """
        if self.output_format == "markdown":
            content = self._to_markdown(molecule)
        elif self.output_format == "jsonl":
            content = json.dumps(molecule, indent=None)
        else:
            raise ValueError(f"Unknown format: {self.output_format}")

        # Write to file if path provided
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                f.write(content)

        return content

    def convert_batch(
        self,
        molecules: List[Dict[str, Any]],
        output_dir: Path,
    ) -> List[Path]:
        """
        Convert multiple molecules to KB format.

        Args:
            molecules: List of molecule dictionaries
            output_dir: Output directory

        Returns:
            List of created file paths
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        created_files = []

        for i, mol in enumerate(molecules):
            # Generate filename from SMILES or ID
            filename = self._generate_filename(mol, i)
            output_path = output_dir / filename

            self.convert(mol, output_path)
            created_files.append(output_path)

        return created_files

    def _to_markdown(self, molecule: Dict[str, Any]) -> str:
        """Convert molecule to Markdown document."""
        lines = []

        # Title (compound name or SMILES hash)
        title = (
            molecule.get("compound_name") or f"Molecule {molecule.get('id', 'unknown')}"
        )
        lines.append(f"# {title}")
        lines.append("")

        # SMILES (prominent display)
        smiles = molecule.get("smiles") or molecule.get("canonical_smiles")
        lines.append(f"**SMILES**: `{smiles}`")
        lines.append("")

        # Source information
        if "source" in molecule:
            lines.append("## Source")
            lines.append(f"- **Paper**: {molecule['source'].get('paper', 'Unknown')}")
            lines.append(f"- **Page**: {molecule['source'].get('page', 'N/A')}")
            lines.append(
                f"- **Image**: `{molecule['source'].get('image_path', 'N/A')}`"
            )
            lines.append("")

        # Properties table
        if self.include_properties and "properties" in molecule:
            lines.append("## Properties")
            lines.append("")
            lines.append("| Property | Value |")
            lines.append("|----------|-------|")

            props = molecule["properties"]
            if "molecular_weight" in props:
                lines.append(
                    f"| Molecular Weight | {props['molecular_weight']:.2f} Da |"
                )
            if "logp" in props:
                lines.append(f"| LogP | {props['logp']:.2f} |")
            if "tpsa" in props:
                lines.append(f"| TPSA | {props['tpsa']:.2f} Å² |")
            if "num_hba" in props:
                lines.append(f"| H-Bond Acceptors | {props['num_hba']} |")
            if "num_hbd" in props:
                lines.append(f"| H-Bond Donors | {props['num_hbd']} |")
            if "num_rotatable_bonds" in props:
                lines.append(f"| Rotatable Bonds | {props['num_rotatable_bonds']} |")
            if "molecular_formula" in props:
                lines.append(f"| Formula | {props['molecular_formula']} |")

            lines.append("")

        # Lipinski analysis
        if "lipinski" in molecule:
            lines.append("## Drug-Likeness")
            lines.append("")
            lip = molecule["lipinski"]
            lines.append(f"- **Lipinski Violations**: {lip.get('violations', 0)}")
            lines.append(f"- **Drug-Likeness**: {lip.get('drug_likeness', 'unknown')}")
            if lip.get("violation_details"):
                lines.append(f"- **Details**: {', '.join(lip['violation_details'])}")
            lines.append("")

        # Validation results
        if "validation" in molecule:
            lines.append("## Validation")
            lines.append("")
            val = molecule["validation"]

            if "level_1_syntax" in val:
                status = "✓" if val["level_1_syntax"].get("is_valid") else "✗"
                lines.append(f"- **Syntax (Indigo)**: {status}")

            if "level_2_chemical" in val:
                status = "✓" if val["level_2_chemical"].get("is_valid") else "✗"
                lines.append(f"- **Chemical Plausibility**: {status}")

            if "level_3_domain" in val:
                status = "✓" if val["level_3_domain"].get("is_valid") else "✗"
                lines.append(f"- **Domain (Alkaloid)**: {status}")

                if val["level_3_domain"].get("compound_class"):
                    lines.append(
                        f"- **Class**: {val['level_3_domain']['compound_class']}"
                    )

            lines.append("")

        # Gold standard matches
        if molecule.get("gold_standard_matches"):
            lines.append("## Similar Known Alkaloids")
            lines.append("")

            for match in molecule["gold_standard_matches"][:5]:  # Top 5
                lines.append(
                    f"- **{match['id']}**: {match['match_type']} "
                    f"(Tanimoto: {match['similarity']:.2f})"
                )
            lines.append("")

        # Confidence score
        if "confidence_score" in molecule:
            lines.append("## Confidence")
            lines.append("")
            score = molecule["confidence_score"]
            status = "HIGH" if score >= 70 else "MEDIUM" if score >= 50 else "LOW"
            lines.append(f"**Score**: {score}/100 ({status})")
            lines.append("")

        # Fingerprints (optional)
        if self.include_fingerprints and molecule.get("fingerprints"):
            lines.append("## Fingerprints")
            lines.append("")
            fps = molecule["fingerprints"]
            if "ecfp4" in fps:
                fp_hex = fps["ecfp4"].get("hex", "")
                lines.append(f"**ECFP4** ({fps['ecfp4'].get('nbits', 2048)} bits):")
                lines.append(f"```")
                lines.append(fp_hex[:200] + "..." if len(fp_hex) > 200 else fp_hex)
                lines.append(f"```")
            lines.append("")

        # Metadata footer
        lines.append("---")
        lines.append(f"*Generated: {datetime.now().isoformat()}*")
        lines.append(f"*Pipeline version: 1.0.0*")

        return "\n".join(lines)

    def _generate_filename(
        self,
        molecule: Dict[str, Any],
        index: int,
    ) -> str:
        """Generate filename for molecule document."""
        # Try to use compound name or gold standard match
        if molecule.get("gold_standard_matches"):
            match = molecule["gold_standard_matches"][0]
            name = match["id"]
            return f"{name}_{index:04d}.md"

        # Fallback to hash
        smiles = molecule.get("canonical_smiles", str(index))
        import hashlib

        hash_id = hashlib.md5(smiles.encode()).hexdigest()[:8]
        return f"molecule_{hash_id}_{index:04d}.md"

    def create_index(
        self,
        molecules: List[Dict[str, Any]],
        output_path: Path,
    ) -> Path:
        """
        Create index file for all molecules.

        Args:
            molecules: List of molecules
            output_path: Path to index file

        Returns:
            Path to created index
        """
        lines = [
            "# Molecule Knowledge Base Index",
            "",
            f"Total molecules: {len(molecules)}",
            f"Generated: {datetime.now().isoformat()}",
            "",
            "## Contents",
            "",
        ]

        # Group by compound class
        by_class = {}
        for mol in molecules:
            compound_class = mol.get("compound_class", "unknown")
            if compound_class not in by_class:
                by_class[compound_class] = []
            by_class[compound_class].append(mol)

        for compound_class, mols in sorted(by_class.items()):
            lines.append(f"### {compound_class.title()} ({len(mols)} molecules)")
            lines.append("")

            for mol in mols[:20]:  # Show first 20 per class
                smiles = mol.get("smiles", mol.get("canonical_smiles", "N/A"))[:50]
                score = mol.get("confidence_score", 0)
                lines.append(f"- {mol.get('id', 'Unknown')} (Score: {score})")
                lines.append(f"  - SMILES: `{smiles}...`")

            if len(mols) > 20:
                lines.append(f"- ... and {len(mols) - 20} more")

            lines.append("")

        # Write index
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write("\n".join(lines))

        return output_path


def convert_to_kb(
    molecule: Dict[str, Any],
    output_format: str = "markdown",
) -> str:
    """
    Convenience function for quick conversion.

    Args:
        molecule: Molecule dictionary
        output_format: "markdown" or "jsonl"

    Returns:
        Formatted document
    """
    converter = KBConverter(output_format)
    return converter.convert(molecule)


if __name__ == "__main__":
    # Test example molecule
    test_molecule = {
        "id": "test_001",
        "smiles": "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC",
        "canonical_smiles": "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC",
        "compound_name": "Mesembrine",
        "properties": {
            "molecular_weight": 289.41,
            "logp": 2.0,
            "tpsa": 32.2,
            "num_hba": 4,
            "num_hbd": 0,
            "molecular_formula": "C17H23NO3",
        },
        "lipinski": {
            "violations": 0,
            "drug_likeness": "good",
        },
        "validation": {
            "level_1_syntax": {"is_valid": True},
            "level_2_chemical": {"is_valid": True},
            "level_3_domain": {
                "is_valid": True,
                "compound_class": "alkaloid",
            },
        },
        "gold_standard_matches": [
            {
                "id": "mesembrine",
                "similarity": 1.0,
                "match_type": "exact",
            }
        ],
        "confidence_score": 95,
        "source": {
            "paper": "2008_-_Gericke_-_Sceletium_a_review_update",
            "page": 5,
        },
    }

    converter = KBConverter()
    markdown_doc = converter.convert(test_molecule)

    print("Generated Markdown Document")
    print("=" * 60)
    print(markdown_doc)
