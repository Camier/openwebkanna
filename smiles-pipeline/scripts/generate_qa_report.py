#!/usr/bin/env python3
"""
SMILES Pipeline - QA Report Generator.

Generates comprehensive QA report including:
- Per-paper extraction statistics
- Backend performance breakdown
- Property distributions
- Gold standard validation
- Rejection analysis

Usage:
    python generate_qa_report.py \\
        --input-dir smiles-pipeline/data/validated \\
        --output-file smiles-pipeline/data/summary/qa_report.json
"""

import argparse
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class QAReportGenerator:
    """Generate comprehensive QA report for extraction results."""

    def __init__(self, input_dir: Path):
        """
        Initialize report generator.

        Args:
            input_dir: Directory with molecules.jsonl
        """
        self.input_dir = input_dir
        self.molecules: List[Dict[str, Any]] = []
        self.stats = {
            "summary": {},
            "per_paper": {},
            "backend_breakdown": {},
            "property_distributions": {},
            "gold_standard_matches": {},
            "rejection_reasons": {},
            "confidence_distribution": {},
        }

    def load_molecules(self) -> int:
        """Load molecules from JSONL file."""
        jsonl_file = self.input_dir / "molecules.jsonl"

        if not jsonl_file.exists():
            raise FileNotFoundError(f"Not found: {jsonl_file}")

        with open(jsonl_file) as f:
            for line in f:
                if line.strip():
                    self.molecules.append(json.loads(line))

        return len(self.molecules)

    def compute_summary_stats(self):
        """Compute overall summary statistics."""
        if not self.molecules:
            return

        # Count by validation status
        syntax_valid = sum(
            1
            for m in self.molecules
            if m.get("validation", {}).get("level_1_syntax", {}).get("is_valid")
        )
        chemically_valid = sum(
            1
            for m in self.molecules
            if m.get("validation", {}).get("level_2_chemical", {}).get("is_valid")
        )
        domain_valid = sum(
            1
            for m in self.molecules
            if m.get("validation", {}).get("level_3_domain", {}).get("is_valid")
        )
        high_confidence = sum(
            1 for m in self.molecules if m.get("confidence_score", 0) >= 70
        )

        # Count by compound class
        by_class = defaultdict(int)
        for mol in self.molecules:
            compound_class = mol.get("compound_class", "unknown")
            by_class[compound_class] += 1

        self.stats["summary"] = {
            "total_molecules": len(self.molecules),
            "syntax_valid": syntax_valid,
            "chemically_valid": chemically_valid,
            "domain_valid": domain_valid,
            "high_confidence": high_confidence,
            "validation_rates": {
                "syntax_rate": syntax_valid / len(self.molecules)
                if self.molecules
                else 0,
                "chemical_rate": chemically_valid / len(self.molecules)
                if self.molecules
                else 0,
                "domain_rate": domain_valid / len(self.molecules)
                if self.molecules
                else 0,
            },
            "compound_classes": dict(by_class),
            "generated_at": datetime.now().isoformat(),
        }

    def compute_per_paper_stats(self):
        """Compute per-paper breakdown."""
        by_paper = defaultdict(list)

        for mol in self.molecules:
            source = mol.get("source", {})
            paper = source.get("paper", source.get("image_path", "unknown"))
            # Extract paper ID from path
            paper_id = paper.split("/")[0] if "/" in paper else paper
            by_paper[paper_id].append(mol)

        for paper_id, mols in by_paper.items():
            high_conf = sum(1 for m in mols if m.get("confidence_score", 0) >= 70)
            self.stats["per_paper"][paper_id] = {
                "total": len(mols),
                "high_confidence": high_conf,
                "avg_confidence": sum(m.get("confidence_score", 0) for m in mols)
                / len(mols)
                if mols
                else 0,
            }

    def compute_backend_breakdown(self):
        """Compute performance by OCSR backend."""
        by_backend = defaultdict(list)

        for mol in self.molecules:
            backend = mol.get("backend_used", "unknown")
            by_backend[backend].append(mol)

        for backend, mols in by_backend.items():
            confidences = [m.get("extraction_confidence", 0) for m in mols]
            self.stats["backend_breakdown"][backend] = {
                "total": len(mols),
                "avg_confidence": sum(confidences) / len(mols) if mols else 0,
                "success_rate": len(mols) / len(self.molecules)
                if self.molecules
                else 0,
            }

    def compute_property_distributions(self):
        """Compute property distributions."""
        if not self.molecules:
            return

        properties = []
        for mol in self.molecules:
            if "properties" in mol:
                props = mol["properties"]
                properties.append(
                    {
                        "mw": props.get("molecular_weight"),
                        "logp": props.get("logp"),
                        "tpsa": props.get("tpsa"),
                        "hba": props.get("num_hba"),
                        "hbd": props.get("num_hbd"),
                    }
                )

        if properties:
            # Filter None values
            mw_values = [p["mw"] for p in properties if p.get("mw")]
            logp_values = [p["logp"] for p in properties if p.get("logp")]
            tpsa_values = [p["tpsa"] for p in properties if p.get("tpsa")]

            self.stats["property_distributions"] = {
                "molecular_weight": {
                    "mean": sum(mw_values) / len(mw_values) if mw_values else 0,
                    "median": sorted(mw_values)[len(mw_values) // 2]
                    if mw_values
                    else 0,
                    "min": min(mw_values) if mw_values else 0,
                    "max": max(mw_values) if mw_values else 0,
                    "std": self._compute_std(mw_values),
                },
                "logp": {
                    "mean": sum(logp_values) / len(logp_values) if logp_values else 0,
                    "median": sorted(logp_values)[len(logp_values) // 2]
                    if logp_values
                    else 0,
                    "min": min(logp_values) if logp_values else 0,
                    "max": max(logp_values) if logp_values else 0,
                },
                "tpsa": {
                    "mean": sum(tpsa_values) / len(tpsa_values) if tpsa_values else 0,
                    "median": sorted(tpsa_values)[len(tpsa_values) // 2]
                    if tpsa_values
                    else 0,
                    "min": min(tpsa_values) if tpsa_values else 0,
                    "max": max(tpsa_values) if tpsa_values else 0,
                },
            }

    def _compute_std(self, values: List[float]) -> float:
        """Compute standard deviation."""
        if not values:
            return 0.0

        mean = sum(values) / len(values)
        variance = sum((x - mean) ** 2 for x in values) / len(values)
        return variance**0.5

    def compute_gold_standard_matches(self):
        """Analyze gold standard matches."""
        matches = defaultdict(list)

        for mol in self.molecules:
            mol_matches = mol.get("gold_standard_matches", [])
            for match in mol_matches:
                alkaloid_id = match.get("id", "unknown")
                similarity = match.get("similarity", 0)
                matches[alkaloid_id].append(similarity)

        for alkaloid_id, similarities in matches.items():
            self.stats["gold_standard_matches"][alkaloid_id] = {
                "total_matches": len(similarities),
                "exact_matches": sum(1 for s in similarities if s >= 0.99),
                "similar_matches": sum(1 for s in similarities if 0.7 <= s < 0.99),
                "max_similarity": max(similarities) if similarities else 0,
                "avg_similarity": sum(similarities) / len(similarities)
                if similarities
                else 0,
            }

    def compute_confidence_distribution(self):
        """Compute confidence score distribution."""
        if not self.molecules:
            return

        scores = [m.get("confidence_score", 0) for m in self.molecules]

        # Bucket by ranges
        buckets = {
            "0-20": 0,
            "21-40": 0,
            "41-60": 0,
            "61-80": 0,
            "81-100": 0,
        }

        for score in scores:
            if score <= 20:
                buckets["0-20"] += 1
            elif score <= 40:
                buckets["21-40"] += 1
            elif score <= 60:
                buckets["41-60"] += 1
            elif score <= 80:
                buckets["61-80"] += 1
            else:
                buckets["81-100"] += 1

        self.stats["confidence_distribution"] = buckets

    def generate_report(self, output_file: Path):
        """Generate and save complete QA report."""
        print(f"Loading molecules from {self.input_dir}...")
        count = self.load_molecules()
        print(f"Loaded {count} molecules")

        print("Computing summary statistics...")
        self.compute_summary_stats()

        print("Computing per-paper breakdown...")
        self.compute_per_paper_stats()

        print("Computing backend performance...")
        self.compute_backend_breakdown()

        print("Computing property distributions...")
        self.compute_property_distributions()

        print("Analyzing gold standard matches...")
        self.compute_gold_standard_matches()

        print("Computing confidence distribution...")
        self.compute_confidence_distribution()

        # Write report
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "w") as f:
            json.dump(self.stats, f, indent=2)

        print(f"\nQA Report saved to: {output_file}")

        # Print summary
        self._print_summary()

        return self.stats

    def _print_summary(self):
        """Print summary to console."""
        print("\n" + "=" * 60)
        print("QA Report Summary")
        print("=" * 60)

        summary = self.stats.get("summary", {})
        print(f"Total molecules: {summary.get('total_molecules', 0)}")
        print(
            f"Syntax valid: {summary.get('syntax_valid', 0)} ({summary.get('validation_rates', {}).get('syntax_rate', 0) * 100:.1f}%)"
        )
        print(
            f"Chemically valid: {summary.get('chemically_valid', 0)} ({summary.get('validation_rates', {}).get('chemical_rate', 0) * 100:.1f}%)"
        )
        print(
            f"Domain valid: {summary.get('domain_valid', 0)} ({summary.get('validation_rates', {}).get('domain_rate', 0) * 100:.1f}%)"
        )
        print(f"High confidence (≥70): {summary.get('high_confidence', 0)}")

        print("\nGold Standard Matches:")
        for alkaloid, data in self.stats.get("gold_standard_matches", {}).items():
            print(
                f"  {alkaloid}: {data['total_matches']} (max: {data['max_similarity']:.2f})"
            )

        mw = self.stats.get("property_distributions", {}).get("molecular_weight", {})
        if mw:
            print(
                f"\nMolecular Weight: {mw.get('mean', 0):.1f} ± {mw.get('std', 0):.1f} Da"
            )


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate QA report for SMILES extraction"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory with molecules.jsonl",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        required=True,
        help="Output JSON report file",
    )

    args = parser.parse_args()

    generator = QAReportGenerator(args.input_dir)
    generator.generate_report(args.output_file)

    return 0


if __name__ == "__main__":
    sys.exit(main())
