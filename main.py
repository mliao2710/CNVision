"""
main.py — Command-line interface for CNVision

Allows users to run CNV mapping in two modes:
1) Genomic coordinate mode
2) Transcript exon-number mode

Loads MANE data, allows transcript selection, and prints CNV consequences.
"""

import os
import sys

# Ensure the project root is on sys.path so the cnvision package is discoverable
repo_root = os.path.dirname(os.path.abspath(__file__))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

# Import from the cnvision package
from cnvision.mane_loader import load_mane_exons
from cnvision.coordinate_mapper import map_cnv_to_exons, map_exon_numbers_to_regions
from cnvision.functional_predictor import predict_cnv_effect


def select_transcript(mane_data, gene):
    """Allow user to pick transcript if multiple exist."""
    transcripts = list(mane_data[gene].keys())
    if len(transcripts) == 1:
        return transcripts[0]

    print(f"\nMultiple MANE transcripts found for {gene}:")
    for i, t in enumerate(transcripts, 1):
        print(f"{i}. {t}")

    while True:
        choice = input("Select transcript number: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(transcripts):
            return transcripts[int(choice) - 1]
        print("Invalid selection. Try again.")


def run_coordinate_mode(mane_data):
    """User enters genomic coordinates."""
    gene = input("Enter gene name: ").strip()
    start = int(input("Enter genomic start coordinate: "))
    end = int(input("Enter genomic end coordinate: "))
    cnv_type = input("Enter CNV type (deletion/duplication): ").strip().lower()

    if gene not in mane_data:
        print(f"Gene {gene} not found in MANE data.")
        return

    transcript = select_transcript(mane_data, gene)
    cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type}
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    effects = predict_cnv_effect(hits)

    print("\n=== Results ===")
    for ef in effects:
        print(
            f"Gene: {gene}, Transcript: {ef['transcript']}, "
            f"Exons hit: {ef['hit_exons']}, Consequence: {ef['predicted_consequence']}"
        )


def run_exon_mode(mane_data):
    """User enters exon numbers instead of genomic coordinates."""
    gene = input("Enter gene name: ").strip()
    first_exon = int(input("Enter first exon number: "))
    last_exon = int(input("Enter last exon number: "))
    cnv_type = input("Enter CNV type (deletion/duplication): ").strip().lower()

    if gene not in mane_data:
        print(f"Gene {gene} not found in MANE data.")
        return

    transcript = select_transcript(mane_data, gene)
    cnv_region = map_exon_numbers_to_regions(gene, transcript, first_exon, last_exon, mane_data)

    if cnv_region is None:
        print("Invalid exon numbers.")
        return

    hits = map_cnv_to_exons(
        {"gene": gene, "start": cnv_region[0], "end": cnv_region[1], "type": cnv_type},
        mane_data,
        transcript
    )
    effects = predict_cnv_effect(hits)

    print("\n=== Results ===")
    for ef in effects:
        print(
            f"Gene: {gene}, Transcript: {ef['transcript']}, "
            f"Exons hit: {ef['hit_exons']}, Consequence: {ef['predicted_consequence']}"
        )


def main():
    """Main CLI entry point for CNVision."""
    mane_path = os.path.join(repo_root, "data", "MANE.GRCh38.v1.4.refseq_genomic.gff.gz")

    print("Loading MANE data...")
    mane_data = load_mane_exons(mane_path)
    print(f"Loaded {len(mane_data)} genes.\n")

    print("Select input mode:")
    print("1. Genomic coordinate mode")
    print("2. Transcript exon-number mode")
    print("3. Show example inputs")

    choice = input("Choice: ").strip()

    if choice == "1":
        run_coordinate_mode(mane_data)
    elif choice == "2":
        run_exon_mode(mane_data)
    elif choice == "3":
        print("\nExamples:\n")
        print("Genomic mode example:")
        print("  Gene: BRCA1")
        print("  Start: 43044295")
        print("  End: 43045805")
        print("  Type: deletion\n")

        print("Exon mode example:")
        print("  Gene: BRCA1")
        print("  First exon: 2")
        print("  Last exon: 5")
        print("  Type: duplication")
    else:
        print("Invalid choice.")


if __name__ == "__main__":
    main()
