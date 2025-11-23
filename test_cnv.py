# cnvision/examples/test_cnv.py

import os
import sys

# Ensure cnvision is importable regardless of where this script is run
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.insert(0, repo_root)

from cnvision.mane_loader import load_mane_exons
from cnvision.coordinate_mapper import (
    map_cnv_to_exons,
    map_exon_numbers_to_coordinates,
)
from cnvision.functional_predictor import predict_cnv_effect


def get_user_cnv_input():
    """
    Allow flexible user input:
    - Gene + start/end coordinates
    - Gene + exon numbers
    """

    gene = input("Enter gene symbol (e.g., BRCA1): ").strip().upper()

    print("\nChoose input type:")
    print("1 → Genomic coordinates")
    print("2 → Exon numbers")
    choice = input("Enter 1 or 2: ").strip()

    if choice == "1":
        # Coordinate-based input
        start = int(input("Enter CNV start genomic coordinate: ").strip())
        end = int(input("Enter CNV end genomic coordinate: ").strip())
        cnv_type = input("Type (deletion/duplication): ").strip().lower()

        return {
            "gene": gene,
            "start": start,
            "end": end,
            "type": cnv_type
        }

    elif choice == "2":
        # Exon-based input
        exon_range = input("Enter exon(s), e.g., '3' or '3-7': ").strip()
        cnv_type = input("Type (deletion/duplication): ").strip().lower()

        if "-" in exon_range:
            e1, e2 = exon_range.split("-")
            exon_numbers = list(range(int(e1), int(e2) + 1))
        else:
            exon_numbers = [int(exon_range)]

        return {
            "gene": gene,
            "exons": exon_numbers,
            "type": cnv_type
        }

    else:
        raise ValueError("Invalid selection. Choose 1 or 2.")


def main():
    # Load MANE exons
    mane_path = os.path.join(repo_root, "cnvision", "data",
                             "MANE.GRCh38.v1.4.refseq_genomic.gff.gz")
    mane_data = load_mane_exons(mane_path)

    print(f"\nLoaded {len(mane_data)} genes from MANE.\n")

    # Get flexible CNV input
    cnv = get_user_cnv_input()

    # If using exon numbers → convert to genomic coordinates
    if "exons" in cnv:
        print(f"Mapping exons {cnv['exons']} to genomic coordinates...")
        cnv["start"], cnv["end"] = map_exon_numbers_to_coordinates(
            cnv["gene"], cnv["exons"], mane_data
        )
        print(f"Mapped range: {cnv['start']} – {cnv['end']}")

    # Map the CNV
    hits = map_cnv_to_exons(cnv, mane_data)

    # Predict functional consequences
    effects = predict_cnv_effect(hits, cnv_type=cnv["type"])

    # Display results
    print("\n=== CNV Consequence Analysis ===")
    for eff in effects:
        print(
            f"\nGene: {cnv['gene']}\n"
            f"Transcript: {eff['transcript']}\n"
            f"Exons hit: {eff['hit_exons']}\n"
            f"Consequence: {eff['predicted_consequence']}"
        )


if __name__ == "__main__":
    main()
