"""
CNVision CLI - Command-Line Interface for Terminal-Based Analysis

Technologies: Python 3.10+ (core), typing module (type hints), interactive I/O

Core Functions:
- Interactive menu system and user input handling
- Input mode coordination (exon, transcript, genomic coordinate modes)
- Analysis loop and result display in terminal
- Programmatic integration with bioinformatics pipelines
"""

import os
import sys
from typing import Dict, Any, Optional

repo_root = os.path.dirname(os.path.abspath(__file__))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from mane_loader import load_mane_exons
from coordinate_mapper import map_cnv_to_exons, map_exon_numbers_to_regions
from functional_predictor import predict_cnv_effect
from utils.ncbi_fetcher import fetch_refseq_info


def select_transcript(mane_data: Dict[str, Dict[str, Any]], gene: str) -> str:
    """
    Interactive transcript selection when multiple transcripts exist.
    
    Args:
        mane_data: MANE gene/exon data
        gene: Gene symbol
    
    Returns:
        Selected transcript ID
    """
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


def run_coordinate_mode(mane_data: Dict[str, Any]) -> None:
    """
    Handle genomic coordinate input mode.
    
    User provides: gene name, genomic start, genomic end, CNV type
    
    Example:
        Gene: BRCA1, Start: 43044295, End: 43045805, Type: deletion
    """
    gene = input("Enter gene name: ").strip()
    
    try:
        start = int(input("Enter genomic start coordinate: ").strip())
        end = int(input("Enter genomic end coordinate: ").strip())
    except ValueError:
        print("Start and end must be integers.")
        return

    cnv_type = input("Enter CNV type (deletion/duplication): ").strip().lower()
    if gene not in mane_data:
        print(f"Gene {gene} not found in MANE data.")
        return

    # ========================================================================
    # STEP 3: Select transcript and process CNV
    # ========================================================================
    transcript = select_transcript(mane_data, gene)
    cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type}
    
    # Map CNV to exons and predict consequences
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    effects = predict_cnv_effect(hits, mane_data)  # Note: pass mane_data for exon length lookup

    # ========================================================================
    # STEP 4: Display results
    # ========================================================================
    print_results(gene, effects)


def run_exon_mode(mane_data: Dict[str, Any]) -> None:
    """
    Handle exon number input mode.
    
    User provides:
    - Gene name
    - First exon number
    - Last exon number
    - CNV type
    
    This mode is easier for users who know which exons are affected
    but don't know the exact genomic coordinates.
    
    Example:
        Gene: BRCA1
        First exon: 2
        Last exon: 5
        Type: duplication
    """
    # ========================================================================
    # STEP 1: Get user input
    # ========================================================================
    gene = input("Enter gene name: ").strip()
    
    try:
        first_exon = int(input("Enter first exon number: ").strip())
        last_exon = int(input("Enter last exon number: ").strip())
    except ValueError:
        print("Exon numbers must be integers.")
        return

    cnv_type = input("Enter CNV type (deletion/duplication): ").strip().lower()

    # ========================================================================
    # STEP 2: Validate gene exists
    # ========================================================================
    if gene not in mane_data:
        print(f"Gene {gene} not found in MANE data.")
        return

    # ========================================================================
    # STEP 3: Convert exon numbers to genomic coordinates
    # ========================================================================
    transcript = select_transcript(mane_data, gene)
    cnv_region = map_exon_numbers_to_regions(gene, transcript, first_exon, last_exon, mane_data)
    
    if cnv_region is None:
        print("Invalid exon numbers or transcript.")
        return

    # ========================================================================
    # STEP 4: Process CNV (same as coordinate mode from here)
    # ========================================================================
    cnv = {"gene": gene, "start": cnv_region[0], "end": cnv_region[1], "type": cnv_type}
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    effects = predict_cnv_effect(hits, mane_data)

    print_results(gene, effects)


def run_nm_mode(mane_data: Dict[str, Any]) -> None:
    """
    Handle RefSeq transcript ID (NM_XXXXX) input mode.
    
    User provides:
    - RefSeq transcript ID (e.g., "NM_007294.4")
    - Genomic coordinates
    - CNV type
    
    This mode is useful when:
    - Gene isn't in MANE database
    - User wants to analyze a specific transcript variant
    - Working with alternative transcripts
    
    The function fetches gene information from NCBI if not in MANE.
    """
    # ========================================================================
    # STEP 1: Get RefSeq ID and fetch gene information
    # ========================================================================
    nm_id = input("Enter NM_XXXXX RefSeq ID: ").strip().upper()
    refseq_info = fetch_refseq_info(nm_id)

    if not refseq_info:
        print(f"Could not fetch information for {nm_id}.")
        return

    # ========================================================================
    # STEP 2: Extract gene symbol and transcript
    # ========================================================================
    gene = refseq_info["gene_symbol"]
    # Use fetched transcript, or let user select if gene is in MANE
    transcript = refseq_info.get("transcript_id") or select_transcript(mane_data, gene)

    # ========================================================================
    # STEP 3: Get genomic coordinates
    # ========================================================================
    try:
        start = int(input("Enter genomic start coordinate: ").strip())
        end = int(input("Enter genomic end coordinate: ").strip())
    except ValueError:
        print("Start and end must be integers.")
        return

    cnv_type = input("Enter CNV type (deletion/duplication): ").strip().lower()

    # ========================================================================
    # STEP 4: Process CNV (same pipeline as coordinate mode)
    # ========================================================================
    cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type}
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    effects = predict_cnv_effect(hits, mane_data)

    print_results(gene, effects)


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def print_results(gene: str, effects: list) -> None:
    """
    Display CNV analysis results to the user.
    
    Formats and prints:
    - Gene name
    - Transcript ID
    - List of affected exons
    - Predicted consequence (frameshift/in-frame)
    
    Args:
        gene: Gene symbol
        effects: List of result dictionaries from predict_cnv_effect()
    """
    print("\n=== Results ===")
    
    # Handle case where no exons were affected
    if not effects:
        print(f"No exons hit for {gene}.")
        return

    # Print each result (usually just one, but function supports multiple)
    for ef in effects:
        print(
            f"Gene: {gene}, Transcript: {ef['transcript']}, "
            f"Exons hit: {ef['hit_exons']}, Consequence: {ef['predicted_consequence']}"
        )


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main() -> None:
    """
    Main entry point for CNVision CLI application.
    
    This function:
    1. Loads MANE gene/exon data
    2. Displays menu of input modes
    3. Routes to appropriate handler based on user choice
    4. Processes CNV and displays results
    
    Called when script is run directly: python main.py
    """
    # ========================================================================
    # STEP 1: Load MANE data (this takes a few seconds)
    # ========================================================================
    mane_path = os.path.join(repo_root, "data", "MANE.GRCh38.v1.4.refseq_genomic.gff.gz")

    print("Loading MANE data...")
    mane_data = load_mane_exons(mane_path)
    print(f"Loaded {len(mane_data)} genes.\n")

    # ========================================================================
    # STEP 2: Display menu and get user choice
    # ========================================================================
    print("Select input mode:")
    print("1. Genomic coordinate mode")
    print("2. Transcript exon-number mode")
    print("3. NM_XXXXX RefSeq ID mode")
    print("4. Show example inputs")

    choice = input("Choice: ").strip()

    # ========================================================================
    # STEP 3: Route to appropriate handler function
    # ========================================================================
    if choice == "1":
        run_coordinate_mode(mane_data)
    elif choice == "2":
        run_exon_mode(mane_data)
    elif choice == "3":
        run_nm_mode(mane_data)
    elif choice == "4":
        show_examples()
    else:
        print("Invalid choice.")


def show_examples() -> None:
    """
    Display example inputs for each mode to help users understand the format.
    
    Shows realistic examples that users can copy/paste or use as templates.
    """
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
    print("  Type: duplication\n")

    print("NM_XXXXX mode example:")
    print("  NM_ID: NM_007294")
    print("  Start: 43044295")
    print("  End: 43045805")
    print("  Type: deletion\n")


if __name__ == "__main__":
    main()
