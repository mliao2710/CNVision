"""
================================================================================
COORDINATE_MAPPER.PY - Map CNVs to Exons and Convert Between Coordinate Systems
================================================================================

PURPOSE:
--------
This module handles the conversion between different ways of describing CNVs:
1. Genomic coordinates (e.g., chromosome position 43000000-43002000)
2. Exon numbers (e.g., exons 2-4 of a gene)

It also identifies which exons are affected by a CNV, which is crucial for
predicting the biological impact.

BIOLOGICAL IMPACT:
------------------
Genes are made of exons (coding regions) separated by introns (non-coding).
When a CNV occurs, we need to know:
- Which exons are affected? (This determines what part of the protein is affected)
- How many exons? (More exons = potentially more severe impact)
- Exact boundaries? (Needed to calculate if reading frame is maintained)

This module provides the foundation for all downstream analysis.

CODING CONCEPTS & PRINCIPLES:
------------------------------
1. COORDINATE SYSTEMS:
   - Genomic coordinates: Absolute positions on chromosome (e.g., 43000000)
   - Exon numbers: Relative positions within gene (e.g., exon 1, 2, 3...)
   - Conversion: Exon numbers → genomic coordinates (and vice versa)

2. OVERLAP DETECTION:
   - Two ranges overlap if: range1_end >= range2_start AND range1_start <= range2_end
   - This is a fundamental algorithm for interval/range problems

3. DATA STRUCTURES:
   - Nested dictionaries: {gene: {transcript: [exon_list]}}
   - Lists of dictionaries: [{exon: 1, start: 100, end: 200}, ...]

4. OPTIONAL VALUES:
   - Uses Optional[Type] for values that might be None
   - Handles missing data gracefully (returns None or empty list)

5. FUNCTIONAL PROGRAMMING:
   - Pure functions: Input → Output, no side effects
   - Easy to test and reason about

================================================================================
"""

from typing import Dict, List, Optional, Tuple, Any

def map_cnv_to_exons(
    cnv: Dict[str, Any],
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]],
    transcript: Optional[str] = None
) -> List[Dict[str, Any]]:
    """
    Find which exons are affected by a CNV given genomic coordinates.
    
    This is the core mapping function: given a CNV region (start-end coordinates),
    it identifies all exons that overlap with that region.
    
    Example:
        CNV: start=43000000, end=43002000
        Gene BRCA1 has exons at positions:
          - Exon 22: 43000000-43000150 (OVERLAPS!)
          - Exon 23: 43000151-43000300 (OVERLAPS!)
          - Exon 24: 43000301-43000500 (doesn't overlap)
        Result: [22, 23]

    Args:
        cnv: Dictionary with keys:
            - 'gene': Gene symbol (e.g., "BRCA1")
            - 'start': Genomic start coordinate
            - 'end': Genomic end coordinate
            - 'type': "deletion" or "duplication"
        mane_data: Nested dictionary structure: {gene: {transcript: [exons]}}
        transcript: Optional transcript ID. If None, uses first available transcript.

    Returns:
        List containing one dictionary with:
            - 'transcript': Transcript ID used
            - 'hit_exons': List of exon numbers that overlap CNV
            - 'cnv_type': Type of CNV
            - 'predicted_consequence': None (filled by functional_predictor)
    """
    # ========================================================================
    # STEP 1: Extract CNV information
    # ========================================================================
    gene = cnv.get("gene")
    start = cnv.get("start")
    end = cnv.get("end")
    cnv_type = cnv.get("type")

    # ========================================================================
    # STEP 2: Validate gene exists in MANE data
    # ========================================================================
    if gene not in mane_data:
        return []  # Return empty list if gene not found

    # ========================================================================
    # STEP 3: Select transcript
    # ========================================================================
    transcripts = mane_data[gene]
    if not transcript:
        # If no transcript specified, use the first one available
        # next(iter(...)) gets first key from dictionary
        transcript = next(iter(transcripts.keys()))

    # Validate transcript exists
    if transcript not in transcripts:
        return []

    # ========================================================================
    # STEP 4: Get exon list for this transcript
    # ========================================================================
    exon_list = transcripts[transcript]
    hit_exons = []  # Will store exon numbers that overlap with CNV

    # ========================================================================
    # STEP 5: Check each exon for overlap with CNV region
    # ========================================================================
    for exon in exon_list:
        exon_start = exon.get("start")
        exon_end = exon.get("end")
        exon_number = exon.get("exon")

        # Skip exons with missing coordinate data
        if exon_start is None or exon_end is None:
            continue

        # OVERLAP ALGORITHM:
        # Two ranges overlap if:
        #   range1_end >= range2_start AND range1_start <= range2_end
        # 
        # Visual example:
        #   CNV:     [========]
        #   Exon:        [========]
        #   Result: OVERLAP! (CNV end >= Exon start AND CNV start <= Exon end)
        #
        #   CNV:     [====]
        #   Exon:              [====]
        #   Result: NO OVERLAP (CNV end < Exon start)
        if end >= exon_start and start <= exon_end:
            hit_exons.append(exon_number)

    # ========================================================================
    # STEP 6: Return results in standardized format
    # ========================================================================
    return [{
        "transcript": transcript,
        "hit_exons": hit_exons,  # List of affected exon numbers
        "cnv_type": cnv_type,
        "predicted_consequence": None  # Will be filled by functional_predictor
    }]


def map_exon_numbers_to_regions(
    gene: str,
    transcript: str,
    first_exon: int,
    last_exon: int,
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]]
) -> Optional[Tuple[int, int]]:
    """
    Convert exon numbers to genomic coordinates.
    
    This function does the reverse of map_cnv_to_exons:
    - Input: Exon numbers (e.g., exons 2-4)
    - Output: Genomic coordinates (e.g., start=1000000, end=1005000)
    
    Example:
        User says: "TP53 exons 2-4 are deleted"
        This function finds:
          - Exon 2: genomic positions 7565097-7565200
          - Exon 3: genomic positions 7565201-7565300
          - Exon 4: genomic positions 7565301-7565400
        Returns: (7565097, 7565400) - the span covering all affected exons

    Args:
        gene: Gene symbol (e.g., "TP53")
        transcript: Transcript ID (e.g., "NM_000546.6")
        first_exon: First exon number in range
        last_exon: Last exon number in range
        mane_data: Nested dictionary with gene/exon data

    Returns:
        Tuple of (start_coordinate, end_coordinate) if valid, None otherwise
    """
    # ========================================================================
    # STEP 1: Validate gene and transcript exist
    # ========================================================================
    if gene not in mane_data:
        return None
    if transcript not in mane_data[gene]:
        return None

    # ========================================================================
    # STEP 2: Get list of exons for this transcript
    # ========================================================================
    exon_list = mane_data[gene][transcript]

    # ========================================================================
    # STEP 3: Filter exons within the specified range
    # ========================================================================
    # List comprehension: creates new list with only exons in range
    # Condition: exon number must be between first_exon and last_exon (inclusive)
    selected_exons = [
        e for e in exon_list
        if e.get("exon") is not None and first_exon <= e["exon"] <= last_exon
    ]

    # If no exons found in range, return None
    if not selected_exons:
        return None

    # ========================================================================
    # STEP 4: Calculate genomic span
    # ========================================================================
    # Find the minimum start position and maximum end position
    # This gives us the full genomic region covering all selected exons
    #
    # Example:
    #   Exon 2: start=1000, end=1200
    #   Exon 3: start=1500, end=1700
    #   Exon 4: start=2000, end=2200
    #   Result: start=1000 (min), end=2200 (max)
    genomic_start = min(e["start"] for e in selected_exons if e.get("start") is not None)
    genomic_end = max(e["end"] for e in selected_exons if e.get("end") is not None)

    return genomic_start, genomic_end