"""
CNVision Coordinate Mapper - Maps CNVs to Gene Exons

Technologies: Python 3.10+ (core), typing module (type hints), interval overlap algorithms

Core Functions:
- map_cnv_to_exons(): identifies which exons overlap with CNV region
- map_exon_numbers_to_regions(): converts exon numbers to genomic coordinates
- Coordinate system conversion (exon numbers ↔ genomic positions)
"""

from typing import Dict, List, Optional, Tuple, Any

def map_cnv_to_exons(
    cnv: Dict[str, Any],
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]],
    transcript: Optional[str] = None
) -> List[Dict[str, Any]]:
    """
    Find which exons are affected by a CNV given genomic coordinates.
    
    Args:
        cnv: Dict with 'gene', 'start', 'end', 'type', 'chromosome'
        mane_data: Nested dict structure {gene: {transcript: [exons]}}
        transcript: Optional transcript ID
    
    Returns:
        List containing one dict with 'transcript', 'hit_exons', 'cnv_type'
    """
    gene = cnv.get("gene")
    start = cnv.get("start")
    end = cnv.get("end")
    cnv_type = cnv.get("type")

    if gene not in mane_data:
        return []

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
    # Normalize chromosome for comparison (if provided).
    # Some CLI paths do not include chromosome in CNV input.
    cnv_chr = cnv.get("chromosome")  # e.g., "chrX"
    if cnv_chr is not None:
        cnv_chr = str(cnv_chr)
        if cnv_chr.startswith("chr"):
            cnv_chr = cnv_chr[3:]  # remove "chr" prefix

    # Check each exon for overlap
    for exon in exon_list:
        exon_start = exon.get("start")
        exon_end = exon.get("end")
        exon_number = exon.get("exon")
        exon_chr = exon.get("chromosome")

        # Normalize exon chromosome
        if exon_chr is not None:
            exon_chr = str(exon_chr)
            if exon_chr.startswith("chr"):
                exon_chr = exon_chr[3:]

        # Skip exons on other chromosomes when chromosome is known.
        if cnv_chr is not None and exon_chr != cnv_chr:
            continue

        # Skip exons with missing coordinates
        if exon_start is None or exon_end is None:
            continue

        # Overlap check - if CNV overlaps with exon, add exon number to results
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
            overlap_start = max(start, exon_start)
            overlap_end = min(end, exon_end)
            overlap_length = overlap_end - overlap_start + 1
            
            # Clinical-grade frame calls should be based on coding sequence only.
            coding_overlap = 0
            for seg in exon.get("cds_segments", []):
                cds_overlap_start = max(overlap_start, seg["start"])
                cds_overlap_end = min(overlap_end, seg["end"])
                if cds_overlap_start <= cds_overlap_end:
                    coding_overlap += cds_overlap_end - cds_overlap_start + 1
            
            hit_exons.append({
                "exon": exon_number,
                "start": overlap_start,
                "end": overlap_end,
                "length": coding_overlap,
                "genomic_length": overlap_length
            })

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
