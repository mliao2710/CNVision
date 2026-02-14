"""
CNVision Coordinate Mapper - Exon/Coordinate Conversion Utilities

Purpose:
- Translate genomic CNV intervals into affected exon hits.
- Translate exon-number ranges into genomic coordinate spans.

Core Functions:
- `map_cnv_to_exons`: interval-overlap mapping from CNV coordinates to exon hits.
- `map_exon_numbers_to_regions`: reverse mapping from exon range to genomic span.
- Annotates overlap detail: coding overlap length and full vs partial exon coverage.
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

    # STEP 1: Select transcript.
    transcripts = mane_data[gene]
    if not transcript:
        transcript = next(iter(transcripts.keys()))

    # Validate transcript exists
    if transcript not in transcripts:
        return []

    # STEP 2: Prepare exon list.
    exon_list = transcripts[transcript]
    hit_exons = []

    # STEP 3: Normalize chromosome labels.
    cnv_chr = cnv.get("chromosome")
    if cnv_chr is not None:
        cnv_chr = str(cnv_chr)
        if cnv_chr.startswith("chr"):
            cnv_chr = cnv_chr[3:]

    # STEP 4: Scan exons and compute overlap payload.
    for exon in exon_list:
        exon_start = exon.get("start")
        exon_end = exon.get("end")
        exon_number = exon.get("exon")
        exon_chr = exon.get("chromosome")

        # Normalize exon chromosome label.
        if exon_chr is not None:
            exon_chr = str(exon_chr)
            if exon_chr.startswith("chr"):
                exon_chr = exon_chr[3:]

        # Skip exons on other chromosomes when CNV chromosome is known.
        if cnv_chr is not None and exon_chr != cnv_chr:
            continue

        # Skip exons with missing coordinates.
        if exon_start is None or exon_end is None:
            continue

        # Interval overlap test.
        if end >= exon_start and start <= exon_end:
            overlap_start = max(start, exon_start)
            overlap_end = min(end, exon_end)
            overlap_length = overlap_end - overlap_start + 1
            full_exon_coverage = overlap_start == exon_start and overlap_end == exon_end
            
            # Frame logic uses CDS overlap length, not full genomic overlap length.
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
                "genomic_length": overlap_length,
                "full_exon_coverage": full_exon_coverage
            })

    # STEP 5: Return normalized hit structure.
    return [{
        "transcript": transcript,
        "hit_exons": hit_exons,
        "cnv_type": cnv_type,
        "predicted_consequence": None
    }]


def map_exon_numbers_to_regions(
    gene: str,
    transcript: str,
    first_exon: int,
    last_exon: int,
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]]
) -> Optional[Tuple[int, int]]:
    """Convert an exon-number range to genomic start/end coordinates."""
    # STEP 1: Validate gene/transcript.
    if gene not in mane_data:
        return None
    if transcript not in mane_data[gene]:
        return None

    # STEP 2: Select exons in requested range.
    exon_list = mane_data[gene][transcript]

    selected_exons = [
        e for e in exon_list
        if e.get("exon") is not None and first_exon <= e["exon"] <= last_exon
    ]

    if not selected_exons:
        return None

    # STEP 3: Compute min/max genomic span.
    genomic_start = min(e["start"] for e in selected_exons if e.get("start") is not None)
    genomic_end = max(e["end"] for e in selected_exons if e.get("end") is not None)

    return genomic_start, genomic_end
