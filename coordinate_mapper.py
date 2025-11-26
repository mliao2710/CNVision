"""
coordinate_mapper.py — Functions for mapping CNVs to exons and regions
"""

from typing import Dict, List, Optional, Tuple, Any

def map_cnv_to_exons(
    cnv: Dict[str, Any],
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]],
    transcript: Optional[str] = None
) -> List[Dict[str, Any]]:
    """
    Map a CNV (genomic start/end) to MANE exons for a given transcript.

    Args:
        cnv: dict with keys 'gene', 'start', 'end', 'type'
        mane_data: dict returned by load_mane_exons, structured as gene -> transcript -> exon list
        transcript: optional transcript id (if None, use first MANE transcript)

    Returns:
        list of dicts with keys:
            - transcript
            - hit_exons (list of overlapping exon numbers)
            - cnv_type
            - predicted_consequence (currently None, placeholder for functional predictor)
    """
    gene = cnv.get("gene")
    start = cnv.get("start")
    end = cnv.get("end")
    cnv_type = cnv.get("type")

    if gene not in mane_data:
        return []

    transcripts = mane_data[gene]
    if not transcript:
        transcript = next(iter(transcripts.keys()))  # Take first transcript

    if transcript not in transcripts:
        return []

    exon_list = transcripts[transcript]
    hit_exons = []

    for exon in exon_list:
        exon_start = exon.get("start")
        exon_end = exon.get("end")
        exon_number = exon.get("exon")

        if exon_start is None or exon_end is None:
            continue

        # Check overlap
        if end >= exon_start and start <= exon_end:
            hit_exons.append(exon_number)

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
    """
    Convert a range of exon numbers to genomic coordinates.

    Args:
        gene: gene symbol
        transcript: transcript id
        first_exon: first exon number
        last_exon: last exon number
        mane_data: dict from load_mane_exons

    Returns:
        tuple (start, end) genomic coordinates, or None if invalid exon numbers
    """
    if gene not in mane_data:
        return None
    if transcript not in mane_data[gene]:
        return None

    exon_list = mane_data[gene][transcript]

    # Select exons within the specified range
    selected_exons = [
        e for e in exon_list
        if e.get("exon") is not None and first_exon <= e["exon"] <= last_exon
    ]

    if not selected_exons:
        return None

    genomic_start = min(e["start"] for e in selected_exons if e.get("start") is not None)
    genomic_end = max(e["end"] for e in selected_exons if e.get("end") is not None)

    return genomic_start, genomic_end
