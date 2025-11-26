"""
coordinate_mapper.py — Functions for mapping CNVs to exons and regions
"""

def map_cnv_to_exons(cnv, mane_data, transcript=None):
    """
    Map a CNV (genomic start/end) to MANE exons for a given transcript.

    Args:
        cnv: dict with keys 'gene', 'start', 'end', 'type'
        mane_data: dict returned by load_mane_exons
        transcript: optional transcript id (if None, take first MANE transcript)

    Returns:
        list of dicts with keys: transcript, hit_exons, cnv_type
    """
    gene = cnv["gene"]
    start = cnv["start"]
    end = cnv["end"]
    cnv_type = cnv["type"]

    if gene not in mane_data:
        return []

    gene_transcripts = mane_data[gene]
    if transcript is None:
        transcript = list(gene_transcripts.keys())[0]

    exon_list = gene_transcripts[transcript]
    hit_exons = []

    for exon in exon_list:
        exon_start = exon["start"]
        exon_end = exon["end"]

        # Check for overlap
        if end < exon_start or start > exon_end:
            continue
        hit_exons.append(exon["exon"])

    return [
        {
            "transcript": transcript,
            "hit_exons": hit_exons,
            "cnv_type": cnv_type,
            "predicted_consequence": None,  # To be filled by functional_predictor
        }
    ]


def map_exon_numbers_to_regions(gene, transcript, first_exon, last_exon, mane_data):
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
    if gene not in mane_data or transcript not in mane_data[gene]:
        return None

    exon_list = mane_data[gene][transcript]

    # Filter only exons in range
    selected = [e for e in exon_list if e["exon"] is not None and first_exon <= e["exon"] <= last_exon]

    if not selected:
        return None

    # Genomic start is the smallest start, end is largest end
    genomic_start = min(e["start"] for e in selected)
    genomic_end = max(e["end"] for e in selected)

    return genomic_start, genomic_end
