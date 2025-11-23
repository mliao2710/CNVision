# cnvision/coordinate_mapper.py

def map_cnv_to_exons(cnv, mane_data):
    """
    Map a CNV (with keys: gene, start, end, type) to MANE exons for each transcript.
    Returns a list of dicts; one per transcript:
      {
        "transcript": transcript_id,
        "exons_hit": [ exon_numbers ],
        "overlaps": [
            {"exon": exon_number, "exon_start": x, "exon_end": y,
             "overlap_start": s, "overlap_end": e, "overlap_len": L}
        ],
        "total_overlap_len": N  # sum of overlap_len for this transcript
      }
    """
    gene = cnv.get("gene")
    start = cnv.get("start")
    end = cnv.get("end")

    if gene is None or start is None or end is None:
        return None

    if gene not in mane_data:
        return None

    results = []
    for transcript, exons in mane_data[gene].items():
        overlaps = []
        exons_hit = []
        total_overlap = 0
        for exon in exons:
            exon_num = exon.get("exon")
            exon_start = exon.get("start")
            exon_end = exon.get("end")

            # check overlap
            ov_start = max(start, exon_start)
            ov_end = min(end, exon_end)
            if ov_start <= ov_end:
                overlap_len = ov_end - ov_start + 1  # inclusive coordinates
                overlaps.append({
                    "exon": exon_num,
                    "exon_start": exon_start,
                    "exon_end": exon_end,
                    "overlap_start": ov_start,
                    "overlap_end": ov_end,
                    "overlap_len": overlap_len
                })
                exons_hit.append(exon_num)
                total_overlap += overlap_len

        results.append({
            "transcript": transcript,
            "exons_hit": exons_hit,
            "overlaps": overlaps,
            "total_overlap_len": total_overlap
        })
        

    return results

def map_exon_numbers_to_coordinates(gene, exon_list, mane_data):
    """
    Takes gene + exon numbers and returns genomic (start, end) coords.
    """

    if gene not in mane_data:
        raise ValueError(f"Gene {gene} not found in MANE")

    gene_transcripts = mane_data[gene]

    # Use MANE Select transcript if available
    transcript = sorted(gene_transcripts.keys())[0]
    exons = gene_transcripts[transcript]

    starts = []
    ends = []

    for exon_num in exon_list:
        if exon_num not in exons:
            raise ValueError(f"Exon {exon_num} not found in MANE for {gene}")

        starts.append(exons[exon_num]["genomic_start"])
        ends.append(exons[exon_num]["genomic_end"])

    return min(starts), max(ends)
