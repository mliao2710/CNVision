"""
CNVision Functional Predictor - Core CNV Impact Analysis

Technologies: Python 3.10+ (core), typing module (type hints), modulo arithmetic

Core Functions:
- predict_cnv_effect(): determines frameshift vs in-frame effect
- Core logic: sum exon lengths, check if divisible by 3
- Consequence prediction (frameshift deletion/duplication, in-frame, etc.)
"""

from typing import List, Dict, Any

def predict_cnv_effect(
    exon_hits: List[Dict[str, Any]],
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]] = None
) -> List[Dict[str, Any]]:
    """
    Predict functional consequence of a CNV based on affected exons.
    
    Core logic: sum lengths of affected exons. If divisible by 3 → in-frame,
    otherwise → frameshift.
    
    Args:
        exon_hits: List of dicts with 'gene', 'transcript', 'hit_exons', 'cnv_type'
        mane_data: MANE gene/exon data structure
    
    Returns:
        Modified exon_hits list with 'predicted_consequence' field added
    Returns:
        List of result dictionaries with predicted consequences
    """
    results = []

    # ========================================================================
    # Process each CNV hit
    # ========================================================================
    for hit in exon_hits:
        # Extract basic information
        gene = hit.get("gene")
        transcript = hit.get("transcript")
        cnv_type = hit.get("cnv_type", "").lower()
        raw_exons = hit.get("hit_exons", [])

        # ====================================================================
        # STEP 1: Normalize exon data to include lengths
        # ====================================================================
        # Exons can come in different formats:
        # - List of numbers: [2, 3, 4] (need to look up lengths)
        # - List of dicts: [{"exon": 2, "length": 150}, ...] (already have lengths)
        
        consequence = "unknown coding impact (NM_XXXXX or non-MANE gene)"
        hit_exons = []  # Will store normalized exon data with lengths

        # Lookup exon lengths from MANE data if needed
        if raw_exons and mane_data:
            mane_gene_data = mane_data.get(gene)
            if mane_gene_data:
                mane_tx_data = mane_gene_data.get(transcript)
                if mane_tx_data:
                    # Process each exon
                    for e in raw_exons:
                        if isinstance(e, int):
                            # Exon is just a number - need to look up its length
                            # Find the exon info in MANE data
                            exon_info = next((ex for ex in mane_tx_data if ex["exon"] == e), None)
                            if exon_info:
                                # Calculate length: end - start + 1
                                # +1 because both start and end positions are inclusive
                                hit_exons.append({
                                    "exon": e,
                                    "length": exon_info.get("cds_length", exon_info["end"] - exon_info["start"] + 1)
                                })
                        elif isinstance(e, dict) and "length" in e:
                            # Exon already has length information
                            hit_exons.append(e)

        # Fallback: if we couldn't get lengths, try to use raw_exons as-is
        if not hit_exons:
            if raw_exons and isinstance(raw_exons[0], dict) and "length" in raw_exons[0]:
                hit_exons = raw_exons

        # ====================================================================
        # STEP 2: Calculate frameshift vs in-frame
        # ====================================================================
        if hit_exons:
            # Sum up all exon lengths
            # Example: exon1=150bp + exon2=123bp = 273bp total
            total_length = sum(exon["length"] for exon in hit_exons)
            if total_length == 0:
                consequence = f"no CDS impact ({len(hit_exons)} exon{'s' if len(hit_exons)>1 else ''}, UTR/non-coding region)"
                results.append({
                    "gene": gene,
                    "transcript": transcript,
                    "hit_exons": hit_exons if hit_exons else raw_exons,
                    "predicted_consequence": consequence
                })
                continue
            
            # THE KEY CALCULATION:
            # Modulo operator (%) gives remainder after division
            # If remainder is 0, number is divisible by 3 → in-frame
            # If remainder is not 0, number is not divisible → frameshift
            #
            # Examples:
            #   273 % 3 = 0 → in-frame (273 ÷ 3 = 91 exactly)
            #   274 % 3 = 1 → frameshift (274 ÷ 3 = 91 remainder 1)
            #   275 % 3 = 2 → frameshift (275 ÷ 3 = 91 remainder 2)
            frame_status = "in-frame" if total_length % 3 == 0 else "frameshift"
            
            num_exons = len(hit_exons)
            
            # ================================================================
            # STEP 3: Create human-readable consequence description
            # ================================================================
            # Format: "{frame_status} {cnv_type} ({num_exons} exon[s])"
            # Examples:
            #   "frameshift deletion (2 exons)"
            #   "in-frame duplication (1 exon)"
            if cnv_type == "deletion":
                consequence = f"{frame_status} deletion ({num_exons} exon{'s' if num_exons>1 else ''})"
            elif cnv_type == "duplication":
                consequence = f"{frame_status} duplication ({num_exons} exon{'s' if num_exons>1 else ''})"
            else:
                consequence = f"{frame_status} {cnv_type} ({num_exons} exon{'s' if num_exons>1 else ''})"

        # ====================================================================
        # STEP 4: Store results
        # ====================================================================
        results.append({
            "gene": gene,
            "transcript": transcript,
            "hit_exons": hit_exons if hit_exons else raw_exons,
            "predicted_consequence": consequence
        })

    return results


# -----------------------------
# Example usage / testing
# -----------------------------
if __name__ == "__main__":
    # Example MANE data stub
    mane_data = {
        "SPOUT1": {
            "NM_XXXXXX": [
                {"exon": 8, "start": 1000, "end": 1119, "strand": "+"},
                {"exon": 9, "start": 1120, "end": 1242, "strand": "+"},
                {"exon": 10, "start": 1243, "end": 1365, "strand": "+"},
                {"exon": 11, "start": 1366, "end": 1491, "strand": "+"},
            ]
        }
    }

    test_hits = [
        {"gene": "BRCA1", "transcript": "NM_007294.4",
         "hit_exons":[{"exon":23,"length":150}], "cnv_type":"deletion"},
        {"gene": "DMD", "transcript": "NM_004006.3",
         "hit_exons":[{"exon":45,"length":123},{"exon":46,"length":129}], "cnv_type":"duplication"},
        {"gene": "SPOUT1", "transcript": "NM_XXXXXX",
         "hit_exons":[8,9,10,11], "cnv_type":"duplication"},
        {"gene": "TP53", "transcript": "NM_000546.6", "hit_exons": [], "cnv_type": "deletion"},
    ]

    results = predict_cnv_effect(test_hits, mane_data)
    for r in results:
        print(
            f"Gene: {r['gene']}, Transcript: {r['transcript']}, "
            f"Exons hit: {r['hit_exons']}, Consequence: {r['predicted_consequence']}"
        )
