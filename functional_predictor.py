"""
functional_predictor.py — Predicts functional consequences of intragenic CNVs

Features:
- Handles deletions and duplications
- Supports multi-exon CNVs
- Works with MANE genes or NM_XXXXX-derived genes
- Determines frameshift vs in-frame by summing exon lengths
- Properly looks up exons by gene+transcript in MANE data
"""

from typing import List, Dict, Any

def predict_cnv_effect(
    exon_hits: List[Dict[str, Any]],
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]] = None
) -> List[Dict[str, Any]]:
    """
    Predict functional impact of CNV based on exons hit.
    exon_hits: list of dicts with:
        - 'gene': str
        - 'transcript': str
        - 'hit_exons': list[int] or list[dict] with 'exon' and 'length'
        - 'cnv_type': 'deletion' or 'duplication'
    mane_data: dict from mane_loader.load_mane_exons()
    """
    results = []

    for hit in exon_hits:
        gene = hit.get("gene")
        transcript = hit.get("transcript")
        cnv_type = hit.get("cnv_type", "").lower()
        raw_exons = hit.get("hit_exons", [])

        # Default fallback
        consequence = "unknown coding impact (NM_XXXXX or non-MANE gene)"
        hit_exons = []

        # Lookup in MANE if needed
        if raw_exons and mane_data:
            mane_gene_data = mane_data.get(gene)
            if mane_gene_data:
                mane_tx_data = mane_gene_data.get(transcript)
                if mane_tx_data:
                    for e in raw_exons:
                        # If numeric, convert to dict with length
                        if isinstance(e, int):
                            exon_info = next((ex for ex in mane_tx_data if ex["exon"] == e), None)
                            if exon_info:
                                hit_exons.append({
                                    "exon": e,
                                    "length": exon_info["end"] - exon_info["start"] + 1
                                })
                        # Already a dict with length
                        elif isinstance(e, dict) and "length" in e:
                            hit_exons.append(e)

        # If hit_exons still empty, fall back to raw_exons as-is
        if not hit_exons:
            # If raw_exons are dicts with lengths, use them
            if raw_exons and isinstance(raw_exons[0], dict) and "length" in raw_exons[0]:
                hit_exons = raw_exons

        # Calculate consequence if we have lengths
        if hit_exons:
            total_length = sum(exon["length"] for exon in hit_exons)
            frame_status = "in-frame" if total_length % 3 == 0 else "frameshift"
            num_exons = len(hit_exons)
            if cnv_type == "deletion":
                consequence = f"{frame_status} deletion ({num_exons} exon{'s' if num_exons>1 else ''})"
            elif cnv_type == "duplication":
                consequence = f"{frame_status} duplication ({num_exons} exon{'s' if num_exons>1 else ''})"
            else:
                consequence = f"{frame_status} {cnv_type} ({num_exons} exon{'s' if num_exons>1 else ''})"

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
