"""
functional_predictor.py — Predicts functional consequences of intragenic CNVs

Current features:
- Handles deletions and duplications
- Supports multi-exon CNVs
- Returns detailed consequence strings
- Placeholder for future CDS-aware analysis
"""

def predict_cnv_effect(exon_hits):
    """
    Predict functional impact of CNV based on exons hit.

    Parameters:
    ----------
    exon_hits : list of dict
        Each dict contains:
            - 'gene': str
            - 'transcript': str
            - 'hit_exons': list of int
            - 'cnv_type': 'deletion' or 'duplication'

    Returns:
    -------
    list of dict
        Each dict contains:
            - 'gene': str
            - 'transcript': str
            - 'hit_exons': list of int
            - 'predicted_consequence': str
    """
    results = []

    for hit in exon_hits:
        gene = hit.get("gene")
        transcript = hit.get("transcript")
        exons = hit.get("hit_exons", [])
        cnv_type = hit.get("cnv_type", "").lower()

        # Determine consequence
        if not exons:
            consequence = "no coding impact"
        else:
            # Simple rule: 1 exon affected -> frameshift or in-frame
            # Multi-exon: could be larger frameshift, depends on exon count (simplified)
            if cnv_type == "deletion":
                if len(exons) == 1:
                    consequence = "frameshift deletion"
                else:
                    consequence = f"frameshift deletion ({len(exons)} exons)"
            elif cnv_type == "duplication":
                if len(exons) == 1:
                    consequence = "frameshift duplication"
                else:
                    consequence = f"frameshift duplication ({len(exons)} exons)"
            else:
                consequence = f"{cnv_type} affecting {len(exons)} exon(s)"

            # Placeholder for CDS-aware logic
            # TODO: Check if exons contain CDS. If only UTR, adjust consequence to 'non-coding impact'
            # TODO: If CNV preserves reading frame exactly, mark as 'in-frame' instead of frameshift

        results.append({
            "gene": gene,
            "transcript": transcript,
            "hit_exons": exons,
            "predicted_consequence": consequence
        })

    return results


# Example usage for testing
if __name__ == "__main__":
    test_hits = [
        {"gene": "BRCA1", "transcript": "NM_007294.4", "hit_exons": [23], "cnv_type": "deletion"},
        {"gene": "DMD", "transcript": "NM_004006.3", "hit_exons": [45, 46], "cnv_type": "duplication"},
        {"gene": "TP53", "transcript": "NM_000546.6", "hit_exons": [], "cnv_type": "deletion"}
    ]

    results = predict_cnv_effect(test_hits)
    for r in results:
        print(
            f"Gene: {r['gene']}, Transcript: {r['transcript']}, "
            f"Exons hit: {r['hit_exons']}, Consequence: {r['predicted_consequence']}"
        )
