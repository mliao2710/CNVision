"""
CNVision Functional Predictor - Consequence Classification Logic

Purpose:
- Convert mapped exon-hit payloads into interpretable functional consequences.

Core Functions:
- Normalize exon-hit structures into a consistent `{exon, length}` format.
- Compute total coding length impact across all hit exons.
- Apply modulo-3 frame rule to classify in-frame vs frameshift outcomes.
- Return standardized consequence labels for deletions/duplications, including
  no-CDS-impact cases.
"""

from typing import List, Dict, Any

def predict_cnv_effect(
    exon_hits: List[Dict[str, Any]],
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]] = None
) -> List[Dict[str, Any]]:
    """Predict consequence labels for mapped exon hits."""
    results = []

    # STEP 1: Process each CNV hit.
    for hit in exon_hits:
        gene = hit.get("gene")
        transcript = hit.get("transcript")
        cnv_type = hit.get("cnv_type", "").lower()
        raw_exons = hit.get("hit_exons", [])

        # STEP 2: Normalize exon payload to {exon, length}.
        consequence = "unknown coding impact (NM_XXXXX or non-MANE gene)"
        hit_exons = []

        if raw_exons and mane_data:
            mane_gene_data = mane_data.get(gene)
            if mane_gene_data:
                mane_tx_data = mane_gene_data.get(transcript)
                if mane_tx_data:
                    for e in raw_exons:
                        if isinstance(e, int):
                            exon_info = next((ex for ex in mane_tx_data if ex["exon"] == e), None)
                            if exon_info:
                                hit_exons.append({
                                    "exon": e,
                                    "length": exon_info.get("cds_length", exon_info["end"] - exon_info["start"] + 1)
                                })
                        elif isinstance(e, dict) and "length" in e:
                            hit_exons.append(e)

        # Fallback when hit payload already contains lengths.
        if not hit_exons:
            if raw_exons and isinstance(raw_exons[0], dict) and "length" in raw_exons[0]:
                hit_exons = raw_exons

        # STEP 3: Classify frame effect.
        if hit_exons:
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
            
            frame_status = "in-frame" if total_length % 3 == 0 else "frameshift"
            num_exons = len(hit_exons)

            # STEP 4: Build human-readable consequence label.
            if cnv_type == "deletion":
                consequence = f"{frame_status} deletion ({num_exons} exon{'s' if num_exons>1 else ''})"
            elif cnv_type == "duplication":
                consequence = f"{frame_status} duplication ({num_exons} exon{'s' if num_exons>1 else ''})"
            else:
                consequence = f"{frame_status} {cnv_type} ({num_exons} exon{'s' if num_exons>1 else ''})"

        # STEP 5: Append normalized result.
        results.append({
            "gene": gene,
            "transcript": transcript,
            "hit_exons": hit_exons if hit_exons else raw_exons,
            "predicted_consequence": consequence
        })

    return results
