def predict_cnv_effect(exon_hits, cnv_type="deletion"):
    effects = []

    for transcript, hit_info in exon_hits.items():
        hit_exons = hit_info["hit_exons"]

        if not hit_exons:
            consequence = "no coding impact"
        else:
            # frameshift if coding length disrupted
            if hit_info["coding_bases_lost"] % 3 != 0:
                consequence = f"frameshift {cnv_type}"
            else:
                consequence = f"in-frame {cnv_type}"

        effects.append({
            "transcript": transcript,
            "hit_exons": hit_exons,
            "predicted_consequence": consequence
        })

    return effects
