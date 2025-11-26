"""
mane_loader.py — Load MANE exon GFF and fetch NM_XXXXX transcripts outside MANE
"""

import gzip
import requests
from xml.etree import ElementTree as ET

def load_mane_exons(gff_path):
    """
    Load MANE exons from a GFF.gz file.

    Returns:
        dict: {gene: {transcript: [ {exon, start, end, strand, length}, ... ] } }
    """
    exon_data = {}

    with gzip.open(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields

            if feature_type != "exon":
                continue

            # Parse attributes
            attr_dict = {}
            for item in attrs.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attr_dict[key] = value

            gene = (
                attr_dict.get("gene")
                or attr_dict.get("gene_name")
                or attr_dict.get("gene_id")
                or attr_dict.get("Name")
            )

            transcript = attr_dict.get("transcript_id")
            if not transcript:
                parent = attr_dict.get("Parent") or attr_dict.get("parent")
                if parent:
                    transcript = parent.split(",")[0]
                    for prefix in ("rna-", "transcript-", "rna:", "transcript:"):
                        if transcript.startswith(prefix):
                            transcript = transcript[len(prefix):]
                            break

            exon_number = attr_dict.get("exon_number")
            if exon_number is None:
                idv = attr_dict.get("ID") or attr_dict.get("id")
                if idv and "-" in idv:
                    last = idv.split("-")[-1]
                    if last.isdigit():
                        exon_number = int(last)
                    else:
                        try:
                            exon_number = int(last.split("_")[-1])
                        except Exception:
                            exon_number = None
            else:
                try:
                    exon_number = int(exon_number)
                except Exception:
                    exon_number = None

            if not gene or not transcript:
                continue

            start = int(start)
            end = int(end)
            length = end - start + 1  # <-- NEW: compute exon length

            exon_data.setdefault(gene, {})
            exon_data[gene].setdefault(transcript, [])
            exon_data[gene][transcript].append({
                "exon": exon_number,
                "start": start,
                "end": end,
                "strand": strand,
                "length": length,   # <-- NEW FIELD
            })

    # Sort exons
    for gene in exon_data:
        for transcript in exon_data[gene]:
            exon_data[gene][transcript] = sorted(
                exon_data[gene][transcript],
                key=lambda x: (x["exon"] if x["exon"] is not None else x["start"])
            )

    return exon_data


def fetch_refseq_exons(nm_id):
    """
    Fetch exon coordinates for a RefSeq NM_XXXXX transcript from NCBI.

    Returns:
        list of dicts: [{'exon': 1, 'start': ..., 'end': ..., 'length': ...}, ...]
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": nm_id,
        "rettype": "gb",
        "retmode": "xml"
    }
    try:
        r = requests.get(url, params=params, timeout=10)
        r.raise_for_status()
    except Exception as e:
        print(f"Error fetching {nm_id}: {e}")
        return None

    exons = []
    try:
        tree = ET.fromstring(r.content)
        for feature in tree.findall(".//GBFeature"):
            if feature.findtext("GBFeature_key") == "exon":
                loc = feature.findtext("GBFeature_location")
                if not loc:
                    continue
                # Handle ranges like "123..456" or "<123..>456"
                start_str, end_str = loc.replace("<", "").replace(">", "").split("..")
                start = int(start_str)
                end = int(end_str)
                length = end - start + 1
                exons.append({
                    "exon": len(exons) + 1,
                    "start": start,
                    "end": end,
                    "strand": None,
                    "length": length  # <-- NEW FIELD
                })
    except Exception as e:
        print(f"Error parsing NCBI response for {nm_id}: {e}")
        return None

    return exons if exons else None


if __name__ == "__main__":
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(description="Load MANE exon GFF and summarize contents")
    parser.add_argument(
        "gff",
        nargs="?",
        default="data/MANE.GRCh38.v1.4.refseq_genomic.gff.gz",
        help="Path to MANE GFF.gz file"
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    gff_path = args.gff
    if not os.path.isabs(gff_path):
        gff_path = os.path.join(script_dir, gff_path)
    if not os.path.exists(gff_path):
        print(f"Error: GFF file not found: {gff_path}")
        sys.exit(2)

    data = load_mane_exons(gff_path)
    print(f"Loaded {len(data)} genes")

    # Example NM fetch
    nm_test = "NM_000546.6"
    exons = fetch_refseq_exons(nm_test)
    if exons:
        print(f"{nm_test} fetched {len(exons)} exons")
