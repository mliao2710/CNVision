import gzip

def load_mane_exons(gff_path):
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

            # Flexible attribute extraction: MANE GFF uses `gene`, `transcript_id`, and
            # exon lines often lack `exon_number`. We'll accept multiple key names
            # and fall back to parsing the ID/Parent when needed.
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
                    # Parent may be like 'rna-<id>' or contain multiple parents separated by commas
                    transcript = parent.split(",")[0]
                    for prefix in ("rna-", "transcript-", "rna:", "transcript:"):
                        if transcript.startswith(prefix):
                            transcript = transcript[len(prefix):]
                            break

            exon_number = attr_dict.get("exon_number")
            if exon_number is None:
                # Try to parse exon number from ID (e.g. ID=exon-NM_001005484.2-1)
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

            # Need at least gene and transcript to proceed
            if not gene or not transcript:
                continue

            # Convert coordinates
            start = int(start)
            end = int(end)

            # Build dict structure
            exon_data.setdefault(gene, {})
            exon_data[gene].setdefault(transcript, [])
            exon_data[gene][transcript].append({
                "exon": exon_number,
                "start": start,
                "end": end,
                "strand": strand,
            })

    # Sort exon lists
    for gene in exon_data:
        for transcript in exon_data[gene]:
            # Sort exon lists: prefer explicit exon number when available,
            # otherwise sort by genomic start coordinate.
            exon_data[gene][transcript] = sorted(
                exon_data[gene][transcript],
                key=lambda x: (x["exon"] if x["exon"] is not None else x["start"]),
            )

    return exon_data


if __name__ == "__main__":
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(description="Load MANE exon GFF and summarize contents")
    parser.add_argument(
        "gff",
        nargs="?",
        default="data/MANE.GRCh38.v1.4.refseq_genomic.gff.gz",
        help="Path to MANE GFF.gz file (default: data/MANE.GRCh38.v1.4.refseq_genomic.gff.gz)",
    )
    args = parser.parse_args()

    # Resolve default/relative paths relative to this script's directory so the
    # loader works regardless of the current working directory when invoked.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gff_path = args.gff
    if not os.path.isabs(gff_path):
        gff_path = os.path.join(script_dir, gff_path)
    if not os.path.exists(gff_path):
        # Try common subdirectory `data/mane/` (some users unpack into a subfolder)
        basename = os.path.basename(gff_path)
        alt = os.path.join(os.path.dirname(gff_path), "mane", basename)
        if os.path.exists(alt):
            print(f"Using GFF found at: {alt}")
            gff_path = alt
        else:
            # Search under the data/ tree for a matching basename
            found = None
            root_search = os.path.dirname(gff_path) or "."
            for root, _, files in os.walk(root_search):
                if basename in files:
                    found = os.path.join(root, basename)
                    break

            if found:
                print(f"Using GFF found at: {found}")
                gff_path = found
            else:
                tried = os.path.abspath(gff_path)
                print(f"Error: GFF file not found: {gff_path}\nTried absolute path: {tried}")
                print("Either place the file at the default path or run:\n  python3 mane_loader.py /full/path/to/MANE.GRCh38.v1.4.refseq_genomic.gff.gz")
                sys.exit(2)

    data = load_mane_exons(gff_path)
    print("Loaded genes:", len(data))
