from flask import Flask, render_template, request
from coordinate_mapper import map_cnv_to_exons, map_exon_numbers_to_regions
from functional_predictor import predict_cnv_effect
from utils.ncbi_fetcher import fetch_refseq_info
from mane_loader import load_mane_exons
import os

app = Flask(__name__)

# Load MANE once
repo_root = os.path.dirname(os.path.abspath(__file__))
mane_path = os.path.join(repo_root, "data", "MANE.GRCh38.v1.4.refseq_genomic.gff.gz")
# Load MANE once
mane_data = load_mane_exons(mane_path)
# Build a case-insensitive index for gene lookup so web input isn't sensitive to
# letter case. Maps lower-case gene symbol -> canonical gene key in mane_data.
_mane_gene_index = {g.lower(): g for g in mane_data}

# Examples for template
EXAMPLES = [
    {"gene": "BRCA1", "start": 43000000, "end": 43002000, "type": "deletion"},
    {"gene": "TP53", "first_exon": 2, "last_exon": 4, "type": "duplication"},
    {"gene": "NM_000546", "start": 7565097, "end": 7590856, "type": "deletion"}
]


def process_gene_input(gene_input, cnv_type, mode, start=None, end=None, first_exon=None, last_exon=None):
    """
    Handles either gene symbol or NM_XXXXX input
    """
    if gene_input.upper().startswith("NM_"):
        refseq_info = fetch_refseq_info(gene_input.upper())
        if not refseq_info:
            return [{"error": f"Could not fetch information for {gene_input}"}]
        gene = refseq_info["gene_symbol"]
        transcript = refseq_info.get("transcript_id")
    else:
        gene = gene_input
        # Resolve case-insensitively
        canonical = _mane_gene_index.get(gene.lower())
        if canonical:
            gene = canonical
        else:
            return [{"error": f"Gene {gene_input} not found in MANE data."}]
        transcript = list(mane_data[gene].keys())[0]

    if mode == "coordinate":
        cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type}
    elif mode == "exon":
        cnv_region = map_exon_numbers_to_regions(gene, transcript, first_exon, last_exon, mane_data)
        if cnv_region is None:
            return [{"error": "Invalid exon numbers"}]
        cnv = {"gene": gene, "start": cnv_region[0], "end": cnv_region[1], "type": cnv_type}
    else:
        return [{"error": "Invalid mode"}]

    # If gene is None (e.g., NM_ lookup returned no gene symbol), try to
    # resolve the gene by searching for the transcript in MANE data.
    if not gene and transcript:
        for g, txs in mane_data.items():
            if transcript in txs:
                gene = g
                # update cnv with resolved gene
                cnv["gene"] = gene
                break

    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    # Ensure each hit dict contains the gene (predictor expects it)
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene

    # Pass mane_data into predictor so it can look up exon lengths and compute frames
    effects = predict_cnv_effect(hits, mane_data)
    return effects

@app.route("/", methods=["GET", "POST"])
def index():
    result = None
    error = None
    if request.method == "POST":
        gene_input = request.form.get("gene_input", "").strip()
        cnv_type = request.form.get("cnv_type", "").strip().lower()
        mode = request.form.get("mode", "coordinate")

        try:
            if mode == "coordinate":
                start = int(request.form.get("start", 0))
                end = int(request.form.get("end", 0))
                result = process_gene_input(gene_input, cnv_type, mode, start=start, end=end)
            elif mode == "exon":
                first_exon = int(request.form.get("first_exon", 0))
                last_exon = int(request.form.get("last_exon", 0))
                result = process_gene_input(gene_input, cnv_type, mode, first_exon=first_exon, last_exon=last_exon)
            else:
                error = "Invalid mode selected"
        except ValueError:
            error = "Invalid numeric input"
        except Exception as e:
            error = f"Unexpected error: {str(e)}"

    return render_template("index.html", result=result, error=error, examples=EXAMPLES)

if __name__ == "__main__":
    app.run(debug=True)
