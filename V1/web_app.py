"""
================================================================================
WEB_APP.PY - CNVision Web Application Main Entry Point
================================================================================

PURPOSE:
--------
This is the main Flask web application that provides a user-friendly web interface
for CNVision, a computational framework to predict the functional impact of 
intragenic Copy Number Variants (CNVs).

BIOLOGICAL IMPACT:
------------------
CNVs (deletions or duplications of DNA segments) can cause genetic diseases. This
tool helps researchers and clinicians understand whether a CNV will:
- Cause a "frameshift" (disrupts protein reading frame, usually severe)
- Be "in-frame" (maintains reading frame, may be less severe)
- Affect specific exons (protein-coding regions) of genes

By predicting these effects, clinicians can better interpret genetic test results
and make informed decisions about patient care.

CODING CONCEPTS & PRINCIPLES:
------------------------------
1. WEB FRAMEWORK (Flask):
   - Flask is a lightweight Python web framework
   - Routes handle HTTP requests (GET for viewing, POST for submitting forms)
   - Templates render HTML with dynamic data

2. MODULAR DESIGN:
   - Separates concerns: web logic here, data processing in other modules
   - Each module has a single responsibility (separation of concerns)

3. DATA FLOW:
   User Input → Form Processing → Gene Lookup → Coordinate Mapping → 
   Exon Analysis → Functional Prediction → Results Display

4. ERROR HANDLING:
   - Try/except blocks catch errors gracefully
   - Returns user-friendly error messages instead of crashing

5. DATA STRUCTURES:
   - Dictionaries store gene/exon data: {gene: {transcript: [exons]}}
   - Lists store results: [{gene, transcript, consequence}, ...]

6. CASE-INSENSITIVE LOOKUP:
   - Builds index mapping lowercase → canonical gene names
   - Makes the tool more user-friendly (BRCA1 = brca1 = Brca1)

================================================================================
"""

from flask import Flask, render_template, request
from coordinate_mapper import map_cnv_to_exons, map_exon_numbers_to_regions
from functional_predictor import predict_cnv_effect
from utils.ncbi_fetcher import fetch_refseq_info
from mane_loader import load_mane_exons
from iscn_parser import find_genes_in_region, format_iscn_output
import os

# Initialize Flask application
# Flask is a web framework that handles HTTP requests and responses
app = Flask(__name__)

# ============================================================================
# INITIALIZATION: Load gene/exon data once at startup (not on every request)
# ============================================================================

# Get the directory where this script is located
# This allows us to find the data file regardless of where Python is run from
repo_root = os.path.dirname(os.path.abspath(__file__))
mane_path = os.path.join(repo_root, "data", "MANE.GRCh38.v1.4.refseq_genomic.gff.gz")

# Load MANE (Matched Annotation from NCBI and EMBL-EBI) data
# MANE provides high-quality, manually curated gene annotations
# This creates a nested dictionary: {gene_name: {transcript_id: [exon_list]}}
MANE_LOAD_ERROR = None
try:
    mane_data = load_mane_exons(mane_path)
    # Build a case-insensitive lookup index
    _mane_gene_index = {g.lower(): g for g in mane_data}
except FileNotFoundError:
    MANE_LOAD_ERROR = (
        f"MANE file not found at {mane_path}.\n"
        "Mount your local `data` folder to /app/data when running the container,\n"
        "for example: `-v $(pwd)/data:/app/data`"
    )
    mane_data = {}
    _mane_gene_index = {}
except Exception as e:
    MANE_LOAD_ERROR = f"Error loading MANE data: {str(e)}"
    mane_data = {}
    _mane_gene_index = {}

# ============================================================================
# EXAMPLE DATA: Pre-defined examples to show users how to use each input mode
# ============================================================================
# These examples demonstrate the two input modes:
# 1. Genomic coordinates (start/end positions with chromosome and build)
# 2. Exon numbers (which exons are affected)
EXAMPLES = [
    {"gene": "BRCA1", "chromosome": "17", "start": 43000000, "end": 43002000, "build": "GRCh38", "type": "deletion"},
    {"gene": "TP53", "first_exon": 2, "last_exon": 4, "type": "duplication"}
]


# ============================================================================
# CORE PROCESSING FUNCTIONS
# ============================================================================

def process_gene_input(gene_input, cnv_type, mode, start=None, end=None, first_exon=None, last_exon=None, chromosome=None, build=None):
    """
    Main processing function for standard input modes (genomic coordinates or exon numbers).
    
    This function orchestrates the entire analysis pipeline:
    1. Resolve gene name (handle gene symbols or RefSeq IDs)
    2. Convert input format to genomic coordinates
    3. Map coordinates to exons
    4. Predict functional consequences
    
    Args:
        gene_input: Gene symbol (e.g., "BRCA1") or RefSeq ID (e.g., "NM_007294.4")
        cnv_type: "deletion" or "duplication"
        mode: "coordinate" or "exon"
        start/end: Genomic coordinates (for coordinate mode)
        first_exon/last_exon: Exon numbers (for exon mode)
    
    Returns:
        List of dictionaries with analysis results, or error dicts
    """
    # ========================================================================
    # STEP 1: Resolve gene name and transcript
    # ========================================================================
    # Users can input either:
    # - Gene symbol: "BRCA1", "TP53"
    # - RefSeq transcript ID: "NM_007294.4" (starts with "NM_")
    
    if gene_input.upper().startswith("NM_"):
        # User provided a RefSeq transcript ID (e.g., "NM_007294.4")
        # Fetch gene information from NCBI database
        refseq_info = fetch_refseq_info(gene_input.upper())
        if not refseq_info:
            return [{"error": f"Could not fetch information for {gene_input}"}]
        gene = refseq_info["gene_symbol"]
        transcript = refseq_info.get("transcript_id")
    else:
        # User provided a gene symbol (e.g., "BRCA1")
        gene = gene_input
        # Resolve case-insensitively using our lookup index
        # Example: "brca1" -> "BRCA1"
        canonical = _mane_gene_index.get(gene.lower())
        if canonical:
            gene = canonical
        else:
            # More helpful guidance when a gene is not present in MANE
            # Common reasons:
            # - The gene is not covered in the MANE GFF file used by this tool
            # - User provided GRCh37 coordinates while MANE uses GRCh38
            # - The user supplied a gene symbol instead of a RefSeq (NM_) transcript
            # Provide actionable next steps instead of a terse error.
            return [{
                "error": (
                    f"Gene '{gene_input}' not found in MANE data. "
                    "This prevents mapping genomic coordinates to MANE exons, so the tool cannot compute coding impact. "
                    "Possible fixes: provide a RefSeq transcript ID (e.g., NM_XXXXX), use GRCh38 coordinates, "
                    "or enable/implement an NCBI fallback to fetch transcript exon coordinates."
                )
            }]
        # Get the first (primary) transcript for this gene
        transcript = list(mane_data[gene].keys())[0]

    # ========================================================================
    # STEP 2: Convert input to standardized CNV format
    # ========================================================================
    # We need genomic coordinates (start/end) regardless of input mode
    # This standardizes the data for downstream processing
    
    if mode == "coordinate":
        # User already provided genomic coordinates
        cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type}
    elif mode == "exon":
        # User provided exon numbers - convert to genomic coordinates
        # Example: exons 2-4 -> genomic coordinates spanning those exons
        cnv_region = map_exon_numbers_to_regions(gene, transcript, first_exon, last_exon, mane_data)
        if cnv_region is None:
            return [{"error": "Invalid exon numbers"}]
        cnv = {"gene": gene, "start": cnv_region[0], "end": cnv_region[1], "type": cnv_type}
    else:
        return [{"error": "Invalid mode"}]

    # ========================================================================
    # STEP 3: Fallback gene resolution (if needed)
    # ========================================================================
    # Sometimes RefSeq lookup doesn't return a gene symbol
    # Search MANE data to find which gene has this transcript
    if not gene and transcript:
        for g, txs in mane_data.items():
            if transcript in txs:
                gene = g
                cnv["gene"] = gene
                break
    # If gene isn't present in MANE but we have refseq info (NM_) fetched earlier,
    # build a temporary fallback mane-like structure so we can map exons and predict
    # consequences without mutating the global `mane_data`.
    used_mane_data = mane_data
    if gene not in mane_data:
        # try to use fetched refseq_info if available in this scope
        # Note: refseq_info is defined only when user provided NM_ id
        try:
            if refseq_info and isinstance(refseq_info.get("exons"), list):
                used_mane_data = {gene: {transcript: refseq_info["exons"]}}
        except NameError:
            # refseq_info not defined in this scope; keep using mane_data
            used_mane_data = mane_data

    # ========================================================================
    # STEP 4: Map CNV coordinates to affected exons
    # ========================================================================
    # Find which exons overlap with the CNV region
    # Returns list of exons that are hit by this CNV
    hits = map_cnv_to_exons(cnv, used_mane_data, transcript)
    
    # Ensure each hit dict contains the gene (predictor expects it)
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene

    # ========================================================================
    # STEP 5: Predict functional consequences
    # ========================================================================
    # Analyze the exons to determine if CNV causes frameshift or in-frame effect
    # This is the core biological prediction
    effects = predict_cnv_effect(hits, used_mane_data)
    
    # Add ISCN metadata ONLY for coordinate mode (not exon mode)
    if mode == "coordinate" and chromosome and build:
        for effect in effects:
            effect["iscn_chromosome"] = chromosome
            effect["iscn_build"] = build
            effect["iscn_start"] = start
            effect["iscn_end"] = end
            effect["iscn_notation"] = format_iscn_output(chromosome, start, end, build, cnv_type)
    
    return effects


def process_iscn_input(gene_input, chromosome, start, end, build, cnv_type):
    """
    Process ISCN (International System for Human Cytogenomic Nomenclature) format input.
    
    ISCN is a standardized notation used in clinical genetics labs. This function
    handles ISCN-specific input while using the same analysis pipeline as other modes.
    
    Key difference: ISCN includes chromosome and genome build information, which
    provides more context for clinical interpretation.
    
    Args:
        gene_input: Gene symbol or NM_XXXXX ID
        chromosome: Chromosome number (as string, e.g., "17", "X", "Y")
        start: Genomic start coordinate
        end: Genomic end coordinate
        build: Genomic build version ("GRCh37" or "GRCh38")
        cnv_type: "deletion" or "duplication"
    
    Returns:
        List of effect dictionaries (same format as process_gene_input)
        with additional ISCN metadata (chromosome, build)
    """
    # ========================================================================
    # STEP 1: Validate genome build
    # ========================================================================
    # MANE data uses GRCh38 coordinates
    # GRCh37 coordinates would need conversion (not implemented yet)
    if build not in ["GRCh37", "GRCh38"]:
        return [{"error": f"Invalid genomic build: {build}. Supported: GRCh37, GRCh38"}]
    
    # Note: GRCh37 coordinates might not align perfectly with GRCh38 MANE data
    # This is a limitation that could be addressed with coordinate conversion
    
    # ========================================================================
    # STEP 2: Resolve gene (same as process_gene_input)
    # ========================================================================
    # Process gene input same as coordinate mode
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
            # Helpful guidance when ISCN-provided gene is not in MANE
            return [{
                "error": (
                    f"Gene '{gene_input}' not found in MANE data. "
                    "For ISCN/genomic coordinate analysis this tool relies on MANE exon coordinates (GRCh38). "
                    "Try supplying a RefSeq transcript (NM_...), use GRCh38 coordinates, or enable NCBI lookup as a fallback."
                )
            }]
        transcript = list(mane_data[gene].keys())[0]
    
    # Create CNV dict with ISCN coordinates
    cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type}
    
    # If gene is None, try to resolve
    if not gene and transcript:
        for g, txs in mane_data.items():
            if transcript in txs:
                gene = g
                cnv["gene"] = gene
                break
    
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    # Ensure each hit dict contains the gene
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene
    
    # ========================================================================
    # STEP 3-5: Same analysis pipeline as process_gene_input
    # ========================================================================
    # Map to exons and predict consequences (same logic as coordinate mode)
    # If gene isn't in MANE, but we fetched refseq exons above, use them as a
    # temporary fallback dataset for mapping/prediction (do not mutate global data).
    used_mane_data = mane_data
    if gene not in mane_data:
        try:
            if refseq_info and isinstance(refseq_info.get("exons"), list):
                used_mane_data = {gene: {transcript: refseq_info["exons"]}}
        except NameError:
            used_mane_data = mane_data

    hits = map_cnv_to_exons(cnv, used_mane_data, transcript)
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene

    effects = predict_cnv_effect(hits, used_mane_data)
    
    # ========================================================================
    # STEP 6: Add ISCN-specific metadata to results
    # ========================================================================
    # This helps users track which chromosome and build were used
    for effect in effects:
        effect["iscn_chromosome"] = chromosome
        effect["iscn_build"] = build
    
    return effects

# ============================================================================
# WEB ROUTE: Handle HTTP requests
# ============================================================================

@app.route("/", methods=["GET", "POST"])
def index():
    """
    Main web route handler for the CNVision application.
    
    Handles two types of requests:
    - GET: User visits the page (show form)
    - POST: User submits CNV data (process and show results)
    
    Flask automatically calls this function when user visits the root URL ("/")
    """
    # Initialize variables with safe defaults
    # Empty list allows template to iterate without errors if no results yet
    result = []
    error = None
    
    # ========================================================================
    # PROCESS FORM SUBMISSION (POST request)
    # ========================================================================
    if request.method == "POST":
        # Extract form data from HTTP POST request
        # .get() returns default value if field is missing (prevents crashes)
        # .strip() removes whitespace (makes input more robust)
        gene_input = request.form.get("gene_input", "").strip()
        cnv_type = request.form.get("cnv_type", "").strip().lower()
        mode = request.form.get("mode", "exon")  # Default to "exon" mode

        # ====================================================================
        # Route to appropriate processing function based on input mode
        # ====================================================================
        try:
            if mode == "coordinate":
                # Genomic coordinate mode: user provides chromosome, build, start/end positions
                chromosome = request.form.get("chromosome", "").strip()
                start = int(request.form.get("start", 0))
                end = int(request.form.get("end", 0))
                build = request.form.get("build", "GRCh38").strip()
                result = process_gene_input(gene_input, cnv_type, mode, start=start, end=end, chromosome=chromosome, build=build)
                
            elif mode == "exon":
                # Exon number mode: user provides exon numbers
                first_exon_val = request.form.get("first_exon", "")
                last_exon_val = request.form.get("last_exon", "")

                try:
                    first_exon = int(first_exon_val)
                    last_exon = int(last_exon_val)
                except ValueError:
                    error = "Please provide valid first and last exon numbers."
                    result = []
                else:
                    if not first_exon or not last_exon:
                        error = "Please provide valid first and last exon numbers."
                        result = []
                    else:
                        result = process_gene_input(gene_input, cnv_type, mode, first_exon=first_exon, last_exon=last_exon)
            else:
                error = "Invalid mode selected"
                
        except ValueError:
            # User entered non-numeric value where number expected
            error = "Invalid numeric input"
        except Exception as e:
            # Catch any unexpected errors and show friendly message
            error = f"Unexpected error: {str(e)}"

    # ========================================================================
    # RENDER TEMPLATE: Send data to HTML template for display
    # ========================================================================
    # Flask's render_template combines HTML template with Python data
    # Template can access: result (analysis results), error (if any), examples
    return render_template("index.html", result=result, error=error, examples=EXAMPLES)


# ============================================================================
# APPLICATION ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    # Only run when executed directly. Bind to 0.0.0.0 so the container port is reachable.
    app.run(host="0.0.0.0", port=5000, debug=True)
