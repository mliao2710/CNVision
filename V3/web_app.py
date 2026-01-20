from flask import Flask, render_template, request
from coordinate_mapper import map_cnv_to_exons, map_exon_numbers_to_regions
from functional_predictor import predict_cnv_effect
from mane_loader import load_mane_exons
from iscn_parser import find_genes_in_region, format_iscn_output
import os
import re

# Try to import ncbi_fetcher if it exists, otherwise define a stub
try:
    from utils.ncbi_fetcher import fetch_refseq_info
except ImportError:
    def fetch_refseq_info(transcript_id):
        """Stub function if ncbi_fetcher doesn't exist."""
        return None

app = Flask(__name__)

# Load BOTH GRCh37 and GRCh38 data
repo_root = os.path.dirname(os.path.abspath(__file__))
mane_path_grch38 = os.path.join(repo_root, "data", "MANE.GRCh38.v1.4.refseq_genomic.gff.gz")
mane_path_grch37 = os.path.join(repo_root, "data", "GRCh37.ref_seq_select.gz")

MANE_LOAD_ERROR = None
mane_data_grch38 = {}
mane_data_grch37 = {}
_mane_gene_index_grch38 = {}
_mane_gene_index_grch37 = {}

try:
    # Load GRCh38 data
    mane_data_grch38 = load_mane_exons(mane_path_grch38)
    _mane_gene_index_grch38 = {g.lower(): g for g in mane_data_grch38}
    print(f"✓ Loaded GRCh38 MANE data: {len(mane_data_grch38)} genes")
except FileNotFoundError:
    MANE_LOAD_ERROR = f"GRCh38 MANE file not found at {mane_path_grch38}."
except Exception as e:
    MANE_LOAD_ERROR = f"Error loading GRCh38 MANE data: {str(e)}"

try:
    # Load GRCh37 data
    mane_data_grch37 = load_mane_exons(mane_path_grch37)
    _mane_gene_index_grch37 = {g.lower(): g for g in mane_data_grch37}
    print(f"✓ Loaded GRCh37 MANE data: {len(mane_data_grch37)} genes")
except FileNotFoundError:
    if not MANE_LOAD_ERROR:
        MANE_LOAD_ERROR = f"GRCh37 MANE file not found at {mane_path_grch37}."
except Exception as e:
    if not MANE_LOAD_ERROR:
        MANE_LOAD_ERROR = f"Error loading GRCh37 MANE data: {str(e)}"

# Helper function to get the correct dataset based on build
def get_mane_data(build):
    """Get the appropriate MANE dataset based on genomic build."""
    if build == "GRCh37":
        return mane_data_grch37, _mane_gene_index_grch37
    else:  # Default to GRCh38
        return mane_data_grch38, _mane_gene_index_grch38


# ============================================================================
# NEW: HGVS cDNA Parser
# ============================================================================

def parse_cdna_notation(cdna_string):
    """
    Parse HGVS cDNA notation into gene/transcript and coordinates.
    
    Supported formats:
    - FBN1:c.5065_5546del
    - c.(5065+1_5066-1)_(5545+1_5546-1)
    - NM_001406716.1:c.5065_5546del
    
    Returns:
        dict with keys: gene, transcript, start, end, variant_type
        or None if parsing fails
    """
    # Remove whitespace
    cdna_string = cdna_string.strip()
    
    # Pattern 1: Gene:c.start_end[del/dup]
    # Example: FBN1:c.5065_5546del
    pattern1 = r'([A-Z0-9]+):c\.(\d+)_(\d+)(del|dup)'
    match1 = re.match(pattern1, cdna_string, re.IGNORECASE)
    
    if match1:
        gene_or_nm = match1.group(1)
        start = int(match1.group(2))
        end = int(match1.group(3))
        var_type = match1.group(4).lower()
        
        # Check if it's a transcript ID or gene symbol
        if gene_or_nm.upper().startswith('NM_'):
            return {
                'transcript': gene_or_nm.upper(),
                'gene': None,  # Will be resolved later
                'start': start,
                'end': end,
                'variant_type': 'deletion' if var_type == 'del' else 'duplication'
            }
        else:
            return {
                'gene': gene_or_nm.upper(),
                'transcript': None,
                'start': start,
                'end': end,
                'variant_type': 'deletion' if var_type == 'del' else 'duplication'
            }
    
    # Pattern 2: c.(start+1_end-1)_(start+1_end-1) - intronic notation
    # Example: c.(5065+1_5066-1)_(5545+1_5546-1)
    # This means the breakpoints are in introns, so we use the exon boundaries
    pattern2 = r'c\.\((\d+)\+\d+_(\d+)-\d+\)_\((\d+)\+\d+_(\d+)-\d+\)'
    match2 = re.match(pattern2, cdna_string, re.IGNORECASE)
    
    if match2:
        # For intronic breakpoints, use the exon end/start as boundaries
        start = int(match2.group(1))
        end = int(match2.group(4))
        return {
            'gene': None,
            'transcript': None,
            'start': start,
            'end': end,
            'variant_type': 'deletion',  # Default, should be specified separately
            'intronic': True
        }
    
    return None


# ============================================================================
# NEW: Genomic Coordinate Parser
# ============================================================================

def parse_genomic_input(genomic_string):
    """
    Parse genomic coordinate notation.
    
    Supported formats:
    - DEL chr15:48,741,090-48,756,096
    - chr15:48,741,090-48,756,096
    - arr[GRCh37]15q21.1(48,741,090_48,756,096)x1
    - 15:48741090-48756096
    
    Returns:
        dict with keys: chromosome, start, end, build (if specified)
        or None if parsing fails
    """
    genomic_string = genomic_string.strip()
    
    # Pattern 1: DEL chr15:48,741,090-48,756,096 or chr15:48,741,090-48,756,096
    pattern1 = r'(?:DEL|DUP)?\s*chr([0-9XY]+):([0-9,]+)-([0-9,]+)'
    match1 = re.search(pattern1, genomic_string, re.IGNORECASE)
    
    if match1:
        chromosome = match1.group(1)
        start = int(match1.group(2).replace(',', ''))
        end = int(match1.group(3).replace(',', ''))
        return {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'build': None  # Will use form input
        }
    
    # Pattern 2: Just the coordinates without 'chr' prefix - 15:48741090-48756096
    pattern1b = r'(?:DEL|DUP)?\s*([0-9XY]+):([0-9,]+)-([0-9,]+)'
    match1b = re.search(pattern1b, genomic_string, re.IGNORECASE)
    
    if match1b:
        chromosome = match1b.group(1)
        start = int(match1b.group(2).replace(',', ''))
        end = int(match1b.group(3).replace(',', ''))
        return {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'build': None  # Will use form input
        }
    
    # Pattern 3: arr[GRCh37]15q21.1(48,741,090_48,756,096)x1 or arr[GRCh37]15q21.1( 48,741,090_48,756,096)x1
    pattern2 = r'arr\[([^\]]+)\]([0-9XY]+)[pq][0-9.]+\(\s*([0-9,_\s]+)\)x([0-9])'
    match2 = re.search(pattern2, genomic_string, re.IGNORECASE)
    
    if match2:
        build = match2.group(1)
        chromosome = match2.group(2)
        coords = match2.group(3).replace(',', '').replace(' ', '').split('_')
        copy_number = match2.group(4)
        
        if len(coords) == 2:
            start = int(coords[0])
            end = int(coords[1])
            return {
                'chromosome': chromosome,
                'start': start,
                'end': end,
                'build': build,
                'copy_number': int(copy_number)
            }
    
    return None


# ============================================================================
# Processing Functions for Each Mode
# ============================================================================

def process_exon_mode(gene_input, cnv_type, first_exon, last_exon, start_pos=None, end_pos=None, build="GRCh38"):
    """
    Process exon mode input with optional positions within exons.
    
    Args:
        gene_input: Gene symbol
        cnv_type: deletion or duplication
        first_exon: First exon number
        last_exon: Last exon number
        start_pos: Optional position in first exon (1 = first base, -1 = last base)
        end_pos: Optional position in last exon (1 = first base, -1 = last base)
        build: GRCh37 or GRCh38
    """
    # Get correct dataset for build
    mane_data, _mane_gene_index = get_mane_data(build)
    
    # Resolve gene name
    gene = gene_input.upper()
    canonical = _mane_gene_index.get(gene.lower())
    if canonical:
        gene = canonical
    else:
        return [{"error": f"Gene '{gene_input}' not found in MANE {build} data"}]
    
    transcript = list(mane_data[gene].keys())[0]
    
    # Convert exon numbers to genomic coordinates
    cnv_region = map_exon_numbers_to_regions(gene, transcript, first_exon, last_exon, mane_data)
    if cnv_region is None:
        return [{"error": "Invalid exon numbers"}]
    
    # Get chromosome from first exon
    exons = mane_data[gene][transcript]
    chromosome = exons[0].get("chromosome", "1") if exons else "1"
    
    cnv = {"gene": gene, "start": cnv_region[0], "end": cnv_region[1], "type": cnv_type, "chromosome": chromosome}
    
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene
    
    effects = predict_cnv_effect(hits, mane_data)
    return effects


def process_transcript_mode(cdna_notation, cnv_type):
    """
    Process HGVS cDNA notation.
    
    Args:
        cdna_notation: HGVS c. notation string
        cnv_type: deletion or duplication (may be overridden by notation)
    """
    parsed = parse_cdna_notation(cdna_notation)
    
    if not parsed:
        return [{"error": "Could not parse cDNA notation. Please use format like: FBN1:c.5065_5546del"}]
    
    # Extract variant type from notation if present
    if 'variant_type' in parsed:
        cnv_type = parsed['variant_type']
    
    # Resolve gene and transcript
    gene = parsed.get('gene')
    transcript = parsed.get('transcript')
    
    if transcript and not gene:
        # Look up gene from transcript ID
        refseq_info = fetch_refseq_info(transcript)
        if refseq_info:
            gene = refseq_info['gene_symbol']
        else:
            # Search in MANE data
            for g, txs in mane_data.items():
                if transcript in txs:
                    gene = g
                    break
    
    if not gene:
        return [{"error": "Could not determine gene from cDNA notation"}]
    
    # Resolve canonical gene name
    canonical = _mane_gene_index.get(gene.lower())
    if canonical:
        gene = canonical
    else:
        return [{"error": f"Gene '{gene}' not found in MANE data"}]
    
    if not transcript:
        transcript = list(mane_data[gene].keys())[0]
    
    # TODO: Convert cDNA coordinates to genomic coordinates
    # For now, create CNV dict with cDNA positions
    # This requires mapping c. positions to genomic positions
    
    return [{"error": "Transcript mode not fully implemented yet. Coming soon!"}]


def process_genomic_mode(genomic_input, build, cnv_type, gene_hint=None):
    """
    Process genomic coordinate input.
    
    Args:
        genomic_input: Genomic coordinate string
        build: GRCh37 or GRCh38
        cnv_type: deletion or duplication
        gene_hint: Optional gene symbol to help identify target
    """
    parsed = parse_genomic_input(genomic_input)
    
    if not parsed:
        return [{"error": "Could not parse genomic coordinates. Use format like: DEL chr15:48,741,090-48,756,096"}]
    
    chromosome = parsed['chromosome']
    start = parsed['start']
    end = parsed['end']
    
    # Use build from ISCN notation if present, otherwise use form input
    if parsed.get('build'):
        build = parsed['build']
    
    # Get correct dataset for build
    mane_data, _mane_gene_index = get_mane_data(build)
    
    # Determine CNV type from copy number if present in ISCN
    if parsed.get('copy_number'):
        copy_num = parsed['copy_number']
        if copy_num > 2:
            cnv_type = 'duplication'
        elif copy_num < 2:
            cnv_type = 'deletion'
    
    # Try to find gene - either from hint or by searching region
    gene = None
    
    if gene_hint:
        # User provided gene hint
        gene = gene_hint.upper()
        canonical = _mane_gene_index.get(gene.lower())
        if canonical:
            gene = canonical
    else:
        # Try to find genes in this region automatically
        genes_in_region = find_genes_in_region(chromosome, start, end, mane_data)
        if genes_in_region:
            # Use the first gene found
            gene = genes_in_region[0]
        else:
            return [{"error": f"No genes found in region chr{chromosome}:{start:,}-{end:,} using {build}. Please provide a gene symbol as a hint."}]
    
    # Validate gene exists in MANE
    if gene not in mane_data:
        return [{"error": f"Gene '{gene}' not found in MANE {build} data"}]
    
    # Get transcript and analyze
    transcript = list(mane_data[gene].keys())[0]
    
    cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type, "chromosome": f"chr{chromosome}"}
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene
    
    effects = predict_cnv_effect(hits, mane_data)
    
    return effects


# ============================================================================
# Web Route
# ============================================================================

@app.route("/", methods=["GET", "POST"])
def index():
    result = []
    error = None
    
    if request.method == "POST":
        mode = request.form.get("mode", "exon")
        cnv_type = request.form.get("cnv_type", "deletion").strip().lower()
        build = request.form.get("build", "GRCh38").strip()
        
        try:
            if mode == "exon":
                gene_input = request.form.get("gene_input", "").strip()
                first_exon = int(request.form.get("first_exon", 0))
                last_exon = int(request.form.get("last_exon", 0))
                start_pos = request.form.get("exon_start_pos", "").strip()
                end_pos = request.form.get("exon_end_pos", "").strip()
                
                result = process_exon_mode(
                    gene_input, cnv_type, first_exon, last_exon,
                    start_pos if start_pos else None,
                    end_pos if end_pos else None,
                    build
                )
                
            elif mode == "transcript":
                cdna_notation = request.form.get("cdna_notation", "").strip()
                result = process_transcript_mode(cdna_notation, cnv_type, build)
                
            elif mode == "coordinate":
                genomic_input = request.form.get("genomic_input", "").strip()
                gene_hint = request.form.get("gene_hint", "").strip()
                
                result = process_genomic_mode(
                    genomic_input, build, cnv_type,
                    gene_hint if gene_hint else None
                )
            
        except ValueError as e:
            error = f"Invalid input: {str(e)}"
        except Exception as e:
            error = f"Error: {str(e)}"
            import traceback
            traceback.print_exc()
    
    # Update examples for new modes
    examples = [
        {"mode": "exon", "description": "FBN1 exons 42-45 deletion"},
        {"mode": "transcript", "description": "FBN1:c.5065_5546del"},
        {"mode": "coordinate", "description": "DEL chr15:48,741,090-48,756,096"}
    ]
    
    return render_template("index.html", result=result, error=error, examples=examples)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
