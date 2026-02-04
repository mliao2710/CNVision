"""
CNVision Web Application - Flask REST API and Web Interface

Technologies: Python 3.10+, Flask (web framework), Jinja2 (template rendering),
HTML5/CSS3/JavaScript (frontend), Bootstrap (UI framework)

Core Functions:
- HTTP request routing and form processing (3 input modes)
- CNV analysis coordination (parsing, mapping exons, predicting effects)
- Result rendering and display via web browser
- NCBI API integration for transcript data fetching
"""

from flask import Flask, render_template, request
from coordinate_mapper import map_cnv_to_exons, map_exon_numbers_to_regions
from functional_predictor import predict_cnv_effect
from mane_loader import load_mane_exons
from iscn_parser import find_genes_in_region, format_iscn_output
import os
import re
import requests
import xml.etree.ElementTree as ET
import time
import warnings
warnings.filterwarnings("ignore", message="Unverified HTTPS request")

_ncbi_cache = {}

def fetch_refseq_info(transcript_id):
    """Fetch gene symbol and genomic exon coordinates from NCBI."""
    base_id = transcript_id.split('.')[0]
    
    if base_id in _ncbi_cache:
        return _ncbi_cache[base_id]
    
    NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    try:
        # ============================================================
        # STEP 1: Find the Gene ID from the transcript accession
        # ============================================================
        r = requests.get(f"{NCBI_BASE}/esearch.fcgi", params={
            "db": "gene",
            "term": f"{base_id}[Gene RefSeq]",
            "retmode": "xml"
        }, timeout=10, verify=False)
        r.raise_for_status()
        
        root = ET.fromstring(r.content)
        id_list = root.findall(".//IdList/Id")
        
        # Fallback search if first didn't find anything
        if not id_list:
            r = requests.get(f"{NCBI_BASE}/esearch.fcgi", params={
                "db": "gene",
                "term": f"{base_id}[RefSeq]",
                "retmode": "xml"
            }, timeout=10, verify=False)
            r.raise_for_status()
            root = ET.fromstring(r.content)
            id_list = root.findall(".//IdList/Id")
        
        if not id_list:
            return None
        
        gene_id = id_list[0].text
        time.sleep(0.4)
        
        # ============================================================
        # STEP 2: Fetch gene_table to get genomic exon coordinates
        # ============================================================
        r = requests.get(f"{NCBI_BASE}/efetch.fcgi", params={
            "db": "gene",
            "id": gene_id,
            "rettype": "gene_table",
            "retmode": "text"
        }, timeout=10, verify=False)
        r.raise_for_status()
        
        gene_table = r.text
        gene_symbol = None
        chromosome = None
        strand = None
        exons = []
        
        lines = gene_table.split('\n')
        in_our_exon_table = False
        past_header = False  # skip the dashed header line
        
        for line in lines:
            stripped = line.strip()
            
            # Gene symbol is first word on first non-empty line
            if not gene_symbol and stripped and not stripped.startswith("Gene ID"):
                gene_symbol = stripped.split()[0]
            
            # Extract strand and chromosome from reference line
            if "Primary Assembly" in stripped:
                if "minus strand" in stripped:
                    strand = "-"
                elif "plus strand" in stripped:
                    strand = "+"
                # Extract chromosome from NC_000017.11 -> chr17
                if "NC_" in stripped:
                    nc_part = stripped.split("NC_")[1]
                    nc_num = nc_part.split(".")[0]  # "000017"
                    chrom_num = str(int(nc_num))    # "17"
                    chromosome = f"chr{chrom_num}"
            
            # Detect our transcript's exon table
            if f"Exon table for" in stripped and base_id in stripped:
                in_our_exon_table = True
                past_header = False
                exons = []
                continue
            
            if in_our_exon_table:
                # Skip the column header and dashed line
                if "---" in stripped:
                    past_header = True
                    continue
                if not past_header:
                    continue
                
                # Empty line = end of this exon table
                if stripped == "":
                    if exons:
                        break
                    continue
                
                # Parse exon line - first field is "7687490-7687377" (genomic interval)
                # Fields are whitespace-separated
                parts = stripped.split()
                if len(parts) >= 1 and "-" in parts[0]:
                    try:
                        coord1, coord2 = parts[0].split("-")
                        start = min(int(coord1), int(coord2))
                        end = max(int(coord1), int(coord2))
                        
                        exons.append({
                            "exon": len(exons) + 1,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "length": end - start + 1,
                            "chromosome": chromosome
                        })
                    except (ValueError, IndexError):
                        continue
        
        if not gene_symbol:
            return None
        
        result = {
            "gene_symbol": gene_symbol,
            "transcript_id": base_id,
            "exons": exons,
            "chromosome": chromosome,
            "strand": strand
        }
        _ncbi_cache[base_id] = result
        return result
            
    except Exception:
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
_transcript_to_gene_grch38 = {}
_transcript_to_gene_grch37 = {}

try:
    mane_data_grch38 = load_mane_exons(mane_path_grch38)
    _mane_gene_index_grch38 = {g.lower(): g for g in mane_data_grch38}
    _transcript_to_gene_grch38 = {}
    for gene, transcripts in mane_data_grch38.items():
        for transcript in transcripts.keys():
            base_transcript = transcript.split('.')[0]
            _transcript_to_gene_grch38[transcript] = gene
            _transcript_to_gene_grch38[base_transcript] = gene
except FileNotFoundError:
    MANE_LOAD_ERROR = f"GRCh38 MANE file not found at {mane_path_grch38}."
    _transcript_to_gene_grch38 = {}
except Exception as e:
    MANE_LOAD_ERROR = f"Error loading GRCh38 MANE data: {str(e)}"
    _transcript_to_gene_grch38 = {}

try:
    mane_data_grch37 = load_mane_exons(mane_path_grch37)
    _mane_gene_index_grch37 = {g.lower(): g for g in mane_data_grch37}
    _transcript_to_gene_grch37 = {}
    for gene, transcripts in mane_data_grch37.items():
        for transcript in transcripts.keys():
            base_transcript = transcript.split('.')[0]
            _transcript_to_gene_grch37[transcript] = gene
            _transcript_to_gene_grch37[base_transcript] = gene
except FileNotFoundError:
    if not MANE_LOAD_ERROR:
        MANE_LOAD_ERROR = f"GRCh37 MANE file not found at {mane_path_grch37}."
        _transcript_to_gene_grch37 = {}
except Exception as e:
    if not MANE_LOAD_ERROR:
        MANE_LOAD_ERROR = f"Error loading GRCh37 MANE data: {str(e)}"
        _transcript_to_gene_grch37 = {}

def get_mane_data(build):
    """Get the appropriate MANE dataset based on genomic build."""
    if build == "GRCh37":
        return mane_data_grch37, _mane_gene_index_grch37, _transcript_to_gene_grch37
    else:
        return mane_data_grch38, _mane_gene_index_grch38, _transcript_to_gene_grch38


def parse_cdna_notation(cdna_string):
    """Parse HGVS cDNA notation."""
    cdna_string = cdna_string.strip()
    
    # Updated pattern to accept transcript IDs with underscores, periods, and version numbers
    # Matches: NM_001406716.1, NM_000138, or gene symbols like FBN1
    pattern1 = r'([A-Z0-9_\.]+):c\.(\d+)_(\d+)(del|dup)'
    match1 = re.match(pattern1, cdna_string, re.IGNORECASE)
    
    if match1:
        gene_or_nm = match1.group(1)
        start = int(match1.group(2))
        end = int(match1.group(3))
        var_type = match1.group(4).lower()
        
        if gene_or_nm.upper().startswith('NM_'):
            return {
                'transcript': gene_or_nm.upper(),
                'gene': None,
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
    
    pattern2 = r'c\.\((\d+)\+\d+_(\d+)-\d+\)_\((\d+)\+\d+_(\d+)-\d+\)'
    match2 = re.match(pattern2, cdna_string, re.IGNORECASE)
    
    if match2:
        start = int(match2.group(1))
        end = int(match2.group(4))
        return {
            'gene': None,
            'transcript': None,
            'start': start,
            'end': end,
            'variant_type': 'deletion',
            'intronic': True
        }
    
    return None


def parse_genomic_input(genomic_string):
    """Parse genomic coordinate notation with DEL/DUP detection."""
    genomic_string = genomic_string.strip().replace(' ', '')
    
    # Detect DEL or DUP prefix
    cnv_type_hint = None
    if genomic_string.upper().startswith('DEL'):
        cnv_type_hint = 'deletion'
    elif genomic_string.upper().startswith('DUP'):
        cnv_type_hint = 'duplication'
    
    # Helper: extract trailing build [GRCh38] or [GRCh37] and copy number x1/x3
    # These can appear at the end in any order: x1[GRCh38] or [GRCh38]x1
    def extract_suffix(s):
        build = None
        copy_number = None
        build_match = re.search(r'\[(GRCh3[78])\]', s, re.IGNORECASE)
        if build_match:
            build = build_match.group(1)
        copy_match = re.search(r'x([0-9])', s, re.IGNORECASE)
        if copy_match:
            copy_number = int(copy_match.group(1))
            if cnv_type_hint is None:
                # x1 = deletion, x3 = duplication (same logic as ISCN)
                pass  # handled below
        return build, copy_number
    
    # Determine cnv_type_hint from copy number if not already set by DEL/DUP prefix
    def resolve_cnv_hint(existing_hint, copy_number):
        if existing_hint:
            return existing_hint
        if copy_number is not None:
            if copy_number < 2:
                return 'deletion'
            elif copy_number > 2:
                return 'duplication'
        return None
    
    # Pattern 1: chr15:48,741,090-48,756,096 (with optional DEL/DUP prefix and x1[GRCh38] suffix)
    pattern1 = r'(?:DEL|DUP)?\s*chr([0-9XY]+):([0-9,]+)-([0-9,]+)'
    match1 = re.search(pattern1, genomic_string, re.IGNORECASE)
    
    if match1:
        build, copy_number = extract_suffix(genomic_string)
        return {
            'chromosome': match1.group(1),
            'start': int(match1.group(2).replace(',', '')),
            'end': int(match1.group(3).replace(',', '')),
            'build': build,
            'cnv_type_hint': resolve_cnv_hint(cnv_type_hint, copy_number)
        }
    
    # Pattern 2: Without 'chr' prefix: 15:48,741,090-48,756,096
    pattern1b = r'(?:DEL|DUP)?\s*([0-9XY]+):([0-9,]+)-([0-9,]+)'
    match1b = re.search(pattern1b, genomic_string, re.IGNORECASE)
    
    if match1b:
        build, copy_number = extract_suffix(genomic_string)
        return {
            'chromosome': match1b.group(1),
            'start': int(match1b.group(2).replace(',', '')),
            'end': int(match1b.group(3).replace(',', '')),
            'build': build,
            'cnv_type_hint': resolve_cnv_hint(cnv_type_hint, copy_number)
        }
    
    # Pattern 3: ISCN format arr[GRCh37]15q21.1(48,741,090_48,756,096)x1
    pattern2 = r'arr\[([^\]]+)\]([0-9XY]+)[pq][0-9.]+\(\s*([0-9,_\s]+)\)x([0-9])'
    match2 = re.search(pattern2, genomic_string, re.IGNORECASE)
    
    if match2:
        coords = match2.group(3).replace(',', '').replace(' ', '').split('_')
        copy_number = int(match2.group(4))
        if len(coords) == 2:
            return {
                'chromosome': match2.group(2),
                'start': int(coords[0]),
                'end': int(coords[1]),
                'build': match2.group(1),
                'copy_number': copy_number,
                'cnv_type_hint': resolve_cnv_hint(cnv_type_hint, copy_number)
            }
    
    return None


def process_exon_mode(gene_input, cnv_type, first_exon, last_exon, build="GRCh38"):
    """Process exon mode input - analyzes full exons only.
    
    Accepts either gene symbol (e.g., SPOUT1) or RefSeq transcript ID (e.g., NM_016390.1)
    If transcript is not in MANE database, fetches exon data from NCBI.
    """
    mane_data, _mane_gene_index, _transcript_to_gene = get_mane_data(build)
    
    gene_input = gene_input.strip().upper()
    gene = None
    transcript = None
    use_ncbi_data = False
    ncbi_exons = None
    
    # Check if input is a transcript ID (starts with NM_)
    if gene_input.startswith('NM_'):
        # It's a transcript ID
        transcript_base = gene_input.split('.')[0]  # Remove version if present
        
        # First, try to find in MANE database
        gene = _transcript_to_gene.get(gene_input) or _transcript_to_gene.get(transcript_base)
        
        if gene:
            # Found in MANE database
            transcript = list(mane_data[gene].keys())[0]
        else:
            ncbi_data = fetch_refseq_info(transcript_base)
            
            if not ncbi_data or not ncbi_data.get('gene_symbol'):
                return [{"error": f"Could not fetch transcript '{gene_input}' from NCBI. Transcript may not exist or NCBI may be unavailable."}]
            
            gene = ncbi_data['gene_symbol']
            transcript = gene_input
            ncbi_exons = ncbi_data.get('exons', [])
            use_ncbi_data = True
            
            if not ncbi_exons:
                return [{"error": f"Transcript '{gene_input}' was found (gene: {gene}) but has no exon data in NCBI."}]
    else:
        # It's a gene symbol
        gene = gene_input
        canonical = _mane_gene_index.get(gene.lower())
        if canonical:
            gene = canonical
            transcript = list(mane_data[gene].keys())[0]
        else:
            return [{"error": f"Gene '{gene_input}' not found in MANE {build} data"}]
    
    # Now analyze based on data source
    if use_ncbi_data:
        # Validate exon range is within transcript bounds
        max_exon = len(ncbi_exons)
        if first_exon < 1 or last_exon > max_exon:
            return [{"error": f"Exon range {first_exon}-{last_exon} exceeds transcript bounds. {transcript} has {max_exon} exons."}]
        
        exons_in_range = [e for e in ncbi_exons if first_exon <= e['exon'] <= last_exon]
        
        if not exons_in_range:
            return [{"error": f"Exons {first_exon}-{last_exon} not found in transcript {transcript}. This transcript has {len(ncbi_exons)} exons."}]
        
        genomic_start = min(e['start'] for e in exons_in_range)
        genomic_end = max(e['end'] for e in exons_in_range)
        
        # Create hits list - these need to match the format expected by predict_cnv_effect
        hits = []
        for exon in exons_in_range:
            hit = {
                'gene': gene,
                'transcript': transcript,
                'exon': exon['exon'],
                'start': exon['start'],
                'end': exon['end'],
                'strand': exon.get('strand') or ncbi_data.get('strand') or '+',
                'length': exon['end'] - exon['start'] + 1,
                'chromosome': exon.get('chromosome') or ncbi_data.get('chromosome') or 'unknown'
            }
            hits.append(hit)
        
        total_length = sum(h['length'] for h in hits)
        is_in_frame = (total_length % 3) == 0
        
        if cnv_type == "deletion":
            if is_in_frame:
                consequence = f"in-frame deletion ({len(hits)} exon{'s' if len(hits) > 1 else ''})"
            else:
                consequence = f"frameshift deletion ({len(hits)} exon{'s' if len(hits) > 1 else ''})"
        else:  # duplication
            if is_in_frame:
                consequence = f"in-frame duplication ({len(hits)} exon{'s' if len(hits) > 1 else ''})"
            else:
                consequence = f"frameshift duplication ({len(hits)} exon{'s' if len(hits) > 1 else ''})"
        
        result = {
            "gene": gene,
            "transcript": transcript,
            "predicted_consequence": consequence,
            "hit_exons": hits,
            "cnv_type": cnv_type,
            "data_source": "NCBI (non-MANE transcript)"
        }
        
        return [result]
        
    else:
        # Use MANE data (original logic)
        
        # Validate exon range is within transcript bounds
        exons = mane_data[gene][transcript]
        max_exon = max(e['exon'] for e in exons)
        if first_exon < 1 or last_exon > max_exon:
            return [{"error": f"Exon range {first_exon}-{last_exon} exceeds transcript bounds. {transcript} has {max_exon} exons."}]
        
        cnv_region = map_exon_numbers_to_regions(gene, transcript, first_exon, last_exon, mane_data)
        if cnv_region is None:
            return [{"error": "Invalid exon numbers"}]
        
        chromosome = exons[0].get("chromosome", "1") if exons else "1"
        
        cnv = {"gene": gene, "start": cnv_region[0], "end": cnv_region[1], "type": cnv_type, "chromosome": chromosome}
        hits = map_cnv_to_exons(cnv, mane_data, transcript)
        
        for h in hits:
            if "gene" not in h or not h.get("gene"):
                h["gene"] = gene
        
        effects = predict_cnv_effect(hits, mane_data)
        
        # Add CNV type
        for effect in effects:
            effect['cnv_type'] = cnv_type
            effect['data_source'] = 'MANE'
        
        return effects


def cdna_to_genomic(gene, transcript, c_position, mane_data):
    """
    Convert cDNA position (c.123) to genomic position.
    
    Args:
        gene: Gene symbol
        transcript: Transcript ID
        c_position: cDNA position (integer)
        mane_data: MANE exon data
    
    Returns:
        genomic_position (integer) or None
    """
    if gene not in mane_data or transcript not in mane_data[gene]:
        return None
    
    exons = mane_data[gene][transcript]
    
    # Sort exons by exon number (already sorted, but be safe)
    exons_sorted = sorted(exons, key=lambda x: x.get('exon', 0))
    
    # Walk through exons adding up lengths
    cumulative = 0
    for exon in exons_sorted:
        exon_start = exon.get('start')
        exon_end = exon.get('end')
        if exon_start is None or exon_end is None:
            continue
        
        exon_length = exon_end - exon_start + 1
        
        # Check if position is in this exon
        if cumulative + exon_length >= c_position:
            # Position is in this exon
            offset = c_position - cumulative - 1  # -1 for 0-based offset
            
            # Handle strand direction
            strand = exon.get('strand', '+')
            if strand == '-':
                # Minus strand: count from end
                genomic_pos = exon_end - offset
            else:
                # Plus strand: count from start
                genomic_pos = exon_start + offset
            
            return genomic_pos
        
        cumulative += exon_length
    
    return None


def process_transcript_mode(cdna_notation, cnv_type, build="GRCh38"):
    """
    Process HGVS cDNA notation.
    
    Args:
        cdna_notation: HGVS c. notation string
        cnv_type: deletion or duplication (may be overridden by notation)
        build: GRCh37 or GRCh38
    """
    parsed = parse_cdna_notation(cdna_notation)
    
    if not parsed:
        return [{"error": "Could not parse cDNA notation. Please use format like: FBN1:c.5065_5546del or NM_001406716.1:c.5065_5546del"}]
    
    warnings = []
    
    # Extract variant type from notation if present and check for mismatch
    if 'variant_type' in parsed:
        notation_type = parsed['variant_type']
        if notation_type != cnv_type:
            warnings.append(f"⚠️ cDNA notation contains '{notation_type}' ('{notation_type[:3]}') but CNV Type dropdown is set to '{cnv_type}'. Using '{notation_type}' from cDNA notation.")
        cnv_type = notation_type
    
    # Get correct dataset
    mane_data, _mane_gene_index, _transcript_to_gene = get_mane_data(build)
    
    # Resolve gene from input
    gene = parsed.get('gene')
    transcript_input = parsed.get('transcript')
    
    # If transcript was provided, look up which gene it belongs to
    if transcript_input and not gene:
        # Strip version number from transcript if present (e.g., NM_001406716.1 -> NM_001406716)
        transcript_base = transcript_input.split('.')[0]
        
        # Look up in the transcript-to-gene index
        gene = _transcript_to_gene.get(transcript_input) or _transcript_to_gene.get(transcript_base)
        
        # If not found in MANE database
        if not gene:
            refseq_info = fetch_refseq_info(transcript_base)
            if refseq_info:
                gene = refseq_info['gene_symbol']
                
                canonical = _mane_gene_index.get(gene.lower())
                if canonical:
                    gene_in_mane = canonical
                    mane_transcript = list(mane_data[gene_in_mane].keys())[0]
                    
                    # This is the key issue: we can't accurately map cDNA coordinates from one isoform to another
                    return [{"error": f"Transcript '{transcript_input}' is not the MANE Select transcript for {gene}. The MANE Select is '{mane_transcript}'. Since different transcript isoforms have different exon structures, cDNA coordinates (c.{parsed['start']}_{parsed['end']}) cannot be accurately transferred between isoforms. Please re-annotate your variant using the MANE Select transcript '{mane_transcript}' to ensure accurate analysis."}]
                else:
                    return [{"error": f"Transcript '{transcript_input}' is not in MANE {build} database. NCBI reports it belongs to gene '{gene}', but '{gene}' is also not in MANE {build} database. Cannot analyze this variant."}]
            else:
                return [{"error": f"Transcript '{transcript_input}' not found in MANE {build} database and could not determine gene via NCBI lookup."}]
    
    # Validate we found a gene
    if not gene:
        if transcript_input:
            return [{"error": f"Could not determine gene for transcript '{transcript_input}'. The transcript may not be in the MANE database for {build}."}]
        else:
            return [{"error": "Could not determine gene from cDNA notation. Please include gene symbol (e.g., FBN1:c.5065_5546del)"}]
    
    # Resolve canonical gene name
    canonical = _mane_gene_index.get(gene.lower())
    if canonical:
        gene = canonical
    else:
        return [{"error": f"Gene '{gene}' not found in MANE {build} data"}]
    
    # Get the MANE Select transcript for this gene (the only one in MANE data)
    mane_transcript = list(mane_data[gene].keys())[0]
    
    # Use the MANE Select transcript for analysis
    transcript = mane_transcript
    
    # Validate cDNA position is within transcript length
    exons = sorted(mane_data[gene][transcript], key=lambda x: x.get('exon', 0))
    cdna_length = sum(e['end'] - e['start'] + 1 for e in exons)
    
    if parsed['start'] > cdna_length or parsed['end'] > cdna_length:
        return [{"error": f"cDNA position c.{parsed['start']}_{parsed['end']} exceeds transcript length. {transcript} has {cdna_length}bp of coding sequence."}]
    
    # Convert cDNA coordinates to genomic coordinates
    c_start = parsed['start']
    c_end = parsed['end']
    
    genomic_start = cdna_to_genomic(gene, transcript, c_start, mane_data)
    genomic_end = cdna_to_genomic(gene, transcript, c_end, mane_data)
    
    if genomic_start is None or genomic_end is None:
        return [{"error": f"Could not convert cDNA positions c.{c_start}_c.{c_end} to genomic coordinates. Position may be outside transcript boundaries."}]
    
    # Ensure start < end (for minus strand genes, might be reversed)
    if genomic_start > genomic_end:
        genomic_start, genomic_end = genomic_end, genomic_start
    
    # Get chromosome
    chromosome = exons[0].get("chromosome", "1") if exons else "1"
    
    # Create CNV and analyze
    cnv = {"gene": gene, "start": genomic_start, "end": genomic_end, "type": cnv_type, "chromosome": chromosome}
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene
    
    effects = predict_cnv_effect(hits, mane_data)
    
    # Add CNV type, cDNA info, and warnings to results
    for effect in effects:
        effect['cnv_type'] = cnv_type
        effect['cdna_notation'] = cdna_notation
        effect['cdna_start'] = c_start
        effect['cdna_end'] = c_end
        if warnings:
            effect['warnings'] = warnings
    
    return effects


def process_genomic_mode(genomic_input, build, cnv_type, gene_hint=None):
    """Process genomic coordinate input with DEL/DUP auto-detection."""
    parsed = parse_genomic_input(genomic_input)
    
    if not parsed:
        return [{"error": "Could not parse genomic coordinates. Use format like: DEL chr15:48,741,090-48,756,096"}]
    
    chromosome = parsed['chromosome']
    start = parsed['start']
    end = parsed['end']
    warnings = []
    
    # Auto-detect from DEL/DUP prefix
    if parsed.get('cnv_type_hint'):
        detected_type = parsed['cnv_type_hint']
        if detected_type != cnv_type:
            warnings.append(f"⚠️ Input contains '{detected_type.upper()}' but CNV Type dropdown is set to '{cnv_type}'. Using '{detected_type}' from input.")
            cnv_type = detected_type
    
    # Use build from input if present, override dropdown
    original_build = build
    if parsed.get('build'):
        build = parsed['build']
        if build != original_build:
            warnings.append(f"⚠️ Input specifies {build} but the dropdown has {original_build} selected. Using {build} from the input.")
    
    mane_data, _mane_gene_index, _transcript_to_gene = get_mane_data(build)
    
    # Auto-detect from ISCN copy number
    if parsed.get('copy_number'):
        copy_num = parsed['copy_number']
        if copy_num > 2:
            auto_cnv_type = 'duplication'
        elif copy_num < 2:
            auto_cnv_type = 'deletion'
        else:
            auto_cnv_type = None
        
        if auto_cnv_type and auto_cnv_type != cnv_type:
            warnings.append(f"⚠️ ISCN copy number (x{copy_num}) indicates {auto_cnv_type} but CNV Type is set to {cnv_type}. Using {auto_cnv_type} from ISCN.")
            cnv_type = auto_cnv_type
    
    # Find gene
    if gene_hint:
        gene = gene_hint.upper().strip()
        canonical = _mane_gene_index.get(gene.lower())
        if canonical:
            gene = canonical
    else:
        genes_in_region = find_genes_in_region(chromosome, start, end, mane_data)
        if not genes_in_region:
            return [{"error": f"No genes found in region chr{chromosome}:{start:,}-{end:,} using {build}. Please provide a gene symbol as a hint."}]
        # Check for multiple genes in region
        if len(genes_in_region) > 1:
            gene_list = ", ".join(genes_in_region)
            return [{"error": f"Multiple genes found in region chr{chromosome}:{start:,}-{end:,}: {gene_list}. This version only supports single-gene analysis. Please provide a specific gene symbol as a hint."}]
        gene = genes_in_region[0]
    
    if gene not in mane_data:
        return [{"error": f"Gene '{gene}' not found in MANE {build} data"}]
    
    transcript = list(mane_data[gene].keys())[0]
    exons = mane_data[gene][transcript]
    
    # Check if CNV affects transcription start or stop sites
    gene_start = min(e['start'] for e in exons)
    gene_end = max(e['end'] for e in exons)
    
    affects_start = start <= gene_start
    affects_end = end >= gene_end
    
    if affects_start and affects_end:
        return [{"error": f"CNV encompasses entire gene {gene} (transcription start site to stop codon). This would result in complete gene deletion/duplication. Region: chr{chromosome}:{start:,}-{end:,}, Gene span: {gene_start:,}-{gene_end:,}."}]
    if affects_start:
        warnings.append(f"⚠️ CNV affects transcription start site of {gene}. This may prevent transcription initiation.")
    if affects_end:
        warnings.append(f"⚠️ CNV affects stop codon of {gene}. This may result in a truncated or extended protein product.")
    
    cnv = {"gene": gene, "start": start, "end": end, "type": cnv_type, "chromosome": f"chr{chromosome}"}
    hits = map_cnv_to_exons(cnv, mane_data, transcript)
    
    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene
    
    effects = predict_cnv_effect(hits, mane_data)
    
    for effect in effects:
        effect['cnv_type'] = cnv_type
        if warnings:
            effect['warnings'] = warnings
    
    return effects


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
                
                result = process_exon_mode(
                    gene_input, cnv_type, first_exon, last_exon, build
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
    
    return render_template("index.html", result=result, error=error)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=False)
