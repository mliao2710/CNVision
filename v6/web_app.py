"""
CNVision Web Application - Flask UI and Analysis Orchestration

Purpose:
- Main runtime entry point for browser-based CNV analysis.
- Accepts three input modes: exon numbers, genomic coordinates, and HGVS genomic (g.) notation.
- Coordinates parsing, mapping, consequence prediction, and result rendering.

Core Functions:
- Parse and normalize user inputs (build aliases, CNV type hints, coordinate formats).
- Resolve genes/transcripts from MANE indices and NCBI fallback when needed.
- Run exon mapping + frame consequence prediction pipeline.
- Return structured outputs to the Jinja template for display.
"""

from flask import Flask, render_template, request
from coordinate_mapper import map_cnv_to_exons, map_exon_numbers_to_regions
from functional_predictor import predict_cnv_effect
from mane_loader import load_mane_exons
from iscn_parser import find_genes_in_region
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
        # STEP 1: Resolve transcript accession to Gene ID.
        r = requests.get(f"{NCBI_BASE}/esearch.fcgi", params={
            "db": "gene",
            "term": f"{base_id}[Gene RefSeq]",
            "retmode": "xml"
        }, timeout=10, verify=False)
        r.raise_for_status()
        
        root = ET.fromstring(r.content)
        id_list = root.findall(".//IdList/Id")
        
        # Fallback query when primary term has no hits.
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
        
        # STEP 2: Parse `gene_table` for strand/chromosome/exon rows.
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
        past_header = False
        
        for line in lines:
            stripped = line.strip()
            
            # Gene symbol is the first token on the first non-empty data line.
            if not gene_symbol and stripped and not stripped.startswith("Gene ID"):
                gene_symbol = stripped.split()[0]
            
            # Parse reference line for strand and chromosome.
            if "Primary Assembly" in stripped:
                if "minus strand" in stripped:
                    strand = "-"
                elif "plus strand" in stripped:
                    strand = "+"
                if "NC_" in stripped:
                    nc_part = stripped.split("NC_")[1]
                    nc_num = nc_part.split(".")[0]
                    chrom_num = str(int(nc_num))
                    chromosome = f"chr{chrom_num}"
            
            # Start parsing once we hit the transcript-specific exon table.
            if f"Exon table for" in stripped and base_id in stripped:
                in_our_exon_table = True
                past_header = False
                exons = []
                continue
            
            if in_our_exon_table:
                if "---" in stripped:
                    past_header = True
                    continue
                if not past_header:
                    continue
                
                # Empty line ends this exon section.
                if stripped == "":
                    if exons:
                        break
                    continue
                
                # Parse genomic interval in first column (e.g. 7687490-7687377).
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

# STEP: Load both reference builds at startup.
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
    normalized_build = normalize_build(build) or "GRCh38"
    if normalized_build == "GRCh37":
        return mane_data_grch37, _mane_gene_index_grch37, _transcript_to_gene_grch37
    else:
        return mane_data_grch38, _mane_gene_index_grch38, _transcript_to_gene_grch38


def normalize_build(build):
    """Normalize build aliases (e.g., hg19/hg38) to GRCh labels."""
    if not build:
        return None
    value = str(build).strip()
    if value.startswith("[") and value.endswith("]"):
        value = value[1:-1]
    value = value.replace("_", "").replace(" ", "")
    alias = {
        "grch37": "GRCh37",
        "hg19": "GRCh37",
        "grch38": "GRCh38",
        "hg38": "GRCh38",
    }
    return alias.get(value.lower())


def _build_from_nc_version(chromosome, version):
    """Infer genome build from NC accession chromosome + version."""
    grch37_versions = {
        "1": 10, "2": 11, "3": 11, "4": 11, "5": 9, "6": 11, "7": 13, "8": 10,
        "9": 11, "10": 10, "11": 9, "12": 11, "13": 10, "14": 8, "15": 9, "16": 9,
        "17": 10, "18": 9, "19": 9, "20": 10, "21": 8, "22": 10, "X": 10, "Y": 9
    }
    grch38_versions = {
        "1": 11, "2": 12, "3": 12, "4": 12, "5": 10, "6": 12, "7": 14, "8": 11,
        "9": 12, "10": 11, "11": 10, "12": 12, "13": 11, "14": 9, "15": 10, "16": 10,
        "17": 11, "18": 10, "19": 10, "20": 11, "21": 9, "22": 11, "X": 11, "Y": 10
    }

    if chromosome in grch37_versions and grch37_versions[chromosome] == version:
        return "GRCh37"
    if chromosome in grch38_versions and grch38_versions[chromosome] == version:
        return "GRCh38"
    return None


def parse_hgvs_genomic_notation(hgvs_string):
    """Parse HGVS genomic notation (g.) into chromosome and interval."""
    text = hgvs_string.strip()
    if not text:
        return None

    # STEP 1: Parse NC accession and map chromosome.
    acc_match = re.search(r"NC_(\d{6})\.(\d+)", text, re.IGNORECASE)
    if not acc_match:
        return None

    chrom_num = int(acc_match.group(1))
    version = int(acc_match.group(2))
    if 1 <= chrom_num <= 22:
        chromosome = str(chrom_num)
    elif chrom_num == 23:
        chromosome = "X"
    elif chrom_num == 24:
        chromosome = "Y"
    else:
        return None

    # STEP 2: Require HGVS genomic scope.
    if ":g." not in text.lower():
        return None

    # STEP 3: Detect CNV type hint from HGVS suffix.
    cnv_type_hint = None
    suffix_match = re.search(r"(del|dup)\s*$", text, re.IGNORECASE)
    if suffix_match:
        token = suffix_match.group(1).lower()
        cnv_type_hint = "deletion" if token == "del" else "duplication"

    # STEP 4: Parse HGVS breakpoint form (a_b)_(c_d) and use b..c.
    pair_match = re.search(
        r"\(\s*(\d[\d,]*)\s*[_-]\s*(\d[\d,]*)\s*\)\s*[_-]\s*\(\s*(\d[\d,]*)\s*[_-]\s*(\d[\d,]*)\s*\)",
        text,
        re.IGNORECASE
    )
    if pair_match:
        first = int(pair_match.group(2).replace(",", ""))
        second = int(pair_match.group(3).replace(",", ""))
    else:
        # Fallback: use first two large integers after g.
        tail = text.split(":g.", 1)[1]
        coords = [int(m.group(0).replace(",", "")) for m in re.finditer(r"\d[\d,]{4,}", tail)]
        if len(coords) < 2:
            return None
        first, second = coords[0], coords[1]

    start, end = (first, second) if first <= second else (second, first)
    inferred_build = _build_from_nc_version(chromosome, version)

    return {
        "chromosome": chromosome,
        "start": start,
        "end": end,
        "build": inferred_build,
        "cnv_type_hint": cnv_type_hint,
        "nc_accession": f"NC_{acc_match.group(1)}.{version}",
        "nc_version": version
    }


def parse_genomic_input(genomic_string):
    """Parse genomic coordinate notation with robust format tolerance."""
    raw = genomic_string.strip()
    if not raw:
        return None

    text = raw
    text_upper = text.upper()

    # STEP 1: Read CNV type hint.
    cnv_type_hint = None
    if re.search(r'^\s*DEL\b', text_upper):
        cnv_type_hint = 'deletion'
    elif re.search(r'^\s*DUP\b', text_upper):
        cnv_type_hint = 'duplication'

    # STEP 2: Detect build aliases.
    build = None
    build_match = re.search(r'\[?\s*(grch\s*3[78]|hg\s*19|hg\s*38)\s*\]?', text, re.IGNORECASE)
    if build_match:
        build = normalize_build(build_match.group(1))

    # STEP 3: Detect copy number from xN / ×N.
    copy_number = None
    copy_match = re.search(r'[x×]\s*([0-9]+)', text, re.IGNORECASE)
    if copy_match:
        copy_number = int(copy_match.group(1))

    def resolve_cnv_hint(existing_hint, copy_num):
        if existing_hint:
            return existing_hint
        if copy_num is None:
            return None
        if copy_num < 2:
            return 'deletion'
        if copy_num > 2:
            return 'duplication'
        return None

    # STEP 4: Detect chromosome.
    # 4a) explicit chr prefix.
    chr_match = re.search(r'chr\s*[:\-]?\s*([0-9]{1,2}|X|Y)\b', text, re.IGNORECASE)
    chromosome = chr_match.group(1) if chr_match else None

    # 4b) ISCN-like lead token (ignore p/q cytoband suffixes).
    if not chromosome:
        lead = text
        lead = re.sub(r'^\s*(DEL|DUP)\b', '', lead, flags=re.IGNORECASE)
        lead = re.sub(r'^\s*(arr|seq|sseq)\b', '', lead, flags=re.IGNORECASE)
        lead = re.sub(r'\[\s*(grch\s*3[78]|hg\s*19|hg\s*38)\s*\]', '', lead, flags=re.IGNORECASE)
        lead = lead.strip()
        lead_match = re.match(r'^([0-9]{1,2}|X|Y)(?=[pq]|[:\(\s_\-]|$)', lead, re.IGNORECASE)
        if lead_match:
            chromosome = lead_match.group(1)

    # STEP 5: Use first two large integers as genomic coordinates.
    coordinate_values = []
    for m in re.finditer(r'\d[\d,]*', text):
        token = m.group(0).replace(',', '')
        if len(token) >= 5:
            try:
                coordinate_values.append(int(token))
            except ValueError:
                continue

    if chromosome and len(coordinate_values) >= 2:
        start = coordinate_values[0]
        end = coordinate_values[1]
        if start > end:
            start, end = end, start
        return {
            'chromosome': chromosome.upper(),
            'start': start,
            'end': end,
            'build': build,
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
    
    # STEP 1: Resolve user input to gene/transcript.
    if gene_input.startswith('NM_'):
        transcript_base = gene_input.split('.')[0]
        
        gene = _transcript_to_gene.get(gene_input) or _transcript_to_gene.get(transcript_base)
        
        if gene:
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
        gene = gene_input
        canonical = _mane_gene_index.get(gene.lower())
        if canonical:
            gene = canonical
            transcript = list(mane_data[gene].keys())[0]
        else:
            return [{"error": f"Gene '{gene_input}' not found in MANE {build} data"}]
    
    # STEP 2: Run exon-range analysis with the resolved data source.
    if use_ncbi_data:
        max_exon = len(ncbi_exons)
        if first_exon < 1 or last_exon > max_exon:
            return [{"error": f"Exon range {first_exon}-{last_exon} exceeds transcript bounds. {transcript} has {max_exon} exons."}]
        
        exons_in_range = [e for e in ncbi_exons if first_exon <= e['exon'] <= last_exon]
        
        if not exons_in_range:
            return [{"error": f"Exons {first_exon}-{last_exon} not found in transcript {transcript}. This transcript has {len(ncbi_exons)} exons."}]
        
        genomic_start = min(e['start'] for e in exons_in_range)
        genomic_end = max(e['end'] for e in exons_in_range)
        
        # Build normalized exon-hit payload.
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
        else:
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
        # MANE-backed path.
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
        
        for effect in effects:
            effect['cnv_type'] = cnv_type
            effect['data_source'] = 'MANE'
        
        return effects


def _analyze_genomic_interval(chromosome, start, end, build, cnv_type, gene_hint=None, warnings_list=None):
    """Shared genomic interval analysis for coordinate and HGVS g. modes."""
    warnings = list(warnings_list or [])
    mane_data, _mane_gene_index, _transcript_to_gene = get_mane_data(build)

    # STEP 1: Resolve target gene (hint or overlap search).
    if gene_hint:
        gene = gene_hint.upper().strip()
        canonical = _mane_gene_index.get(gene.lower())
        if canonical:
            gene = canonical
    else:
        genes_in_region = find_genes_in_region(chromosome, start, end, mane_data)
        if not genes_in_region:
            return [{"error": f"No genes found in region chr{chromosome}:{start:,}-{end:,} using {build}. Please provide a gene symbol as a hint."}]
        if len(genes_in_region) > 1:
            gene_list = ", ".join(genes_in_region)
            return [{"error": f"Multiple genes found in region chr{chromosome}:{start:,}-{end:,}: {gene_list}. This version only supports single-gene analysis. Please provide a specific gene symbol as a hint."}]
        gene = genes_in_region[0]

    if gene not in mane_data:
        return [{"error": f"Gene '{gene}' not found in MANE {build} data"}]

    transcript = list(mane_data[gene].keys())[0]
    exons = mane_data[gene][transcript]

    # STEP 2: Add boundary warnings for start/stop involvement.
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

    # STEP 3: Return intronic result when no exon overlap is found.
    if not hits or not hits[0].get("hit_exons"):
        result = {
            "gene": gene,
            "transcript": transcript,
            "hit_exons": [],
            "predicted_consequence": "intronic deletion (no exon overlap in MANE Select transcript)"
        }
        result["cnv_type"] = cnv_type
        if warnings:
            result["warnings"] = warnings
        return [result]

    for h in hits:
        if "gene" not in h or not h.get("gene"):
            h["gene"] = gene

    # STEP 4: Predict frame consequence from mapped hits.
    effects = predict_cnv_effect(hits, mane_data)
    for effect in effects:
        effect["cnv_type"] = cnv_type
        if warnings:
            effect["warnings"] = warnings
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
    
    # STEP 1: Apply parser hints (type/build/copy number).
    if parsed.get('cnv_type_hint'):
        detected_type = parsed['cnv_type_hint']
        if detected_type != cnv_type:
            warnings.append(f"⚠️ Input contains '{detected_type.upper()}' but CNV Type dropdown is set to '{cnv_type}'. Using '{detected_type}' from input.")
            cnv_type = detected_type
    
    original_build = build
    if parsed.get('build'):
        build = parsed['build']
        if build != original_build:
            warnings.append(f"⚠️ Input specifies {build} but the dropdown has {original_build} selected. Using {build} from the input.")
    
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

    return _analyze_genomic_interval(chromosome, start, end, build, cnv_type, gene_hint, warnings)


def process_hgvs_mode(hgvs_notation, build, cnv_type, gene_hint=None):
    """Process HGVS genomic notation (g.) as genomic coordinate input."""
    parsed = parse_hgvs_genomic_notation(hgvs_notation)
    if not parsed:
        return [{"error": "Could not parse HGVS g. notation. Use format like: NC_000019.9:g.(13325430_13335371)_(13346648_13352244)del"}]

    warnings = []

    # STEP 1: Apply parser CNV type hint.
    if parsed.get("cnv_type_hint") and parsed["cnv_type_hint"] != cnv_type:
        detected_type = parsed["cnv_type_hint"]
        warnings.append(f"⚠️ HGVS notation indicates '{detected_type}' but CNV Type dropdown is set to '{cnv_type}'. Using '{detected_type}' from HGVS.")
        cnv_type = detected_type

    # STEP 2: Apply build inferred from NC accession version.
    original_build = build
    if parsed.get("build"):
        build = parsed["build"]
        if build != original_build:
            warnings.append(f"⚠️ HGVS accession version indicates {build} but the dropdown has {original_build} selected. Using {build} from HGVS accession.")
    else:
        warnings.append("⚠️ Could not infer build from NC accession version. Using the selected build from the dropdown.")

    effects = _analyze_genomic_interval(
        parsed["chromosome"],
        parsed["start"],
        parsed["end"],
        build,
        cnv_type,
        gene_hint,
        warnings
    )

    for effect in effects:
        if not effect.get("error"):
            effect["hgvs_notation"] = hgvs_notation
            effect["nc_accession"] = parsed.get("nc_accession")
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
                
            elif mode == "coordinate":
                genomic_input = request.form.get("genomic_input", "").strip()
                gene_hint = request.form.get("gene_hint", "").strip()
                
                result = process_genomic_mode(
                    genomic_input, build, cnv_type,
                    gene_hint if gene_hint else None
                )
            elif mode == "hgvs":
                hgvs_notation = request.form.get("hgvs_notation", "").strip()
                build_hgvs = request.form.get("build_hgvs", "").strip()
                gene_hint_hgvs = request.form.get("gene_hint_hgvs", "").strip()
                selected_build = build_hgvs if build_hgvs else build
                result = process_hgvs_mode(
                    hgvs_notation, selected_build, cnv_type,
                    gene_hint_hgvs if gene_hint_hgvs else None
                )

        except ValueError as e:
            error = f"Invalid input: {str(e)}"
        except Exception as e:
            error = f"Error: {str(e)}"
    
    return render_template("index.html", result=result, error=error)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=False)
