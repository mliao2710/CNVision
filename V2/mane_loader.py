"""
================================================================================
MANE_LOADER.PY - Load Gene/Exon Data from MANE Database
================================================================================

PURPOSE:
--------
This module loads and parses gene/exon annotation data from the MANE database.
MANE (Matched Annotation from NCBI and EMBL-EBI) provides high-quality, manually
curated gene annotations that are the gold standard for clinical genetics.

BIOLOGICAL IMPACT:
------------------
Accurate gene/exon coordinates are ESSENTIAL for CNV analysis. Errors here would
lead to incorrect predictions about disease severity. MANE data ensures we're
using the most reliable, clinically validated gene annotations available.

The data structure created here is used throughout the entire application:
- Web app uses it to find genes
- Coordinate mapper uses it to find exons
- Functional predictor uses it to calculate exon lengths

CODING CONCEPTS & PRINCIPLES:
------------------------------
1. FILE PARSING:
   - GFF (General Feature Format) is a standard genomics file format
   - Each line represents a genomic feature (gene, exon, etc.)
   - Fields are tab-separated: chromosome, start, end, feature_type, attributes...

2. GZIP COMPRESSION:
   - MANE files are compressed (.gz) to save space
   - gzip.open() automatically decompresses while reading
   - "rt" mode = read text (not binary)

3. STRING PARSING:
   - Split strings by delimiters (tabs, semicolons, equals signs)
   - Extract key-value pairs from attribute strings
   - Handle missing or malformed data gracefully

4. DATA STRUCTURE BUILDING:
   - Build nested dictionaries incrementally
   - setdefault() creates nested structure if it doesn't exist
   - Sort data for consistent ordering

5. EXTERNAL API CALLS:
   - fetch_refseq_exons() calls NCBI API to get transcript data
   - Uses requests library for HTTP calls
   - XML parsing to extract structured data

6. ERROR HANDLING:
   - Try/except blocks prevent crashes on bad data
   - Returns None on failure (caller can check)

================================================================================
"""

import gzip
import requests
from xml.etree import ElementTree as ET

def load_mane_exons(gff_path):
    """
    Load and parse MANE exon data from a GFF (General Feature Format) file.
    
    GFF Format (tab-separated):
        chromosome  source  feature_type  start  end  score  strand  phase  attributes
    
    Example line:
        chr17  MANE  exon  7565097  7565200  .  +  .  gene=TP53;transcript_id=NM_000546.6;exon_number=1
    
    This function:
    1. Reads compressed GFF file line by line
    2. Filters for "exon" features only
    3. Extracts gene name, transcript ID, exon number, coordinates
    4. Builds nested dictionary structure: {gene: {transcript: [exons]}}
    5. Calculates exon lengths (needed for frameshift prediction)
    
    Returns:
        Nested dictionary structure:
        {
            "BRCA1": {
                "NM_007294.4": [
                    {"exon": 1, "start": 43000000, "end": 43000150, "strand": "+", "length": 151},
                    {"exon": 2, "start": 43010000, "end": 43010100, "strand": "+", "length": 101},
                    ...
                ]
            },
            ...
        }
    """
    # Initialize empty dictionary to store parsed data
    exon_data = {}

    # ========================================================================
    # STEP 1: Open and read GFF file
    # ========================================================================
    # gzip.open() automatically decompresses .gz files
    # "rt" mode = read as text (not binary)
    # "with" statement ensures file is closed automatically
    with gzip.open(gff_path, "rt") as f:
        for line in f:
            # Skip comment lines (start with #)
            if line.startswith("#"):
                continue

            # ================================================================
            # STEP 2: Parse tab-separated fields
            # ================================================================
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # Skip malformed lines

            # GFF format: seqid, source, feature_type, start, end, score, strand, phase, attributes
            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields

            # ================================================================
            # STEP 3: Filter for exon features only
            # ================================================================
            # GFF files contain many feature types (gene, transcript, exon, CDS, etc.)
            # We only care about exons
            if feature_type != "exon":
                continue

            # ================================================================
            # STEP 4: Parse attributes field
            # ================================================================
            # Attributes are semicolon-separated key=value pairs
            # Example: "gene=TP53;transcript_id=NM_000546.6;exon_number=1"
            attr_dict = {}
            for item in attrs.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attr_dict[key] = value

            # ================================================================
            # STEP 5: Extract gene name (try multiple possible attribute names)
            # ================================================================
            # Different GFF versions use different attribute names
            # Try all common variations to be robust
            gene = (
                attr_dict.get("gene")
                or attr_dict.get("gene_name")
                or attr_dict.get("gene_id")
                or attr_dict.get("Name")
            )

            # ================================================================
            # STEP 6: Extract transcript ID
            # ================================================================
            transcript = attr_dict.get("transcript_id")
            if not transcript:
                # Sometimes transcript is in "Parent" attribute
                parent = attr_dict.get("Parent") or attr_dict.get("parent")
                if parent:
                    transcript = parent.split(",")[0]  # Take first if multiple
                    # Remove common prefixes if present
                    for prefix in ("rna-", "transcript-", "rna:", "transcript:"):
                        if transcript.startswith(prefix):
                            transcript = transcript[len(prefix):]
                            break

            # ================================================================
            # STEP 7: Extract exon number
            # ================================================================
            exon_number = attr_dict.get("exon_number")
            if exon_number is None:
                # Try to extract from ID attribute (e.g., "exon-1" → 1)
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

            # ================================================================
            # STEP 8: Validate required data
            # ================================================================
            if not gene or not transcript:
                continue  # Skip if missing critical information

            # ================================================================
            # STEP 9: Convert coordinates and calculate length
            # ================================================================
            start = int(start)
            end = int(end)
            # Length = end - start + 1 (both positions are inclusive)
            # Example: start=100, end=150 → length = 150-100+1 = 51 bp
            length = end - start + 1

            # ================================================================
            # STEP 10: Store in nested dictionary structure
            # ================================================================
            # setdefault() creates nested structure if it doesn't exist
            # This avoids KeyError if gene/transcript seen for first time
            exon_data.setdefault(gene, {})
            exon_data[gene].setdefault(transcript, [])
            exon_data[gene][transcript].append({
                "exon": exon_number,
                "start": start,
                "end": end,
                "strand": strand,
                "length": length,
                "chromosome": seqid,  # <--- ADD THIS LINE to save the chromosome (e.g., 'chr15')
            })
    # ========================================================================
    # STEP 11: Sort exons by exon number (or start position if number missing)
    # ========================================================================
    # Ensures exons are in correct order for downstream processing
    for gene in exon_data:
        for transcript in exon_data[gene]:
            exon_data[gene][transcript] = sorted(
                exon_data[gene][transcript],
                key=lambda x: (x["exon"] if x["exon"] is not None else x["start"])
            )

    return exon_data


def fetch_refseq_exons(nm_id):
    """
    Fetch exon coordinates for a RefSeq transcript from NCBI database.
    
    This function is used when a user provides an NM_XXXXX transcript ID that
    isn't in the MANE database. It queries NCBI's public API to get the data.
    
    NCBI API:
    - Uses Entrez eUtils (public API for accessing NCBI databases)
    - Returns data in XML format
    - Requires internet connection
    
    Args:
        nm_id: RefSeq transcript ID (e.g., "NM_000546.6")
    
    Returns:
        List of exon dictionaries with coordinates, or None if fetch fails
    """
    # ========================================================================
    # STEP 1: Build API request URL
    # ========================================================================
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",      # Database: nucleotide sequences
        "id": nm_id,              # Transcript ID to fetch
        "rettype": "gb",          # Return type: GenBank format
        "retmode": "xml"          # Return mode: XML (structured data)
    }
    
    # ========================================================================
    # STEP 2: Make HTTP request to NCBI
    # ========================================================================
    try:
        # requests.get() sends HTTP GET request
        # timeout=10 means give up after 10 seconds (prevents hanging)
        r = requests.get(url, params=params, timeout=10)
        r.raise_for_status()  # Raise exception if HTTP error (404, 500, etc.)
    except Exception as e:
        print(f"Error fetching {nm_id}: {e}")
        return None

    # ========================================================================
    # STEP 3: Parse XML response
    # ========================================================================
    exons = []
    try:
        # Parse XML string into tree structure
        tree = ET.fromstring(r.content)
        
        # Find all "exon" features in the XML
        # ".//GBFeature" means "find GBFeature anywhere in tree"
        for feature in tree.findall(".//GBFeature"):
            if feature.findtext("GBFeature_key") == "exon":
                # Extract location string (e.g., "123..456" or "<123..>456")
                loc = feature.findtext("GBFeature_location")
                if not loc:
                    continue
                
                # ============================================================
                # STEP 4: Parse location string
                # ============================================================
                # Handle different formats:
                #   "123..456" → start=123, end=456
                #   "<123..>456" → start=123, end=456 (remove < >)
                start_str, end_str = loc.replace("<", "").replace(">", "").split("..")
                start = int(start_str)
                end = int(end_str)
                length = end - start + 1
                
                # Store exon data
                # exon number = position in list + 1 (1-indexed)
                exons.append({
                    "exon": len(exons) + 1,
                    "start": start,
                    "end": end,
                    "strand": None,  # Not always available in NCBI data
                    "length": length
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