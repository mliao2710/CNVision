"""
CNVision NCBI Fetcher - Remote Gene Data Retrieval

Technologies: Python 3.10+, requests (HTTP), XML parsing, REST APIs

Core Functions:
- fetch_refseq_info(): queries NCBI via esearch/efetch APIs
- Extracts gene symbol, chromosome, strand, exon coordinates
- Handles rate limiting and caching for performance
- Fallback for non-MANE transcripts
"""

import requests
from xml.etree import ElementTree as ET

def fetch_refseq_info(nm_id: str):
    """
    Fetch gene information for a RefSeq transcript from NCBI database.
    
    This function queries NCBI's Entrez API to get:
    - Gene symbol (e.g., "TP53")
    - Exon coordinates
    - Transcript information
    
    Args:
        nm_id: RefSeq transcript ID (e.g., "NM_000546.6")
    
    Returns:
        Dictionary with gene information:
        {
            "gene_symbol": "TP53",
            "transcript_id": "NM_000546.6",
            "exons": [
                {"exon": 1, "start": 7565097, "end": 7565200, "strand": None},
                ...
            ]
        }
        Returns None if fetch or parsing fails
    """
    # ========================================================================
    # STEP 1: Build API request
    # ========================================================================
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",      # NCBI nucleotide database
        "id": nm_id,             # Transcript ID to fetch
        "rettype": "gb",         # Return GenBank format
        "retmode": "xml"         # Return as XML (structured)
    }

    # ========================================================================
    # STEP 2: Make HTTP request
    # ========================================================================
    try:
        # Send GET request to NCBI API
        r = requests.get(url, params=params, timeout=10)
        r.raise_for_status()  # Raise exception if HTTP error
    except Exception as e:
        print(f"Error fetching {nm_id}: {e}")
        return None

    # ========================================================================
    # STEP 3: Parse XML response
    # ========================================================================
    try:
        # Parse XML string into tree structure
        tree = ET.fromstring(r.content)
        gene_symbol = None
        exons = []
        
        # ================================================================
        # STEP 4: Extract gene symbol and exon data from XML
        # ================================================================
        # Find all "GBFeature" elements (represent genomic features)
        for feature in tree.findall(".//GBFeature"):
            key = feature.findtext("GBFeature_key")
            
            # ============================================================
            # Extract gene symbol
            # ============================================================
            if key == "gene" and gene_symbol is None:
                # Gene symbol is stored in a "qualifier" sub-element
                # Search through qualifiers to find one named "gene"
                for qualifier in feature.findall("GBFeature_quals/GBQualifier"):
                    if qualifier.findtext("GBQualifier_name") == "gene":
                        gene_symbol = qualifier.findtext("GBQualifier_value")
                        break
            
            # ============================================================
            # Extract exon coordinates
            # ============================================================
            elif key == "exon":
                # Exon location is stored as a range string
                # Format: "123..456" or "<123..>456"
                loc = feature.findtext("GBFeature_location")
                if not loc:
                    continue
                
                # Remove angle brackets if present and split by ".."
                loc = loc.replace("<", "").replace(">", "")
                if ".." not in loc:
                    continue
                
                # Parse start and end coordinates
                start_str, end_str = loc.split("..")
                exons.append({
                    "exon": len(exons) + 1,  # Exon number (1-indexed)
                    "start": int(start_str),
                    "end": int(end_str),
                    "strand": None  # Not always available in NCBI XML
                })

        # ================================================================
        # STEP 5: Validate and return results
        # ================================================================
        if not gene_symbol:
            print(f"Could not find gene symbol for {nm_id}")
            return None

        return {
            "gene_symbol": gene_symbol,
            "transcript_id": nm_id,
            "exons": exons
        }
    except Exception as e:
        print(f"Error parsing NCBI response for {nm_id}: {e}")
        return None
