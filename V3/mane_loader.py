"""
================================================================================
MANE_LOADER.PY - Load Gene/Exon Data from MANE or RefSeq Select Database
================================================================================
Updated to support both MANE GFF and UCSC RefSeq Select GTF formats
Uses pre-loaded mapping file for fast gene symbol lookup
================================================================================
"""

import gzip
import requests
from xml.etree import ElementTree as ET
import re
import os

# Global cache for transcript -> gene mapping
TRANSCRIPT_TO_GENE_CACHE = {}

def load_gene_mapping(mapping_file="data/gene2refseq"):
    """
    Load transcript -> gene mapping from NCBI gene2refseq file.
    This avoids having to fetch each gene individually from NCBI API.
    
    Format of gene2refseq:
    #tax_id GeneID  status  RNA_nucleotide_accession.version  ...  Symbol  ...
    9606    1       REVIEWED        NM_130786.4     ...     A1BG    ...
    """
    mapping = {}
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(script_dir, mapping_file)
    
    if not os.path.exists(full_path):
        print(f"Warning: Gene mapping file not found at {full_path}")
        print("Will use transcript IDs as gene names. Download with:")
        print("  cd data && curl -O https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz && gunzip gene2refseq.gz")
        return mapping
    
    print(f"Loading gene mappings from {mapping_file}...")
    
    with open(full_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 16:
                continue
            
            tax_id = fields[0]
            # Only process human genes (tax_id 9606)
            if tax_id != '9606':
                continue
            
            rna_acc = fields[3]  # RNA_nucleotide_accession.version (NM_016390.4)
            gene_symbol = fields[15]  # Symbol is in column 16 (0-indexed = 15) - SPOUT1
            
            if rna_acc and rna_acc != '-' and gene_symbol and gene_symbol != '-':
                # Store both with and without version
                mapping[rna_acc] = gene_symbol
                base_id = rna_acc.split('.')[0]
                mapping[base_id] = gene_symbol
    
    print(f"Loaded {len(mapping)} transcript → gene mappings")
    return mapping


def fetch_gene_from_transcript(transcript_id):
    """
    Get gene symbol from transcript ID using pre-loaded mapping.
    Falls back to transcript ID if not found.
    """
    # Check cache first
    if transcript_id in TRANSCRIPT_TO_GENE_CACHE:
        return TRANSCRIPT_TO_GENE_CACHE[transcript_id]
    
    # Try without version
    base_id = transcript_id.split('.')[0]
    if base_id in TRANSCRIPT_TO_GENE_CACHE:
        return TRANSCRIPT_TO_GENE_CACHE[base_id]
    
    # Fallback: use transcript ID as gene name
    return base_id


def load_mane_exons(gff_path):
    """
    Load and parse exon data from MANE GFF or RefSeq Select GTF file.
    
    Automatically detects file format and parses accordingly.
    """
    global TRANSCRIPT_TO_GENE_CACHE
    
    # Load gene mapping file once if not already loaded
    if not TRANSCRIPT_TO_GENE_CACHE:
        TRANSCRIPT_TO_GENE_CACHE = load_gene_mapping()
    
    exon_data = {}
    transcript_to_gene = {}  # Local cache for this file
    
    with gzip.open(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields

            # Only process exon features
            if feature_type != "exon":
                continue

            # Parse attributes
            attr_dict = {}
            
            # Check if it's GTF format (UCSC RefSeq Select) or GFF format (MANE)
            if 'gene_id "' in attrs:
                # GTF format (UCSC RefSeq Select)
                # Example: gene_id "NM_032291.4"; transcript_id "NM_032291.4";
                for match in re.finditer(r'(\w+)\s+"([^"]+)"', attrs):
                    key = match.group(1)
                    value = match.group(2)
                    attr_dict[key] = value
            else:
                # GFF format (MANE)
                # Example: gene=TP53;transcript_id=NM_000546.6;exon_number=1
                for item in attrs.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attr_dict[key] = value

            # Extract transcript ID
            transcript = attr_dict.get("transcript_id") or attr_dict.get("gene_id")
            if not transcript:
                continue

            # Extract or fetch gene symbol
            gene = (
                attr_dict.get("gene")
                or attr_dict.get("gene_name")
                or attr_dict.get("gene_id")
                or attr_dict.get("Name")
            )
            
            # If no gene symbol in file, try to fetch it from NCBI
            if not gene or gene.startswith("NM_"):
                # Check cache first
                if transcript in transcript_to_gene:
                    gene = transcript_to_gene[transcript]
                else:
                    # For RefSeq Select files, we need to fetch gene symbols
                    # Only fetch for unique transcripts to avoid hammering NCBI
                    fetched_gene = fetch_gene_from_transcript(transcript)
                    if fetched_gene:
                        gene = fetched_gene
                        transcript_to_gene[transcript] = gene
                    else:
                        # Skip this exon if we can't get gene name
                        continue

            # Extract exon number
            exon_number = attr_dict.get("exon_number")
            if exon_number is None:
                # Try to extract from ID or count exons
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

            if not gene or not transcript:
                continue

            # Convert coordinates
            start = int(start)
            end = int(end)
            length = end - start + 1

            # Store in nested dictionary
            exon_data.setdefault(gene, {})
            exon_data[gene].setdefault(transcript, [])
            exon_data[gene][transcript].append({
                "exon": exon_number,
                "start": start,
                "end": end,
                "strand": strand,
                "length": length,
                "chromosome": seqid,
            })

    # Sort exons and assign exon numbers based on strand
    for gene in exon_data:
        for transcript in exon_data[gene]:
            exons = exon_data[gene][transcript]
            
            # Determine strand (all exons should have same strand)
            strand = exons[0]["strand"] if exons else "+"
            
            # Sort by start position (genomic order)
            exons_sorted = sorted(exons, key=lambda x: x["start"])
            
            # Assign exon numbers based on strand direction
            if strand == "-":
                # Minus strand: highest position = exon 1
                # Reverse the list so highest position gets lowest exon number
                exons_sorted.reverse()
            
            # Assign exon numbers
            for idx, exon in enumerate(exons_sorted, start=1):
                if exon["exon"] is None:
                    exon["exon"] = idx
            
            # Sort back by genomic position for storage
            exon_data[gene][transcript] = sorted(exons_sorted, key=lambda x: x["start"])

    return exon_data


def fetch_refseq_exons(nm_id):
    """
    Fetch exon coordinates for a RefSeq transcript from NCBI database.
    
    Args:
        nm_id: RefSeq transcript ID (e.g., "NM_000546.6")
    
    Returns:
        List of exon dictionaries with coordinates, or None if fetch fails
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
                
                start_str, end_str = loc.replace("<", "").replace(">", "").split("..")
                start = int(start_str)
                end = int(end_str)
                length = end - start + 1
                
                exons.append({
                    "exon": len(exons) + 1,
                    "start": start,
                    "end": end,
                    "strand": None,
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

    parser = argparse.ArgumentParser(description="Load MANE/RefSeq exon data and summarize contents")
    parser.add_argument(
        "gff",
        nargs="?",
        default="data/MANE.GRCh38.v1.4.refseq_genomic.gff.gz",
        help="Path to MANE GFF.gz or RefSeq Select GTF.gz file"
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    gff_path = args.gff
    if not os.path.isabs(gff_path):
        gff_path = os.path.join(script_dir, gff_path)
    if not os.path.exists(gff_path):
        print(f"Error: File not found: {gff_path}")
        sys.exit(2)

    print(f"Loading {gff_path}...")
    data = load_mane_exons(gff_path)
    print(f"Loaded {len(data)} genes")
    
    # Show sample
    if data:
        sample_gene = list(data.keys())[0]
        sample_transcript = list(data[sample_gene].keys())[0]
        sample_exons = data[sample_gene][sample_transcript]
        print(f"\nSample: {sample_gene} ({sample_transcript}) has {len(sample_exons)} exons")