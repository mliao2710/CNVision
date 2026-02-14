"""
CNVision MANE Loader - Reference Annotation Parser

Purpose:
- Load MANE/RefSeq annotation files into CNVision in-memory data structures.

Core Functions:
- Parse compressed GFF/GTF records into `gene -> transcript -> exon` maps.
- Resolve transcript IDs to gene symbols via local gene2refseq mapping.
- Attach CDS segment overlap and CDS length to each exon entry.
- Normalize exon numbering and strand-aware ordering for downstream analysis.
"""

import gzip
import requests
from xml.etree import ElementTree as ET
import re
import os

# Cache for transcript -> gene symbols from gene2refseq.
TRANSCRIPT_TO_GENE_CACHE = {}

def load_gene_mapping(mapping_file="data/gene2refseq"):
    """Load human transcript->gene mapping from `gene2refseq`."""
    mapping = {}
    
    # STEP 1: Resolve mapping file path.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(script_dir, mapping_file)
    
    if not os.path.exists(full_path):
        print(f"Warning: Gene mapping file not found at {full_path}")
        print("Will use transcript IDs as gene names. Download with:")
        print("  cd data && curl -O https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz && gunzip gene2refseq.gz")
        return mapping
    
    print(f"Loading gene mappings from {mapping_file}...")
    
    # STEP 2: Parse mapping rows.
    with open(full_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 16:
                continue
            
            tax_id = fields[0]
            # Keep human genes only.
            if tax_id != '9606':
                continue
            
            rna_acc = fields[3]
            gene_symbol = fields[15]
            
            if rna_acc and rna_acc != '-' and gene_symbol and gene_symbol != '-':
                # Store both versioned and base accession forms.
                mapping[rna_acc] = gene_symbol
                base_id = rna_acc.split('.')[0]
                mapping[base_id] = gene_symbol
    
    print(f"Loaded {len(mapping)} transcript → gene mappings")
    return mapping


def fetch_gene_from_transcript(transcript_id):
    """Resolve a transcript accession to gene symbol, with fallback."""
    # STEP 1: Cache lookup with version and base accession.
    if transcript_id in TRANSCRIPT_TO_GENE_CACHE:
        return TRANSCRIPT_TO_GENE_CACHE[transcript_id]
    
    # Try without version
    base_id = transcript_id.split('.')[0]
    if base_id in TRANSCRIPT_TO_GENE_CACHE:
        return TRANSCRIPT_TO_GENE_CACHE[base_id]
    
    # STEP 2: Fallback to base transcript accession.
    return base_id


def load_mane_exons(gff_path):
    """Load exon/CDS annotations into `{gene: {transcript: [exons]}}`."""
    global TRANSCRIPT_TO_GENE_CACHE
    
    # STEP 1: Ensure transcript->gene mapping cache is populated.
    if not TRANSCRIPT_TO_GENE_CACHE:
        TRANSCRIPT_TO_GENE_CACHE = load_gene_mapping()
    
    exon_data = {}
    cds_data = {}
    transcript_to_gene = {}  # Local cache for this file
    
    # STEP 2: Parse feature rows.
    with gzip.open(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields

            # Keep exon and CDS features only.
            if feature_type not in {"exon", "CDS"}:
                continue

            # Parse attributes for GTF or GFF syntax.
            attr_dict = {}
            
            # Check if it's GTF format (UCSC RefSeq Select) or GFF format (MANE)
            if 'gene_id "' in attrs:
                for match in re.finditer(r'(\w+)\s+"([^"]+)"', attrs):
                    key = match.group(1)
                    value = match.group(2)
                    attr_dict[key] = value
            else:
                for item in attrs.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attr_dict[key] = value

            # Extract transcript ID; CDS rows may use Parent=rna-NM_xxx.
            transcript = attr_dict.get("transcript_id") or attr_dict.get("gene_id")
            if not transcript:
                parent = attr_dict.get("Parent")
                if parent:
                    transcript = parent.split(",")[0]
                    if transcript.startswith("rna-"):
                        transcript = transcript[4:]
            if not transcript:
                continue

            # Extract or resolve gene symbol.
            gene = (
                attr_dict.get("gene")
                or attr_dict.get("gene_name")
                or attr_dict.get("gene_id")
                or attr_dict.get("Name")
            )
            
            # Resolve missing symbols via transcript mapping.
            if not gene or gene.startswith("NM_"):
                if transcript in transcript_to_gene:
                    gene = transcript_to_gene[transcript]
                else:
                    fetched_gene = fetch_gene_from_transcript(transcript)
                    if fetched_gene:
                        gene = fetched_gene
                        transcript_to_gene[transcript] = gene
                    else:
                        continue

            # Extract exon number if available.
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

            if not gene or not transcript:
                continue

            # Parse coordinate fields.
            start = int(start)
            end = int(end)
            length = end - start + 1

            # Store exons and CDS segments separately; merge later.
            if feature_type == "exon":
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
            elif feature_type == "CDS":
                cds_data.setdefault(gene, {})
                cds_data[gene].setdefault(transcript, [])
                cds_data[gene][transcript].append({
                    "start": start,
                    "end": end,
                    "phase": phase
                })

    # STEP 3: Normalize exon order and attach CDS overlap per exon.
    for gene in exon_data:
        for transcript in exon_data[gene]:
            exons = exon_data[gene][transcript]
            
            # Determine transcript strand.
            strand = exons[0]["strand"] if exons else "+"
            
            # Start with genomic ordering.
            exons_sorted = sorted(exons, key=lambda x: x["start"])
            
            if strand == "-":
                # On minus strand, biological exon order is reverse genomic order.
                exons_sorted.reverse()
            
            # Fill missing exon indices.
            for idx, exon in enumerate(exons_sorted, start=1):
                if exon["exon"] is None:
                    exon["exon"] = idx
            
            # Attach CDS segments and CDS lengths per exon.
            tx_cds_segments = cds_data.get(gene, {}).get(transcript, [])
            for exon in exons_sorted:
                cds_segments = []
                for cds in tx_cds_segments:
                    overlap_start = max(exon["start"], cds["start"])
                    overlap_end = min(exon["end"], cds["end"])
                    if overlap_start <= overlap_end:
                        cds_segments.append({
                            "start": overlap_start,
                            "end": overlap_end,
                            "phase": cds.get("phase")
                        })

                exon["cds_segments"] = cds_segments
                exon["cds_length"] = sum(seg["end"] - seg["start"] + 1 for seg in cds_segments)

            # Store in genomic order for downstream interval operations.
            exon_data[gene][transcript] = sorted(exons_sorted, key=lambda x: x["start"])

    return exon_data


def fetch_refseq_exons(nm_id):
    """Fetch exon coordinates for a RefSeq transcript from NCBI."""
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
    
    # Show one sample transcript for a quick sanity check.
    if data:
        sample_gene = list(data.keys())[0]
        sample_transcript = list(data[sample_gene].keys())[0]
        sample_exons = data[sample_gene][sample_transcript]
        print(f"\nSample: {sample_gene} ({sample_transcript}) has {len(sample_exons)} exons")
