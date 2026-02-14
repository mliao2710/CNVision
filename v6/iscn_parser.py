"""
CNVision ISCN Parser - Clinical String Parsing Helpers

Purpose:
- Support ISCN-style genomic notation parsing for coordinate-mode workflows.

Core Functions:
- Parse ISCN text into normalized fields (build, chromosome, coordinates, copy number).
- Derive CNV type from copy number for deletion/duplication hints.
- Find genes with exon overlap in a queried genomic interval.
- Provide optional ISCN-style output formatting for display/use downstream.
"""

import re
from typing import Dict, Optional, List, Any

def parse_iscn_notation(iscn_string: str) -> Optional[Dict[str, Any]]:
    """Parse ISCN text into normalized genomic fields."""
    iscn_string = iscn_string.strip()
    
    build_match = re.search(r'\[(GRCh\d+)\]', iscn_string)
    build = build_match.group(1) if build_match else "GRCh38"
    
    chrom_match = re.search(r'(?:arr|sseq)\s*(?:\[GRCh\d+\])?\s*\(?(\d+|X|Y|\d+-\d+)\)?', iscn_string)
    if not chrom_match:
        return None
    chromosome = chrom_match.group(1)
    
    copy_match = re.search(r'×(\d+)', iscn_string)
    if not copy_match:
        return None
    copy_number = int(copy_match.group(1))
    
    if copy_number > 2:
        cnv_type = "duplication"
    elif copy_number < 2:
        cnv_type = "deletion"
    else:
        cnv_type = "normal"
    
    coord_match = re.search(r'\(([\d,]+)_([\d,]+)\)', iscn_string)
    if coord_match:
        start_str = coord_match.group(1).replace(',', '')
        end_str = coord_match.group(2).replace(',', '')
        try:
            start = int(start_str)
            end = int(end_str)
        except ValueError:
            return None
    else:
        start = None
        end = None
    
    return {
        "chromosome": chromosome,
        "start": start,
        "end": end,
        "build": build,
        "cnv_type": cnv_type,
        "copy_number": copy_number
    }


def find_genes_in_region(
    chromosome: str,
    start: int,
    end: int,
    mane_data: Dict[str, Dict[str, List[Dict[str, Any]]]]
) -> List[str]:
    """Return genes with any exon overlap in a genomic interval."""
    overlapping_genes = []
    
    # STEP 1: Normalize query chromosome.
    chromosome = str(chromosome).replace('chr', '').replace('Chr', '').upper()
    
    # STEP 2: Scan MANE exons for overlap.
    for gene, transcripts in mane_data.items():
        for transcript, exons in transcripts.items():
            for exon in exons:
                exon_start = exon.get("start")
                exon_end = exon.get("end")
                exon_chr = exon.get("chromosome")
                
                # Skip exons with incomplete coordinates.
                if not exon_chr or exon_start is None or exon_end is None:
                    continue
                
                # Compare normalized chromosome labels.
                exon_chr_normalized = str(exon_chr).replace('chr', '').replace('Chr', '').upper()
                if exon_chr_normalized != chromosome:
                    continue 
                
                # Record gene on first overlapping exon.
                if end >= exon_start and start <= exon_end:
                    if gene not in overlapping_genes:
                        overlapping_genes.append(gene)
                    break
    
    return overlapping_genes


def format_iscn_output(
    chromosome: str,
    start: int,
    end: int,
    build: str,
    cnv_type: str
) -> str:
    """Format normalized values as an ISCN-style string."""
    copy_number = "×1" if cnv_type.lower() == "deletion" else "×3"
    start_formatted = f"{start:,}"
    end_formatted = f"{end:,}"
    iscn_string = f"arr[{build}] {chromosome}({start_formatted}_{end_formatted}){copy_number}"
    return iscn_string
