"""
================================================================================
ISCN_PARSER.PY - Parse ISCN (International System for Human Cytogenomic 
                 Nomenclature) Notation
================================================================================

PURPOSE:
--------
ISCN is the standardized notation system used by clinical genetics laboratories
worldwide to describe chromosomal abnormalities. This module parses ISCN strings
to extract genomic information needed for CNV analysis.

BIOLOGICAL IMPACT:
------------------
ISCN notation is the "language" of clinical genetics. When a lab reports a CNV,
they use ISCN format. By supporting ISCN input, CNVision can directly analyze
clinical lab reports without manual conversion, reducing errors and saving time.

Example ISCN strings:
- arr[GRCh38] 18p11.32q23(102,328_79,093,443)×3
  Translation: Array analysis, GRCh38 build, chromosome 18, duplication (×3),
               coordinates 102,328 to 79,093,443

CODING CONCEPTS & PRINCIPLES:
------------------------------
1. REGULAR EXPRESSIONS (regex):
   - Pattern matching for extracting structured data from text
   - Example: r'×(\d+)' finds copy number after × symbol
   - Powerful but can be complex - used here for parsing structured notation

2. STRING MANIPULATION:
   - .strip() removes whitespace
   - .replace() removes formatting characters (commas, underscores)
   - .split() breaks strings into parts

3. PATTERN MATCHING:
   - re.search() finds first match of pattern
   - Groups in parentheses capture specific parts
   - Example: r'\[(GRCh\d+)\]' captures "GRCh38" from "[GRCh38]"

4. DEFAULT VALUES:
   - If pattern not found, use sensible defaults
   - Example: Default to GRCh38 if build not specified

5. DATA VALIDATION:
   - Check if required patterns found
   - Return None if parsing fails (caller can handle error)

================================================================================
"""

import re
from typing import Dict, Optional, Tuple, List, Any

def parse_iscn_notation(iscn_string: str) -> Optional[Dict[str, Any]]:
    """
    Parse ISCN notation string to extract genomic information.
    
    ISCN Format Examples:
        arr[GRCh38] 18p11.32q23(102,328_79,093,443)×3
        arr[GRCh37] 21q11.2q22.3(13,531,865_46,914,745)×3
        arr (18)×3
    
    Parsing Strategy:
    1. Extract genome build (GRCh37 or GRCh38)
    2. Extract chromosome number
    3. Extract copy number (determines deletion vs duplication)
    4. Extract genomic coordinates (if present)
    
    Args:
        iscn_string: ISCN notation string from clinical lab report
    
    Returns:
        Dictionary with parsed information:
        {
            "chromosome": "18",
            "start": 102328,
            "end": 79093443,
            "build": "GRCh38",
            "cnv_type": "duplication",
            "copy_number": 3
        }
        Returns None if parsing fails
    """
    # Remove leading/trailing whitespace
    iscn_string = iscn_string.strip()
    
    # ========================================================================
    # STEP 1: Extract genome build version
    # ========================================================================
    # Pattern: [GRCh38] or [GRCh37]
    # Regex explanation:
    #   \[      - Literal [ character (escaped because [ is special in regex)
    #   (GRCh\d+) - Capture group: "GRCh" followed by digits
    #   \]      - Literal ] character
    build_match = re.search(r'\[(GRCh\d+)\]', iscn_string)
    build = build_match.group(1) if build_match else "GRCh38"  # Default to GRCh38
    
    # ========================================================================
    # STEP 2: Extract chromosome number
    # ========================================================================
    # Pattern matches: "arr 18" or "arr[GRCh38] 18" or "arr (18)"
    # Can be: 1-22, X, Y, or range like "1-22"
    chrom_match = re.search(r'(?:arr|sseq)\s*(?:\[GRCh\d+\])?\s*\(?(\d+|X|Y|\d+-\d+)\)?', iscn_string)
    if not chrom_match:
        return None  # Can't parse without chromosome
    chromosome = chrom_match.group(1)
    
    # ========================================================================
    # STEP 3: Extract copy number
    # ========================================================================
    # Pattern: ×3, ×2, ×1, ×0
    # × symbol (multiplication sign) is used in ISCN notation
    # Copy number determines if it's deletion or duplication
    copy_match = re.search(r'×(\d+)', iscn_string)
    if not copy_match:
        return None  # Can't parse without copy number
    copy_number = int(copy_match.group(1))
    
    # ========================================================================
    # STEP 4: Determine CNV type from copy number
    # ========================================================================
    # Biological logic:
    #   Normal = 2 copies (diploid)
    #   ×3 or ×2 = duplication (extra copies)
    #   ×1 or ×0 = deletion (missing copies)
    # Biological interpretation:
    #   Normal = 2 copies (diploid)
    #   >2     = duplication (extra copies)
    #   <2     = deletion (missing copies)
    if copy_number > 2:
        cnv_type = "duplication"
    elif copy_number < 2:
        cnv_type = "deletion"
    else:
        # copy_number == 2 represents normal diploid state; treat as no CNV
        cnv_type = "normal"
    
    # ========================================================================
    # STEP 5: Extract genomic coordinates
    # ========================================================================
    # Pattern: (102,328_79,093,443)
    # Format: (start_end) with commas as thousands separators
    coord_match = re.search(r'\(([\d,]+)_([\d,]+)\)', iscn_string)
    if coord_match:
        # Remove commas and convert to integers
        start_str = coord_match.group(1).replace(',', '')
        end_str = coord_match.group(2).replace(',', '')
        try:
            start = int(start_str)
            end = int(end_str)
        except ValueError:
            return None  # Invalid coordinate format
    else:
        # Some ISCN strings don't include coordinates
        # User will need to provide them separately
        start = None
        end = None
    
    # ========================================================================
    # STEP 6: Return parsed data
    # ========================================================================
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
    """
    Find all genes that have exons overlapping a given genomic region.
    
    This is useful for ISCN input where you know coordinates but not which
    genes are affected. It searches through all genes to find overlaps.
    
    Note: Currently doesn't filter by chromosome (searches all genes).
    This could be optimized by adding chromosome information to MANE data.
    
    Args:
        chromosome: Chromosome number (currently not used for filtering)
        start: Genomic start coordinate
        end: Genomic end coordinate
        mane_data: MANE data structure from load_mane_exons
    
    Returns:
        List of gene symbols that have at least one exon overlapping the region
    
    Example:
        Region: start=43000000, end=43002000
        Finds: ["BRCA1"] (if BRCA1 has exons in that range)
    """
    overlapping_genes = []
    
    # ========================================================================
    # Search through all genes in MANE data
    # ========================================================================
    for gene, transcripts in mane_data.items():
        for transcript, exons in transcripts.items():
            # ================================================================
            # Check each exon for overlap with the region
            # ================================================================
            for exon in exons:
                exon_start = exon.get("start")
                exon_end = exon.get("end")
                
                # Skip exons with missing coordinate data
                if exon_start is None or exon_end is None:
                    continue
                
                # ============================================================
                # OVERLAP CHECK: Same algorithm as map_cnv_to_exons()
                # ============================================================
                # Two ranges overlap if:
                #   region_end >= exon_start AND region_start <= exon_end
                if end >= exon_start and start <= exon_end:
                    # Found overlap! Add gene to list (if not already there)
                    if gene not in overlapping_genes:
                        overlapping_genes.append(gene)
                    # Break: found overlap for this gene, no need to check more exons
                    break
    
    return overlapping_genes


def format_iscn_output(
    chromosome: str,
    start: int,
    end: int,
    build: str,
    cnv_type: str
) -> str:
    """
    Format CNV results in ISCN (International System for Human Cytogenomic 
    Nomenclature) notation.
    
    ISCN Format:
        arr[GRCh38] 18p11.32q23(102,328_79,093,443)×3
    
    Where:
        - arr = array analysis method
        - [GRCh38] = genome build
        - 18 = chromosome number
        - p11.32q23 = cytogenetic bands (simplified, we'll use coordinates)
        - (102,328_79,093,443) = genomic coordinates with commas
        - ×3 = copy number (×3 for duplication, ×1 for deletion)
    
    Args:
        chromosome: Chromosome number (e.g., "17", "X", "Y")
        start: Genomic start coordinate
        end: Genomic end coordinate
        build: Genome build (GRCh37 or GRCh38)
        cnv_type: "deletion" or "duplication"
    
    Returns:
        ISCN formatted string
    """
    # Determine copy number from CNV type
    # Normal = 2 copies
    # Deletion = 1 or 0 copies (we'll use ×1)
    # Duplication = 3 or more copies (we'll use ×3)
    copy_number = "×1" if cnv_type.lower() == "deletion" else "×3"
    
    # Format coordinates with commas as thousands separators
    # Example: 102328 -> "102,328"
    start_formatted = f"{start:,}"
    end_formatted = f"{end:,}"
    
    # Build ISCN string
    # Format: arr[GRCh38] 18(102,328_79,093,443)×3
    iscn_string = f"arr[{build}] {chromosome}({start_formatted}_{end_formatted}){copy_number}"
    
    return iscn_string

