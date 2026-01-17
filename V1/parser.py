"""
================================================================================
PARSER.PY - Input Validation and Normalization Module
================================================================================

PURPOSE:
--------
This module validates and normalizes user input before it reaches the analysis
pipeline. It ensures data integrity by:
- Checking required fields are present
- Validating data types and ranges
- Converting input to standardized format
- Providing clear error messages for invalid input

BIOLOGICAL IMPACT:
------------------
Input validation prevents errors that could lead to incorrect CNV predictions.
For example, catching invalid exon numbers or coordinates that don't make sense
prevents downstream analysis from producing misleading results.

CODING CONCEPTS & PRINCIPLES:
------------------------------
1. INPUT VALIDATION:
   - Check data types (int, str, list)
   - Validate ranges (start < end, positive integers)
   - Verify required fields present
   - Raise descriptive errors

2. DATA NORMALIZATION:
   - Convert different input formats to standard structure
   - Ensures consistent data format throughout pipeline
   - Makes downstream code simpler (one format to handle)

3. ERROR HANDLING:
   - Raise ValueError with descriptive messages
   - Fail fast (catch errors early, before processing)
   - Helpful error messages guide users to fix input

4. TYPE SAFETY:
   - Type hints in function signatures
   - Explicit type conversions (int(), str())
   - Validates types before use

5. DESIGN PATTERN:
   - "Parser" pattern: separates input validation from business logic
   - Single Responsibility: only validates/normalizes, doesn't analyze
   - Easy to test: can test validation logic independently

================================================================================
"""


# ============================================================================
# PARSER FUNCTIONS: Validate and normalize different input types
# ============================================================================

def parse_genomic_input(chrom: str, start: int, end: int, cnv_type: str):
    """
    Validate and normalize genomic coordinate input.
    
    Validates:
    - CNV type is valid ("deletion" or "duplication")
    - Start coordinate < end coordinate (logical requirement)
    - Coordinates are integers (enforced by type hints)
    
    Args:
        chrom: Chromosome number (as string, e.g., "17", "X")
        start: Genomic start coordinate
        end: Genomic end coordinate
        cnv_type: "deletion" or "duplication"
    
    Returns:
        Standardized CNV dictionary
    
    Raises:
        ValueError: If input is invalid
    """
    # ========================================================================
    # VALIDATION: Check CNV type
    # ========================================================================
    if cnv_type not in {"deletion", "duplication"}:
        raise ValueError("cnv_type must be 'deletion' or 'duplication'.")

    # ========================================================================
    # VALIDATION: Check coordinate logic
    # ========================================================================
    # Start must be less than end (can't have zero or negative length region)
    if start >= end:
        raise ValueError("Start coordinate must be less than end coordinate.")

    # ========================================================================
    # NORMALIZATION: Convert to standard format
    # ========================================================================
    return {
        "mode": "genomic",
        "chrom": str(chrom),      # Ensure string type
        "start": int(start),      # Ensure integer type
        "end": int(end),          # Ensure integer type
        "type": cnv_type,
    }


def parse_transcript_input(gene: str, transcript: str, exons: list, cnv_type: str):
    """
    Validate and normalize transcript-based exon input.
    
    Validates:
    - Gene symbol is provided (non-empty)
    - Transcript ID is provided (non-empty)
    - Exons is a non-empty list/tuple
    - All exon numbers are positive integers
    - CNV type is valid
    
    Args:
        gene: Gene symbol (e.g., "BRCA1")
        transcript: Transcript ID (e.g., "NM_007294.4")
        exons: List of exon numbers (e.g., [2, 3, 4])
        cnv_type: "deletion" or "duplication"
    
    Returns:
        Standardized CNV dictionary
    
    Raises:
        ValueError: If input is invalid
    """
    # ========================================================================
    # VALIDATION: Check required fields
    # ========================================================================
    if not gene:
        raise ValueError("A gene symbol must be supplied.")

    if not transcript:
        raise ValueError("A transcript ID (e.g., NM_xxxxx) must be supplied.")

    # ========================================================================
    # VALIDATION: Check exons format
    # ========================================================================
    # Must be a list or tuple (not a string, int, etc.)
    if not isinstance(exons, (list, tuple)) or not exons:
        raise ValueError("Exons must be a non-empty list, e.g., [3, 4].")

    # ========================================================================
    # VALIDATION: Check each exon number
    # ========================================================================
    # Exon numbers must be positive integers (exons are numbered 1, 2, 3...)
    for e in exons:
        if not isinstance(e, int) or e <= 0:
            raise ValueError("All exon numbers must be positive integers.")

    # ========================================================================
    # VALIDATION: Check CNV type
    # ========================================================================
    if cnv_type not in {"deletion", "duplication"}:
        raise ValueError("cnv_type must be 'deletion' or 'duplication'.")

    # ========================================================================
    # NORMALIZATION: Convert to standard format
    # ========================================================================
    return {
        "mode": "transcript",
        "gene": gene,
        "transcript": transcript,
        "exons": list(exons),  # Ensure list type (convert tuple if needed)
        "type": cnv_type,
    }


def parse_user_input(input_dict: dict):
    """
    Master parser function - routes input to appropriate parser based on type.
    
    This is the main entry point for parsing. It:
    1. Determines input type ("genomic" or "transcript")
    2. Validates required fields are present
    3. Routes to appropriate specialized parser
    4. Returns standardized CNV dictionary
    
    Args:
        input_dict: Raw user input dictionary with "type" key
    
    Returns:
        Standardized CNV dictionary ready for analysis pipeline
    
    Raises:
        ValueError: If input format is invalid or missing required fields
    
    Example Usage:
        input_dict = {
            "type": "genomic",
            "chrom": "17",
            "start": 43044295,
            "end": 43045805,
            "cnv_type": "deletion"
        }
        cnv = parse_user_input(input_dict)
    """
    # ========================================================================
    # STEP 1: Validate input type is specified
    # ========================================================================
    if "type" not in input_dict:
        raise ValueError("Input must contain key 'type' = 'genomic' or 'transcript'.")

    mode = input_dict["type"]

    # ========================================================================
    # STEP 2: Route to appropriate parser based on type
    # ========================================================================
    if mode == "genomic":
        # ====================================================================
        # Genomic mode: requires chromosome, start, end, cnv_type
        # ====================================================================
        required = ["chrom", "start", "end", "cnv_type"]
        for r in required:
            if r not in input_dict:
                raise ValueError(f"Missing field for genomic CNV: {r}")

        return parse_genomic_input(
            chrom=input_dict["chrom"],
            start=input_dict["start"],
            end=input_dict["end"],
            cnv_type=input_dict["cnv_type"],
        )

    elif mode == "transcript":
        # ====================================================================
        # Transcript mode: requires gene, transcript, exons, cnv_type
        # ====================================================================
        required = ["gene", "transcript", "exons", "cnv_type"]
        for r in required:
            if r not in input_dict:
                raise ValueError(f"Missing field for transcript CNV: {r}")

        return parse_transcript_input(
            gene=input_dict["gene"],
            transcript=input_dict["transcript"],
            exons=input_dict["exons"],
            cnv_type=input_dict["cnv_type"],
        )

    else:
        # Invalid type specified
        raise ValueError("type must be 'genomic' or 'transcript'.")
