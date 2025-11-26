"""
parser.py
---------

Handles all user input parsing for CNVision.

This module converts raw user input (genomic coordinates or transcript-based
coordinates) into standardized CNV dictionaries that the rest of the pipeline
can understand.

Supported Input Types
---------------------
1. Genomic CNV Input:
    {
        "type": "genomic",
        "chrom": "17",
        "start": 43044295,
        "end": 43045805,
        "cnv_type": "deletion"
    }

2. Transcript-Based CNV Input:
    {
        "type": "transcript",
        "gene": "BRCA1",
        "transcript": "NM_007294.4",
        "exons": [22, 23],
        "cnv_type": "duplication"
    }

The parser ensures all inputs are valid before downstream processing.
"""


def parse_genomic_input(chrom: str, start: int, end: int, cnv_type: str):
    """Parse user genomic coordinate input into internal CNV dict format."""
    if cnv_type not in {"deletion", "duplication"}:
        raise ValueError("cnv_type must be 'deletion' or 'duplication'.")

    if start >= end:
        raise ValueError("Start coordinate must be less than end coordinate.")

    return {
        "mode": "genomic",
        "chrom": str(chrom),
        "start": int(start),
        "end": int(end),
        "type": cnv_type,
    }


def parse_transcript_input(gene: str, transcript: str, exons: list, cnv_type: str):
    """Parse user transcript-based exon CNV input."""
    if not gene:
        raise ValueError("A gene symbol must be supplied.")

    if not transcript:
        raise ValueError("A transcript ID (e.g., NM_xxxxx) must be supplied.")

    if not isinstance(exons, (list, tuple)) or not exons:
        raise ValueError("Exons must be a non-empty list, e.g., [3, 4].")

    for e in exons:
        if not isinstance(e, int) or e <= 0:
            raise ValueError("All exon numbers must be positive integers.")

    if cnv_type not in {"deletion", "duplication"}:
        raise ValueError("cnv_type must be 'deletion' or 'duplication'.")

    return {
        "mode": "transcript",
        "gene": gene,
        "transcript": transcript,
        "exons": list(exons),
        "type": cnv_type,
    }


def parse_user_input(input_dict: dict):
    """
    Master parser – determines input type and routes accordingly.

    Parameters
    ----------
    input_dict : dict
        Raw user input in one of two forms (see file header).

    Returns
    -------
    dict :
        Standardized CNV dictionary ready for coordinate mapping.

    Raises
    ------
    ValueError :
        If input format is invalid or missing fields.
    """
    if "type" not in input_dict:
        raise ValueError("Input must contain key 'type' = 'genomic' or 'transcript'.")

    mode = input_dict["type"]

    if mode == "genomic":
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
        raise ValueError("type must be 'genomic' or 'transcript'.")
