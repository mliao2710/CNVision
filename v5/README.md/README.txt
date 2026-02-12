================================================================================
                          CNVision - CNV Impact Predictor
================================================================================

OVERVIEW
========

CNVision is a bioinformatics tool that analyzes Copy Number Variations (CNVs)
in human genes to predict their functional impact on protein function.

WHAT IS A CNV?
A Copy Number Variation is a section of DNA that appears in different numbers
of copies. For example:
  - Normal: 2 copies (everyone has 2 copies of most genes)
  - Deletion: 1 copy (one copy missing)
  - Duplication: 3 copies (extra copy)

WHY DOES THIS MATTER?
CNVs that affect exons (protein-coding regions) can cause disease by either:
  1. Deleting important amino acids → protein doesn't work
  2. Causing a "frameshift" → all downstream amino acids become wrong (severe)
  3. Maintaining reading frame → partial function preserved (less severe)

WHAT CNVISION DOES
Given a CNV location, CNVision:
  1. Identifies which exons are affected
  2. Calculates total length of affected exons
  3. Determines if it causes FRAMESHIFT or IN-FRAME effect
  4. Predicts functional consequence

Example:
  Input: "Gene BRCA1, exons 5-7 deletion"
  Analysis: 670 bp total → 670 % 3 = 1 → NOT divisible by 3
  Output: "FRAMESHIFT DELETION (likely severe impact)"


PROJECT STRUCTURE
=================

  /Users/milesliao/Desktop/cnvision/
  ├── web_app.py                                    [Flask web application]
  ├── main.py                                       [Command-line interface]
  ├── coordinate_mapper.py                          [Maps genomic coords to exons]
  ├── functional_predictor.py                       [Predicts frameshift/in-frame]
  ├── iscn_parser.py                                [Parses ISCN notation]
  ├── mane_loader.py                                [Loads gene/exon reference data]
  ├── parser.py                                     [Additional parsing utilities]
  ├── templates/
  │   └── index.html                                [Web UI - HTML/CSS/JavaScript]
  ├── static/
  │   └── logo.png                                  [Logo image]
  ├── utils/
  │   ├── ncbi_fetcher.py                           [Fetches gene data from NCBI]
  │   ├── gene_cache.py                             [Caches gene lookups]
  │   └── __init__.py
  ├── data/
  │   └── gene2refseq                               [Gene reference database]
  |   ├── GRCh37.ref_seq_select.gz                  [GRCh3 refseq select database from UCSB browser]
  │   └── MANE.GRCh38.v1.4.refseq_genomic.gff.gz    [MANE GRCh38 refseq database]
  ├── dockerfile                                    [Docker configuration]
  ├── pyproject.toml                                [Python project metadata]
  └── README.md                                     [Project documentation]


KEY CONCEPTS
============

1. FRAMESHIFT vs IN-FRAME
   The core concept: DNA is read in groups of 3 nucleotides (codons).
   
   - IN-FRAME: Total deletion length divisible by 3 (e.g., 270 bp)
     → Reading frame maintained
     → Only deleted region affected
     → Less severe
   
   - FRAMESHIFT: Total deletion length NOT divisible by 3 (e.g., 275 bp)
     → Reading frame shifts
     → All downstream amino acids wrong
     → Usually more severe

2. MANE DATABASE
   MANE = "Matched Annotation from NCBI and EBI"
   - Official reference database of human genes
   - One "best" transcript per gene (MANE Select)
   - ~20,000 human genes with exon coordinates
   - Available in GRCh37 and GRCh38 builds

3. GENOME BUILDS
   GRCh37 (hg19): Older reference, still used
   GRCh38 (hg38): Current reference
   → Same gene has different coordinates in each build!

4. INPUT FORMATS SUPPORTED
   Three ways to specify a CNV:
   
   a) Exon Numbers: "Gene FBN1, exons 42-45"
   b) cDNA Notation: "FBN1:c.5065_5546del"
   c) Genomic Coordinates: "DEL chr15:48,741,090-48,756,096"
      or ISCN: "arr[GRCh38]15q21.31(48,741,090_48,756,096)x1"


KEY FILES EXPLAINED
===================

web_app.py (778 lines)
  Flask web application providing browser-based interface
  - Handles HTTP requests from web UI
  - Coordinates analysis workflow
  - Returns formatted results to browser
  - Supports 3 input modes (exon, transcript, genomic)

main.py (321 lines)
  Command-line interface for terminal-based analysis
  - Interactive menu system
  - Same analysis as web interface
  - Useful for scripting and batch processing

coordinate_mapper.py (190 lines)
  Maps CNVs to affected exons
  - Takes genomic coordinates, identifies overlapping exons
  - Converts exon numbers ↔ genomic coordinates
  - Core algorithm: interval overlap detection

functional_predictor.py (163 lines)
  Predicts biological impact of CNVs
  - Core logic: sum exon lengths, check if divisible by 3
  - Determines frameshift vs in-frame
  - Returns consequence prediction

mane_loader.py (313 lines)
  Loads reference gene database
  - Parses GFF/GTF files (MANE databases)
  - Creates in-memory lookup structures
  - Fast exon coordinate retrieval

iscn_parser.py (103 lines)
  Parses clinical ISCN notation
  - Extracts genomic coordinates from microarray reports
  - Detects genome build and CNV type automatically
  - Finds genes in genomic regions

utils/ncbi_fetcher.py (128 lines)
  Fetches gene data from NCBI online
  - Queries NCBI APIs for non-MANE transcripts
  - Extracts gene symbol, exon coordinates
  - Fallback for transcripts not in MANE database

utils/gene_cache.py (36 lines)
  Local caching layer
  - Stores gene lookups to disk (JSON)
  - Reduces NCBI API calls
  - Improves performance


TECHNOLOGIES USED
=================

Languages:
  - Python 3.10+ (main programming language)
  - HTML5 / CSS3 / JavaScript (web interface)

Frameworks & Libraries:
  - Flask (web framework)
  - Jinja2 (HTML templating)
  - Bootstrap (CSS framework)
  - requests (HTTP client)
  - ElementTree (XML parsing)

Data Formats:
  - GFF/GTF (gene annotations)
  - JSON (caching)
  - XML (NCBI API responses)

Deployment:
  - Docker (containerization)
  - Local Flask development server


INSTALLATION & USAGE
====================

LOCAL DEPLOYMENT:
  1. Install Python 3.10+
  2. Install dependencies: pip install flask requests
  3. Download MANE data to data/ folder
  4. Run: python web_app.py
  5. Open browser: http://localhost:5000

DOCKER DEPLOYMENT:
  1. docker build -t cnvision .
  2. docker run -p 5000:5000 -v /path/to/data:/app/data cnvision
  3. Open browser: http://localhost:5000

CLI USAGE:
  1. Run: python main.py
  2. Follow interactive menu
  3. Results displayed in terminal


COMMON QUESTIONS
================

Q: What's the difference between GRCh37 and GRCh38?
A: Two different reference assemblies of the human genome. Same genes have
   different coordinates. Always use the correct build for your data!

Q: What if my transcript isn't in MANE?
A: CNVision automatically fetches from NCBI online database for non-MANE
   transcripts (requires internet connection).

Q: Can I use this for real clinical diagnosis?
A: This is a research/interpretation tool. For clinical use, consult genetic
   counselors and use certified clinical tools. This helps understand biology
   but shouldn't be sole basis for clinical decisions.

Q: What if a CNV partially affects an exon?
A: CNVision only counts exons that are FULLY within the CNV region. Partial
   exon overlaps require manual calculation.

Q: How accurate is the prediction?
A: The frameshift calculation is mathematically accurate (length % 3).
   However, biological impact also depends on: gene location, protein function,
   splice sites, etc. This gives one piece of the interpretation puzzle.


EXAMPLE ANALYSIS
================

User provides: "BRCA1, exons 10-12, deletion"

Step 1: Look up BRCA1 in MANE database
  → Found: NM_007294.4 (MANE Select transcript)

Step 2: Get exon coordinates
  → Exon 10: 150 bp
  → Exon 11: 180 bp
  → Exon 12: 165 bp
  → Total: 495 bp

Step 3: Check frameshift
  → 495 % 3 = 0 (divisible by 3)
  → IN-FRAME deletion

Step 4: Return result
  → "IN-FRAME DELETION (3 exons, 495 bp)"
  → "Lower severity: reading frame maintained"


TROUBLESHOOTING
===============

"Gene not found in MANE database"
  → Check spelling. Try NCBI lookup. Gene may not be in MANE.

"Could not fetch from NCBI"
  → Check internet connection. NCBI APIs may be down. Try again later.

"Exon numbers out of range"
  → Gene doesn't have that many exons. Use valid exon numbers.

"Could not parse coordinates"
  → Check format: chr17:43000000-43002000 or arr[GRCh38]17q21.1(...)x1

Flask server won't start
  → Check Python 3.10+: python --version
  → Install dependencies: pip install flask requests
  → Port 5000 may be in use: lsof -i :5000


FUTURE ENHANCEMENTS
===================

- Support for partial exon deletions
- Intron-spanning deletions
- Copy number > 3
- Disease database integration (ClinVar)
- Protein structure impact visualization
- Batch processing
- RESTful API endpoints
- Mobile app


STRUCTURE OF ANALYSIS
=====================

Input → Parse → Lookup → Map → Predict → Output
  ↓       ↓       ↓      ↓       ↓        ↓
User    Extract  Query  Find   Sum     Display
input   format   MANE   exons  lengths  result


For more detailed technical information, see explanation.txt


================================================================================
                         Questions? See explanation.txt
                        or visit GitHub repository
================================================================================

  Prints results cleanly

Example output:

Loaded 19,388 genes from MANE
Gene: BRCA1, Transcript: NM_007294.4, Exons hit: [23], Consequence: frameshift deletion
Gene: DMD, Transcript: NM_004006.3, Exons hit: [12, 13], Consequence: multi-exon in-frame deletion


This file is meant to be modified by you for testing CNVision on any samples or CNVs you want.

Installation

Navigate to the folder containing the cnvision/ directory and run:

pip install -e.


If you don’t have a setup.py, replace this with:

export PYTHONPATH="$PYTHONPATH:/path/to/your/project"


Or place your repo into a virtual environment.

How to Run the Example
python cnvision/examples/test_cnv.py

Usage (Programmatically)
from cnvision.mane_loader import load_mane_exons
from cnvision.coordinate_mapper import map_cnv_to_exons
from cnvision.functional_predictor import predict_cnv_effect
