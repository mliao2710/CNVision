CNVision: A Computational Framework to Predict Functional Impact of Intragenic Copy Number Variants

Intragenic copy number variants (CNVs), also known as exonic or monogenic CNVs, are deletions or duplications of genomic material occurring within a single gene. These alterations can range from a single exon to multiple exons or an entire gene. Intragenic CNVs can lead to diverse molecular and functional consequences, including shifts in the reading frame that introduce premature stop codons and trigger nonsense-mediated decay, or in-frame losses or gains that alter protein length and function. With the expanding clinical use of high-resolution microarrays and next-generation sequencing, intragenic CNVs are now identified at an unprecedented rate. Recent studies estimate that intragenic CNVs account for 4.7% to 35% of pathogenic variants in Mendelian disease genes. Despite their clinical significance, existing annotation pipelines often lack the precision to interpret these variants accurately. To address this gap, we developed CNVision, a program designed to analyze the functional impact of intragenic CNVs using standardized genomic resources and scalable computational infrastructure.

This tool is designed for:
  Clinical genomics pipelines
  Variant interpretation
  Research involving exon-level CNV annotation
  Efficient universal CNV identification and classification

CNVision supports two CNV input modes:
  Genomic coordinates (start, end)  
  Exon numbers (exons_hit: [3, 4])
And automatically determines which transcripts are affected and what the predicted functional impact is.

Features
✅ Load MANE v1.4 exon annotations from GFF
✅ Map CNVs to specific MANE transcripts
✅ Accept CNVs as coordinates or exon numbers
✅ Predict consequences:
  frameshift deletion
  frameshift duplication
  in-frame deletion
  exon loss
  whole-gene duplication
  no coding impact
✅ Clean, modular design
✅ Example script showing user-input CNVs
Project Structure
cnvision/
│
├── __init__.py
├── mane_loader.py
├── coordinate_mapper.py
├── functional_predictor.py
│
├── data/
│   └── MANE.GRCh38.v1.4.refseq_genomic.gff.gz
│
├── examples/
│   └── test_cnv.py
│
└── README.md

File-by-File Explanation
mane_loader.py
  Responsible for loading MANE exon-level genomic annotations.
Purpose
  Parse the MANE RefSeq GFF file
  Extract gene → transcript → ordered exon coordinates
  Return structured dictionaries used for CNV mapping
  Key Output
  {
    "BRCA1": {
        "NM_007294.4": [
            {"exon": 1, "start": ..., "end": ...},
            {"exon": 2, "start": ..., "end": ...},
            ...
        ]
    }
  }

coordinate_mapper.py
  Maps CNVs to exons.
Purpose
Given either:
  genomic coordinates, or
  exon numbers
It returns the list of exons hit per transcript.
Key Responsibilities
  Overlap checking
  Range comparisons
  Flexible input type handling
  Returning consistent structured output

functional_predictor.py
  Predicts biological consequences of the CNV.
Purpose
Decides whether a CNV causes:
  frameshift deletion
  frameshift duplication
  exon loss
  partial exon truncation
  intronic/no coding impact
  Logic Considerations
  Coding frame phase
  Exon boundaries
  Deletion vs duplication
  Multi-exon involvement

data/MANE.GRCh38.v1.4.refseq_genomic.gff.gz
  The official MANE v1.4 GFF file used for exon coordinates.

examples/test_cnv.py
  An example script demonstrating how the package works.
Purpose
  Loads MANE annotations
  Accepts both coordinate-based and exon-based CNVs
  Maps them to transcripts
  Predicts consequences
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
