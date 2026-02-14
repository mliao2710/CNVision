================================================================================
                      CNVISION - FULL TECHNICAL PROJECT EXPLANATION
================================================================================

Audience:
- This document is written for readers with minimal biology and minimal
  programming background.
- Goal: after reading, you should understand what CNVision does, why it exists,
  how each file contributes, and how data flows through the system.


TABLE OF CONTENTS
================================================================================
1. Project Mission (Plain Language)
2. Biology Basics You Need
3. Why CNVision Exists
4. What the Tool Predicts vs What It Cannot Predict
5. End-to-End System Architecture
6. Input Modes and What They Mean
7. Deep File-by-File Explanation
8. How Files Connect to Each Other (Dependency Map)
9. Data Structures Used Internally
10. Walkthrough of a Real Analysis Request
10A. Literal Execution Trace (Function-by-Function)
11. Output Interpretation Guide
12. Common Pitfalls and Why Results Sometimes Differ
13. Local Runtime and Docker Runtime
14. Practical Validation Strategy
15. Current Scope and Future Improvements


1) PROJECT MISSION (PLAIN LANGUAGE)
================================================================================
CNVision is a CNV interpretation tool.

Given a copy-number variant (CNV), it answers a focused question:

  "Does this coding change preserve reading frame (in-frame)
   or disrupt reading frame (frameshift/out-of-frame)?"

CNVision is not trying to solve all clinical genetics questions. It is focused
on one core interpretation component that is frequently needed and often
repeated manually.


2) BIOLOGY BASICS YOU NEED
================================================================================
2.1 What is a CNV?
- CNV = Copy Number Variant.
- A region of DNA has fewer or more copies than expected.

Typical copy states:
- x2 = typical diploid state
- x1 = deletion (one copy missing)
- x3 = duplication (one extra copy)

2.2 Why exons matter
- Genes contain exons and introns.
- Exons are the transcribed segments retained in mature RNA.
- Protein-coding sequence (CDS) is the part of exon sequence that is translated.
- Not all exon bases are coding (UTRs exist at transcript ends).

2.3 Reading frame
- Proteins are read in codons: groups of 3 nucleotides.
- If coding-length change is divisible by 3 -> in-frame.
- If not divisible by 3 -> frameshift/out-of-frame.

2.4 Why this is clinically relevant
- Frameshift variants often have stronger functional disruption.
- In-frame variants may preserve some downstream protein structure.
- This is a general rule, not a complete pathogenicity verdict.


3) WHY CNVISION EXISTS
================================================================================
Many workflows produce raw coordinates first, not immediate biological meaning.
Examples:
- Microarray reports
- Structural variant pipelines
- Copy-number calling software outputs

People then need to manually:
- Parse weird coordinate strings
- Resolve genome build
- Find affected exons
- Compute coding-length effect
- Decide in-frame vs out-of-frame

CNVision automates these repetitive steps and standardizes the result format.


4) WHAT THE TOOL PREDICTS VS WHAT IT CANNOT PREDICT
================================================================================
4.1 What it predicts reliably
- In-frame deletion/duplication
- Frameshift deletion/duplication
- No CDS impact (for overlaps outside coding sequence)

4.2 What it does not claim to predict
- Full ACMG/clinical pathogenicity class
- Exact protein structural damage
- Splicing regulatory consequences in depth
- Tissue-specific expression consequences
- All transcript-isoform-specific outcomes unless transcript is represented by
  available references


5) END-TO-END SYSTEM ARCHITECTURE
================================================================================
High-level architecture:

User (browser) -> Flask web server -> parser/mapping/predictor modules
               -> rendered HTML result

Pipeline stages:
1. Input collection (web form)
2. Input parsing and normalization
3. Build/transcript/gene resolution
4. Coordinate overlap mapping to exons/CDS
5. Frame consequence classification
6. Result rendering with explanatory context

Core principle:
- Web layer orchestrates.
- Mapping layer computes overlaps.
- Predictor layer classifies consequence.
- Loader layer provides reference annotations.


6) INPUT MODES AND WHAT THEY MEAN
================================================================================
6.1 Exon Number Mode
You provide:
- Gene (or transcript)
- First exon, last exon
- CNV type

Used when reports already say "exons X-Y affected".

6.2 Genomic Coordinate Mode
You provide:
- Genomic coordinate expression (multiple formats accepted)

Parser supports flexible forms including:
- chr-based coordinates
- arr/seq/sseq-like strings
- build aliases (hg19/hg38)
- copy-number hints (x1/x3)
- varied separators (comma, underscore, parentheses, etc.)

6.3 HGVS (g.) Mode
You provide:
- HGVS genomic notation (example:
  NC_000019.9:g.(13325430_13335371)_(13346648_13352244)del)

System treats HGVS g. as genomic coordinates by:
- extracting chromosome from NC accession (NC_000019 -> chr19)
- using the middle breakpoints as analysis interval (13335371-13346648)
- inferring CNV type from del/dup suffix
- inferring GRCh37/GRCh38 from NC accession version when recognized


7) DEEP FILE-BY-FILE EXPLANATION
================================================================================
This section explains not just "what each file does", but why it exists and
what part of the pipeline it owns.

-------------------------------------------------------------------------------
7.1 web_app.py
-------------------------------------------------------------------------------
What it is:
- The main application entry point.
- Runs the Flask server and owns all web request handling.

Why it exists:
- It is the coordinator that connects all other modules together.
- Without it, there is no user-facing application workflow.

Main responsibilities:
- Load reference data at startup (GRCh37 and GRCh38).
- Parse and normalize incoming form inputs.
- Dispatch to mode-specific analysis functions.
- Combine warnings, hints, and results into render-ready objects.
- Return results to `templates/index.html`.

Key internal functions:
- `fetch_refseq_info(...)`
  - NCBI fallback for transcript metadata/exons when needed.
- `normalize_build(...)`
  - Converts aliases like hg19/hg38 to GRCh labels.
- `parse_genomic_input(...)`
  - Robust parser for many genomic string formats.
- `parse_hgvs_genomic_notation(...)`
  - Parses HGVS g. strings into chromosome/start/end/build/type hints.
- `process_exon_mode(...)`
  - Exon-range analysis pipeline.
- `process_genomic_mode(...)`
  - Genomic string to consequence pipeline.
- `process_hgvs_mode(...)`
  - HGVS g. to genomic consequence pipeline.
- `index()`
  - Flask route for form display and submission handling.

How it connects to other files:
- Calls `load_mane_exons` from `mane_loader.py`.
- Calls mapping functions from `coordinate_mapper.py`.
- Calls consequence classification from `functional_predictor.py`.
- Calls overlap-based gene detection from `iscn_parser.py`.

-------------------------------------------------------------------------------
7.2 mane_loader.py
-------------------------------------------------------------------------------
What it is:
- Annotation ingestion and normalization layer.

Why it exists:
- Reference files are large and raw (GFF/GTF format).
- The app needs a consistent in-memory structure for fast lookup.

Main responsibilities:
- Read compressed annotation files from `data/`.
- Parse transcript/exon/CDS records.
- Resolve transcript->gene mapping via local `gene2refseq`.
- Normalize exon numbering and strand-aware order.
- Attach CDS segments and CDS length to exon entries.

Critical output:
- A nested dictionary:
  `mane_data[gene][transcript] = [exon_dict, exon_dict, ...]`

Why CDS attachment matters:
- CNVision frame logic is coding-aware.
- It needs exon-level CDS overlap info to compute coding-length impact.

How it connects:
- `web_app.py` calls it during startup.
- `coordinate_mapper.py` and `functional_predictor.py` consume its structures.

-------------------------------------------------------------------------------
7.3 coordinate_mapper.py
-------------------------------------------------------------------------------
What it is:
- Interval math engine.

Why it exists:
- Mapping coordinates to biologically meaningful exon hits is a separate
  concern from parsing strings or rendering UI.

Main responsibilities:
- `map_cnv_to_exons(...)`
  - Given genomic CNV interval + transcript exon list:
    - determine overlap per exon
    - compute overlap start/end
    - compute coding overlap length
    - flag full exon vs partial exon coverage
- `map_exon_numbers_to_regions(...)`
  - Given exon range:
    - return genomic min/max span for those exons

Important design choice:
- Outputs include both:
  - `genomic_length` (full interval overlap on genome)
  - `length` (coding overlap length used for frame logic)

How it connects:
- Called by `web_app.py` mode handlers.
- Its output is passed directly to `functional_predictor.py`.

-------------------------------------------------------------------------------
7.4 functional_predictor.py
-------------------------------------------------------------------------------
What it is:
- Consequence classification layer.

Why it exists:
- Frame consequence logic should be isolated from parser/web concerns.

Main responsibilities:
- Normalize hit payloads into consistent exon-length data.
- Sum coding lengths from mapped hits.
- Apply modulo-3 rule for frame classification.
- Return standard human-readable consequence labels.

Outcome patterns:
- `in-frame deletion (N exons)`
- `frameshift deletion (N exons)`
- `in-frame duplication ...`
- `frameshift duplication ...`
- `no CDS impact (...)`

How it connects:
- Always downstream of `coordinate_mapper.py` in the analysis flow.

-------------------------------------------------------------------------------
7.5 iscn_parser.py
-------------------------------------------------------------------------------
What it is:
- ISCN helper module and region overlap utility.

Why it exists:
- Some genomic inputs are expressed in clinical ISCN-like notation.
- Gene detection in a genomic interval is needed when no gene hint is provided.

Main responsibilities:
- Parse core ISCN fields where needed.
- Find genes whose exons overlap a provided region.
- (Optional utility) format ISCN-style output string.

How it connects:
- `web_app.py` uses `find_genes_in_region(...)` in genomic mode when gene hint
  is absent.

-------------------------------------------------------------------------------
7.6 templates/index.html
-------------------------------------------------------------------------------
What it is:
- Presentation layer (HTML/CSS + Jinja templating).
- Uses Bootstrap components plus custom CSS for layout and visual consistency.

Why it exists:
- Separates result rendering from computational logic.
- Ensures technical results are formatted so users can read and compare outputs
  quickly without misinterpreting dense data.

Main responsibilities:
- Render input forms for three modes.
- Render errors/warnings/results.
- Display summary values and detailed exon table.
- Show concise exon range, length checks, and overlap notes.
- Apply consistent formatting hierarchy (labels, values, tables, notes) so
  results are understandable at a glance.

How it connects:
- Receives `result` and `error` objects from `web_app.py`.

-------------------------------------------------------------------------------
7.7 static/logo.png
-------------------------------------------------------------------------------
What it is:
- Brand/UI asset.

Why it exists:
- Visual identity for web interface.

-------------------------------------------------------------------------------
7.8 data/ directory
-------------------------------------------------------------------------------
What it is:
- Reference annotations and transcript mapping files.

Why it exists:
- CNVision inference quality depends on annotation quality and build matching.

Required files:
- `MANE.GRCh38.v1.4.refseq_genomic.gff.gz`
- `GRCh37.ref_seq_select.gz`
- `gene2refseq`


8) HOW FILES CONNECT TO EACH OTHER (DEPENDENCY MAP)
================================================================================
Dependency direction (simplified):

`web_app.py`
  -> `mane_loader.py`
  -> `coordinate_mapper.py`
  -> `functional_predictor.py`
  -> `iscn_parser.py`
  -> `templates/index.html`

Key idea:
- `web_app.py` orchestrates, other modules are specialized workers.

Data-flow example:
1. `web_app.py` parses input
2. gets annotation dataset from `mane_loader.py`
3. maps overlaps via `coordinate_mapper.py`
4. classifies consequence via `functional_predictor.py`
5. renders in `templates/index.html`


9) DATA STRUCTURES USED INTERNALLY
================================================================================
9.1 Annotation store (`mane_data`)

Shape:
- `mane_data[gene][transcript] -> list of exons`

Exon entry contains fields such as:
- exon number
- start/end genomic coordinates
- strand
- chromosome
- CDS segments
- CDS length

9.2 CNV mapping payload
- transcript
- hit_exons list where each item includes:
  - exon id
  - overlap start/end
  - coding length (`length`)
  - genomic overlap length (`genomic_length`)
  - full/partial coverage flag

9.3 Prediction result payload
- gene
- transcript
- hit_exons
- predicted_consequence
- optional warnings/cnv_type/HGVS metadata


10) WALKTHROUGH OF A REAL ANALYSIS REQUEST
================================================================================
Example request:
- `del hg38:chr17:56,594,272-56,594,337` with gene hint

Step-by-step:
1. Parser extracts:
   - build GRCh38
   - chromosome 17
   - start/end coordinates
   - deletion hint
2. Build dataset selected (GRCh38 annotation)
3. Gene/transcript resolved
4. Exon overlap computed
5. Coding overlap lengths summed
6. `% 3` computed
7. Result rendered with exon table + summary lines + coverage note

Why this is useful:
- Users see both the verdict and the evidence behind the verdict.


10A) LITERAL EXECUTION TRACE (FUNCTION-BY-FUNCTION)
================================================================================
This section shows the exact function path for one real request.

Sample request (genomic mode):
- Input: `del hg38:chr17:56,594,272-56,594,337`
- Gene hint: `NOG`
- CNV type dropdown: `deletion`

Trace:

1. Browser sends POST `/`
- Function called: `web_app.index()`
- Reads form values:
  - `mode = "coordinate"`
  - `genomic_input = "del hg38:chr17:56,594,272-56,594,337"`
  - `gene_hint = "NOG"`
  - `cnv_type = "deletion"`
  - `build = "GRCh38"`

2. Mode dispatch
- `index()` calls:
  - `process_genomic_mode(genomic_input, build, cnv_type, gene_hint)`

3. Input parsing
- `process_genomic_mode(...)` calls:
  - `parse_genomic_input(...)`
- Parsed object:
  {
    "chromosome": "17",
    "start": 56594272,
    "end": 56594337,
    "build": "GRCh38",
    "copy_number": null,
    "cnv_type_hint": "deletion"
  }

4. Build and annotation resolution
- `get_mane_data("GRCh38")` returns:
  - GRCh38 MANE dataset
  - gene-name index
  - transcript->gene index

5. Gene/transcript resolution
- Because `gene_hint` exists, gene is set to `NOG` (canonicalized).
- Transcript selected from MANE for NOG (example: `NM_005450.6`).

6. CNV payload creation
- `process_genomic_mode(...)` builds:
  {
    "gene": "NOG",
    "start": 56594272,
    "end": 56594337,
    "type": "deletion",
    "chromosome": "chr17"
  }

7. Exon/CDS overlap mapping
- Calls:
  - `coordinate_mapper.map_cnv_to_exons(cnv, mane_data, transcript)`
- Mapper does:
  - interval overlap test exon-by-exon
  - overlap start/end extraction
  - coding overlap length calculation from CDS segments
  - full/partial exon coverage flag
- Returned hit object (shape):
  [
    {
      "transcript": "NM_005450.6",
      "hit_exons": [
        {
          "exon": 1,
          "start": 56594272,
          "end": 56594337,
          "length": 66,
          "genomic_length": 66,
          "full_exon_coverage": false
        }
      ],
      "cnv_type": "deletion",
      "predicted_consequence": null
    }
  ]

8. Consequence classification
- Calls:
  - `functional_predictor.predict_cnv_effect(hits, mane_data)`
- Predictor computes:
  - total coding length = 66
  - 66 % 3 = 0
  - classification = in-frame
- Result label:
  - `in-frame deletion (1 exon)`

9. Result decoration
- `process_genomic_mode(...)` adds warnings/CNV type metadata as needed.

10. HTML rendering
- `index()` returns:
  - `render_template(\"index.html\", result=result, error=None)`
- Template displays:
  - gene, transcript, consequence
  - exon table + concise range
  - total length + amino-acid check
  - coverage summary and partial-overlap note (when applicable)


11) OUTPUT INTERPRETATION GUIDE
================================================================================
11.1 Consequence line
- Primary classification output.

11.2 Exons hit table
- Shows exactly which exons were affected and coding lengths used.

11.3 Total length / amino-acid check
- Quick validation of frame math (`total / 3`).

11.4 Exon coverage summary
- Indicates whole-exon vs partial overlap context.
- Helps explain why genomic mode can differ from exon-number mode.


12) COMMON PITFALLS AND WHY RESULTS SOMETIMES DIFFER
================================================================================
12.1 Build mismatch (GRCh37 vs GRCh38)
- Same gene has different coordinates across builds.
- Wrong build can map to wrong exons or no gene.

12.2 Transcript mismatch
- Different transcripts of same gene may have different exon numbering and
  boundaries.
- Validation tables based on different transcript models can disagree.

12.3 Partial overlap vs full exon assumptions
- Exon-number mode assumes full exons by definition.
- Genomic/HGVS modes can hit partial exon segments.
- This can change coding-length totals and therefore frame classification.

12.4 Clinical context beyond frame logic
- Frame consequence is one important signal, not full diagnosis.


13) LOCAL RUNTIME AND DOCKER RUNTIME
================================================================================
Local run:
1. Install Python 3.10+
2. `pip install flask requests`
3. Ensure `data/` files exist
4. `python3 web_app.py`
5. Open http://localhost:5000

Docker run:
1. `docker build -t cnvision:latest .`
2. `docker run -it --rm -p 5000:5000 cnvision:latest`
3. Open http://localhost:5000

Optional mounted data:
- `docker run -it --rm -v $(pwd)/data:/app/data -p 5000:5000 cnvision:latest`


14) PRACTICAL VALIDATION STRATEGY
================================================================================
For robust validation:
- Validate a mixed panel of:
  - whole-exon deletions
  - partial-exon deletions
  - intronic/no-CDS cases
  - boundary-sensitive 1-bp-shift pairs
- Always record transcript and genome build used.
- Compare consequences only when transcript/build are matched.


15) CURRENT SCOPE AND FUTURE IMPROVEMENTS
================================================================================
Current scope:
- Web-based frame consequence interpretation with CDS-aware overlap logic.

Potential future extensions:
- Transcript override in genomic mode for non-MANE validation alignment
- Batch input processing and downloadable reports
- Additional interpretation annotations (domain/context overlays)
- Optional API interface for external integration


================================================================================
END OF DOCUMENT
================================================================================
