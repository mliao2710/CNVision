# CNVision

CNVision is a web-based bioinformatics tool for interpreting intragenic copy
number variants (CNVs). It predicts whether coding impact is **in-frame**,
**frameshift/out-of-frame**, or **no CDS impact** using exon/CDS overlap and
reading-frame logic.

This README is an operational + technical overview for users and reviewers.
For a beginner-heavy deep dive, see `explanation.txt`.

## What CNVision Solves

CNV reports often provide coordinates, but not immediate frame interpretation.
CNVision automates the repeated interpretation steps:

- parse coordinate/cDNA inputs
- resolve build and gene/transcript context
- map CNV interval to exon/CDS overlap
- compute coding-length impact
- classify in-frame vs frameshift

## Scope

CNVision is optimized for frame-consequence interpretation, not full clinical
classification.

### What it provides

- In-frame deletion/duplication calls
- Frameshift deletion/duplication calls
- No-CDS-impact output when overlap is outside coding sequence
- Exon coverage context (whole vs partial overlap)

### What it does not fully provide

- Full ACMG/clinical pathogenicity assertion
- Protein-structure impact modeling
- Comprehensive splice/regulatory interpretation
- Transcript-agnostic concordance across all isoforms by default

## Input Modes

### 1) Exon Number Mode

Use when you already know affected exon numbers.

Inputs:
- Gene symbol or transcript ID
- First exon
- Last exon
- CNV type

### 2) Transcript/cDNA Mode

Use HGVS-like cDNA notation.

Examples:
- `FBN1:c.5065_5546del`
- `NM_001406716.1:c.5065_5546del`

### 3) Genomic Coordinate Mode

Use raw genomic coordinates / ISCN-like strings.

Parser is tolerant to varied formatting and supports:
- `hg19` -> `GRCh37`, `hg38` -> `GRCh38`
- copy-number hints (`x1` deletion, `x3` duplication)
- prefixes like `arr`, `seq`, `sseq`
- mixed separators (commas, underscores, parentheses, hyphen ranges)

Examples:
- `DEL chr15:48,741,090-48,756,096`
- `arr[GRCh37]15q21.1(48,741,090_48,756,096)x1`
- `seq[GRCh37]19p13.2p13.2(13335371_13346648)x1`

## Core Architecture

Pipeline:

`Input -> Parse -> Resolve -> Map -> Predict -> Render`

Detailed stages:
1. Parse/normalize input string(s)
2. Resolve build + gene/transcript context
3. Map genomic interval to transcript exon/CDS overlap
4. Sum coding overlap length
5. Apply modulo-3 rule for frame consequence
6. Render result summary + supporting details in UI

## Project Structure (Runtime Core)

- `web_app.py`
  - Flask entry point and orchestration layer
  - mode routing, parser dispatch, result assembly
- `mane_loader.py`
  - loads GRCh37/GRCh38 references into in-memory structures
  - attaches CDS segments and CDS lengths to exons
- `coordinate_mapper.py`
  - overlap engine: CNV -> exon hits
  - computes coding overlap + partial/full exon flags
- `functional_predictor.py`
  - consequence classifier (in-frame/frameshift/no-CDS)
- `iscn_parser.py`
  - ISCN helper parsing + genes-in-region lookup
- `templates/index.html`
  - UI rendering, result tables, coverage note display
- `static/logo.png`
  - UI branding asset
- `data/`
  - reference datasets

## Required Data Files

Place these files in `data/`:

- `MANE.GRCh38.v1.4.refseq_genomic.gff.gz`
- `GRCh37.ref_seq_select.gz`
- `gene2refseq`

## Literal Trace (One Request)

Example request:
- Mode: genomic coordinates
- Input: `del hg38:chr17:56,594,272-56,594,337`
- Gene hint: `NOG`

Execution path:
1. `web_app.index()` reads form values.
2. `process_genomic_mode(...)` handles coordinate-mode workflow.
3. `parse_genomic_input(...)` extracts chromosome/start/end/build/type hints.
4. `get_mane_data(...)` selects GRCh38 annotation set.
5. `coordinate_mapper.map_cnv_to_exons(...)` computes exon/CDS overlap.
6. `functional_predictor.predict_cnv_effect(...)` computes frame consequence.
7. `templates/index.html` renders result + evidence context.

Typical parsed object:
```json
{
  "chromosome": "17",
  "start": 56594272,
  "end": 56594337,
  "build": "GRCh38",
  "cnv_type_hint": "deletion"
}
```

Typical consequence:
- `in-frame deletion (1 exon)` when total coding overlap is divisible by 3.

## Output Interpretation

Each result card provides:
- gene + transcript used
- consequence label
- exon table (with coding lengths)
- concise exon range
- total coding length and amino-acid check (`length / 3`)
- overlap context (whole/partial exon)

Important interpretation note:
- Exon-number mode assumes selected full exons.
- Genomic/cDNA inputs can produce partial exon overlap.
- So exon-mode and genomic-mode can disagree if overlap boundaries differ.

## Build and Transcript Caveats

### Build mismatch
Using the wrong build is a common source of incorrect mapping.

### Transcript mismatch
Different transcripts of the same gene can change:
- exon numbering
- exon boundaries
- coding-length totals
- final frame consequence

If comparing to an external validation table, transcript/build must match.

## Local Run

1. Install Python 3.10+
2. Install dependencies:
   ```bash
   pip install flask requests
   ```
3. Ensure `data/` contains required files
4. Start app:
   ```bash
   python3 web_app.py
   ```
5. Open:
   `http://localhost:5000`

## Docker Run

Build image:
```bash
docker build -t cnvision:latest .
```

Run container:
```bash
docker run -it --rm -p 5000:5000 cnvision:latest
```

Optional local data mount:
```bash
docker run -it --rm -v $(pwd)/data:/app/data -p 5000:5000 cnvision:latest
```

## Validation Guidance

For robust checks, include mixed cases:
- whole-exon events
- partial-exon events
- intronic/no-CDS events
- 1-bp boundary-shift pairs

Always record:
- genome build used
- transcript used
- expected consequence definition source

## Current Status

Current implementation is a focused, CDS-aware CNV frame interpreter with a
web-first workflow and flexible input parsing.

For full conceptual and educational explanation, see `explanation.txt`.
