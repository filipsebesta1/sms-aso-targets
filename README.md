# RAI1 ASO Target Prioritization Portal

**Live site: [sms-aso-targets.web.app](https://sms-aso-targets.web.app)**

A computational evidence map for antisense oligonucleotide (ASO) target selection in Smith-Magenis Syndrome (SMS). This portal evaluates all computationally tractable ASO strategies for upregulating RAI1 expression from the remaining allele in 17p11.2 deletion SMS.

## What this does

The pipeline queries public genomic databases (Ensembl, UCSC, GWIPS-viz, HPA, GTEx, TargetScan) and evaluates candidate therapeutic targets across 8 mechanisms:

- uORF blocking (AUG and near-cognate start codons)
- 5'UTR and 3'UTR secondary structure disruption
- miRNA binding site blocking
- Splice switching / poison exon inclusion
- G-quadruplex resolution
- Natural antisense transcript (NAT) targeting

All candidates are ranked on a unified 0-10 composite score incorporating canonical transcript membership, Ribo-seq evidence, start codon strength, evolutionary conservation, and ASO designability. See the [live dashboard](https://sms-aso-targets.web.app) for current results.

## This is a research planning tool, not a drug design tool

All evidence is computational. Experimental validation (reporter assays, ASO walks, off-target analysis) is required before any therapeutic use. See the [methodology page](https://sms-aso-targets.web.app/methodology.html) for scoring details and [limitations](https://sms-aso-targets.web.app/methodology.html#limitations).

## Running the pipeline

```bash
# Set up
python3 -m venv .venv
source .venv/bin/activate
pip install requests pandas pyBigWig jinja2 ViennaRNA

# Run the data pipeline (queries public APIs, takes ~2 minutes)
python src/pipeline.py

# Build the static site
python src/build_site.py

# Output in site/
open site/index.html
```

## Project structure

```
src/
  pipeline.py      - Data pipeline (16 steps, queries Ensembl/UCSC/GTEx/HPA/GWIPS-viz)
  build_site.py    - Static site generator (TSV -> JSON -> HTML via Jinja2)
templates/         - Jinja2 HTML templates
site/              - Generated static site (deployed to Firebase)
data/
  derived/         - Pipeline output TSV files (18 files)
  raw/             - Raw API responses (Ensembl JSON, TargetScan)
reports/           - Decision memo
```

## Data sources

| Source | URL | Data used |
|--------|-----|-----------|
| Ensembl REST API | rest.ensembl.org | Gene/transcript coordinates, exon structure |
| UCSC phyloP100way | hgdownload.cse.ucsc.edu | Conservation scores |
| GWIPS-viz (via UCSC) | genome.ucsc.edu | Aggregate ribosome profiling |
| Human Protein Atlas | proteinatlas.org | Tissue expression |
| GTEx Portal | gtexportal.org | eQTL/sQTL absence |
| TargetScan 8.0 | targetscan.org | miRNA binding site predictions |
| ViennaRNA | tbi.univie.ac.at/RNA | RNA secondary structure (MFE) |
| ENCODE (via UCSC) | encodeproject.org | TF/RBP binding sites |

All coordinates are GRCh38/hg38.

## License

This project is released under the [MIT License](LICENSE). You are free to use, modify, and distribute this software for any purpose.

## Citation

If you use data or code from this project, please cite:

```
SMS ASO Target Prioritization Portal.
Available at: https://sms-aso-targets.web.app
Source code: https://github.com/filipsebesta1/sms-aso-targets
```
