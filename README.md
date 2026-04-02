# RAI1 ASO Target Prioritization Portal

**Live site: [sms-aso-targets.web.app](https://sms-aso-targets.web.app)**

A computational evidence map for antisense oligonucleotide (ASO) target selection in Smith-Magenis Syndrome (SMS). This portal evaluates all computationally tractable ASO strategies for upregulating RAI1 expression from the remaining allele in 17p11.2 deletion SMS.

## What this does

The pipeline queries public genomic databases (Ensembl, UCSC, GWIPS-viz, HPA, GTEx, TargetScan) and scores candidate therapeutic targets across 8 mechanisms:

| Mechanism | Targets found | Verdict |
|-----------|--------------|---------|
| AUG uORF blocking | 5 candidates | Primary strategy |
| Near-cognate uORF blocking | 7 validated | Complementary |
| 5'UTR structure disruption | 5 hotspots | Complementary |
| 3'UTR structure disruption | 5 hotspots | Secondary |
| miRNA site blocking | 6 families | Weak |
| Splice switching / NMD | 0 poison exons | Ruled out |
| G-quadruplex resolution | 0 motifs | Ruled out |
| Natural antisense (RAI1-AS1) | 0.5/4.5 score | Ruled out |

All targets are ranked on a unified 0-10 composite score incorporating canonical transcript membership, Ribo-seq evidence, start codon strength, evolutionary conservation, and ASO designability.

### Ranked Targets

| Rank | Target | Type | Score | Canonical | Ribo-seq | Conservation | Clean ASOs |
|------|--------|------|-------|-----------|----------|-------------|------------|
| 1 | ENST00000353383:262 | ATG | 6.4 | Yes (NM_030665.4) | sum=976 (68% of CDS) | phyloP 0.32 | 0/15 |
| 2 | CTG@358 | near-cognate | 6.3 | canonical UTR | 2.57x enrichment | phyloP 1.58 | 11/15 |
| 3 | CTG@344 | near-cognate | 6.3 | canonical UTR | 2.13x enrichment | phyloP 1.45 | 12/15 |
| 4 | GTG@114 | near-cognate | 5.8 | canonical UTR | 3.57x enrichment | phyloP 0.29 | 15/15 |
| 5 | ATC@342 | near-cognate | 5.8 | canonical UTR | 2.62x enrichment | phyloP 1.49 | 12/15 |
| 6 | ACG@110 | near-cognate | 5.4 | canonical UTR | 3.06x enrichment | phyloP 0.28 | 15/15 |
| 7 | CTG@96 | near-cognate | 5.2 | canonical UTR | 2.13x enrichment | phyloP 0.39 | 8/11 |
| 8 | ATT@2 | near-cognate | 4.2 | canonical UTR | 2.60x enrichment | phyloP 0.65 | 0/13 |
| 9 | ENST00000918590:144 | ATG | 3.2 | non-canonical | sum=29 | phyloP 2.82 | 10/16 |
| 10 | ENST00000471135:217 | ATG | 2.3 | non-canonical | sum=34 | phyloP 0.08 | 12/15 |
| 11 | ENST00000395774:18 | ATG | 1.5 | non-canonical | sum=19 | phyloP 0.72 | 11/13 |
| 12 | ENST00000918590:89 | ATG | 1.4 | non-canonical | sum=11 | phyloP 1.85 | 0/11 |

Key finding: the canonical AUG uORF (#1) has the strongest biological evidence but zero clean ASO candidates due to extreme GC content (78%). The near-cognate targets (#2-7) in the same canonical UTR offer better ASO designability and strong conservation.

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
