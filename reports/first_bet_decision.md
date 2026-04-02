    # RAI1 Upregulation: First-Bet Decision Memo

    **Date:** 2026-03-30
    **Gene:** RAI1 (ENSG00000108557), chr17p11.2
    **Context:** Classic 17p11.2 deletion Smith–Magenis Syndrome — goal is to upregulate the remaining RAI1 allele.

    ---

    ## Executive Summary

    **Recommended first bet: 5′UTR / uORF lever**

    uORF candidates present with conservation and/or strong Kozak context

    ---

    ## Evidence Scorecard

    | Lever | Score | Key Evidence |
    |-------|-------|-------------|
    | **A) 5′UTR / uORF** | 6.5/7 | 5 uORFs found; multiple uORFs suggest regulatory complexity; strong Kozak context on top uORFs; moderate canonical uORF conservation (phyloP=0.32); zero eQTLs supports post-transcriptional bottleneck; Ribo-seq confirms ribosome occupancy at canonical uORF (sum=976 vs CDS=1436) |
    | **B) NAT (RAI1-AS1)** | 0.5/5.5 | overlap quality score=0.5 (large_overlap+0.5) |
    | **C) Splicing/NMD** | — | Sanity check only; no sQTL signal in GTEx |

    ---

    ## A) 5′UTR / uORF Lever — Detailed Evidence

    ### Transcript Inventory
    - RAI1 protein-coding transcripts queried from Ensembl REST.
    - Canonical transcript: **ENST00000353383**
    - See: `data/derived/rai1_transcripts.tsv`

    ### uORF Candidates
    - **5 unique uORFs** identified across supported transcripts.
    - Top candidates (by Kozak score and length):

    | Transcript | uAUG pos | Stop pos | AA length | Kozak | Notes |
    |-----------|----------|----------|-----------|-------|-------|
    | ENST00000918590 | 144.0 | 177.0 | 9.0 | 3.0 |  |
| ENST00000471135 | 217.0 | 375.0 | 51.0 | 2.0 | overlaps_mainAUG_region |
| ENST00000395774 | 18.0 | 100.0 | 26.0 | 2.0 | overlaps_mainAUG_region |

    - Full data: `data/derived/rai1_utr_candidates.tsv`

    ### Triage: Canonical Transcript Membership
    Each uORF is classified by whether it falls in the canonical/RefSeq-supported transcript or only in less-supported isoforms.

    | Candidate | Verdict | Source Quality | uAUG Genomic | RefSeq |
    |-----------|---------|---------------|-------------|--------|
    | ENST00000353383:262 | canonical | high (ensembl_havana) | 17681720 | NM_030665.4 |
| ENST00000471135:217 | non-canonical_only | medium (havana) | 17683694 |  |
| ENST00000395774:18 | non-canonical_only | medium (havana) | 17782150 |  |
| ENST00000918590:89 | non-canonical_only | low (havana_tagene) | 17723811 |  |
| ENST00000918590:144 | non-canonical_only | low (havana_tagene) | 17723866 |  |

    - **canonical** = uORF discovered directly in the canonical transcript (MANE/RefSeq).
    - **shared_with_canonical** = uORF is in a non-canonical isoform, but the uAUG genomic position falls within canonical 5′UTR exons (i.e. the uORF is present in the canonical mRNA).
    - **non-canonical_only** = uORF exists only in alternative isoforms; not therapeutically actionable unless isoform is expressed.
    - Full data: `data/derived/rai1_uorf_triage.tsv`

    ### Ribo-seq Reality Check (GWIPS-viz)
    Ribosome profiling signal queried from GWIPS-viz aggregate track (hosted on UCSC as gwipsvizRiboseq bigWig). Values represent aggregate ribosome footprint counts in a ±50 bp window around each uAUG.

    | Candidate | Sum | Max | Mean | Verdict |
    |-----------|-----|-----|------|---------|
    | ENST00000353383:262 | 976 | 115 | 12.0 | canonical |
| ENST00000471135:217 | 34 | 4 | 1.8 | non-canonical_only |
| ENST00000395774:18 | 19 | 4 | 1.7 | non-canonical_only |
| ENST00000918590:89 | 11 | 4 | 1.8 | non-canonical_only |
| ENST00000918590:144 | 29 | 17 | 3.2 | non-canonical_only |
| CDS_start_control | 1436 | 151 | 15.1 | positive_control |
| intron_control | 0 | 0 | 0.0 | negative_control |

    - **Interpretation:** CDS_start_control validates signal presence; intron_control confirms specificity. uORF candidates with sum >> intron background have detectable ribosome occupancy.
    - Full data: `data/derived/rai1_riboseq_check.tsv`
    - Source: [GWIPS-viz Ribo-seq (via UCSC)](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=gwipsvizRiboseq)

    ### Conservation
    Conservation scores retrieved via UCSC phyloP100way (hg38) bigWig.

    | Candidate | Region | Coordinates | Mean phyloP |
    |----------|--------|-------------|-------------|
    | ENST00000918590:144 | uAUG_window | chr17:17723851-17723881 | 2.8244 |
| ENST00000918590:144 | uORF_body | chr17:17723866-17723899 | 2.3714 |
| ENST00000918590:144 | main_start_neighborhood | chr17:17792933-17792949 | 2.7872 |
| ENST00000471135:217 | uAUG_window | chr17:17683679-17683709 | 0.0778 |
| ENST00000471135:217 | uORF_body | chr17:17683694-17683852 | -0.1129 |
| ENST00000471135:217 | main_start_neighborhood | chr17:17792933-17792949 | 2.7872 |
| ENST00000395774:18 | uAUG_window | chr17:17782135-17782165 | 0.7199 |
| ENST00000395774:18 | uORF_body | chr17:17782150-17782232 | 1.0992 |
| ENST00000395774:18 | main_start_neighborhood | chr17:17792933-17792949 | 2.7872 |
| ENST00000353383:262 | uAUG_window | chr17:17681705-17681735 | 0.3202 |
| ENST00000353383:262 | uORF_body | chr17:17681720-17681780 | 0.3409 |
| ENST00000353383:262 | main_start_neighborhood | chr17:17792933-17792949 | 2.7872 |
| ENST00000918590:89 | uAUG_window | chr17:17723796-17723826 | 1.8545 |
| ENST00000918590:89 | uORF_body | chr17:17723811-17723835 | 2.167 |
| ENST00000918590:89 | main_start_neighborhood | chr17:17792933-17792949 | 2.7872 |

    - Full data: `data/derived/rai1_utr_conservation.tsv`
    - Source: [UCSC phyloP100way](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/)

    ---

    ## B) NAT Lever (RAI1-AS1) — Detailed Evidence

    ### Overlap Analysis
    - **RAI1:** chr17:17681458-17811453 (strand 1)
    - **RAI1-AS1:** chr17:17758962-17774684 (strand -1)
    - **Overlap:** 15,722 bp, category: fully_inside_gene_body
    - **Promoter overlap:** NO
    - **Exon overlap:** NO
    - **5′UTR exon overlap:** NO
    - **Quality score:** 0.5/4.5 (large_overlap+0.5)
    - Full data: `data/derived/rai1_as1_overlap.tsv`

    ### Expression in Brain
    **RAI1 brain expression:** Cerebellum: 18.2 nTPM; Cerebral cortex: 14.2 nTPM; Hypothalamus: 9.2 nTPM

    **RAI1-AS1:** expression data not programmatically retrieved. Check [GTEx Portal](https://gtexportal.org/) and [HPA](https://www.proteinatlas.org/) manually.

    - Full data: `data/derived/rai1_as1_expression.tsv`
    - Source: [Human Protein Atlas](https://www.proteinatlas.org/ENSG00000108557-RAI1/tissue)

    ---

    ## eQTL / sQTL Context

    **RAI1 has zero significant eQTLs and zero sQTLs** across all GTEx tissues (confirmed via API v8 and v10; verified with SORT1 as positive control).

    **Interpretation:** RAI1 expression is tightly controlled and not naturally variable via common genetic variants. This makes post-transcriptional levers (uORF disruption) more attractive, since the bottleneck may be at translation rather than transcription. However, the lack of natural tunability may reflect selection pressure against RAI1 dosage changes — upregulation needs careful titration.
    - Full data: `data/derived/rai1_eqtl_sqtl_summary.tsv`
    - Source: [GTEx Portal API](https://gtexportal.org/api/v2/)

    ---

    ## Recommendation & Rationale

    **First bet: 5′UTR / uORF lever**

    The 5′UTR of RAI1 contains 5 uORFs that could serve as translation-brake targets. ASO-mediated disruption of uORF start codons or their Kozak contexts could increase translation from the main ORF without altering mRNA levels — a post-transcriptional upregulation strategy with precedent (e.g., PCSK9 uORF targeting).
    Key supporting evidence:
    - **Conservation:** uORF regions show phyloP conservation (mean uAUG window=1.16, body=1.17), suggesting functional constraint.
    - **eQTL absence:** RAI1 has zero eQTLs across all GTEx tissues, indicating tight transcriptional control. This makes post-transcriptional levers more attractive — the bottleneck may be at translation.
    - **Ribo-seq:** GWIPS-viz aggregate track shows ribosome occupancy at the canonical uORF (sum=976), confirming active translation in the 5′UTR.
    - **NAT weakness:** RAI1-AS1 sits entirely within gene body introns (no promoter/exon overlap, quality score 0.5/4.5), making it a weak NAT target.



    ### Top Target Windows (if UTR lever)
    Prioritized by: (1) canonical transcript membership, (2) Ribo-seq evidence, (3) Kozak score, (4) conservation.

1. **uORF at UTR position 262** (ENST00000353383, Kozak=1, Ribo-seq sum=976) **[CANONICAL — primary target]**: ASO window spanning uAUG ± 15 nt

2. **uORF at UTR position 144** (ENST00000918590, Kozak=3, Ribo-seq sum=29) ⚠ [non-canonical isoform only]: ASO window spanning uAUG ± 15 nt

3. **uORF at UTR position 217** (ENST00000471135, Kozak=2, Ribo-seq sum=34) ⚠ [non-canonical isoform only]: ASO window spanning uAUG ± 15 nt

4. **uORF at UTR position 18** (ENST00000395774, Kozak=2, Ribo-seq sum=19) ⚠ [non-canonical isoform only]: ASO window spanning uAUG ± 15 nt

---

## Next Cheapest Experiments

1. **Dual-luciferase reporter assay** — Clone the RAI1 5′UTR upstream of firefly luciferase. Mutate each uAUG (ATG→AAG) individually and in combination. Measure translation fold-change in HEK293 or neuronal cell line (SH-SY5Y). Cost: ~2 weeks, ~$2K reagents.

2. **ASO walk on top uORFs** — Design 3–5 LNA gapmers (16–20 nt) targeting each top uAUG Kozak region. Transfect into SH-SY5Y, measure RAI1 protein by Western blot at 48h. Cost: ~$3K for oligos + 2 weeks.

3. **Polysome profiling (if reporter confirms)** — Validate that uORF disruption shifts RAI1 mRNA from light to heavy polysomes, confirming translational upregulation mechanism.

---

## Data Sources

| Source | URL |
|--------|-----|
| Ensembl REST API (GRCh38) | https://rest.ensembl.org |
| UCSC phyloP100way | https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/ |
| GTEx Portal API | https://gtexportal.org/api/v2/ |
| Human Protein Atlas | https://www.proteinatlas.org/ |
| GWIPS-viz Ribo-seq (via UCSC) | https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=gwipsvizRiboseq |
| UCSC Genome Browser | https://genome.ucsc.edu |

---

*Generated by automated pipeline: `src/pipeline.py`*
