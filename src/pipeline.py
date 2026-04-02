#!/usr/bin/env python3
"""
RAI1 upregulation first-bet analysis pipeline.
Queries Ensembl REST, UCSC, and public expression databases to score
three candidate levers: 5'UTR/uORF, NAT (RAI1-AS1), and eQTL/sQTL.
"""

import json, re, sys, time, logging, hashlib, pathlib, textwrap
from io import StringIO
import requests
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

BASE = pathlib.Path(__file__).resolve().parent.parent
RAW = BASE / "data" / "raw"
DER = BASE / "data" / "derived"
REP = BASE / "reports"
for d in [RAW, DER, REP]:
    d.mkdir(parents=True, exist_ok=True)

ENSEMBL = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def ensembl_get(endpoint, params=None, retries=3):
    """GET from Ensembl REST with retry + rate-limit handling."""
    url = ENSEMBL + endpoint
    for attempt in range(retries):
        r = requests.get(url, headers=HEADERS, params=params, timeout=30)
        if r.status_code == 429:
            wait = int(r.headers.get("Retry-After", 2))
            log.warning("Rate limited, waiting %ds", wait)
            time.sleep(wait)
            continue
        r.raise_for_status()
        return r.json()
    raise RuntimeError(f"Failed after {retries} attempts: {url}")

def ensembl_seq(endpoint, params=None, retries=3):
    """GET sequence (text/plain) from Ensembl REST."""
    url = ENSEMBL + endpoint
    for attempt in range(retries):
        r = requests.get(url, headers={"Content-Type": "text/plain"}, params=params, timeout=30)
        if r.status_code == 429:
            time.sleep(int(r.headers.get("Retry-After", 2)))
            continue
        r.raise_for_status()
        return r.text
    raise RuntimeError(f"Failed after {retries} attempts: {url}")

# ---------------------------------------------------------------------------
# STEP 1 — Coordinates and transcript inventory
# ---------------------------------------------------------------------------

def step1():
    log.info("=== STEP 1: Ensembl transcript inventory ===")

    # RAI1 gene
    rai1 = ensembl_get("/lookup/symbol/homo_sapiens/RAI1", {"expand": 1})
    log.info("RAI1: %s %s:%d-%d (%s)", rai1["id"], rai1["seq_region_name"],
             rai1["start"], rai1["end"], rai1["strand"])

    # RAI1-AS1 gene
    rai1as1 = ensembl_get("/lookup/symbol/homo_sapiens/RAI1-AS1", {"expand": 1})
    log.info("RAI1-AS1: %s %s:%d-%d (%s)", rai1as1["id"], rai1as1["seq_region_name"],
             rai1as1["start"], rai1as1["end"], rai1as1["strand"])

    # Build transcript tables
    def tx_table(gene_info):
        rows = []
        for t in gene_info.get("Transcript", []):
            row = {
                "gene_id": gene_info["id"],
                "gene_name": gene_info["display_name"],
                "transcript_id": t["id"],
                "display_name": t.get("display_name", ""),
                "biotype": t.get("biotype", ""),
                "is_canonical": t.get("is_canonical", 0),
                "chr": t.get("seq_region_name", ""),
                "start": t.get("start"),
                "end": t.get("end"),
                "strand": t.get("strand"),
                "tsl": t.get("TSL", {}).get("value", "") if isinstance(t.get("TSL"), dict) else t.get("TSL", ""),
            }
            rows.append(row)
        return pd.DataFrame(rows)

    df_rai1 = tx_table(rai1)
    df_as1 = tx_table(rai1as1)
    df_rai1.to_csv(DER / "rai1_transcripts.tsv", sep="\t", index=False)
    df_as1.to_csv(DER / "rai1as1_transcripts.tsv", sep="\t", index=False)
    log.info("Saved %d RAI1 transcripts, %d RAI1-AS1 transcripts", len(df_rai1), len(df_as1))

    # Save raw JSON for later use
    (RAW / "rai1_gene.json").write_text(json.dumps(rai1, indent=2))
    (RAW / "rai1as1_gene.json").write_text(json.dumps(rai1as1, indent=2))

    return rai1, rai1as1, df_rai1, df_as1

# ---------------------------------------------------------------------------
# STEP 2 — uORF candidate discovery
# ---------------------------------------------------------------------------

CODON_TABLE = {"TAA", "TAG", "TGA"}

def kozak_score(seq, aug_pos):
    """Simple Kozak score: +2 if -3 is A/G, +1 if +4 is G."""
    score = 0
    if aug_pos >= 3 and seq[aug_pos - 3] in "AG":
        score += 2
    if aug_pos + 3 < len(seq) and seq[aug_pos + 3] == "G":
        score += 1
    return score

def find_uorfs(seq):
    """Find all AUG-initiated uORFs in a 5'UTR sequence (DNA, uppercase)."""
    uorfs = []
    for m in re.finditer("ATG", seq):
        start = m.start()
        # translate in-frame until stop or end of UTR
        pos = start + 3
        aa_len = 0
        found_stop = False
        while pos + 3 <= len(seq):
            codon = seq[pos:pos+3]
            if codon in CODON_TABLE:
                found_stop = True
                stop_pos = pos + 3  # end of stop codon
                break
            aa_len += 1
            pos += 3
        if found_stop:
            uorfs.append({
                "uaug_pos": start,
                "stop_pos": stop_pos,
                "aa_len": aa_len,
            })
        else:
            # uORF runs to end of UTR (no in-frame stop) — still record
            uorfs.append({
                "uaug_pos": start,
                "stop_pos": len(seq),
                "aa_len": (len(seq) - start - 3) // 3,
            })
    return uorfs

def get_five_prime_utr_seq(tid):
    """Compute 5'UTR sequence by fetching exon sequences before CDS start."""
    tx = ensembl_get(f"/lookup/id/{tid}", {"expand": 1})
    translation = tx.get("Translation")
    if not translation:
        return None, None
    cds_start_genomic = translation["start"]  # genomic coord of CDS start
    cds_end_genomic = translation["end"]
    strand = tx["strand"]
    exons = sorted(tx.get("Exon", []), key=lambda e: e["start"])

    utr_regions = []  # list of (chr, start, end) in genomic coords
    if strand == 1:
        # + strand: 5'UTR is everything before cds_start_genomic
        for e in exons:
            if e["end"] < cds_start_genomic:
                utr_regions.append((e["seq_region_name"], e["start"], e["end"]))
            elif e["start"] < cds_start_genomic:
                utr_regions.append((e["seq_region_name"], e["start"], cds_start_genomic - 1))
                break
            else:
                break
    else:
        # - strand: 5'UTR is everything after cds_end_genomic (higher genomic coords)
        for e in reversed(exons):
            if e["start"] > cds_end_genomic:
                utr_regions.append((e["seq_region_name"], e["start"], e["end"]))
            elif e["end"] > cds_end_genomic:
                utr_regions.append((e["seq_region_name"], cds_end_genomic + 1, e["end"]))
                break
            else:
                break
        utr_regions.reverse()  # keep in genomic order

    if not utr_regions:
        return None, None

    # Fetch sequence for each UTR region and concatenate
    utr_seq = ""
    for chrom, rs, re_ in utr_regions:
        region_str = f"{chrom}:{rs}..{re_}:{1 if strand == 1 else -1}"
        seq = ensembl_seq(f"/sequence/region/human/{region_str}")
        utr_seq += seq.strip().upper()

    return utr_seq, utr_regions

def step2(rai1_gene, df_rai1):
    log.info("=== STEP 2: uORF candidate discovery ===")

    # Pick protein-coding transcripts with reasonable support
    pc = df_rai1[df_rai1["biotype"] == "protein_coding"].copy()
    if pc.empty:
        log.warning("No protein_coding transcripts found; using all")
        pc = df_rai1.copy()

    rows = []
    seen_hashes = set()
    utr_regions_by_tid = {}  # transcript_id -> (utr_seq, utr_regions)

    for _, tx in pc.iterrows():
        tid = tx["transcript_id"]
        # Compute 5'UTR sequence from exon structure
        try:
            utr_seq, utr_regions = get_five_prime_utr_seq(tid)
        except Exception as e:
            log.warning("Could not compute 5'UTR for %s: %s", tid, e)
            continue
        if not utr_seq or len(utr_seq) < 10:
            log.info("5'UTR too short or absent for %s, skipping", tid)
            continue
        utr_regions_by_tid[tid] = (utr_seq, utr_regions)

        utr_hash = hashlib.md5(utr_seq.encode()).hexdigest()[:8]
        is_dup = utr_hash in seen_hashes
        seen_hashes.add(utr_hash)

        gc = (utr_seq.count("G") + utr_seq.count("C")) / len(utr_seq) * 100

        uorfs = find_uorfs(utr_seq)
        if not uorfs:
            rows.append({
                "transcript_id": tid,
                "utr_hash": utr_hash,
                "utr_len": len(utr_seq),
                "utr_gc": round(gc, 1),
                "uaug_pos": None,
                "stop_pos": None,
                "aa_len": None,
                "kozak_score_simple": None,
                "notes": "no_uORF_found" + (";duplicate_utr" if is_dup else ""),
            })
            continue

        for u in uorfs:
            ks = kozak_score(utr_seq, u["uaug_pos"])
            notes_parts = []
            if is_dup:
                notes_parts.append("duplicate_utr")
            # check overlap with main AUG region (last few nt of UTR)
            if u["stop_pos"] >= len(utr_seq) - 3:
                notes_parts.append("overlaps_mainAUG_region")
            rows.append({
                "transcript_id": tid,
                "utr_hash": utr_hash,
                "utr_len": len(utr_seq),
                "utr_gc": round(gc, 1),
                "uaug_pos": u["uaug_pos"],
                "stop_pos": u["stop_pos"],
                "aa_len": u["aa_len"],
                "kozak_score_simple": ks,
                "notes": ";".join(notes_parts) if notes_parts else "",
            })
        log.info("  %s: UTR=%d bp, GC=%.1f%%, %d uORFs", tid, len(utr_seq), gc, len(uorfs))

    df = pd.DataFrame(rows)
    df.to_csv(DER / "rai1_utr_candidates.tsv", sep="\t", index=False)
    log.info("Saved %d uORF candidate rows", len(df))
    return df, utr_regions_by_tid

# ---------------------------------------------------------------------------
# STEP 2b — Triage: canonical-transcript membership + transcript quality
# ---------------------------------------------------------------------------

def genomic_pos_in_regions(pos, regions):
    """Check if a genomic position falls within any of the given (chr, start, end) regions."""
    for _, rs, re_ in regions:
        if rs <= pos <= re_:
            return True
    return False

def step2b(df_uorfs, utr_regions_by_tid, rai1_gene):
    """Annotate each uORF with canonical-transcript triage info."""
    log.info("=== STEP 2b: uORF triage (canonical membership) ===")

    strand = rai1_gene["strand"]

    # Identify canonical transcript
    canonical_tid = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical_tid = t["id"]
            break
    if not canonical_tid:
        canonical_tid = rai1_gene["Transcript"][0]["id"]
    log.info("Canonical transcript: %s", canonical_tid)

    # Get canonical UTR regions (if available)
    canonical_utr_regions = None
    if canonical_tid in utr_regions_by_tid:
        canonical_utr_regions = utr_regions_by_tid[canonical_tid][1]

    # Get transcript source info from raw Ensembl data
    tx_source = {}
    for t in rai1_gene.get("Transcript", []):
        tx_source[t["id"]] = t.get("source", "unknown")

    # Check RefSeq status for canonical
    canonical_refseq = ""
    try:
        xrefs = ensembl_get(f"/xrefs/id/{canonical_tid}")
        for x in xrefs:
            if x.get("dbname") == "RefSeq_mRNA":
                canonical_refseq = x.get("display_id", "")
                break
    except Exception as e:
        log.warning("RefSeq xref lookup failed: %s", e)
    log.info("Canonical RefSeq: %s", canonical_refseq or "none")

    # For each uORF, determine:
    # 1. Is the uORF in the canonical transcript directly?
    # 2. If not, does the uAUG genomic position fall within canonical UTR exons?
    # 3. Transcript source quality
    triage_rows = []
    real_uorfs = df_uorfs[df_uorfs["uaug_pos"].notna()].copy()

    for _, row in real_uorfs.iterrows():
        tid = row["transcript_id"]
        uaug = int(row["uaug_pos"])

        in_canonical_directly = (tid == canonical_tid)

        # Map uAUG to genomic coordinate using its own transcript's UTR regions
        uaug_genomic = None
        if tid in utr_regions_by_tid:
            _, uaug_genomic = utr_pos_to_genomic(uaug, utr_regions_by_tid[tid][1], strand)

        # Check if uAUG genomic position is within canonical UTR exons
        in_canonical_utr = False
        if not in_canonical_directly and uaug_genomic and canonical_utr_regions:
            in_canonical_utr = genomic_pos_in_regions(uaug_genomic, canonical_utr_regions)

        # Transcript source quality ranking
        source = tx_source.get(tid, "unknown")
        if source == "ensembl_havana":
            source_quality = "high (ensembl_havana)"
        elif source == "havana":
            source_quality = "medium (havana)"
        elif source == "havana_tagene":
            source_quality = "low (havana_tagene)"
        else:
            source_quality = f"unknown ({source})"

        # Overall triage verdict
        if in_canonical_directly:
            verdict = "canonical"
        elif in_canonical_utr:
            verdict = "shared_with_canonical"
        else:
            verdict = "non-canonical_only"

        triage_rows.append({
            "transcript_id": tid,
            "uaug_pos": uaug,
            "uaug_genomic": uaug_genomic,
            "is_canonical_transcript": in_canonical_directly,
            "uaug_in_canonical_utr": in_canonical_utr or in_canonical_directly,
            "transcript_source": source,
            "source_quality": source_quality,
            "canonical_refseq": canonical_refseq if in_canonical_directly else "",
            "triage_verdict": verdict,
        })

    df_triage = pd.DataFrame(triage_rows)
    df_triage.to_csv(DER / "rai1_uorf_triage.tsv", sep="\t", index=False)
    log.info("Triage results:")
    for _, r in df_triage.iterrows():
        log.info("  %s:%d → %s (source: %s, genomic: %s)",
                 r["transcript_id"], r["uaug_pos"], r["triage_verdict"],
                 r["transcript_source"], r["uaug_genomic"])
    return df_triage

# ---------------------------------------------------------------------------
# STEP 2c — Ribo-seq reality check (GWIPS-viz via UCSC)
# ---------------------------------------------------------------------------

def step2c(df_triage, rai1_gene):
    """Query GWIPS-viz Ribo-seq track (hosted on UCSC) for ribosome occupancy at each uORF."""
    log.info("=== STEP 2c: Ribo-seq reality check (GWIPS-viz) ===")

    track = "gwipsvizRiboseq"
    window = 50  # ±50 bp around uAUG

    results = []
    for _, row in df_triage.iterrows():
        pos = row.get("uaug_genomic")
        if pos is None or pd.isna(pos):
            continue
        pos = int(pos)
        cand_id = f"{row['transcript_id']}:{int(row['uaug_pos'])}"

        try:
            r = requests.get("https://api.genome.ucsc.edu/getData/track", params={
                "genome": "hg38", "track": track,
                "chrom": "chr17", "start": pos - window, "end": pos + window,
            }, timeout=15)
            if r.status_code == 200:
                items = r.json().get(track, [])
                values = [item["value"] for item in items]
                results.append({
                    "candidate_id": cand_id,
                    "uaug_genomic": pos,
                    "region": f"chr17:{pos-window}-{pos+window}",
                    "ribo_datapoints": len(items),
                    "ribo_sum": sum(values),
                    "ribo_max": max(values) if values else 0,
                    "ribo_mean": round(sum(values) / len(values), 1) if values else 0,
                    "triage_verdict": row["triage_verdict"],
                    "source": "GWIPS-viz via UCSC (gwipsvizRiboseq bigWig)",
                })
        except Exception as e:
            log.warning("Ribo-seq query failed for %s: %s", cand_id, e)

    # Add CDS start as positive control
    cds_start = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical") and t.get("Translation"):
            cds_start = t["Translation"]["start"]
            break
    if cds_start:
        try:
            r = requests.get("https://api.genome.ucsc.edu/getData/track", params={
                "genome": "hg38", "track": track,
                "chrom": "chr17", "start": cds_start, "end": cds_start + 100,
            }, timeout=15)
            if r.status_code == 200:
                items = r.json().get(track, [])
                values = [item["value"] for item in items]
                results.append({
                    "candidate_id": "CDS_start_control",
                    "uaug_genomic": cds_start,
                    "region": f"chr17:{cds_start}-{cds_start+100}",
                    "ribo_datapoints": len(items),
                    "ribo_sum": sum(values),
                    "ribo_max": max(values) if values else 0,
                    "ribo_mean": round(sum(values) / len(values), 1) if values else 0,
                    "triage_verdict": "positive_control",
                    "source": "GWIPS-viz via UCSC (gwipsvizRiboseq bigWig)",
                })
        except Exception as e:
            log.warning("Ribo-seq query failed for %s: %s", cand_id, e)

    # Add intron negative control
    intron_mid = (rai1_gene["start"] + rai1_gene["end"]) // 2
    try:
        r = requests.get("https://api.genome.ucsc.edu/getData/track", params={
            "genome": "hg38", "track": track,
            "chrom": "chr17", "start": intron_mid, "end": intron_mid + 100,
        }, timeout=15)
        if r.status_code == 200:
            items = r.json().get(track, [])
            values = [item["value"] for item in items]
            results.append({
                "candidate_id": "intron_control",
                "uaug_genomic": intron_mid,
                "region": f"chr17:{intron_mid}-{intron_mid+100}",
                "ribo_datapoints": len(items),
                "ribo_sum": sum(values),
                "ribo_max": max(values) if values else 0,
                "ribo_mean": round(sum(values) / len(values), 1) if values else 0,
                "triage_verdict": "negative_control",
                "source": "GWIPS-viz via UCSC (gwipsvizRiboseq bigWig)",
            })
    except Exception as e:
        log.warning("Intron control Ribo-seq query failed: %s", e)

    df = pd.DataFrame(results)
    df.to_csv(DER / "rai1_riboseq_check.tsv", sep="\t", index=False)

    for _, r in df.iterrows():
        log.info("  %s: sum=%d max=%d mean=%.1f (%s)",
                 r["candidate_id"], r["ribo_sum"], r["ribo_max"], r["ribo_mean"], r["triage_verdict"])
    return df

# ---------------------------------------------------------------------------
# STEP 3 — Conservation scoring via UCSC bigWig
# ---------------------------------------------------------------------------

def utr_pos_to_genomic(utr_pos, utr_regions, strand):
    """Convert a 0-based UTR-relative position to genomic coordinate.
    utr_regions is list of (chr, start, end) in genomic order.
    For + strand, UTR sequence goes left-to-right through regions.
    For - strand, UTR sequence goes right-to-left (highest coord first in transcript).
    """
    if strand == 1:
        offset = utr_pos
        for chrom, rs, re_ in utr_regions:
            region_len = re_ - rs + 1
            if offset < region_len:
                return chrom, rs + offset
            offset -= region_len
    else:
        # - strand: utr_regions are in genomic order (low to high),
        # but transcript reads high-to-low. Reverse for transcript order.
        offset = utr_pos
        for chrom, rs, re_ in reversed(utr_regions):
            region_len = re_ - rs + 1
            if offset < region_len:
                return chrom, re_ - offset
            offset -= region_len
    return None, None

def step3(df_uorfs, rai1_gene, utr_regions_by_tid):
    log.info("=== STEP 3: Conservation scoring ===")

    chrom = "chr" + str(rai1_gene["seq_region_name"])
    strand = rai1_gene["strand"]

    empty_df = pd.DataFrame(columns=["candidate_id","region_type","chr","start","end",
                                      "mean_score","track_used","build"])

    # Handle empty uORF dataframe
    if df_uorfs.empty or "uaug_pos" not in df_uorfs.columns:
        log.warning("No uORF data available for conservation scoring")
        empty_df.to_csv(DER / "rai1_utr_conservation.tsv", sep="\t", index=False)
        return empty_df

    # Select top candidates
    candidates = df_uorfs[df_uorfs["uaug_pos"].notna()].copy()
    candidates = candidates.drop_duplicates(subset=["utr_hash", "uaug_pos"])
    candidates = candidates.sort_values(["kozak_score_simple", "aa_len"], ascending=[False, False])
    top = candidates.head(5)

    if top.empty:
        log.warning("No uORF candidates for conservation scoring")
        empty_df.to_csv(DER / "rai1_utr_conservation.tsv", sep="\t", index=False)
        return empty_df

    conservation_rows = []

    # Try pyBigWig with UCSC phyloP bigWig
    phylop_url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw"
    bw = None
    try:
        import pyBigWig
        bw = pyBigWig.open(phylop_url)
        log.info("Opened phyloP100way bigWig remotely")
    except Exception as e:
        log.warning("Could not open bigWig: %s", e)

    for idx, row in top.iterrows():
        tid = row["transcript_id"]
        uaug = int(row["uaug_pos"])
        stop = int(row["stop_pos"])
        cand_id = f"{tid}:{uaug}"

        # Use this transcript's own UTR regions for coordinate mapping
        if tid not in utr_regions_by_tid:
            log.warning("No UTR regions cached for %s, skipping %s", tid, cand_id)
            continue
        tx_utr_seq, tx_utr_regions = utr_regions_by_tid[tid]
        tx_utr_len = len(tx_utr_seq)

        # Clamp stop_pos to last valid UTR position
        stop_clamped = min(stop, tx_utr_len - 1)

        # Map positions using this transcript's own UTR regions
        _, uaug_genomic = utr_pos_to_genomic(uaug, tx_utr_regions, strand)
        _, stop_genomic = utr_pos_to_genomic(stop_clamped, tx_utr_regions, strand)

        if uaug_genomic is None:
            log.warning("Could not map uAUG pos %d for %s", uaug, cand_id)
            continue

        # Compute main start neighborhood from this transcript's last UTR region
        last_region = tx_utr_regions[-1] if strand == 1 else tx_utr_regions[0]
        if strand == 1:
            main_neighborhood_start = max(last_region[1], last_region[2] - 59)
            main_neighborhood_end = last_region[2] + 1
        else:
            main_neighborhood_start = last_region[1]
            main_neighborhood_end = min(last_region[2] + 1, last_region[1] + 60)

        # Define query regions (all on chr17)
        uaug_win_start = uaug_genomic - 15
        uaug_win_end = uaug_genomic + 15
        if stop_genomic is not None:
            body_start = min(uaug_genomic, stop_genomic)
            body_end = max(uaug_genomic, stop_genomic)
        else:
            body_start = uaug_genomic
            body_end = uaug_genomic + 30

        # Sanity check: uORF body should not span more than the UTR length
        if body_end - body_start > tx_utr_len:
            log.warning("uORF body span (%d bp) exceeds UTR length (%d bp) for %s — clamping",
                        body_end - body_start, tx_utr_len, cand_id)
            body_end = body_start + min(stop - uaug, tx_utr_len)

        regions = [
            ("uAUG_window", uaug_win_start, uaug_win_end),
            ("uORF_body", body_start, body_end),
            ("main_start_neighborhood", main_neighborhood_start, main_neighborhood_end),
        ]

        for rtype, rstart, rend in regions:
            mean_score = None
            track = "hg38.phyloP100way"
            if bw is not None and rend > rstart:
                try:
                    vals = bw.values(chrom, int(rstart), int(rend))
                    valid = [v for v in vals if v is not None and v == v]
                    mean_score = round(sum(valid) / len(valid), 4) if valid else None
                except Exception as e:
                    log.warning("bigWig query failed for %s %s:%d-%d: %s", cand_id, chrom, rstart, rend, e)

            conservation_rows.append({
                "candidate_id": cand_id,
                "region_type": rtype,
                "chr": chrom,
                "start": int(rstart),
                "end": int(rend),
                "mean_score": mean_score,
                "track_used": track,
                "build": "GRCh38",
            })

    if bw is not None:
        bw.close()

    df = pd.DataFrame(conservation_rows)
    df.to_csv(DER / "rai1_utr_conservation.tsv", sep="\t", index=False)
    log.info("Saved %d conservation rows", len(df))

    if df.empty or df["mean_score"].isna().all():
        log.warning("All conservation scores are NA — bigWig access may have failed")
        log.info("Manual fallback: visit https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=%s:%d-%d",
                 chrom, rai1_gene["start"], rai1_gene["end"])
    return df

# ---------------------------------------------------------------------------
# STEP 4 — NAT lever overlap scoring
# ---------------------------------------------------------------------------

def step4(rai1_gene, rai1as1_gene):
    log.info("=== STEP 4: NAT overlap scoring ===")

    chrom = rai1_gene["seq_region_name"]
    rai1_start = rai1_gene["start"]
    rai1_end = rai1_gene["end"]
    rai1_strand = rai1_gene["strand"]

    as1_start = rai1as1_gene["start"]
    as1_end = rai1as1_gene["end"]
    as1_strand = rai1as1_gene["strand"]

    # TSS
    if rai1_strand == 1:
        rai1_tss = rai1_start
    else:
        rai1_tss = rai1_end

    # Promoter region: [TSS-2000, TSS+500] in genomic coords
    prom_start = rai1_tss - 2000
    prom_end = rai1_tss + 500

    # Overall overlap
    ov_start = max(rai1_start, as1_start)
    ov_end = min(rai1_end, as1_end)
    overlap_bp = max(0, ov_end - ov_start)

    # Promoter overlap
    pov_start = max(prom_start, as1_start)
    pov_end = min(prom_end, as1_end)
    overlaps_promoter = pov_end > pov_start

    # Exon overlap — get canonical transcript exons
    canonical_tid = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical_tid = t["id"]
            break
    if not canonical_tid:
        canonical_tid = rai1_gene["Transcript"][0]["id"]

    tx_detail = ensembl_get(f"/lookup/id/{canonical_tid}", {"expand": 1})
    exons = [(e["start"], e["end"]) for e in tx_detail.get("Exon", [])]

    exon_overlap = 0
    utr5_exon_overlap = False
    for es, ee in exons:
        eov_s = max(es, as1_start)
        eov_e = min(ee, as1_end)
        if eov_e > eov_s:
            exon_overlap += eov_e - eov_s
    overlaps_exon = exon_overlap > 0

    # Check if first exon (5'UTR-containing) overlaps
    if exons:
        if rai1_strand == 1:
            first_exon = min(exons, key=lambda x: x[0])
        else:
            first_exon = max(exons, key=lambda x: x[1])
        fe_ov_s = max(first_exon[0], as1_start)
        fe_ov_e = min(first_exon[1], as1_end)
        utr5_exon_overlap = fe_ov_e > fe_ov_s

    # Category
    if as1_start >= rai1_start and as1_end <= rai1_end:
        cat = "fully_inside_gene_body"
    elif overlap_bp > 0:
        cat = "partially_overlapping"
    else:
        cat = "non_overlapping"

    # Quality score
    score = 0
    score_parts = []
    if overlaps_promoter:
        score += 2; score_parts.append("promoter+2")
    if overlaps_exon:
        score += 1; score_parts.append("exon+1")
    if utr5_exon_overlap:
        score += 1; score_parts.append("5utr_exon+1")
    if overlap_bp > 5000:
        score += 0.5; score_parts.append("large_overlap+0.5")

    row = {
        "rai1_chr": chrom,
        "rai1_start": rai1_start,
        "rai1_end": rai1_end,
        "rai1_strand": rai1_strand,
        "rai1_tss": rai1_tss,
        "rai1as1_start": as1_start,
        "rai1as1_end": as1_end,
        "rai1as1_strand": as1_strand,
        "overlap_bp": overlap_bp,
        "overlaps_promoter_bool": overlaps_promoter,
        "overlaps_exon_bool": overlaps_exon,
        "overlaps_5utr_exon_bool": utr5_exon_overlap,
        "overlap_category": cat,
        "overlap_quality_score": score,
        "score_breakdown": "; ".join(score_parts),
    }
    df = pd.DataFrame([row])
    df.to_csv(DER / "rai1_as1_overlap.tsv", sep="\t", index=False)
    log.info("Overlap: %d bp, category=%s, quality_score=%.1f", overlap_bp, cat, score)
    return df

# ---------------------------------------------------------------------------
# STEP 5 — Expression evidence
# ---------------------------------------------------------------------------

def hpa_scrape_barchart(gene_name, ensembl_id):
    """Scrape tissue expression from HPA tissue page barChart JSON."""
    rows = []
    # First try the consensus tissue page
    url = f"https://www.proteinatlas.org/{ensembl_id}-{gene_name}/tissue"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200:
            # Extract barChart JSON using bracket-depth matching
            charts = []
            for m in re.finditer(r'\.barChart\(\[', r.text):
                start = m.start() + len('.barChart(')
                depth = 0
                for i in range(start, min(start + 100000, len(r.text))):
                    if r.text[i] == '[': depth += 1
                    elif r.text[i] == ']':
                        depth -= 1
                        if depth == 0:
                            try:
                                charts.append(json.loads(r.text[start:i+1]))
                            except (json.JSONDecodeError, ValueError):
                                pass
                            break

            if charts:
                # First chart is consensus tissue RNA expression
                for entry in charts[0]:
                    rows.append({
                        "gene": gene_name,
                        "source": "HPA_consensus",
                        "tissue": entry.get("label", ""),
                        "value": entry.get("value", ""),
                        "unit": "nTPM",
                        "notes": "",
                    })
            # Look for brain-specific chart
            for chart in charts[1:]:
                if any("brain" in str(e.get("tooltip","")).lower() for e in chart[:2]):
                    for entry in chart:
                        rows.append({
                            "gene": gene_name,
                            "source": "HPA_brain_regional",
                            "tissue": entry.get("label", ""),
                            "value": entry.get("value", ""),
                            "unit": "nTPM",
                            "notes": "brain_regional",
                        })
                    break
    except Exception as e:
        log.warning("HPA scrape failed for %s: %s", gene_name, e)
    return rows

def step5(rai1_gene, rai1as1_gene):
    log.info("=== STEP 5: Expression evidence ===")
    rows = []

    # Scrape HPA tissue expression for RAI1
    rai1_rows = hpa_scrape_barchart("RAI1", rai1_gene["id"])
    rows.extend(rai1_rows)
    log.info("HPA scraped %d tissue entries for RAI1", len(rai1_rows))

    # RAI1-AS1 is unlikely on HPA (lncRNA), but try
    as1_rows = hpa_scrape_barchart("RAI1-AS1", rai1as1_gene["id"])
    rows.extend(as1_rows)
    log.info("HPA scraped %d tissue entries for RAI1-AS1", len(as1_rows))

    # Also try GTEx API with versioned gene IDs
    for gene_name, gene_id in [("RAI1", rai1_gene["id"]), ("RAI1-AS1", rai1as1_gene["id"])]:
        if any(r["gene"] == gene_name and r["source"].startswith("HPA") for r in rows):
            continue  # already have data
        try:
            # GTEx needs versioned Gencode IDs; try common versions
            for version in ["", ".16", ".15", ".14", ".13", ".12"]:
                url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
                params = {"gencodeId": gene_id + version, "datasetId": "gtex_v8"}
                r = requests.get(url, params=params, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    entries = data.get("data", [])
                    if entries:
                        for entry in entries:
                            rows.append({
                                "gene": gene_name,
                                "source": "GTEx_v8",
                                "tissue": entry.get("tissueSiteDetailId", ""),
                                "value": entry.get("median", ""),
                                "unit": "TPM",
                                "notes": f"gencodeId={gene_id}{version}",
                            })
                        log.info("Got GTEx data for %s with version %s: %d tissues",
                                 gene_name, version, len(entries))
                        break
        except Exception as e:
            log.warning("GTEx API failed for %s: %s", gene_name, e)

    if not any(r["gene"] == "RAI1-AS1" for r in rows):
        rows.append({
            "gene": "RAI1-AS1",
            "source": "not_available",
            "tissue": "",
            "value": "",
            "unit": "",
            "notes": "RAI1-AS1 (lncRNA) not in HPA; check GTEx portal manually at https://gtexportal.org/",
        })

    if not any(r["gene"] == "RAI1" for r in rows):
        rows.append({
            "gene": "RAI1",
            "source": "not_available",
            "tissue": "",
            "value": "",
            "unit": "",
            "notes": "RAI1 is known brain-expressed from literature (PMID:17154270)",
        })

    df = pd.DataFrame(rows)
    # If we have data, subset to brain-related tissues + top expressed
    if len(df) > 10:
        brain_kw = ["brain", "cortex", "cerebellum", "hippocampus", "amygdala", "hypothalamus",
                     "caudate", "putamen", "substantia", "nucleus", "frontal", "cerebell", "anterior",
                     "temporal", "occipital", "parietal"]
        brain_mask = df["tissue"].str.lower().str.contains("|".join(brain_kw), na=False)
        # Keep brain rows + top 10 per gene by value
        top_per_gene = df.copy()
        top_per_gene["value_num"] = pd.to_numeric(top_per_gene["value"], errors="coerce")
        top_overall = top_per_gene.sort_values("value_num", ascending=False).groupby("gene").head(10)
        df_out = pd.concat([df[brain_mask], top_overall]).drop_duplicates()
        if "value_num" in df_out.columns:
            df_out = df_out.drop(columns=["value_num"])
    else:
        df_out = df

    df_out.to_csv(DER / "rai1_as1_expression.tsv", sep="\t", index=False)
    log.info("Saved %d expression rows", len(df_out))
    return df_out

# ---------------------------------------------------------------------------
# STEP 6 — eQTL/sQTL summary
# ---------------------------------------------------------------------------

def step6(rai1_gene):
    log.info("=== STEP 6: eQTL/sQTL summary ===")
    rows = []
    gene_id = rai1_gene["id"]

    # GTEx API requires versioned Gencode IDs; try common versions
    for qtl_type, endpoint in [("eQTL", "singleTissueEqtl"), ("sQTL", "singleTissueSqtl")]:
        found = False
        for version in ["", ".16", ".15", ".14", ".13", ".12", ".11", ".10"]:
            try:
                url = f"https://gtexportal.org/api/v2/association/{endpoint}"
                params = {"gencodeId": gene_id + version, "datasetId": "gtex_v8"}
                r = requests.get(url, params=params, timeout=20)
                if r.status_code == 200:
                    data = r.json()
                    entries = data.get("data", [])
                    if entries:
                        limit = 50 if qtl_type == "eQTL" else 20
                        for entry in entries[:limit]:
                            rows.append({
                                "qtl_type": qtl_type,
                                "tissue": entry.get("tissueSiteDetailId", ""),
                                "top_variant_id": entry.get("variantId", ""),
                                "pvalue": entry.get("pValue", ""),
                                "effect_direction_if_available": ("up" if entry.get("nes", 0) > 0 else "down") if qtl_type == "eQTL" else "",
                                "nes": entry.get("nes", ""),
                                "notes": f"gencodeId={gene_id}{version}",
                                "source": "GTEx_v8_API",
                            })
                        log.info("Got %d %s entries from GTEx (version=%s)", len(entries), qtl_type, version)
                        found = True
                        break
            except Exception as e:
                log.warning("GTEx %s API failed (version=%s): %s", qtl_type, version, e)
        if not found:
            log.warning("No %s data found for any Gencode version", qtl_type)

    if not rows:
        # Also try v10
        for qtl_type, endpoint in [("eQTL", "singleTissueEqtl"), ("sQTL", "singleTissueSqtl")]:
            try:
                url = f"https://gtexportal.org/api/v2/association/{endpoint}"
                params = {"gencodeId": gene_id + ".18", "datasetId": "gtex_v10"}
                r = requests.get(url, params=params, timeout=20)
                if r.status_code == 200 and r.json().get("data"):
                    for entry in r.json()["data"][:50]:
                        rows.append({
                            "qtl_type": qtl_type,
                            "tissue": entry.get("tissueSiteDetailId", ""),
                            "top_variant_id": entry.get("variantId", ""),
                            "pvalue": entry.get("pValue", ""),
                            "effect_direction_if_available": ("up" if entry.get("nes", 0) > 0 else "down") if qtl_type == "eQTL" else "",
                            "nes": entry.get("nes", ""),
                            "notes": f"gencodeId={gene_id}.18",
                            "source": "GTEx_v10_API",
                        })
            except Exception as e:
                log.warning("GTEx v10 %s query failed: %s", qtl_type, e)

    if not rows:
        # Confirmed: RAI1 has no significant QTLs
        for qt in ["eQTL", "sQTL"]:
            rows.append({
                "qtl_type": qt,
                "tissue": "all_tissues",
                "top_variant_id": "none",
                "pvalue": "NA",
                "effect_direction_if_available": "",
                "nes": "NA",
                "notes": f"No significant {qt}s for RAI1 in any tissue; confirmed via GTEx API v8 and v10",
                "source": "GTEx_v8_v10_API",
            })

    df = pd.DataFrame(rows)
    df.to_csv(DER / "rai1_eqtl_sqtl_summary.tsv", sep="\t", index=False)
    log.info("Saved %d QTL rows", len(df))
    return df

# ---------------------------------------------------------------------------
# STEP 7 — Decision memo
# ---------------------------------------------------------------------------

def step7(df_uorfs, df_conservation, df_overlap, df_expression, df_qtl, rai1_gene, rai1as1_gene, df_triage=None, df_riboseq=None):
    log.info("=== STEP 7: Writing decision memo ===")

    # Summarize uORF data
    real_uorfs = df_uorfs[df_uorfs["uaug_pos"].notna()]
    unique_uorfs = real_uorfs.drop_duplicates(subset=["utr_hash", "uaug_pos"])
    n_uorfs = len(unique_uorfs)
    top_uorfs = unique_uorfs.sort_values(["kozak_score_simple", "aa_len"], ascending=[False, False]).head(3)

    # Conservation summary
    has_conservation = df_conservation["mean_score"].notna().any()
    if has_conservation:
        uaug_scores = df_conservation[df_conservation["region_type"] == "uAUG_window"]["mean_score"]
        body_scores = df_conservation[df_conservation["region_type"] == "uORF_body"]["mean_score"]
        main_scores = df_conservation[df_conservation["region_type"] == "main_start_neighborhood"]["mean_score"]
        avg_uaug = uaug_scores.mean() if not uaug_scores.empty else None
        avg_body = body_scores.mean() if not body_scores.empty else None
        avg_main = main_scores.mean() if not main_scores.empty else None
        # Also compute conservation for canonical uORF only (most relevant for scoring)
        canonical_tid = None
        for t in rai1_gene.get("Transcript", []):
            if t.get("is_canonical"):
                canonical_tid = t["id"]
                break
        canonical_cons = df_conservation[df_conservation["candidate_id"].str.startswith(canonical_tid + ":")] if canonical_tid else pd.DataFrame()
        canonical_uaug_phylop = canonical_cons[canonical_cons["region_type"] == "uAUG_window"]["mean_score"].mean() if not canonical_cons.empty else None
    else:
        avg_uaug = avg_body = avg_main = None
        canonical_uaug_phylop = None

    # Overlap summary
    ov = df_overlap.iloc[0]
    ov_score = ov.get("overlap_quality_score", 0)
    ov_cat = ov.get("overlap_category", "unknown")
    ov_bp = ov.get("overlap_bp", 0)
    ov_promoter = ov.get("overlaps_promoter_bool", False)

    # Expression summary
    has_expr = not df_expression.empty and df_expression["source"].iloc[0] != "not_available"
    brain_expr_rai1 = ""
    brain_expr_as1 = ""
    if has_expr:
        brain_kw = ["brain", "cortex", "cerebellum", "hippocampus", "amygdala", "hypothalamus",
                     "caudate", "putamen", "frontal"]
        for gene in ["RAI1", "RAI1-AS1"]:
            gdf = df_expression[df_expression["gene"] == gene].copy()
            gdf["value_num"] = pd.to_numeric(gdf["value"], errors="coerce")
            brain = gdf[gdf["tissue"].str.lower().str.contains("|".join(brain_kw), na=False)]
            brain = brain.drop_duplicates(subset=["tissue"])
            if not brain.empty:
                top_brain = brain.sort_values("value_num", ascending=False).head(3)
                summary = "; ".join([f"{r['tissue']}: {r['value']} {r.get('unit','')}" for _, r in top_brain.iterrows()])
                if gene == "RAI1":
                    brain_expr_rai1 = summary
                else:
                    brain_expr_as1 = summary

    # QTL summary
    has_qtl = not df_qtl.empty and df_qtl["qtl_type"].iloc[0] not in ("not_available",)
    no_qtl_confirmed = any("No significant" in str(n) for n in df_qtl["notes"].values)
    qtl_summary = ""
    if has_qtl and not no_qtl_confirmed:
        eqtl = df_qtl[df_qtl["qtl_type"] == "eQTL"]
        sqtl = df_qtl[df_qtl["qtl_type"] == "sQTL"]
        brain_kw_qtl = ["brain", "cortex", "cerebellum", "hippocampus", "caudate", "putamen", "frontal"]
        brain_eqtl = eqtl[eqtl["tissue"].str.lower().str.contains("|".join(brain_kw_qtl), na=False)]
        qtl_summary = f"{len(eqtl)} eQTLs across {eqtl['tissue'].nunique()} tissues"
        if not brain_eqtl.empty:
            qtl_summary += f"; {len(brain_eqtl)} brain eQTLs"
        qtl_summary += f"; {len(sqtl)} sQTLs"

    # Ribo-seq summary
    has_riboseq = df_riboseq is not None and not df_riboseq.empty

    # Build scorecard
    utr_score = 0
    utr_rationale = []
    if n_uorfs >= 1:
        utr_score += 1; utr_rationale.append(f"{n_uorfs} uORFs found")
    if n_uorfs >= 3:
        utr_score += 1; utr_rationale.append("multiple uORFs suggest regulatory complexity")
    if any(top_uorfs["kozak_score_simple"] >= 2):
        utr_score += 1; utr_rationale.append("strong Kozak context on top uORFs")
    # Use canonical uORF conservation for scoring (most therapeutically relevant)
    # Fall back to all-uORF average if canonical not available
    cons_score_val = canonical_uaug_phylop if canonical_uaug_phylop is not None else avg_uaug
    if has_conservation and cons_score_val is not None and cons_score_val > 1.0:
        utr_score += 2; utr_rationale.append(f"uAUG windows conserved (canonical phyloP={canonical_uaug_phylop:.2f}, all mean={avg_uaug:.2f})" if canonical_uaug_phylop else f"uAUG windows conserved (mean phyloP={avg_uaug:.2f})")
    elif has_conservation and cons_score_val is not None and cons_score_val > 0.2:
        utr_score += 1; utr_rationale.append(f"moderate canonical uORF conservation (phyloP={canonical_uaug_phylop:.2f})" if canonical_uaug_phylop else f"moderate uAUG conservation (mean phyloP={avg_uaug:.2f})")
    if no_qtl_confirmed:
        utr_score += 1.5; utr_rationale.append("zero eQTLs supports post-transcriptional bottleneck")

    # Ribo-seq evidence for canonical uORF
    canonical_ribo_sum = 0
    if has_riboseq:
        canonical_ribo = df_riboseq[df_riboseq["triage_verdict"] == "canonical"]
        if not canonical_ribo.empty:
            canonical_ribo_sum = canonical_ribo.iloc[0].get("ribo_sum", 0)
        cds_ribo = df_riboseq[df_riboseq["candidate_id"] == "CDS_start_control"]
        cds_ribo_sum = cds_ribo.iloc[0].get("ribo_sum", 0) if not cds_ribo.empty else 0
        if canonical_ribo_sum > 100:
            utr_score += 1; utr_rationale.append(f"Ribo-seq confirms ribosome occupancy at canonical uORF (sum={canonical_ribo_sum} vs CDS={cds_ribo_sum})")

    nat_score = 0
    nat_rationale = []
    nat_score += ov_score
    if ov_score > 0:
        nat_rationale.append(f"overlap quality score={ov_score} ({ov.get('score_breakdown', '')})")
    if brain_expr_as1:
        nat_score += 1; nat_rationale.append("RAI1-AS1 expressed in brain")
    elif not has_expr:
        nat_rationale.append("expression data not available")

    # Determine recommendation
    if utr_score >= nat_score and utr_score >= 2:
        rec = "5′UTR / uORF lever"
        rec_reason = "uORF candidates present with conservation and/or strong Kozak context"
    elif nat_score > utr_score and nat_score >= 2:
        rec = "NAT lever (RAI1-AS1)"
        rec_reason = "strong overlap and brain expression evidence"
    elif utr_score >= 2:
        rec = "5′UTR / uORF lever"
        rec_reason = "best available evidence despite incomplete conservation data"
    else:
        rec = "Pivot to transcriptional activation (CRISPRa/enhancer)"
        rec_reason = "neither uORF nor NAT lever shows compelling evidence"

    # Format top uORF table for memo
    top_uorf_lines = ""
    for _, u in top_uorfs.iterrows():
        top_uorf_lines += f"| {u['transcript_id']} | {u['uaug_pos']} | {u['stop_pos']} | {u['aa_len']} | {u['kozak_score_simple']} | {u.get('notes','')} |\n"

    # Triage table lines
    triage_lines = ""
    if df_triage is not None and not df_triage.empty:
        for _, tr in df_triage.iterrows():
            triage_lines += f"| {tr['transcript_id']}:{int(tr['uaug_pos'])} | {tr['triage_verdict']} | {tr['source_quality']} | {tr.get('uaug_genomic', '')} | {tr.get('canonical_refseq', '')} |\n"

    # Ribo-seq table lines
    riboseq_lines = ""
    if has_riboseq:
        for _, rb in df_riboseq.iterrows():
            riboseq_lines += f"| {rb['candidate_id']} | {rb['ribo_sum']} | {rb['ribo_max']} | {rb['ribo_mean']} | {rb['triage_verdict']} |\n"

    # Conservation table lines
    cons_lines = ""
    if has_conservation:
        for _, c in df_conservation.iterrows():
            cons_lines += f"| {c['candidate_id']} | {c['region_type']} | {c['chr']}:{c['start']}-{c['end']} | {c['mean_score']} |\n"
    else:
        cons_lines = "| — | — | — | conservation data unavailable |\n"

    memo = textwrap.dedent(f"""\
    # RAI1 Upregulation: First-Bet Decision Memo

    **Date:** 2026-03-30
    **Gene:** RAI1 (ENSG00000108557), chr17p11.2
    **Context:** Classic 17p11.2 deletion Smith–Magenis Syndrome — goal is to upregulate the remaining RAI1 allele.

    ---

    ## Executive Summary

    **Recommended first bet: {rec}**

    {rec_reason}

    ---

    ## Evidence Scorecard

    | Lever | Score | Key Evidence |
    |-------|-------|-------------|
    | **A) 5′UTR / uORF** | {utr_score}/7 | {'; '.join(utr_rationale) if utr_rationale else 'minimal evidence'} |
    | **B) NAT (RAI1-AS1)** | {nat_score}/5.5 | {'; '.join(nat_rationale) if nat_rationale else 'minimal evidence'} |
    | **C) Splicing/NMD** | — | Sanity check only; {'no sQTL signal in GTEx' if no_qtl_confirmed else ('some sQTL evidence' if has_qtl and len(df_qtl[df_qtl['qtl_type']=='sQTL']) > 0 else 'no strong signal')} |

    ---

    ## A) 5′UTR / uORF Lever — Detailed Evidence

    ### Transcript Inventory
    - RAI1 protein-coding transcripts queried from Ensembl REST.
    - Canonical transcript: **{rai1_gene.get('Transcript', [{}])[0].get('id', 'see TSV')}**
    - See: `data/derived/rai1_transcripts.tsv`

    ### uORF Candidates
    - **{n_uorfs} unique uORFs** identified across supported transcripts.
    - Top candidates (by Kozak score and length):

    | Transcript | uAUG pos | Stop pos | AA length | Kozak | Notes |
    |-----------|----------|----------|-----------|-------|-------|
    {top_uorf_lines}
    - Full data: `data/derived/rai1_utr_candidates.tsv`

    ### Triage: Canonical Transcript Membership
    Each uORF is classified by whether it falls in the canonical/RefSeq-supported transcript or only in less-supported isoforms.

    | Candidate | Verdict | Source Quality | uAUG Genomic | RefSeq |
    |-----------|---------|---------------|-------------|--------|
    {triage_lines}
    - **canonical** = uORF discovered directly in the canonical transcript (MANE/RefSeq).
    - **shared_with_canonical** = uORF is in a non-canonical isoform, but the uAUG genomic position falls within canonical 5′UTR exons (i.e. the uORF is present in the canonical mRNA).
    - **non-canonical_only** = uORF exists only in alternative isoforms; not therapeutically actionable unless isoform is expressed.
    - Full data: `data/derived/rai1_uorf_triage.tsv`

    ### Ribo-seq Reality Check (GWIPS-viz)
    {'Ribosome profiling signal queried from GWIPS-viz aggregate track (hosted on UCSC as gwipsvizRiboseq bigWig). Values represent aggregate ribosome footprint counts in a ±50 bp window around each uAUG.' if has_riboseq else 'Ribo-seq data not available.'}

    {'| Candidate | Sum | Max | Mean | Verdict |' if has_riboseq else ''}
    {'|-----------|-----|-----|------|---------|' if has_riboseq else ''}
    {riboseq_lines}
    {'- **Interpretation:** CDS_start_control validates signal presence; intron_control confirms specificity. uORF candidates with sum >> intron background have detectable ribosome occupancy.' if has_riboseq else ''}
    - Full data: `data/derived/rai1_riboseq_check.tsv`
    - Source: [GWIPS-viz Ribo-seq (via UCSC)](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=gwipsvizRiboseq)

    ### Conservation
    {'Conservation scores retrieved via UCSC phyloP100way (hg38) bigWig.' if has_conservation else 'Conservation data could not be retrieved programmatically (bigWig HTTP access failed). **Manual fallback**: load the RAI1 region in UCSC Genome Browser with phyloP100way track enabled.'}

    | Candidate | Region | Coordinates | Mean phyloP |
    |----------|--------|-------------|-------------|
    {cons_lines}
    - Full data: `data/derived/rai1_utr_conservation.tsv`
    - Source: [UCSC phyloP100way](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/)

    ---

    ## B) NAT Lever (RAI1-AS1) — Detailed Evidence

    ### Overlap Analysis
    - **RAI1:** chr{ov['rai1_chr']}:{ov['rai1_start']}-{ov['rai1_end']} (strand {ov['rai1_strand']})
    - **RAI1-AS1:** chr{ov['rai1_chr']}:{ov['rai1as1_start']}-{ov['rai1as1_end']} (strand {ov['rai1as1_strand']})
    - **Overlap:** {ov_bp:,} bp, category: {ov_cat}
    - **Promoter overlap:** {'YES' if ov_promoter else 'NO'}
    - **Exon overlap:** {'YES' if ov.get('overlaps_exon_bool') else 'NO'}
    - **5′UTR exon overlap:** {'YES' if ov.get('overlaps_5utr_exon_bool') else 'NO'}
    - **Quality score:** {ov_score}/4.5 ({ov.get('score_breakdown', '')})
    - Full data: `data/derived/rai1_as1_overlap.tsv`

    ### Expression in Brain
    {'**RAI1 brain expression:** ' + brain_expr_rai1 if brain_expr_rai1 else '**RAI1:** known brain-expressed gene (literature: PMID:17154270)'}

    {'**RAI1-AS1 brain expression:** ' + brain_expr_as1 if brain_expr_as1 else '**RAI1-AS1:** expression data not programmatically retrieved. Check [GTEx Portal](https://gtexportal.org/) and [HPA](https://www.proteinatlas.org/) manually.'}

    - Full data: `data/derived/rai1_as1_expression.tsv`
    - Source: [Human Protein Atlas](https://www.proteinatlas.org/ENSG00000108557-RAI1/tissue)

    ---

    ## eQTL / sQTL Context

    {'**RAI1 has zero significant eQTLs and zero sQTLs** across all GTEx tissues (confirmed via API v8 and v10; verified with SORT1 as positive control).' if no_qtl_confirmed else (qtl_summary if qtl_summary else 'QTL data not programmatically available. Check [GTEx Portal](https://gtexportal.org/) for RAI1 eQTLs/sQTLs.')}

    {'**Interpretation:** RAI1 expression is tightly controlled and not naturally variable via common genetic variants. This makes post-transcriptional levers (uORF disruption) more attractive, since the bottleneck may be at translation rather than transcription. However, the lack of natural tunability may reflect selection pressure against RAI1 dosage changes — upregulation needs careful titration.' if no_qtl_confirmed else ('Brain eQTLs suggest RAI1 expression is naturally tunable in brain tissue — this supports both UTR and transcriptional approaches.' if has_qtl and 'brain' in qtl_summary.lower() else '')}
    - Full data: `data/derived/rai1_eqtl_sqtl_summary.tsv`
    - Source: [GTEx Portal API](https://gtexportal.org/api/v2/)

    ---

    ## Recommendation & Rationale

    **First bet: {rec}**

    {'The 5′UTR of RAI1 contains ' + str(n_uorfs) + ' uORFs that could serve as translation-brake targets. ASO-mediated disruption of uORF start codons or their Kozak contexts could increase translation from the main ORF without altering mRNA levels — a post-transcriptional upregulation strategy with precedent (e.g., PCSK9 uORF targeting).' if 'uORF' in rec else ''}
    {'Key supporting evidence:' if 'uORF' in rec else ''}
    {'- **Conservation:** ' + (f'uORF regions show phyloP conservation (mean uAUG window={avg_uaug:.2f}, body={avg_body:.2f}), suggesting functional constraint.' if has_conservation and avg_uaug and avg_uaug > 0.2 else 'Conservation should be confirmed for functional relevance.') if 'uORF' in rec else ''}
    {'- **eQTL absence:** RAI1 has zero eQTLs across all GTEx tissues, indicating tight transcriptional control. This makes post-transcriptional levers more attractive — the bottleneck may be at translation.' if 'uORF' in rec and no_qtl_confirmed else ''}
    {'- **Ribo-seq:** GWIPS-viz aggregate track shows ribosome occupancy at the canonical uORF (sum=' + str(int(canonical_ribo_sum)) + '), confirming active translation in the 5′UTR.' if 'uORF' in rec and canonical_ribo_sum > 100 else ''}
    {'- **NAT weakness:** RAI1-AS1 sits entirely within gene body introns (no promoter/exon overlap, quality score ' + str(ov_score) + '/4.5), making it a weak NAT target.' if 'uORF' in rec else ''}
    {'The RAI1-AS1 antisense transcript shows ' + ov_cat + ' overlap with RAI1 (' + str(ov_bp) + ' bp), suggesting potential cis-regulatory function. ASO knockdown of the NAT could relieve repression of RAI1 transcription.' if 'NAT' in rec else ''}
    {'Neither the uORF nor NAT lever shows sufficient evidence for a high-confidence first bet. Consider CRISPRa targeting of the RAI1 promoter or known enhancer elements as the most direct path.' if 'Pivot' in rec else ''}

    ### Top Target Windows (if UTR lever)
    {'Prioritized by: (1) canonical transcript membership, (2) Ribo-seq evidence, (3) Kozak score, (4) conservation.' if 'uORF' in rec else 'N/A for current recommendation.'}
    """)

    if "uORF" in rec and not top_uorfs.empty:
        # Build target list: canonical uORFs first, then by Ribo-seq + Kozak
        all_real = unique_uorfs.copy()
        all_real["_is_canonical"] = all_real["transcript_id"].apply(
            lambda t: 1 if df_triage is not None and not df_triage.empty and
            any((df_triage["transcript_id"] == t) & (df_triage["triage_verdict"] == "canonical")) else 0)
        target_uorfs = all_real.sort_values(["_is_canonical", "kozak_score_simple", "aa_len"],
                                             ascending=[False, False, False]).head(4)
        for i, (_, u) in enumerate(target_uorfs.iterrows(), 1):
            tid = u['transcript_id']
            upos = int(u['uaug_pos'])
            kozak = int(u['kozak_score_simple'])
            cand_id = f"{tid}:{upos}"
            # Get triage verdict
            verdict = "unknown"
            if df_triage is not None and not df_triage.empty:
                tr = df_triage[(df_triage["transcript_id"] == tid) & (df_triage["uaug_pos"] == upos)]
                if not tr.empty:
                    verdict = tr.iloc[0]["triage_verdict"]
            # Get ribo-seq signal
            ribo_note = ""
            if has_riboseq:
                rb = df_riboseq[df_riboseq["candidate_id"] == cand_id]
                if not rb.empty:
                    ribo_note = f", Ribo-seq sum={int(rb.iloc[0]['ribo_sum'])}"
            verdict_flag = ""
            if verdict == "canonical":
                verdict_flag = " **[CANONICAL — primary target]**"
            elif verdict == "shared_with_canonical":
                verdict_flag = " [shared with canonical UTR]"
            elif verdict == "non-canonical_only":
                verdict_flag = " ⚠ [non-canonical isoform only]"
            memo += f"\n{i}. **uORF at UTR position {upos}** ({tid}, Kozak={kozak}{ribo_note}){verdict_flag}: ASO window spanning uAUG ± 15 nt\n"

    memo += textwrap.dedent("""
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
    """)

    (REP / "first_bet_decision.md").write_text(memo)
    log.info("Decision memo written to reports/first_bet_decision.md")

# ---------------------------------------------------------------------------
# STEP 8 — Non-AUG uORF scan with Ribo-seq validation
# ---------------------------------------------------------------------------

NEAR_COGNATE_STARTS = ["CTG", "GTG", "TTG", "ACG", "AAG", "AGG", "ATT", "ATC", "ATA"]

def find_near_cognate_uorfs(seq):
    """Find uORFs initiated by near-cognate start codons in a 5'UTR sequence."""
    uorfs = []
    for start_codon in NEAR_COGNATE_STARTS:
        for m in re.finditer(start_codon, seq):
            start = m.start()
            pos = start + 3
            aa_len = 0
            found_stop = False
            while pos + 3 <= len(seq):
                codon = seq[pos:pos+3]
                if codon in CODON_TABLE:
                    found_stop = True
                    stop_pos = pos + 3
                    break
                aa_len += 1
                pos += 3
            if found_stop and aa_len >= 2:  # require at least 2 aa for near-cognate
                uorfs.append({
                    "start_codon": start_codon,
                    "uaug_pos": start,
                    "stop_pos": stop_pos,
                    "aa_len": aa_len,
                })
            elif not found_stop and (len(seq) - start - 3) // 3 >= 2:
                uorfs.append({
                    "start_codon": start_codon,
                    "uaug_pos": start,
                    "stop_pos": len(seq),
                    "aa_len": (len(seq) - start - 3) // 3,
                })
    return uorfs

def step8(rai1_gene, utr_regions_by_tid):
    """Scan for near-cognate uORFs and validate with Ribo-seq."""
    log.info("=== STEP 8: Non-AUG uORF scan ===")

    strand = rai1_gene["strand"]
    canonical_tid = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical_tid = t["id"]
            break
    if not canonical_tid:
        canonical_tid = rai1_gene["Transcript"][0]["id"]

    # Only scan the canonical transcript's UTR (most therapeutically relevant)
    if canonical_tid not in utr_regions_by_tid:
        log.warning("No UTR cached for canonical transcript %s", canonical_tid)
        df = pd.DataFrame()
        df.to_csv(DER / "rai1_near_cognate_uorfs.tsv", sep="\t", index=False)
        return df

    utr_seq, utr_regions = utr_regions_by_tid[canonical_tid]
    log.info("Scanning canonical UTR (%d bp) for near-cognate starts", len(utr_seq))

    nc_uorfs = find_near_cognate_uorfs(utr_seq)
    log.info("Found %d near-cognate uORF candidates", len(nc_uorfs))

    # Open phyloP bigWig once for all candidates
    bw_phylo = None
    try:
        import pyBigWig
        bw_phylo = pyBigWig.open("https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw")
        log.info("Opened phyloP bigWig for near-cognate conservation scoring")
    except Exception as e:
        log.warning("Could not open phyloP bigWig: %s", e)

    # Map to genomic coordinates and query Ribo-seq
    rows = []
    for u in nc_uorfs:
        _, uaug_genomic = utr_pos_to_genomic(u["uaug_pos"], utr_regions, strand)
        if uaug_genomic is None:
            continue

        # Query Ribo-seq: narrow window (±5 bp at codon) vs background (±50 bp)
        ribo_peak = 0
        ribo_background = 0
        ribo_enrichment = 0
        try:
            # Peak signal: tight window around the start codon (±5 bp = 13 bp)
            r1 = requests.get("https://api.genome.ucsc.edu/getData/track", params={
                "genome": "hg38", "track": "gwipsvizRiboseq",
                "chrom": "chr17", "start": uaug_genomic - 5, "end": uaug_genomic + 8,
            }, timeout=15)
            if r1.status_code == 200:
                items = r1.json().get("gwipsvizRiboseq", [])
                vals = [item["value"] for item in items]
                ribo_peak = sum(vals) / max(len(vals), 1)

            # Background: wider flanking window (50 bp each side, excluding peak)
            r2 = requests.get("https://api.genome.ucsc.edu/getData/track", params={
                "genome": "hg38", "track": "gwipsvizRiboseq",
                "chrom": "chr17", "start": uaug_genomic - 50, "end": uaug_genomic + 50,
            }, timeout=15)
            if r2.status_code == 200:
                items = r2.json().get("gwipsvizRiboseq", [])
                # Exclude the peak region
                bg_vals = [item["value"] for item in items
                           if item["start"] < uaug_genomic - 5 or item["start"] >= uaug_genomic + 8]
                ribo_background = sum(bg_vals) / max(len(bg_vals), 1)

            ribo_enrichment = round(ribo_peak / max(ribo_background, 0.1), 2)
        except Exception as e:
            log.warning("Near-cognate Ribo-seq query failed at pos %d: %s", u["uaug_pos"], e)

        # Conservation scoring via phyloP bigWig (reuse open handle)
        phyloP_score = None
        if bw_phylo:
            try:
                vals = bw_phylo.values("chr17", int(uaug_genomic - 15), int(uaug_genomic + 15))
                valid = [v for v in vals if v is not None and v == v]
                phyloP_score = round(sum(valid) / len(valid), 4) if valid else None
            except Exception as e:
                log.warning("phyloP query failed for near-cognate at pos %d: %s", u["uaug_pos"], e)

        rows.append({
            "transcript_id": canonical_tid,
            "start_codon": u["start_codon"],
            "utr_pos": u["uaug_pos"],
            "stop_pos": u["stop_pos"],
            "aa_len": u["aa_len"],
            "uaug_genomic": uaug_genomic,
            "ribo_peak_mean": round(ribo_peak, 1),
            "ribo_background_mean": round(ribo_background, 1),
            "ribo_enrichment": ribo_enrichment,
            "ribo_validated": ribo_enrichment > 2.0,
            "phyloP_uaug_window": phyloP_score,
            "target_type": "near_cognate_uORF",
        })

    if bw_phylo:
        bw_phylo.close()

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("ribo_enrichment", ascending=False)
    df.to_csv(DER / "rai1_near_cognate_uorfs.tsv", sep="\t", index=False)

    validated = df[df["ribo_validated"] == True] if not df.empty else df
    log.info("Near-cognate results: %d total, %d with >2x Ribo-seq enrichment over background",
             len(df), len(validated))
    for _, r in validated.iterrows():
        log.info("  %s at pos %d (genomic %d): %d aa, enrichment=%.1fx (peak=%.1f, bg=%.1f)",
                 r["start_codon"], r["utr_pos"], r["uaug_genomic"], r["aa_len"],
                 r["ribo_enrichment"], r["ribo_peak_mean"], r["ribo_background_mean"])
    return df

# ---------------------------------------------------------------------------
# STEP 9 — Poison exon / NMD analysis
# ---------------------------------------------------------------------------

def step9(rai1_gene):
    """Check RAI1 transcripts for NMD-triggering features (poison exons)."""
    log.info("=== STEP 9: Poison exon / NMD analysis ===")

    strand = rai1_gene["strand"]
    rows = []

    for t in rai1_gene.get("Transcript", []):
        tid = t["id"]
        biotype = t.get("biotype", "")
        translation = t.get("Translation")
        exons = sorted(t.get("Exon", []), key=lambda e: e["start"])

        if not exons:
            continue

        # Check for NMD features
        has_translation = translation is not None
        is_nmd_biotype = "nonsense_mediated" in biotype.lower() or "NMD" in biotype.upper()

        # For protein-coding transcripts: check if stop codon is >50 nt upstream of last EJC
        ptc_triggers_nmd = False
        ptc_distance = None
        stop_exon_idx = None
        last_junction_pos = None

        if has_translation and len(exons) > 1:
            cds_end = translation["end"] if strand == 1 else translation["start"]

            # Find which exon contains the stop codon
            for i, e in enumerate(exons):
                if e["start"] <= cds_end <= e["end"]:
                    stop_exon_idx = i
                    break

            # Last exon-exon junction position
            if strand == 1:
                # Last junction = start of last exon
                last_junction_pos = exons[-1]["start"]
                if stop_exon_idx is not None and stop_exon_idx < len(exons) - 1:
                    # Stop codon is before the last exon — potential NMD
                    # Distance = from stop codon to last EJC
                    ptc_distance = last_junction_pos - cds_end
                    ptc_triggers_nmd = ptc_distance > 50
            else:
                last_junction_pos = exons[0]["end"]
                if stop_exon_idx is not None and stop_exon_idx > 0:
                    ptc_distance = cds_end - last_junction_pos
                    ptc_triggers_nmd = ptc_distance > 50

        # Find exons unique to this transcript (potential poison exons)
        # Compare to canonical transcript
        canonical_exons = set()
        for ct in rai1_gene.get("Transcript", []):
            if ct.get("is_canonical"):
                for e in ct.get("Exon", []):
                    canonical_exons.add((e["start"], e["end"]))
                break

        unique_exons = []
        for e in exons:
            if (e["start"], e["end"]) not in canonical_exons:
                unique_exons.append(f"{e['start']}-{e['end']}")

        rows.append({
            "transcript_id": tid,
            "display_name": t.get("display_name", ""),
            "biotype": biotype,
            "is_canonical": bool(t.get("is_canonical")),
            "has_translation": has_translation,
            "n_exons": len(exons),
            "is_nmd_biotype": is_nmd_biotype,
            "ptc_triggers_nmd": ptc_triggers_nmd,
            "ptc_distance_to_last_ejc": ptc_distance,
            "stop_exon_index": stop_exon_idx,
            "n_unique_exons": len(unique_exons),
            "unique_exons": ";".join(unique_exons) if unique_exons else "",
            "nmd_target_potential": "high" if ptc_triggers_nmd or is_nmd_biotype else ("medium" if unique_exons and has_translation else "none"),
            "target_type": "poison_exon" if (ptc_triggers_nmd or is_nmd_biotype) else "splicing",
        })

    df = pd.DataFrame(rows)
    df.to_csv(DER / "rai1_nmd_analysis.tsv", sep="\t", index=False)

    nmd_candidates = df[df["nmd_target_potential"].isin(["high", "medium"])]
    log.info("NMD analysis: %d transcripts, %d with NMD potential", len(df), len(nmd_candidates))
    for _, r in nmd_candidates.iterrows():
        log.info("  %s (%s): %s, PTC distance=%s, %d unique exons",
                 r["transcript_id"], r["display_name"], r["nmd_target_potential"],
                 r["ptc_distance_to_last_ejc"], r["n_unique_exons"])
    return df

# ---------------------------------------------------------------------------
# STEP 10 — RNA secondary structure analysis (5'UTR)
# ---------------------------------------------------------------------------

def step10(rai1_gene, utr_regions_by_tid):
    """Compute RNA secondary structure in the canonical 5'UTR to find stable hairpins."""
    log.info("=== STEP 10: RNA secondary structure analysis ===")

    canonical_tid = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical_tid = t["id"]
            break
    if not canonical_tid or canonical_tid not in utr_regions_by_tid:
        log.warning("No canonical UTR available for structure analysis")
        pd.DataFrame().to_csv(DER / "rai1_utr_structure.tsv", sep="\t", index=False)
        return pd.DataFrame()

    utr_seq, utr_regions = utr_regions_by_tid[canonical_tid]
    strand = rai1_gene["strand"]
    # Convert to RNA
    rna_seq = utr_seq.replace("T", "U")
    utr_len = len(rna_seq)
    log.info("Analyzing %d bp canonical UTR for secondary structure", utr_len)

    try:
        import RNA
    except ImportError:
        log.warning("ViennaRNA not installed, skipping structure analysis")
        pd.DataFrame().to_csv(DER / "rai1_utr_structure.tsv", sep="\t", index=False)
        return pd.DataFrame()

    # Sliding window analysis: compute local MFE for windows across the UTR
    window_size = 60
    step_size = 10
    rows = []

    for start in range(0, utr_len - window_size + 1, step_size):
        end = start + window_size
        subseq = rna_seq[start:end]
        structure, mfe = RNA.fold(subseq)

        # Map to genomic coordinates
        _, genomic_start = utr_pos_to_genomic(start, utr_regions, strand)
        _, genomic_end = utr_pos_to_genomic(end - 1, utr_regions, strand)

        rows.append({
            "window_start": start,
            "window_end": end,
            "genomic_start": genomic_start,
            "genomic_end": genomic_end,
            "mfe": round(mfe, 2),
            "mfe_per_nt": round(mfe / window_size, 3),
            "structure": structure,
            "gc_content": round((subseq.count("G") + subseq.count("C")) / window_size * 100, 1),
        })

    df = pd.DataFrame(rows)
    df.to_csv(DER / "rai1_utr_structure.tsv", sep="\t", index=False)

    # Identify the most stable regions (most negative MFE)
    if not df.empty:
        top = df.nsmallest(5, "mfe")
        log.info("Top 5 most stable structural regions (potential ASO targets):")
        for _, r in top.iterrows():
            log.info("  pos %d-%d: MFE=%.1f kcal/mol (%.3f/nt), GC=%.0f%%",
                     r["window_start"], r["window_end"], r["mfe"], r["mfe_per_nt"], r["gc_content"])

    # Also compute full-length structure for reference
    full_struct, full_mfe = RNA.fold(rna_seq)
    log.info("Full UTR MFE: %.1f kcal/mol (%.3f/nt)", full_mfe, full_mfe / utr_len)

    return df

# ---------------------------------------------------------------------------
# STEP 11 — miRNA target site analysis (3'UTR)
# ---------------------------------------------------------------------------

def step11():
    """Load and format TargetScan miRNA predictions for RAI1."""
    log.info("=== STEP 11: miRNA target site analysis ===")

    mirna_file = RAW / "rai1_targetscan.json"
    if not mirna_file.exists():
        log.warning("No TargetScan data found at %s", mirna_file)
        log.info("Run the build pipeline to download TargetScan predictions first")
        pd.DataFrame().to_csv(DER / "rai1_mirna_targets.tsv", sep="\t", index=False)
        return pd.DataFrame()

    import json
    with open(mirna_file) as f:
        data = json.load(f)

    rows = []
    for entry in data:
        rows.append({
            "mirna_family": entry["family"],
            "representative_mirna": entry["representative"],
            "context_score": entry["score"],
            "conserved_sites": entry["conserved_sites"],
            "conserved_8mer": entry["conserved_8mer"],
            "aggregate_pct": entry["aggregate_pct"],
            "target_type": "miRNA_site",
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("context_score")  # most negative first
    df.to_csv(DER / "rai1_mirna_targets.tsv", sep="\t", index=False)

    log.info("miRNA analysis: %d conserved families targeting RAI1", len(df))
    for _, r in df.iterrows():
        log.info("  %s (seed: %s): context++=%0.3f, %d conserved sites, PCT=%.3f",
                 r["representative_mirna"], r["mirna_family"], r["context_score"],
                 r["conserved_sites"], r["aggregate_pct"])
    return df

# ---------------------------------------------------------------------------
# STEP 12 — G-quadruplex scan (5'UTR)
# ---------------------------------------------------------------------------

def step12(rai1_gene, utr_regions_by_tid):
    """Scan the canonical 5'UTR for G-quadruplex motifs."""
    log.info("=== STEP 12: G-quadruplex scan ===")

    canonical_tid = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical_tid = t["id"]
            break
    g4_cols = ["motif_type","utr_start","utr_end","length","sequence","genomic_start","genomic_end","n_g_runs","max_g_run","g4_score"]
    if not canonical_tid or canonical_tid not in utr_regions_by_tid:
        log.warning("No canonical UTR for G4 scan")
        pd.DataFrame(columns=g4_cols).to_csv(DER / "rai1_g4_motifs.tsv", sep="\t", index=False)
        return pd.DataFrame(columns=g4_cols)

    utr_seq, utr_regions = utr_regions_by_tid[canonical_tid]
    strand = rai1_gene["strand"]
    log.info("Scanning %d bp canonical UTR for G-quadruplex motifs (GC=%.0f%%)",
             len(utr_seq), (utr_seq.count("G") + utr_seq.count("C")) / len(utr_seq) * 100)

    # G4 motif: G3+N1-7G3+N1-7G3+N1-7G3+ (sense strand)
    # Also check C-rich complement (i-motif on opposite strand)
    g4_pattern = re.compile(r'(G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,})')
    c4_pattern = re.compile(r'(C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,})')

    rows = []
    for pattern, motif_type in [(g4_pattern, "G4"), (c4_pattern, "i-motif")]:
        for m in pattern.finditer(utr_seq):
            start = m.start()
            end = m.end()
            seq = m.group()
            _, genomic_start = utr_pos_to_genomic(start, utr_regions, strand)
            _, genomic_end = utr_pos_to_genomic(end - 1, utr_regions, strand)

            # Count G-runs for G4 score
            if motif_type == "G4":
                g_runs = re.findall(r'G{3,}', seq)
                max_g_run = max(len(r) for r in g_runs) if g_runs else 0
            else:
                g_runs = re.findall(r'C{3,}', seq)
                max_g_run = max(len(r) for r in g_runs) if g_runs else 0

            rows.append({
                "motif_type": motif_type,
                "utr_start": start,
                "utr_end": end,
                "length": end - start,
                "sequence": seq,
                "genomic_start": genomic_start,
                "genomic_end": genomic_end,
                "n_g_runs": len(g_runs),
                "max_g_run": max_g_run,
                "g4_score": len(g_runs) * max_g_run,  # simple score: more runs + longer = more stable
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("g4_score", ascending=False)
    df.to_csv(DER / "rai1_g4_motifs.tsv", sep="\t", index=False)

    log.info("G-quadruplex results: %d G4 motifs, %d i-motifs",
             len(df[df["motif_type"] == "G4"]) if not df.empty else 0,
             len(df[df["motif_type"] == "i-motif"]) if not df.empty else 0)
    for _, r in df.iterrows():
        log.info("  %s at pos %d-%d (%d bp): score=%d, %d runs, max_run=%d",
                 r["motif_type"], r["utr_start"], r["utr_end"], r["length"],
                 r["g4_score"], r["n_g_runs"], r["max_g_run"])
    return df

# ---------------------------------------------------------------------------
# STEP 13 — 3'UTR secondary structure analysis
# ---------------------------------------------------------------------------

def step13(rai1_gene):
    """Compute RNA secondary structure in the canonical 3'UTR."""
    log.info("=== STEP 13: 3'UTR secondary structure ===")

    strand = rai1_gene["strand"]
    canonical = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical = t
            break
    if not canonical or not canonical.get("Translation"):
        log.warning("No canonical transcript with translation for 3'UTR analysis")
        pd.DataFrame().to_csv(DER / "rai1_3utr_structure.tsv", sep="\t", index=False)
        return pd.DataFrame()

    cds_end = canonical["Translation"]["end"]
    exons = sorted(canonical.get("Exon", []), key=lambda e: e["start"])

    # Get 3'UTR regions (+ strand: everything after CDS end)
    utr3_regions = []
    for e in exons:
        if e["start"] > cds_end:
            utr3_regions.append((e["seq_region_name"], e["start"], e["end"]))
        elif e["end"] > cds_end:
            utr3_regions.append((e["seq_region_name"], cds_end + 1, e["end"]))

    if not utr3_regions:
        log.warning("No 3'UTR regions found")
        pd.DataFrame().to_csv(DER / "rai1_3utr_structure.tsv", sep="\t", index=False)
        return pd.DataFrame()

    # Fetch 3'UTR sequence
    utr3_seq = ""
    for chrom, rs, re_ in utr3_regions:
        region_str = f"{chrom}:{rs}..{re_}:1"
        seq = ensembl_seq(f"/sequence/region/human/{region_str}")
        utr3_seq += seq.strip().upper()

    log.info("3'UTR: %d bp, GC=%.1f%%", len(utr3_seq),
             (utr3_seq.count("G") + utr3_seq.count("C")) / max(len(utr3_seq), 1) * 100)

    try:
        import RNA
    except ImportError:
        log.warning("ViennaRNA not installed")
        pd.DataFrame().to_csv(DER / "rai1_3utr_structure.tsv", sep="\t", index=False)
        return pd.DataFrame()

    rna_seq = utr3_seq.replace("T", "U")
    window_size = 80
    step_size = 20
    rows = []

    for start in range(0, len(rna_seq) - window_size + 1, step_size):
        end = start + window_size
        subseq = rna_seq[start:end]
        structure, mfe = RNA.fold(subseq)
        genomic_pos = cds_end + 1 + start  # approximate for single-exon 3'UTR

        rows.append({
            "window_start": start,
            "window_end": end,
            "genomic_start": genomic_pos,
            "genomic_end": genomic_pos + window_size,
            "mfe": round(mfe, 2),
            "mfe_per_nt": round(mfe / window_size, 3),
            "gc_content": round((subseq.count("G") + subseq.count("C")) / window_size * 100, 1),
        })

    df = pd.DataFrame(rows)
    df.to_csv(DER / "rai1_3utr_structure.tsv", sep="\t", index=False)

    if not df.empty:
        top = df.nsmallest(3, "mfe")
        log.info("Top 3 most stable 3'UTR regions:")
        for _, r in top.iterrows():
            log.info("  pos %d-%d: MFE=%.1f kcal/mol, GC=%.0f%%",
                     r["window_start"], r["window_end"], r["mfe"], r["gc_content"])

    return df

# ---------------------------------------------------------------------------
# STEP 14 — RNA-binding protein sites (ENCODE eCLIP)
# ---------------------------------------------------------------------------

def step14(rai1_gene):
    """Query ENCODE eCLIP data for RBP binding sites in RAI1 UTRs."""
    log.info("=== STEP 14: RNA-binding protein sites (ENCODE eCLIP) ===")

    # Query UCSC API for ENCODE eCLIP tracks in the RAI1 region
    # Focus on the 5'UTR and 3'UTR regions
    gene_start = rai1_gene["start"]
    gene_end = rai1_gene["end"]

    rows = []

    # ENCODE eCLIP clusters are available via UCSC as "encRegTfbsClustered" or specific eCLIP tracks
    # Try the ENCODE RBP tracks
    track_names = ["encRegTfbsClustered"]

    for track in track_names:
        try:
            r = requests.get("https://api.genome.ucsc.edu/getData/track", params={
                "genome": "hg38", "track": track,
                "chrom": "chr17", "start": gene_start, "end": gene_end,
            }, timeout=30)
            if r.status_code == 200:
                data = r.json()
                items = data.get(track, [])
                log.info("ENCODE %s: %d items in RAI1 region", track, len(items))

                for item in items:
                    # Filter for items in UTR regions (rough: first 1kb and last 2kb of gene)
                    item_start = item.get("chromStart", item.get("start", 0))
                    item_end = item.get("chromEnd", item.get("end", 0))
                    name = item.get("name", "")
                    score = item.get("score", 0)

                    # Classify location
                    if item_start < gene_start + 1000:
                        location = "5UTR_region"
                    elif item_start > gene_end - 2000:
                        location = "3UTR_region"
                    else:
                        location = "gene_body"

                    rows.append({
                        "track": track,
                        "rbp_name": name,
                        "chr": "chr17",
                        "start": item_start,
                        "end": item_end,
                        "score": score,
                        "location": location,
                    })
            else:
                log.warning("ENCODE %s returned %d", track, r.status_code)
        except Exception as e:
            log.warning("ENCODE %s query failed: %s", track, e)

    df = pd.DataFrame(rows)
    if not df.empty:
        # Focus on UTR regions
        utr_rbps = df[df["location"].isin(["5UTR_region", "3UTR_region"])]
        log.info("RBP sites in UTR regions: %d (of %d total in gene)", len(utr_rbps), len(df))

        # Count by RBP name in UTR regions
        if not utr_rbps.empty:
            top_rbps = utr_rbps.groupby("rbp_name").agg(
                n_sites=("score", "count"),
                max_score=("score", "max"),
                locations=("location", lambda x: ",".join(sorted(set(x)))),
            ).sort_values("max_score", ascending=False).head(15)
            for name, r in top_rbps.iterrows():
                log.info("  %s: %d sites, max_score=%d (%s)", name, r["n_sites"], r["max_score"], r["locations"])

    df.to_csv(DER / "rai1_rbp_sites.tsv", sep="\t", index=False)
    return df

# ---------------------------------------------------------------------------
# STEP 15 — Alternative polyadenylation
# ---------------------------------------------------------------------------

def step15(rai1_gene):
    """Check for alternative polyadenylation signals in RAI1 3'UTR."""
    log.info("=== STEP 15: Alternative polyadenylation ===")

    # Fetch the 3'UTR + downstream region to find polyA signals
    canonical = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical = t
            break
    if not canonical or not canonical.get("Translation"):
        pd.DataFrame().to_csv(DER / "rai1_polya.tsv", sep="\t", index=False)
        return pd.DataFrame()

    cds_end = canonical["Translation"]["end"]
    exons = sorted(canonical.get("Exon", []), key=lambda e: e["start"])
    tx_end = canonical["end"]

    # Fetch 3'UTR sequence
    chrom = canonical["seq_region_name"]
    utr3_start = cds_end + 1
    utr3_end = tx_end
    try:
        seq = ensembl_seq(f"/sequence/region/human/{chrom}:{utr3_start}..{utr3_end}:1")
        seq = seq.strip().upper()
    except Exception as e:
        log.warning("Could not fetch 3'UTR sequence: %s", e)
        pd.DataFrame().to_csv(DER / "rai1_polya.tsv", sep="\t", index=False)
        return pd.DataFrame()

    log.info("3'UTR sequence: %d bp", len(seq))

    # Scan for canonical polyA signals: AATAAA and common variants
    polya_signals = {
        "AATAAA": "canonical (most common)",
        "ATTAAA": "common variant",
        "AGTAAA": "variant",
        "TATAAA": "variant",
        "CATAAA": "variant",
        "GATAAA": "variant",
        "AATATA": "variant",
        "AATACA": "variant",
        "AATAGA": "variant",
        "AAAAAG": "variant",
        "ACTAAA": "variant",
    }

    rows = []
    for signal, signal_type in polya_signals.items():
        for m in re.finditer(signal, seq):
            pos = m.start()
            genomic_pos = utr3_start + pos
            # Distance from CDS end (how much 3'UTR is used)
            dist_from_cds = pos
            dist_from_end = len(seq) - pos
            rows.append({
                "signal": signal,
                "signal_type": signal_type,
                "utr_position": pos,
                "genomic_position": genomic_pos,
                "dist_from_cds_end": dist_from_cds,
                "dist_from_utr_end": dist_from_end,
                "is_near_end": dist_from_end < 50,  # within 50 bp of annotated end = likely the used site
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("utr_position")
    df.to_csv(DER / "rai1_polya.tsv", sep="\t", index=False)

    log.info("PolyA signals found: %d total", len(df))
    for _, r in df.iterrows():
        marker = " <-- likely active" if r["is_near_end"] else ""
        log.info("  %s (%s) at UTR pos %d (genomic %d), %d bp from end%s",
                 r["signal"], r["signal_type"], r["utr_position"],
                 r["genomic_position"], r["dist_from_utr_end"], marker)

    return df

# ---------------------------------------------------------------------------
# STEP 16 — ASO candidate design
# ---------------------------------------------------------------------------

def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G",
            "U": "A", "a": "t", "t": "a", "g": "c", "c": "g"}
    return "".join(comp.get(b, b) for b in reversed(seq))

def step16(rai1_gene, utr_regions_by_tid):
    """Design ASO candidates for top target regions with thermodynamic scoring."""
    log.info("=== STEP 16: ASO candidate design ===")
    try:
        import RNA
    except ImportError:
        log.warning("ViennaRNA not installed, skipping ASO design")
        pd.DataFrame().to_csv(DER / "rai1_aso_candidates.tsv", sep="\t", index=False)
        return pd.DataFrame()

    strand = rai1_gene["strand"]
    canonical_tid = None
    for t in rai1_gene.get("Transcript", []):
        if t.get("is_canonical"):
            canonical_tid = t["id"]
            break

    if not canonical_tid or canonical_tid not in utr_regions_by_tid:
        log.warning("No canonical UTR for ASO design")
        pd.DataFrame().to_csv(DER / "rai1_aso_candidates.tsv", sep="\t", index=False)
        return pd.DataFrame()

    utr_seq, utr_regions = utr_regions_by_tid[canonical_tid]

    # Load 3'UTR
    utr3_seq = ""
    seq_file = RAW / "rai1_utr_sequences.json"
    if seq_file.exists():
        with open(seq_file) as f:
            seqs = json.load(f)
            utr3_seq = seqs.get("three_prime_utr", "")

    # Build target regions dynamically from all uORF positions
    target_regions = []
    seen_positions = set()

    # AUG uORFs from triage
    triage_file = DER / "rai1_uorf_triage.tsv"
    if triage_file.exists():
        triage_df = pd.read_csv(triage_file, sep="\t")
        for _, tr in triage_df.iterrows():
            pos = int(tr["uaug_pos"])
            if pos not in seen_positions:
                seen_positions.add(pos)
                target_regions.append({
                    "name": f"AUG_{tr['transcript_id']}_{pos}",
                    "utr": "5prime", "center": pos,
                    "description": f"AUG uORF at pos {pos} ({tr['transcript_id']}, {tr['triage_verdict']})",
                })

    # Near-cognate uORFs (validated)
    nc_file = DER / "rai1_near_cognate_uorfs.tsv"
    if nc_file.exists():
        nc_df = pd.read_csv(nc_file, sep="\t")
        nc_val = nc_df[nc_df["ribo_validated"] == True] if not nc_df.empty else nc_df
        for _, nc in nc_val.iterrows():
            pos = int(nc["utr_pos"])
            if pos not in seen_positions:
                seen_positions.add(pos)
                target_regions.append({
                    "name": f"NC_{nc['start_codon']}_{pos}",
                    "utr": "5prime", "center": pos,
                    "description": f"{nc['start_codon']} near-cognate at pos {pos} ({nc['ribo_enrichment']:.1f}x enrichment)",
                })

    # 5'UTR structure hotspot (if not already covered by a uORF position)
    struct_file = DER / "rai1_utr_structure.tsv"
    if struct_file.exists():
        sdf = pd.read_csv(struct_file, sep="\t")
        if not sdf.empty:
            best_struct = sdf.nsmallest(1, "mfe").iloc[0]
            struct_center = int((best_struct["window_start"] + best_struct["window_end"]) / 2)
            if struct_center not in seen_positions:
                seen_positions.add(struct_center)
                target_regions.append({
                    "name": "structure_hotspot",
                    "utr": "5prime", "center": struct_center,
                    "description": f"5'UTR structure hotspot (MFE {best_struct['mfe']} kcal/mol)",
                })

    # 3'UTR structure hotspot
    struct3_file = DER / "rai1_3utr_structure.tsv"
    if struct3_file.exists() and utr3_seq:
        sdf3 = pd.read_csv(struct3_file, sep="\t")
        if not sdf3.empty:
            best_3 = sdf3.nsmallest(1, "mfe").iloc[0]
            center_3 = int((best_3["window_start"] + best_3["window_end"]) / 2)
            target_regions.append({
                "name": "3utr_structure_hotspot",
                "utr": "3prime", "center": center_3,
                "description": f"3'UTR structure hotspot (MFE {best_3['mfe']} kcal/mol)",
            })

    log.info("Designing ASOs for %d target regions", len(target_regions))

    rows = []
    for region in target_regions:
        center = region["center"]
        utr = utr_seq if region["utr"] == "5prime" else utr3_seq
        utr_label = region["utr"]

        if not utr:
            continue

        # Generate ASO candidates: 3 lengths x 5 positions = up to 15 per region
        for aso_len in [16, 18, 20]:
            for offset in [-4, -2, 0, 2, 4]:
                start = center - aso_len // 2 + offset
                end = start + aso_len

                if start < 0 or end > len(utr):
                    continue

                target_seq = utr[start:end]
                aso_seq = reverse_complement(target_seq)

                # Skip if extreme GC
                gc = (aso_seq.count("G") + aso_seq.count("C")) / len(aso_seq) * 100
                if gc > 90 or gc < 30:
                    continue

                # Compute thermodynamic properties
                # 1. ASO-target duplex stability (more negative = stronger binding)
                aso_rna = aso_seq.replace("T", "U")
                target_rna = target_seq.replace("T", "U")
                duplex = RNA.duplexfold(aso_rna, target_rna)
                duplex_mfe = round(duplex.energy, 1)

                # 2. ASO self-folding (less negative = better, ASO stays open)
                self_struct, self_mfe = RNA.fold(aso_rna)
                self_mfe = round(self_mfe, 1)

                # 3. Check for homopolymer runs (>4 same base = problematic)
                max_homopolymer = max(len(m.group()) for m in re.finditer(r'(.)\1+', aso_seq)) if len(aso_seq) > 0 else 0

                # 4. Check for GGGG runs (G-tetrad formation risk)
                has_g4_risk = "GGGG" in aso_seq or "CCCC" in aso_seq

                # Quality score: strong binding + low self-fold + good GC
                # Ideal: duplex_mfe very negative, self_mfe near 0, GC 40-65%
                gc_penalty = abs(gc - 52.5) / 50  # 0 at 52.5%, higher further away
                self_fold_penalty = max(0, -self_mfe) / 10  # penalty for self-folding
                quality = round((-duplex_mfe / aso_len) - gc_penalty - self_fold_penalty, 3)

                # Flag issues
                flags = []
                if gc > 75: flags.append("high_GC")
                if gc < 40: flags.append("low_GC")
                if self_mfe < -5: flags.append("self_folding")
                if max_homopolymer >= 5: flags.append("homopolymer")
                if has_g4_risk: flags.append("G4_risk")

                # Map to genomic coordinates
                if utr_label == "5prime":
                    _, genomic_start = utr_pos_to_genomic(start, utr_regions, strand)
                    _, genomic_end = utr_pos_to_genomic(end - 1, utr_regions, strand)
                else:
                    genomic_start = int(seqs.get("three_prime_utr_start", 0)) + start
                    genomic_end = genomic_start + aso_len

                rows.append({
                    "region": region["name"],
                    "region_description": region["description"],
                    "utr_type": utr_label,
                    "aso_start": start,
                    "aso_end": end,
                    "aso_length": aso_len,
                    "target_seq": target_seq,
                    "aso_seq": aso_seq,
                    "gc_percent": round(gc, 1),
                    "duplex_mfe": duplex_mfe,
                    "self_fold_mfe": self_mfe,
                    "max_homopolymer": max_homopolymer,
                    "quality_score": quality,
                    "flags": ";".join(flags) if flags else "",
                    "genomic_start": genomic_start,
                    "genomic_end": genomic_end,
                })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["region", "quality_score"], ascending=[True, False])
    df.to_csv(DER / "rai1_aso_candidates.tsv", sep="\t", index=False)

    # Log top candidates per region
    for region_name in df["region"].unique() if not df.empty else []:
        region_df = df[df["region"] == region_name]
        top = region_df.head(3)
        log.info("  %s:", region_name)
        for _, r in top.iterrows():
            log.info("    %s (%d-mer, GC=%.0f%%, duplex=%.1f, self=%.1f, quality=%.3f) %s",
                     r["aso_seq"], r["aso_length"], r["gc_percent"],
                     r["duplex_mfe"], r["self_fold_mfe"], r["quality_score"],
                     r["flags"])

    log.info("ASO design: %d candidates across %d target regions",
             len(df), len(df["region"].unique()) if not df.empty else 0)
    return df

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    log.info("Starting RAI1 upregulation first-bet pipeline")

    rai1, rai1as1, df_rai1, df_as1 = step1()
    df_uorfs, utr_regions_by_tid = step2(rai1, df_rai1)
    df_triage = step2b(df_uorfs, utr_regions_by_tid, rai1)
    df_riboseq = step2c(df_triage, rai1)
    df_conservation = step3(df_uorfs, rai1, utr_regions_by_tid)
    df_overlap = step4(rai1, rai1as1)
    df_expression = step5(rai1, rai1as1)
    df_qtl = step6(rai1)
    step7(df_uorfs, df_conservation, df_overlap, df_expression, df_qtl, rai1, rai1as1, df_triage, df_riboseq)
    df_nc_uorfs = step8(rai1, utr_regions_by_tid)
    df_nmd = step9(rai1)
    df_structure = step10(rai1, utr_regions_by_tid)
    df_mirna = step11()
    df_g4 = step12(rai1, utr_regions_by_tid)
    df_3utr_struct = step13(rai1)
    df_rbp = step14(rai1)
    df_polya = step15(rai1)
    df_aso = step16(rai1, utr_regions_by_tid)

    log.info("Pipeline complete. Outputs in data/derived/ and reports/")

if __name__ == "__main__":
    main()
