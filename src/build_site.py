#!/usr/bin/env python3
"""
Build the SMS ASO Target Prioritization Portal.
Reads pipeline TSV outputs and generates a static website.
"""

import json, sys, shutil, pathlib
import pandas as pd
from jinja2 import Environment, FileSystemLoader

ROOT = pathlib.Path(__file__).resolve().parent.parent
DER = ROOT / "data" / "derived"
RAW = ROOT / "data" / "raw"
SITE = ROOT / "site"
TMPL = ROOT / "templates"

# ---------------------------------------------------------------------------
# Phase 1: Load and validate data
# ---------------------------------------------------------------------------

REQUIRED_FILES = [
    "rai1_transcripts.tsv", "rai1_utr_candidates.tsv", "rai1_uorf_triage.tsv",
    "rai1_riboseq_check.tsv", "rai1_utr_conservation.tsv", "rai1_as1_overlap.tsv",
    "rai1_as1_expression.tsv", "rai1_eqtl_sqtl_summary.tsv", "rai1as1_transcripts.tsv",
]

def load_data():
    for f in REQUIRED_FILES:
        if not (DER / f).exists():
            print(f"ERROR: Missing required file: data/derived/{f}", file=sys.stderr)
            print("Run 'python src/pipeline.py' first.", file=sys.stderr)
            sys.exit(1)

    data = {}
    for f in REQUIRED_FILES:
        key = f.replace(".tsv", "")
        data[key] = pd.read_csv(DER / f, sep="\t")

    # Load raw gene JSON for exon coordinates
    with open(RAW / "rai1_gene.json") as fh:
        data["rai1_gene"] = json.load(fh)

    return data

# ---------------------------------------------------------------------------
# Phase 2: Build unified target records
# ---------------------------------------------------------------------------

def traffic(val, green_thresh, yellow_thresh):
    if val >= green_thresh: return "green"
    if val >= yellow_thresh: return "yellow"
    return "red"

# Start codon strength: AUG > CUG/CTG > GUG/GTG > ACG > others
START_CODON_STRENGTH = {
    "ATG": 1.0, "CTG": 0.5, "CUG": 0.5, "GTG": 0.4, "GUG": 0.4,
    "ACG": 0.3, "TTG": 0.25, "UUG": 0.25,
    "AAG": 0.15, "AGG": 0.15, "ATC": 0.1, "ATT": 0.1, "ATA": 0.1,
}

def unified_score(canonical_weight, ribo_signal, codon_strength, phyloP, cds_control_sum, designability=1.0):
    """Unified composite score for both AUG and near-cognate uORFs (0-10).

    Components (max 8 raw, scaled to 10):
      - Canonical membership:  3 * weight (1.0 = canonical, 0.33 = shared, 0 = non-canonical)
      - Ribo-seq evidence:     2 * min(signal / cds_control, 1.0)
      - Start codon strength:  1 * codon_strength (AUG=1.0, CUG=0.5, GTG=0.4, etc.)
      - Conservation:          1 * min(|phyloP|, 3.0) / 3.0
      - ASO designability:     1 * designability (1.0 = clean candidates, 0 = no clean ASOs)
    """
    return round((
        3.0 * canonical_weight
        + 2.0 * min(ribo_signal / max(cds_control_sum, 1.0), 1.0)
        + 1.0 * codon_strength
        + 1.0 * min(abs(phyloP), 3.0) / 3.0
        + 1.0 * designability
    ) / 8.0 * 10.0, 1)

def compute_designability(utr_pos, aso_df, utr_type="5prime"):
    """Compute ASO designability for a given UTR position.
    Returns (score 0-1, n_clean, n_total, best_candidates)."""
    if aso_df is None or aso_df.empty:
        return 0.5, 0, 0, []  # unknown = neutral

    # Find ASO candidates that cover this UTR position (center within ±12 of target)
    mask = (aso_df["utr_type"] == utr_type) & \
           (aso_df["aso_start"] <= utr_pos + 5) & \
           (aso_df["aso_end"] >= utr_pos - 5)
    nearby = aso_df[mask]

    if nearby.empty:
        return 0.5, 0, 0, []  # no ASOs designed for this position = neutral

    clean = nearby[nearby["flags"].isna() | (nearby["flags"] == "")]
    n_clean = len(clean)
    n_total = len(nearby)

    # Designability: fraction of clean candidates, min 0, max 1
    score = n_clean / n_total if n_total > 0 else 0

    # Get ALL candidates for this target, deduplicated by sequence, sorted by quality
    sorted_nearby = nearby.drop_duplicates(subset=["aso_seq"]).sort_values("quality_score", ascending=False)
    best = []
    for _, r in sorted_nearby.iterrows():
        best.append({
            "seq": r["aso_seq"],
            "length": int(r["aso_length"]),
            "gc": float(r["gc_percent"]),
            "duplex_mfe": float(r["duplex_mfe"]),
            "self_mfe": float(r["self_fold_mfe"]),
            "quality": round(float(r["quality_score"]), 3),
            "flags": str(r["flags"]) if pd.notna(r["flags"]) and r["flags"] else "",
            "target_seq": r["target_seq"],
        })

    return round(score, 2), n_clean, n_total, best

def build_targets(data):
    triage = data["rai1_uorf_triage"]
    candidates = data["rai1_utr_candidates"]
    riboseq = data["rai1_riboseq_check"]
    conservation = data["rai1_utr_conservation"]

    # Load ASO candidates for designability scoring
    aso_file = DER / "rai1_aso_candidates.tsv"
    aso_df = pd.read_csv(aso_file, sep="\t") if aso_file.exists() else None

    # Get CDS control for normalization
    cds_row = riboseq[riboseq["candidate_id"] == "CDS_start_control"]
    cds_control_sum = float(cds_row.iloc[0]["ribo_sum"]) if not cds_row.empty else 1.0
    if cds_control_sum == 0:
        cds_control_sum = 1.0

    targets = []

    # --- AUG uORFs ---
    for _, tr in triage.iterrows():
        tid = tr["transcript_id"]
        upos = int(tr["uaug_pos"])
        cand_id = f"{tid}:{upos}"
        slug = f"{tid}-{upos}"

        cand = candidates[(candidates["transcript_id"] == tid) & (candidates["uaug_pos"] == upos)]
        kozak = float(cand.iloc[0]["kozak_score_simple"]) if not cand.empty else 0
        aa_len = int(cand.iloc[0]["aa_len"]) if not cand.empty and pd.notna(cand.iloc[0]["aa_len"]) else 0
        utr_len = int(cand.iloc[0]["utr_len"]) if not cand.empty else 0
        utr_gc = float(cand.iloc[0]["utr_gc"]) if not cand.empty else 0
        stop_pos = int(cand.iloc[0]["stop_pos"]) if not cand.empty and pd.notna(cand.iloc[0]["stop_pos"]) else 0
        raw_notes = cand.iloc[0].get("notes", "") if not cand.empty else ""
        notes = str(raw_notes) if pd.notna(raw_notes) and raw_notes != "" else ""

        ribo = riboseq[riboseq["candidate_id"] == cand_id]
        ribo_sum = float(ribo.iloc[0]["ribo_sum"]) if not ribo.empty else 0
        ribo_max = float(ribo.iloc[0]["ribo_max"]) if not ribo.empty else 0
        ribo_mean = float(ribo.iloc[0]["ribo_mean"]) if not ribo.empty else 0
        ribo_region = str(ribo.iloc[0]["region"]) if not ribo.empty else ""

        cons = conservation[conservation["candidate_id"] == cand_id]
        cons_dict = {}
        for _, c in cons.iterrows():
            rtype = c["region_type"]
            cons_dict[rtype] = {
                "coords": f"{c['chr']}:{int(c['start'])}-{int(c['end'])}",
                "phyloP": round(float(c["mean_score"]), 4) if pd.notna(c["mean_score"]) else None,
            }

        uaug_phyloP = cons_dict.get("uAUG_window", {}).get("phyloP", 0) or 0

        verdict = tr["triage_verdict"]
        canonical_weight = 1.0 if verdict == "canonical" else (0.33 if verdict == "shared_with_canonical" else 0.0)

        # For AUG uORFs, codon_strength includes Kozak: AUG base (1.0) scaled by Kozak
        # Kozak 3 = full AUG strength, Kozak 0 = reduced
        codon_str = START_CODON_STRENGTH["ATG"] * (0.5 + 0.5 * kozak / 3.0)

        # ASO designability for this target position
        design_score, n_clean, n_total, best_asos = compute_designability(upos, aso_df, "5prime")

        composite = unified_score(canonical_weight, ribo_sum, codon_str, uaug_phyloP, cds_control_sum, design_score)

        lights = {
            "canonical": "green" if verdict == "canonical" else ("yellow" if verdict == "shared_with_canonical" else "red"),
            "riboseq": traffic(ribo_sum, 100, 10),
            "conservation": traffic(abs(uaug_phyloP), 1.0, 0.2),
            "initiation": traffic(codon_str, 0.8, 0.4),
            "designability": traffic(design_score, 0.5, 0.1),
        }

        uaug_genomic = int(tr["uaug_genomic"]) if pd.notna(tr["uaug_genomic"]) else 0
        aso_start = uaug_genomic - 15
        aso_end = uaug_genomic + 15

        tx_df = data["rai1_transcripts"]
        tx_row = tx_df[tx_df["transcript_id"] == tid]
        tx_name = str(tx_row.iloc[0]["display_name"]) if not tx_row.empty else tid

        targets.append({
            "id": cand_id,
            "slug": slug,
            "transcript_id": tid,
            "transcript_name": tx_name,
            "target_type": "AUG_uORF",
            "start_codon": "ATG",
            "codon_strength": round(codon_str, 2),
            "uaug_pos": upos,
            "stop_pos": stop_pos,
            "uaug_genomic": uaug_genomic,
            "aa_len": aa_len,
            "utr_len": utr_len,
            "utr_gc": utr_gc,
            "kozak": kozak,
            "triage_verdict": verdict,
            "source_quality": str(tr.get("source_quality", "")),
            "refseq": str(tr.get("canonical_refseq", "")) if pd.notna(tr.get("canonical_refseq", "")) else "",
            "ribo": {"sum": ribo_sum, "max": ribo_max, "mean": ribo_mean, "region": ribo_region},
            "conservation": cons_dict,
            "composite_score": composite,
            "lights": lights,
            "designability": {"score": design_score, "clean": n_clean, "total": n_total},
            "best_asos": best_asos,
            "aso_window": {"chr": "chr17", "start": aso_start, "end": aso_end},
            "ucsc_link": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17:{aso_start-200}-{aso_end+200}",
            "notes": notes,
            "cds_control_sum": cds_control_sum,
        })

    # --- Near-cognate uORFs (validated only) ---
    nc_file = DER / "rai1_near_cognate_uorfs.tsv"
    if nc_file.exists():
        nc_df = pd.read_csv(nc_file, sep="\t")
        nc_validated = nc_df[nc_df["ribo_validated"] == True] if not nc_df.empty else nc_df

        # Get canonical transcript info
        canonical_tid = None
        for _, tr in triage.iterrows():
            if tr["triage_verdict"] == "canonical":
                canonical_tid = tr["transcript_id"]
                break

        tx_df = data["rai1_transcripts"]
        tx_row = tx_df[tx_df["transcript_id"] == canonical_tid] if canonical_tid else pd.DataFrame()
        tx_name = str(tx_row.iloc[0]["display_name"]) if not tx_row.empty else canonical_tid or ""
        utr_len_canonical = 0
        cand_canonical = candidates[candidates["transcript_id"] == canonical_tid]
        if not cand_canonical.empty:
            utr_len_canonical = int(cand_canonical.iloc[0]["utr_len"])

        for _, nc in nc_validated.iterrows():
            start_codon = nc["start_codon"]
            utr_pos = int(nc["utr_pos"])
            uaug_genomic = int(nc["uaug_genomic"])
            aa_len = int(nc["aa_len"])
            stop_pos = int(nc["stop_pos"])
            ribo_enrichment = float(nc["ribo_enrichment"])
            ribo_peak = float(nc.get("ribo_peak_mean", 0))
            phyloP = float(nc["phyloP_uaug_window"]) if pd.notna(nc.get("phyloP_uaug_window")) else 0

            # Near-cognate uORFs in canonical UTR get canonical weight = 1.0
            canonical_weight = 1.0

            codon_str = START_CODON_STRENGTH.get(start_codon, 0.1)

            ribo_signal = round(ribo_peak * ribo_enrichment, 2)

            # ASO designability
            design_score, n_clean, n_total, best_asos = compute_designability(utr_pos, aso_df, "5prime")

            composite = unified_score(canonical_weight, ribo_signal, codon_str, phyloP, cds_control_sum, design_score)

            slug = f"nc-{start_codon}-{utr_pos}"
            cand_id = f"{start_codon}@{utr_pos}"
            aso_start = uaug_genomic - 15
            aso_end = uaug_genomic + 15

            lights = {
                "canonical": "green",
                "riboseq": traffic(ribo_enrichment, 3.0, 2.0),
                "conservation": traffic(abs(phyloP), 1.0, 0.2),
                "initiation": traffic(codon_str, 0.8, 0.4),
                "designability": traffic(design_score, 0.5, 0.1),
            }

            targets.append({
                "id": cand_id,
                "slug": slug,
                "transcript_id": canonical_tid or "",
                "transcript_name": tx_name,
                "target_type": "near_cognate_uORF",
                "start_codon": start_codon,
                "codon_strength": round(codon_str, 2),
                "uaug_pos": utr_pos,
                "stop_pos": stop_pos,
                "uaug_genomic": uaug_genomic,
                "aa_len": aa_len,
                "utr_len": utr_len_canonical,
                "utr_gc": 0,
                "kozak": 0,
                "triage_verdict": "canonical_utr",
                "source_quality": "canonical UTR scan",
                "refseq": "",
                "ribo": {"sum": ribo_signal, "max": ribo_peak, "mean": ribo_peak,
                         "region": f"chr17:{uaug_genomic-30}-{uaug_genomic+30}",
                         "enrichment": ribo_enrichment},
                "conservation": {"uAUG_window": {"coords": f"chr17:{uaug_genomic-15}-{uaug_genomic+15}", "phyloP": phyloP}},
                "composite_score": composite,
                "lights": lights,
                "designability": {"score": design_score, "clean": n_clean, "total": n_total},
                "best_asos": best_asos,
                "aso_window": {"chr": "chr17", "start": aso_start, "end": aso_end},
                "ucsc_link": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17:{aso_start-200}-{aso_end+200}",
                "notes": f"near-cognate {start_codon}, {ribo_enrichment}x Ribo-seq enrichment",
                "cds_control_sum": cds_control_sum,
            })

    # Sort all targets by composite score descending
    targets.sort(key=lambda t: t["composite_score"], reverse=True)
    for i, t in enumerate(targets):
        t["rank"] = i + 1

    return targets

# ---------------------------------------------------------------------------
# Phase 3: Build expression data
# ---------------------------------------------------------------------------

def build_expression(data):
    df = data["rai1_as1_expression"]
    rai1 = df[df["gene"] == "RAI1"].copy()
    rai1 = rai1[rai1["source"] != "not_available"]
    rai1 = rai1.drop_duplicates(subset=["tissue"])
    rai1["value_num"] = pd.to_numeric(rai1["value"], errors="coerce")
    rai1 = rai1.sort_values("value_num", ascending=False).head(15)

    tissues = []
    for _, r in rai1.iterrows():
        brain_kw = ["brain", "cortex", "cerebellum", "hippocampus", "amygdala",
                     "hypothalamus", "caudate", "putamen", "midbrain", "pituitary"]
        is_brain = any(kw in str(r["tissue"]).lower() for kw in brain_kw)
        tissues.append({
            "tissue": r["tissue"],
            "value": float(r["value_num"]) if pd.notna(r["value_num"]) else 0,
            "is_brain": is_brain,
        })

    return tissues

# ---------------------------------------------------------------------------
# Phase 4: Build transcript structure for genome viewer
# ---------------------------------------------------------------------------

def build_transcripts(data):
    gene = data["rai1_gene"]
    txs = []
    for t in gene.get("Transcript", []):
        if t.get("biotype") != "protein_coding":
            continue
        exons = sorted(t.get("Exon", []), key=lambda e: e["start"])
        translation = t.get("Translation", {})
        txs.append({
            "id": t["id"],
            "name": t.get("display_name", t["id"]),
            "is_canonical": bool(t.get("is_canonical")),
            "start": t["start"],
            "end": t["end"],
            "exons": [{"start": e["start"], "end": e["end"]} for e in exons],
            "cds_start": translation.get("start"),
            "cds_end": translation.get("end"),
        })
    return {
        "gene_id": gene["id"],
        "gene_name": gene["display_name"],
        "chr": gene["seq_region_name"],
        "start": gene["start"],
        "end": gene["end"],
        "strand": gene["strand"],
        "transcripts": txs,
    }

# ---------------------------------------------------------------------------
# Phase 5: Build overlap summary
# ---------------------------------------------------------------------------

def build_overlap(data):
    df = data["rai1_as1_overlap"]
    if df.empty:
        return {"rai1": "", "rai1as1": "", "overlap_bp": 0, "category": "unknown",
                "quality_score": 0, "promoter": False, "exon": False}
    ov = df.iloc[0]
    return {
        "rai1": f"chr{ov['rai1_chr']}:{ov['rai1_start']:,}-{ov['rai1_end']:,}",
        "rai1as1": f"chr{ov['rai1_chr']}:{ov['rai1as1_start']:,}-{ov['rai1as1_end']:,}",
        "overlap_bp": int(ov["overlap_bp"]),
        "category": ov["overlap_category"],
        "quality_score": float(ov["overlap_quality_score"]),
        "promoter": bool(ov["overlaps_promoter_bool"]),
        "exon": bool(ov["overlaps_exon_bool"]),
    }

# ---------------------------------------------------------------------------
# Phase 6: Build near-cognate uORFs
# ---------------------------------------------------------------------------

def build_structure(data):
    f = DER / "rai1_utr_structure.tsv"
    if not f.exists():
        return []
    df = pd.read_csv(f, sep="\t")
    if df.empty:
        return []
    # Return top 5 most stable windows as ASO structure targets
    top = df.nsmallest(5, "mfe")
    results = []
    for _, r in top.iterrows():
        results.append({
            "window_start": int(r["window_start"]),
            "window_end": int(r["window_end"]),
            "genomic_start": int(r["genomic_start"]) if pd.notna(r["genomic_start"]) else 0,
            "genomic_end": int(r["genomic_end"]) if pd.notna(r["genomic_end"]) else 0,
            "mfe": float(r["mfe"]),
            "mfe_per_nt": float(r["mfe_per_nt"]),
            "gc_content": float(r["gc_content"]),
            "structure": str(r["structure"]),
        })
    return results

def build_mirna(data):
    f = DER / "rai1_mirna_targets.tsv"
    if not f.exists():
        return []
    df = pd.read_csv(f, sep="\t")
    results = []
    for _, r in df.iterrows():
        results.append({
            "mirna_family": r["mirna_family"],
            "representative": r["representative_mirna"],
            "context_score": float(r["context_score"]),
            "conserved_sites": int(r["conserved_sites"]),
            "conserved_8mer": int(r["conserved_8mer"]),
            "aggregate_pct": float(r["aggregate_pct"]) if pd.notna(r["aggregate_pct"]) else 0,
        })
    return results

# Also get all structure windows for the profile chart
def build_structure_profile(data):
    f = DER / "rai1_utr_structure.tsv"
    if not f.exists():
        return []
    df = pd.read_csv(f, sep="\t")
    return [{"pos": int(r["window_start"]), "mfe": float(r["mfe"])} for _, r in df.iterrows()]

def build_g4(data):
    f = DER / "rai1_g4_motifs.tsv"
    if not f.exists() or f.stat().st_size <= 1:
        return []
    df = pd.read_csv(f, sep="\t")
    return [row.to_dict() for _, row in df.iterrows()]

def build_3utr_structure(data):
    f = DER / "rai1_3utr_structure.tsv"
    if not f.exists():
        return {"top": [], "profile": []}
    df = pd.read_csv(f, sep="\t")
    top = df.nsmallest(5, "mfe")
    return {
        "top": [{"pos": f"{int(r['window_start'])}-{int(r['window_end'])}", "mfe": float(r["mfe"]),
                 "mfe_per_nt": float(r["mfe_per_nt"]), "gc": float(r["gc_content"]),
                 "genomic_start": int(r["genomic_start"])} for _, r in top.iterrows()],
        "profile": [{"pos": int(r["window_start"]), "mfe": float(r["mfe"])} for _, r in df.iterrows()],
    }

def build_rbp(data):
    f = DER / "rai1_rbp_sites.tsv"
    if not f.exists():
        return {"total": 0, "utr_total": 0, "top_5utr": [], "top_3utr": []}
    df = pd.read_csv(f, sep="\t")
    utr = df[df["location"].isin(["5UTR_region", "3UTR_region"])]
    result = {"total": len(df), "utr_total": len(utr), "top_5utr": [], "top_3utr": []}
    for loc, key in [("5UTR_region", "top_5utr"), ("3UTR_region", "top_3utr")]:
        sub = utr[utr["location"] == loc]
        if not sub.empty:
            top = sub.groupby("rbp_name").agg(count=("score","count"), max_score=("score","max")).reset_index()
            top = top.sort_values("max_score", ascending=False).head(10)
            result[key] = [{"name": r["rbp_name"], "count": int(r["count"]), "score": int(r["max_score"])} for _, r in top.iterrows()]
    return result

def build_polya(data):
    f = DER / "rai1_polya.tsv"
    if not f.exists():
        return []
    df = pd.read_csv(f, sep="\t")
    return [{"signal": r["signal"], "type": r["signal_type"], "utr_pos": int(r["utr_position"]),
             "dist_from_end": int(r["dist_from_utr_end"]), "active": bool(r["is_near_end"])} for _, r in df.iterrows()]

def build_aso_candidates(data):
    f = DER / "rai1_aso_candidates.tsv"
    if not f.exists():
        return {}
    df = pd.read_csv(f, sep="\t")
    if df.empty:
        return {}
    result = {}
    for region in df["region"].unique():
        rdf = df[df["region"] == region].copy().drop_duplicates(subset=["aso_seq"])
        clean = rdf[rdf["flags"].isna() | (rdf["flags"] == "")]
        flagged = rdf[~(rdf["flags"].isna() | (rdf["flags"] == ""))]
        # Return ALL candidates (deduplicated), clean first then flagged, sorted by quality
        all_candidates = []
        for _, r in clean.sort_values("quality_score", ascending=False).iterrows():
            all_candidates.append(_aso_row(r, False))
        for _, r in flagged.sort_values("quality_score", ascending=False).iterrows():
            all_candidates.append(_aso_row(r, True))
        result[region] = {
            "description": rdf.iloc[0]["region_description"],
            "utr_type": rdf.iloc[0]["utr_type"],
            "total": len(rdf),
            "clean": len(clean),
            "flagged": len(flagged),
            "candidates": all_candidates,
        }
    return result

def _aso_row(r, is_flagged):
    return {
        "seq": r["aso_seq"],
        "length": int(r["aso_length"]),
        "gc": float(r["gc_percent"]),
        "duplex_mfe": float(r["duplex_mfe"]),
        "self_mfe": float(r["self_fold_mfe"]),
        "quality": float(r["quality_score"]),
        "flags": str(r["flags"]) if pd.notna(r["flags"]) and r["flags"] else "",
        "is_flagged": is_flagged,
        "target_seq": r["target_seq"],
    }

def build_near_cognate(data):
    f = DER / "rai1_near_cognate_uorfs.tsv"
    if not f.exists():
        return []
    df = pd.read_csv(f, sep="\t")
    if df.empty:
        return []
    # Only return validated ones
    validated = df[df["ribo_validated"] == True].copy()
    results = []
    for _, r in validated.iterrows():
        results.append({
            "start_codon": r["start_codon"],
            "utr_pos": int(r["utr_pos"]),
            "stop_pos": int(r["stop_pos"]),
            "aa_len": int(r["aa_len"]),
            "uaug_genomic": int(r["uaug_genomic"]),
            "ribo_peak": float(r.get("ribo_peak_mean", 0)),
            "ribo_background": float(r.get("ribo_background_mean", 0)),
            "ribo_enrichment": float(r["ribo_enrichment"]),
            "target_type": "near_cognate_uORF",
            "ucsc_link": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17:{int(r['uaug_genomic'])-200}-{int(r['uaug_genomic'])+200}",
        })
    return results

# ---------------------------------------------------------------------------
# Phase 7: Build NMD analysis
# ---------------------------------------------------------------------------

def build_nmd(data):
    f = DER / "rai1_nmd_analysis.tsv"
    if not f.exists():
        return []
    df = pd.read_csv(f, sep="\t")
    results = []
    for _, r in df.iterrows():
        results.append({
            "transcript_id": r["transcript_id"],
            "display_name": r["display_name"],
            "biotype": r["biotype"],
            "is_canonical": bool(r["is_canonical"]),
            "has_translation": bool(r["has_translation"]),
            "n_exons": int(r["n_exons"]),
            "ptc_triggers_nmd": bool(r["ptc_triggers_nmd"]),
            "ptc_distance": r["ptc_distance_to_last_ejc"] if pd.notna(r["ptc_distance_to_last_ejc"]) else None,
            "n_unique_exons": int(r["n_unique_exons"]),
            "unique_exons": str(r["unique_exons"]) if pd.notna(r["unique_exons"]) else "",
            "nmd_potential": r["nmd_target_potential"],
        })
    return results

# ---------------------------------------------------------------------------
# Phase 8: Render site
# ---------------------------------------------------------------------------

def render_site(targets, expression, transcripts, overlap, near_cognate, nmd,
                structure, structure_profile, mirna, g4, utr3_struct, rbp, polya, aso, all_viable):
    env = Environment(loader=FileSystemLoader(str(TMPL)), autoescape=False)

    # Write JSON data files
    (SITE / "data" / "targets.json").write_text(json.dumps(targets, indent=2))
    (SITE / "data" / "expression.json").write_text(json.dumps(expression, indent=2))
    (SITE / "data" / "transcripts.json").write_text(json.dumps(transcripts, indent=2))

    # Write additional data
    (SITE / "data" / "near_cognate.json").write_text(json.dumps(near_cognate, indent=2))
    (SITE / "data" / "nmd.json").write_text(json.dumps(nmd, indent=2))

    # Common context
    ctx = {
        "targets": targets,
        "expression": expression,
        "transcripts": transcripts,
        "overlap": overlap,
        "near_cognate": near_cognate,
        "nmd": nmd,
        "structure": structure,
        "structure_profile": structure_profile,
        "mirna": mirna,
        "structure_json": json.dumps(structure),
        "structure_profile_json": json.dumps(structure_profile),
        "mirna_json": json.dumps(mirna),
        "g4": g4,
        "utr3_struct": utr3_struct,
        "rbp": rbp,
        "polya": polya,
        "utr3_profile_json": json.dumps(utr3_struct.get("profile", [])),
        "aso": aso,
        "all_viable": all_viable,
        "targets_json": json.dumps(targets),
        "expression_json": json.dumps(expression),
        "transcripts_json": json.dumps(transcripts),
        "near_cognate_json": json.dumps(near_cognate),
        "nmd_json": json.dumps(nmd),
    }

    # Render pages
    pages = {
        "index.html": "index.html",
        "methodology.html": "methodology.html",
        "about.html": "about.html",
        "compare.html": "compare.html",
    }

    for out_name, tmpl_name in pages.items():
        tmpl = env.get_template(tmpl_name)
        html = tmpl.render(**ctx)
        (SITE / out_name).write_text(html)
        print(f"  Rendered {out_name}")

    # Per-target pages
    tmpl = env.get_template("target.html")
    for t in targets:
        html = tmpl.render(target=t, **ctx)
        out = SITE / "target" / f"{t['slug']}.html"
        out.write_text(html)
        print(f"  Rendered target/{t['slug']}.html")

    # Copy downloads
    dl = SITE / "downloads"
    dl.mkdir(exist_ok=True)
    for f in DER.glob("*.tsv"):
        shutil.copy2(f, dl / f.name)
    print(f"  Copied {len(list(DER.glob('*.tsv')))} TSV files to downloads/")

def main():
    print("Building SMS ASO Target Portal...")
    data = load_data()
    print("  Data loaded OK")

    targets = build_targets(data)
    print(f"  Built {len(targets)} target records")
    for t in targets:
        print(f"    #{t['rank']} {t['id']} (composite={t['composite_score']}, verdict={t['triage_verdict']})")

    expression = build_expression(data)
    print(f"  Built {len(expression)} expression records")

    transcripts = build_transcripts(data)
    print(f"  Built transcript structure ({len(transcripts['transcripts'])} coding transcripts)")

    overlap = build_overlap(data)
    print(f"  Built overlap summary ({overlap['overlap_bp']:,} bp)")

    near_cognate = build_near_cognate(data)
    print(f"  Built {len(near_cognate)} validated near-cognate uORFs")

    nmd = build_nmd(data)
    nmd_candidates = [n for n in nmd if n["nmd_potential"] in ("high", "medium")]
    print(f"  Built NMD analysis ({len(nmd)} transcripts, {len(nmd_candidates)} with potential)")

    structure = build_structure(data)
    structure_profile = build_structure_profile(data)
    print(f"  Built {len(structure)} top structural targets, {len(structure_profile)} profile windows")

    mirna = build_mirna(data)
    print(f"  Built {len(mirna)} miRNA predictions")

    g4 = build_g4(data)
    print(f"  Built {len(g4)} G-quadruplex motifs")

    utr3_struct = build_3utr_structure(data)
    print(f"  Built 3'UTR structure: {len(utr3_struct['top'])} top windows, {len(utr3_struct['profile'])} profile")

    rbp = build_rbp(data)
    print(f"  Built RBP sites: {rbp['utr_total']} in UTRs ({rbp['total']} total)")

    polya = build_polya(data)
    print(f"  Built {len(polya)} polyA signals")

    aso = build_aso_candidates(data)
    total_aso = sum(r["total"] for r in aso.values())
    total_clean = sum(r["clean"] for r in aso.values())
    print(f"  Built {total_aso} ASO candidates ({total_clean} clean) across {len(aso)} regions")

    # Build unified target list across all mechanisms, sorted by tier
    all_viable = []
    for t in targets:
        tier = 1 if t["composite_score"] >= 5 else (2 if t["composite_score"] >= 3 else 3)
        all_viable.append({
            "tier": tier, "mechanism": f"uORF blocking{' (' + t['start_codon'] + ')' if t['target_type'] == 'near_cognate_uORF' else ''}",
            "target_name": t["id"], "region": f"5'UTR pos {t['uaug_pos']}",
            "key_metric": f"Ribo={int(t['ribo']['sum'])}" if t["target_type"] == "AUG_uORF" else f"{t['ribo'].get('enrichment','')}x enrich",
            "key_metric_2": f"phyloP={t['conservation'].get('uAUG_window',{}).get('phyloP','N/A')}",
            "feasibility_light": t["lights"]["designability"],
            "feasibility_text": f"{t['designability']['clean']}/{t['designability']['total']}",
            "link": f"target/{t['slug']}.html", "link_type": "internal",
            "sort_score": t["composite_score"],
        })
    for s in structure:
        all_viable.append({
            "tier": 1, "mechanism": "5'UTR structure",
            "target_name": f"Hairpin {s['window_start']}-{s['window_end']}",
            "region": f"5'UTR pos {s['window_start']}-{s['window_end']}",
            "key_metric": f"MFE={s['mfe']}", "key_metric_2": f"GC={s['gc_content']}%",
            "feasibility_light": "yellow" if s["gc_content"] > 80 else "green",
            "feasibility_text": "High GC" if s["gc_content"] > 80 else "OK",
            "link": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17:{s['genomic_start']-100}-{s['genomic_end']+100}",
            "link_type": "external",
            "sort_score": 6.0 + abs(s["mfe"]) / 100,  # Tier 1, sorted by MFE strength
        })
    for s in utr3_struct.get("top", []):
        all_viable.append({
            "tier": 2, "mechanism": "3'UTR structure",
            "target_name": f"Hairpin {s['pos']}",
            "region": f"3'UTR pos {s['pos']}",
            "key_metric": f"MFE={s['mfe']}", "key_metric_2": f"GC={s['gc']}%",
            "feasibility_light": "green" if s["gc"] < 75 else "yellow",
            "feasibility_text": "OK" if s["gc"] < 75 else "High GC",
            "link": "", "link_type": "none",
            "sort_score": 3.5 + abs(s["mfe"]) / 100,
        })
    for m in mirna:
        all_viable.append({
            "tier": 3, "mechanism": "miRNA blocking",
            "target_name": m["representative"],
            "region": "3'UTR",
            "key_metric": f"ctx++={m['context_score']}", "key_metric_2": f"{m['conserved_sites']} sites",
            "feasibility_light": "red",
            "feasibility_text": "Weak",
            "link": "", "link_type": "none",
            "sort_score": 1.0 + abs(m["context_score"]),
        })
    all_viable.sort(key=lambda x: (-x["sort_score"]))  # highest score first (tier 1 naturally on top)
    print(f"  Built {len(all_viable)} unified viable targets ({sum(1 for v in all_viable if v['tier']==1)} tier 1, {sum(1 for v in all_viable if v['tier']==2)} tier 2, {sum(1 for v in all_viable if v['tier']==3)} tier 3)")

    render_site(targets, expression, transcripts, overlap, near_cognate, nmd,
                structure, structure_profile, mirna, g4, utr3_struct, rbp, polya, aso, all_viable)
    print("Build complete. Output in site/")

if __name__ == "__main__":
    main()
