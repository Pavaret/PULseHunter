#!/usr/bin/env python3
"""
integrate_dbcan.py - merge a dbCAN3 (run_dbcan + CGCFinder) run into the
PULseHunter results and report where the two approaches agree.

Adds to final/:
  PUL_summary.tsv              + dbcan3_cgcs, dbcan3_overlap_frac, call_support
  gene_annotation_master.tsv   + dbcan3_hmmer, dbcan3_diamond, dbcan3_n_tools,
                                 dbcan3_cgc
  dbcan/dbcan_cgc_regions.tsv  every dbCAN3 CGC as one region
  dbcan/locus_concordance.tsv  every locus from either method, classified as
                               both / PULseHunter_only / dbCAN3_only
  dbcan/cazyme_concordance.tsv per-gene CAZyme agreement
  dbcan/concordance_summary.tsv headline counts
"""

import argparse
from pathlib import Path

import pandas as pd


def read_cgc_standard(path):
    """cgc_standard.out: CGC#, Gene Type, Contig ID, Protein ID, Gene Start,
    Gene Stop, Direction, Protein Family. CGC numbers restart per contig."""
    cols = ["cgc_num", "gene_type", "contig", "gene_id", "start", "end",
            "strand", "family"]
    rows = []
    if not Path(path).exists():
        return pd.DataFrame(columns=cols)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7 or p[0].startswith("CGC#") or "+++++" in line:
                continue
            p += [""] * (8 - len(p))
            rows.append(p[:8])
    df = pd.DataFrame(rows, columns=cols)
    if df.empty:
        return df
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def read_overview(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("-")
    df = df.rename(columns={"Gene ID": "gene_id", "HMMER": "dbcan3_hmmer",
                            "DIAMOND": "dbcan3_diamond",
                            "#ofTools": "dbcan3_n_tools"})
    df["dbcan3_n_tools"] = df["dbcan3_n_tools"].astype(int)
    return df[["gene_id", "dbcan3_hmmer", "dbcan3_diamond", "dbcan3_n_tools"]]


def main():
    ap = argparse.ArgumentParser(
        description="Compare PULseHunter loci with dbCAN3 CGCs")
    ap.add_argument("--final-dir", required=True)
    ap.add_argument("--dbcan-dir", required=True,
                    help="run_dbcan output directory")
    ap.add_argument("--min-overlap-genes", type=int, default=1)
    args = ap.parse_args()

    final = Path(args.final_dir)
    dbdir = Path(args.dbcan_dir)
    out = final / "dbcan"
    out.mkdir(exist_ok=True)

    summary = pd.read_csv(final / "PUL_summary.tsv", sep="\t")
    genes = pd.read_csv(final / "gene_annotation_master.tsv", sep="\t",
                        low_memory=False)
    drop = [c for c in genes.columns if c.startswith("dbcan3_")]
    genes = genes.drop(columns=drop)
    summary = summary.drop(columns=[c for c in ("dbcan3_cgcs",
                                                "dbcan3_overlap_frac",
                                                "call_support")
                                    if c in summary.columns])
    sample = summary["sample"].iloc[0] if len(summary) else ""

    # ── per-gene CAZyme calls ──────────────────────────────────────────
    genes = genes.merge(read_overview(dbdir / "overview.txt"), on="gene_id",
                        how="left")
    genes["dbcan3_n_tools"] = genes["dbcan3_n_tools"].fillna(0).astype(int)
    for c in ("dbcan3_hmmer", "dbcan3_diamond"):
        genes[c] = genes[c].fillna("-")

    # ── CGC regions ────────────────────────────────────────────────────
    cgc_genes = read_cgc_standard(dbdir / "cgc_standard.out")
    regions = []
    gene_to_cgc = {}
    if not cgc_genes.empty:
        grouped = sorted(cgc_genes.groupby(["contig", "cgc_num"], sort=False),
                         key=lambda kv: (kv[0][0], kv[1]["start"].min()))
        for n, ((contig, _), sub) in enumerate(grouped, 1):
            cgc_id = f"CGC_{n:04d}"
            types = sub["gene_type"].value_counts()
            fams = sorted(set(sub.loc[sub["gene_type"] == "CAZyme", "family"]))
            for g in sub["gene_id"]:
                gene_to_cgc[g] = cgc_id
            regions.append({
                "sample": sample, "cgc_id": cgc_id, "contig": contig,
                "start": int(sub["start"].min()), "end": int(sub["end"].max()),
                "length_bp": int(sub["end"].max() - sub["start"].min() + 1),
                "num_genes": len(sub),
                "cazyme_count": int(types.get("CAZyme", 0)),
                "tc_count": int(types.get("TC", 0)),
                "tf_count": int(types.get("TF", 0)),
                "stp_count": int(types.get("STP", 0)),
                "cazyme_families": ",".join(fams),
            })
    regions = pd.DataFrame(regions, columns=[
        "sample", "cgc_id", "contig", "start", "end", "length_bp",
        "num_genes", "cazyme_count", "tc_count", "tf_count", "stp_count",
        "cazyme_families"])
    regions.to_csv(out / "dbcan_cgc_regions.tsv", sep="\t", index=False)
    genes["dbcan3_cgc"] = genes["gene_id"].map(gene_to_cgc).fillna("")

    # ── locus concordance (shared genes) ───────────────────────────────
    in_cand = genes[genes["candidate_id"].notna() & (genes["candidate_id"] != "")]
    cand_to_cgc, cgc_to_cand, cand_frac = {}, {}, {}
    for cid, sub in in_cand.groupby("candidate_id"):
        hits = sub.loc[sub["dbcan3_cgc"] != "", "dbcan3_cgc"]
        counts = hits.value_counts()
        counts = counts[counts >= args.min_overlap_genes]
        cand_to_cgc[cid] = list(counts.index)
        cand_frac[cid] = round(len(hits) / len(sub), 2)
        for c in counts.index:
            cgc_to_cand.setdefault(c, []).append(cid)
    summary["dbcan3_cgcs"] = summary["candidate_id"].map(
        lambda c: ";".join(sorted(cand_to_cgc.get(c, []))))
    summary["dbcan3_overlap_frac"] = summary["candidate_id"].map(
        cand_frac).fillna(0.0)
    summary["call_support"] = summary["dbcan3_cgcs"].map(
        lambda s: "both" if s else "PULseHunter_only")

    conc = []
    for r in summary.itertuples():
        conc.append({
            "sample": sample, "locus": r.candidate_id,
            "category": r.call_support, "contig": r.contig,
            "start": r.start, "end": r.end, "length_bp": r.length_bp,
            "pulsehunter_confidence": r.confidence,
            "pulsehunter_type": r.pul_type,
            "susCD_pairs": r.susCD_pairs, "cazyme_count": r.cazyme_count,
            "cazyme_families": r.cazyme_families,
            "matching_locus": r.dbcan3_cgcs,
        })
    for r in regions.itertuples():
        if r.cgc_id in cgc_to_cand:
            continue
        conc.append({
            "sample": sample, "locus": r.cgc_id, "category": "dbCAN3_only",
            "contig": r.contig, "start": r.start, "end": r.end,
            "length_bp": r.length_bp, "pulsehunter_confidence": "",
            "pulsehunter_type": "", "susCD_pairs": 0,
            "cazyme_count": r.cazyme_count,
            "cazyme_families": r.cazyme_families, "matching_locus": "",
        })
    conc = pd.DataFrame(conc)
    if not conc.empty:
        conc = conc.sort_values(["contig", "start"])
    conc.to_csv(out / "locus_concordance.tsv", sep="\t", index=False)

    # ── CAZyme concordance ─────────────────────────────────────────────
    ph = genes["is_cazyme"].astype(bool)
    d1 = genes["dbcan3_n_tools"] >= 1
    d2 = genes["dbcan3_n_tools"] >= 2
    cz = genes[ph | d1][["gene_id", "contig", "start", "end", "candidate_id",
                         "cazyme_family", "cazyme_evidence", "dbcan3_hmmer",
                         "dbcan3_diamond", "dbcan3_n_tools"]].copy()
    cz["category"] = [
        "both" if a and b else ("PULseHunter_only" if a else "dbCAN3_only")
        for a, b in zip(ph[ph | d1], d1[ph | d1])]
    cz.to_csv(out / "cazyme_concordance.tsv", sep="\t", index=False)

    stats = [
        ("pulsehunter_loci", len(summary)),
        ("pulsehunter_loci_high", int((summary["confidence"] == "high").sum())),
        ("dbcan3_cgcs", len(regions)),
        ("pulsehunter_loci_supported_by_dbcan3",
         int((summary["call_support"] == "both").sum())),
        ("pulsehunter_only_loci",
         int((summary["call_support"] == "PULseHunter_only").sum())),
        ("dbcan3_cgcs_overlapping_pulsehunter", len(cgc_to_cand)),
        ("dbcan3_only_cgcs", len(regions) - len(cgc_to_cand)),
        ("cazymes_pulsehunter", int(ph.sum())),
        ("cazymes_dbcan3_any_tool", int(d1.sum())),
        ("cazymes_dbcan3_both_tools", int(d2.sum())),
        ("cazymes_both_methods", int((ph & d1).sum())),
        ("cazymes_pulsehunter_only", int((ph & ~d1).sum())),
        ("cazymes_dbcan3_only", int((~ph & d1).sum())),
    ]
    pd.DataFrame(stats, columns=["metric", "value"]).assign(
        sample=sample).to_csv(out / "concordance_summary.tsv", sep="\t",
                              index=False)
    for k, v in stats:
        print(f"{k}\t{v}")

    genes.to_csv(final / "gene_annotation_master.tsv", sep="\t", index=False)
    summary.to_csv(final / "PUL_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
