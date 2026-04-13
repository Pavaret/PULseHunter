#!/usr/bin/env python3
"""
call_puls_hmm.py – Call candidate PULs from HMM / BLAST annotations and
export wet-lab-friendly outputs.

Fixed version – key changes from original:
  1. Gene-ID reconciliation between Prodigal GFF (ID=seqnum_genenum) and
     Prodigal FASTA (contigname_genenum).  Both are now mapped to the
     canonical FASTA-header form so that HMM/BLAST joins work.
  2. Pfam accessions are normalised to unversioned form (PF00593, not
     PF00593.31) so matching works with any Pfam release.
  3. classify_suscd uses the module-level constant sets consistently.
  4. Custom-HMM hits are also aggregated (not just best-hit) so
     multi-domain proteins keep all annotations.
  5. Block extension has a hard cap to prevent runaway operon merging.
  6. meta_by_pul picks the most-informative row per PUL instead of
     whichever comes first.
  7. Diagnostic counts are printed to stderr for easier debugging.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
import re
import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Pfam domain sets – UNVERSIONED so they match any Pfam release
# ---------------------------------------------------------------------------
PFAM_SUSC_CORE  = {"PF00593", "PF07715"}
PFAM_SUSC_EXTRA = {"PF13715", "PF07660"}
PFAM_SUSD       = {"PF07980", "PF12741", "PF12771", "PF14322"}


def _strip_pfam_version(acc: str) -> str:
    """PF00593.31 -> PF00593"""
    return acc.rsplit(".", 1)[0] if re.match(r"^PF\d+\.\d+$", acc) else acc


def log(msg: str):
    print(msg, file=sys.stderr, flush=True)


# ── GFF / FASTA reading ────────────────────────────────────────────────────

def parse_attrs(attr_str):
    d = {}
    for item in attr_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def read_prodigal_gff(gff_path):
    """Read Prodigal GFF and build a canonical gene_id that matches the
    FASTA header produced by ``prodigal -a``.

    Prodigal GFF line (example):
        NZ_ABC.1  Prodigal  CDS  337  2799  266  +  0  ID=1_1;partial=10;...
    Prodigal FAA header:
        >NZ_ABC.1_1 # 337 # 2799 # 1 # ID=1_1;partial=10;...

    The GFF ``ID`` field is ``seqnum_genenum`` (e.g. 1_1 for the first
    gene on the first sequence Prodigal encountered).  The FASTA id is
    ``contigname_genenum``.  We need the FASTA form because that is what
    HMMER and BLAST will report as the query/target name.

    Strategy: extract the gene-number suffix from the GFF ID and append
    it to the contig name.  Additionally we build a lookup from
    (contig, start, end) so we can reconcile any edge-case mismatches.
    """
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = parse_attrs(parts[8])
            raw_id = attrs.get("ID", "")
            contig = parts[0]

            # Extract the gene number from the Prodigal ID (e.g. "1_1" -> "1")
            if "_" in raw_id:
                gene_num = raw_id.rsplit("_", 1)[1]
            else:
                gene_num = raw_id

            # Build the canonical FASTA-style gene id
            gene_id = f"{contig}_{gene_num}"

            rows.append({
                "contig": contig,
                "source": parts[1],
                "feature": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": parts[5],
                "strand": parts[6],
                "phase": parts[7],
                "attributes": parts[8],
                "gene_id": gene_id,
                "prodigal_raw_id": raw_id,
            })

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit(f"No CDS entries parsed from {gff_path}")
    df = df.sort_values(["contig", "start", "end"]).reset_index(drop=True)
    df["gene_order"] = df.groupby("contig").cumcount() + 1
    df["gene_len_nt"] = df["end"] - df["start"] + 1
    return df


def read_fasta_dict(path):
    return {rec.id.split()[0]: rec for rec in SeqIO.parse(path, "fasta")}


def _reconcile_gene_ids(genes, proteins):
    """If the GFF-derived gene_ids don't overlap with the FASTA keys,
    fall back to a coordinate-based reconciliation.

    This handles non-standard Prodigal invocations or post-processed
    GFF files where the ID scheme has been altered.
    """
    gff_ids = set(genes["gene_id"])
    fasta_ids = set(proteins.keys())
    overlap = gff_ids & fasta_ids
    if len(overlap) >= 0.5 * len(gff_ids):
        # Good enough – the canonical join already works
        return genes

    log(f"WARNING: only {len(overlap)}/{len(gff_ids)} GFF gene_ids match "
        f"FASTA headers.  Attempting coordinate-based reconciliation …")

    # Build a lookup from Prodigal FASTA headers:
    # >contig_N # start # end # strand # attrs
    coord_to_fasta = {}
    for fid in fasta_ids:
        # Prodigal header after split()[0] is e.g. NZ_ABC.1_3
        # The full description has "# start # end # ..."
        rec = proteins[fid]
        desc_parts = rec.description.split(" # ")
        if len(desc_parts) >= 3:
            try:
                s, e = int(desc_parts[1]), int(desc_parts[2])
                coord_to_fasta[(s, e)] = fid
            except ValueError:
                pass

    new_ids = []
    for _, row in genes.iterrows():
        key = (row["start"], row["end"])
        new_ids.append(coord_to_fasta.get(key, row["gene_id"]))
    genes = genes.copy()
    genes["gene_id"] = new_ids

    overlap2 = set(genes["gene_id"]) & fasta_ids
    log(f"  After reconciliation: {len(overlap2)}/{len(genes)} gene_ids match.")
    return genes


# ── HMMER / BLAST parsing ──────────────────────────────────────────────────

def parse_hmmer_domtblout(path, evalue_cutoff=None, cov_cutoff=None):
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=[
            "target_name", "target_len", "query_name", "query_len",
            "full_evalue", "full_score",
            "i_evalue", "domain_score",
            "hmm_from", "hmm_to", "ali_from", "ali_to", "hmm_cov"
        ])

    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(maxsplit=22)
            if len(parts) < 22:
                continue
            try:
                target_name = parts[0]
                target_len  = int(parts[2])
                query_name  = parts[3]
                query_len   = int(parts[5])
                full_evalue = float(parts[6])
                full_score  = float(parts[7])
                i_evalue    = float(parts[12])
                domain_score = float(parts[13])
                hmm_from    = int(parts[15])
                hmm_to      = int(parts[16])
                ali_from    = int(parts[17])
                ali_to      = int(parts[18])
            except ValueError:
                continue
            hmm_cov = (hmm_to - hmm_from + 1) / query_len if query_len else 0
            rows.append({
                "target_name": target_name,
                "target_len": target_len,
                "query_name": query_name,
                "query_len": query_len,
                "full_evalue": full_evalue,
                "full_score": full_score,
                "i_evalue": i_evalue,
                "domain_score": domain_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "hmm_cov": hmm_cov,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    if evalue_cutoff is not None:
        df = df[df["i_evalue"] <= evalue_cutoff].copy()
    if cov_cutoff is not None:
        df = df[df["hmm_cov"] >= cov_cutoff].copy()
    return df


def best_hits(df):
    if df.empty:
        return df
    return (df.sort_values(
                ["target_name", "i_evalue", "full_score", "hmm_cov"],
                ascending=[True, True, False, False])
              .drop_duplicates(subset=["target_name"], keep="first")
              .reset_index(drop=True))


def aggregate_hits(df, field):
    if df.empty:
        return pd.DataFrame(columns=["target_name", field])
    agg = (df.sort_values(
                ["target_name", "i_evalue", "full_score"],
                ascending=[True, True, False])
             .groupby("target_name")["query_name"]
             .apply(lambda x: ";".join(dict.fromkeys(map(str, x))))
             .reset_index(name=field))
    return agg


# ── SusC / SusD classification (Pfam-based) ────────────────────────────────

def classify_suscd(pfam_hits):
    """Classify target proteins as SusC-like or SusD-like based on Pfam
    domain architecture.

    Domain names are normalised to unversioned accessions before
    matching so the logic works regardless of which Pfam release was
    used to build the HMM library.
    """
    if pfam_hits.empty:
        return pd.DataFrame(columns=[
            "target_name", "susC_like_pfam", "susD_like_pfam", "suscd_domains"])
    out = []
    for gene, sub in pfam_hits.groupby("target_name"):
        # Normalise to unversioned accessions
        raw_domains = set(sub["query_name"].astype(str))
        domains = {_strip_pfam_version(d) for d in raw_domains}

        susC = (
            ("PF07715" in domains
             and bool(domains & (PFAM_SUSC_CORE | PFAM_SUSC_EXTRA)))
            or ({"PF00593", "PF13715"} <= domains)
        )
        susD = bool(domains & PFAM_SUSD)

        out.append({
            "target_name": gene,
            "susC_like_pfam": susC,
            "susD_like_pfam": susD,
            "suscd_domains": ";".join(sorted(raw_domains)),
        })
    return pd.DataFrame(out)


# ── BLAST parsing ───────────────────────────────────────────────────────────

def parse_blast(path):
    cols = [
        "qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"
    ]
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    return pd.read_csv(path, sep="\t", names=cols, header=None)


def tc_annotation_is_susc(x):
    return bool(re.search(r"(^|[^0-9])1\.B\.14(\.|$|[^0-9])", str(x)))


def tc_annotation_is_susd(x):
    return bool(re.search(r"(^|[^0-9])8\.A\.46(\.|$|[^0-9])", str(x)))


# ── Annotation integration ─────────────────────────────────────────────────

def join_gene_annotations(genes, proteins, custom_hits, cluster_manifest,
                          pfam_hits, dbcan_hits, blast_hits, ref_table):
    genes = genes.copy()
    genes["protein_len_aa"] = genes["gene_id"].map(
        lambda x: len(proteins[x].seq) if x in proteins else pd.NA)
    genes["short_gene_lt_60aa"] = genes["protein_len_aa"].fillna(0).astype(float) < 60

    # ── custom HMM hits ────────────────────────────────────────────────
    if not custom_hits.empty:
        cm = pd.read_csv(cluster_manifest, sep="\t", dtype=str).fillna("")
        # Best single hit per protein
        best_custom = (best_hits(custom_hits)
                       .merge(cm, left_on="query_name",
                              right_on="cluster_id", how="left"))
        best_custom = best_custom.rename(columns={
            "target_name":        "gene_id",
            "query_name":         "best_custom_cluster",
            "full_score":         "best_custom_score",
            "i_evalue":           "best_custom_ievalue",
            "hmm_cov":            "best_custom_cov",
            "majority_class":     "best_custom_majority_class",
            "dominant_annotation": "best_custom_annotation",
            "dominant_pul":       "best_custom_pul",
        })
        keep = [c for c in [
            "gene_id", "best_custom_cluster", "best_custom_score",
            "best_custom_ievalue", "best_custom_cov",
            "best_custom_majority_class", "best_custom_annotation",
            "best_custom_pul"
        ] if c in best_custom.columns]
        genes = genes.merge(best_custom[keep], on="gene_id", how="left")

        # Aggregate: all custom families for each protein
        custom_all = (aggregate_hits(custom_hits, "custom_families")
                      .rename(columns={"target_name": "gene_id"}))
        genes = genes.merge(custom_all, on="gene_id", how="left")
    else:
        for c in ["best_custom_cluster", "best_custom_score",
                  "best_custom_ievalue", "best_custom_cov",
                  "best_custom_majority_class", "best_custom_annotation",
                  "best_custom_pul", "custom_families"]:
            genes[c] = pd.NA

    # ── Pfam SusC/SusD ─────────────────────────────────────────────────
    suscd = classify_suscd(pfam_hits)
    if not suscd.empty:
        suscd = suscd.rename(columns={"target_name": "gene_id"})
        genes = genes.merge(suscd, on="gene_id", how="left")
    else:
        genes["susC_like_pfam"] = False
        genes["susD_like_pfam"] = False
        genes["suscd_domains"] = pd.NA
    genes["susC_like_pfam"] = genes["susC_like_pfam"].fillna(False)
    genes["susD_like_pfam"] = genes["susD_like_pfam"].fillna(False)

    # ── dbCAN hits ─────────────────────────────────────────────────────
    if not dbcan_hits.empty:
        dbcan_best = best_hits(dbcan_hits).rename(columns={
            "target_name": "gene_id",
            "query_name":  "best_dbcan_family",
            "full_score":  "best_dbcan_score",
            "i_evalue":    "best_dbcan_ievalue",
            "hmm_cov":     "best_dbcan_cov",
        })
        genes = genes.merge(
            dbcan_best[["gene_id", "best_dbcan_family",
                        "best_dbcan_score", "best_dbcan_ievalue",
                        "best_dbcan_cov"]],
            on="gene_id", how="left")
        dbcan_all = (aggregate_hits(dbcan_hits, "dbcan_families")
                     .rename(columns={"target_name": "gene_id"}))
        genes = genes.merge(dbcan_all, on="gene_id", how="left")
    else:
        for c in ["best_dbcan_family", "best_dbcan_score",
                  "best_dbcan_ievalue", "best_dbcan_cov", "dbcan_families"]:
            genes[c] = pd.NA

    # ── BLAST hits ─────────────────────────────────────────────────────
    if not blast_hits.empty:
        ref = ref_table.copy()
        top_blast = (
            blast_hits
            .sort_values(["qseqid", "evalue", "bitscore"],
                         ascending=[True, True, False])
            .drop_duplicates("qseqid")
            .merge(ref[["full_id", "pul_id", "ref_class", "ref_annotation"]],
                   left_on="sseqid", right_on="full_id", how="left"))
        top_blast = top_blast.rename(columns={
            "qseqid":   "gene_id",
            "sseqid":   "top_blast_subject",
            "pul_id":   "top_blast_pul",
            "ref_class": "top_blast_class",
            "ref_annotation": "top_blast_annotation",
            "bitscore": "top_blast_bitscore",
            "evalue":   "top_blast_evalue",
            "qcovs":    "top_blast_qcovs",
        })
        keep = [c for c in [
            "gene_id", "top_blast_subject", "top_blast_pul",
            "top_blast_class", "top_blast_annotation",
            "top_blast_bitscore", "top_blast_evalue", "top_blast_qcovs"
        ] if c in top_blast.columns]
        genes = genes.merge(top_blast[keep], on="gene_id", how="left")
    else:
        for c in ["top_blast_subject", "top_blast_pul", "top_blast_class",
                  "top_blast_annotation", "top_blast_bitscore",
                  "top_blast_evalue", "top_blast_qcovs"]:
            genes[c] = pd.NA

    # ── Composite boolean flags ────────────────────────────────────────
    genes["custom_susC_like"] = genes["best_custom_majority_class"].eq("SusC_like")
    genes["custom_susD_like"] = genes["best_custom_majority_class"].eq("SusD_like")
    genes["blast_susC_like"]  = (genes["top_blast_class"].eq("TC")
                                 & genes["top_blast_annotation"].map(tc_annotation_is_susc))
    genes["blast_susD_like"]  = (genes["top_blast_class"].eq("TC")
                                 & genes["top_blast_annotation"].map(tc_annotation_is_susd))

    genes["susC_like"] = genes[["susC_like_pfam", "custom_susC_like",
                                "blast_susC_like"]].any(axis=1)
    genes["susD_like"] = genes[["susD_like_pfam", "custom_susD_like",
                                "blast_susD_like"]].any(axis=1)

    genes["is_cazyme"] = (
        genes["best_dbcan_family"].notna()
        | genes["best_custom_majority_class"].eq("CAZyme")
        | genes["top_blast_class"].eq("CAZyme")
    )
    genes["is_tc"] = (
        genes["best_custom_majority_class"].isin(
            ["TC_other", "SusC_like", "SusD_like"])
        | genes["top_blast_class"].eq("TC")
        | genes["susC_like"]
        | genes["susD_like"]
    )
    genes["is_tf"]  = (genes["best_custom_majority_class"].eq("TF")
                       | genes["top_blast_class"].eq("TF"))
    genes["is_stp"] = (genes["best_custom_majority_class"].eq("STP")
                       | genes["top_blast_class"].eq("STP"))

    genes["has_reference_hmm_hit"] = genes["best_custom_cluster"].notna()
    genes["signal_gene"] = genes[[
        "susC_like", "susD_like", "is_cazyme",
        "has_reference_hmm_hit", "is_tf", "is_stp"
    ]].any(axis=1)
    genes["signal_weight"] = (
        genes["susC_like"].astype(int) * 3
        + genes["susD_like"].astype(int) * 3
        + genes["is_cazyme"].astype(int) * 2
        + genes["has_reference_hmm_hit"].astype(int)
        + genes["is_tf"].astype(int)
        + genes["is_stp"].astype(int)
    )
    return genes


# ── Candidate block building ───────────────────────────────────────────────

def build_candidate_blocks(genes, max_signal_gene_gap=5, flank_genes=2,
                           max_intergenic=500, max_extension_genes=10):
    """Cluster signal genes into candidate PUL blocks, extend into flanks,
    and merge overlapping blocks.

    ``max_extension_genes`` caps how many non-signal genes the boundary-
    extension loop can pull in on each side, preventing runaway merging
    on densely packed genomes.
    """
    blocks = []
    for contig, sub in genes.groupby("contig", sort=False):
        sub = sub.sort_values("gene_order").reset_index(drop=True)
        signal_idx = sub.index[sub["signal_gene"]].tolist()
        if not signal_idx:
            continue
        current = [signal_idx[0]]
        for idx in signal_idx[1:]:
            prev = current[-1]
            gap_genes = idx - prev - 1
            bp_gap = max(0, int(sub.loc[idx, "start"])
                            - int(sub.loc[prev, "end"]) - 1)
            if gap_genes <= max_signal_gene_gap or bp_gap <= 12000:
                current.append(idx)
            else:
                blocks.append((contig, min(current), max(current)))
                current = [idx]
        blocks.append((contig, min(current), max(current)))

    extended = []
    for contig, left, right in blocks:
        sub = genes[genes["contig"] == contig].sort_values(
            "gene_order").reset_index(drop=True)
        orig_left, orig_right = left, right
        left  = max(0, left - flank_genes)
        right = min(len(sub) - 1, right + flank_genes)

        changed = True
        while changed:
            changed = False
            # Enforce extension cap
            if (orig_left - left) >= max_extension_genes:
                pass  # stop extending left
            elif left > 0:
                gap = max(0, int(sub.loc[left, "start"])
                             - int(sub.loc[left - 1, "end"]) - 1)
                same_strand = (sub.loc[left, "strand"]
                               == sub.loc[left - 1, "strand"])
                not_short = not bool(sub.loc[left - 1, "short_gene_lt_60aa"])
                if gap <= max_intergenic and same_strand and not_short:
                    left -= 1
                    changed = True

            if (right - orig_right) >= max_extension_genes:
                pass  # stop extending right
            elif right < len(sub) - 1:
                gap = max(0, int(sub.loc[right + 1, "start"])
                             - int(sub.loc[right, "end"]) - 1)
                same_strand = (sub.loc[right, "strand"]
                               == sub.loc[right + 1, "strand"])
                not_short = not bool(
                    sub.loc[right + 1, "short_gene_lt_60aa"])
                if gap <= max_intergenic and same_strand and not_short:
                    right += 1
                    changed = True
        extended.append((contig, left, right))

    # Merge overlapping blocks on the same contig
    merged = []
    by_contig = defaultdict(list)
    for contig, left, right in extended:
        by_contig[contig].append((left, right))
    for contig, intervals in by_contig.items():
        intervals.sort()
        cur_l, cur_r = intervals[0]
        for l, r in intervals[1:]:
            if l <= cur_r + 1:
                cur_r = max(cur_r, r)
            else:
                merged.append((contig, cur_l, cur_r))
                cur_l, cur_r = l, r
        merged.append((contig, cur_l, cur_r))
    return merged


# ── Block scoring ───────────────────────────────────────────────────────────

def score_block(block_df, meta_by_pul):
    counts = {
        "num_genes":          len(block_df),
        "susC_count":         int(block_df["susC_like"].sum()),
        "susD_count":         int(block_df["susD_like"].sum()),
        "cazyme_count":       int(block_df["is_cazyme"].sum()),
        "tc_count":           int(block_df["is_tc"].sum()),
        "tf_count":           int(block_df["is_tf"].sum()),
        "stp_count":          int(block_df["is_stp"].sum()),
        "ref_hmm_gene_count": int(block_df["has_reference_hmm_hit"].sum()),
    }

    # Vote for best reference PUL via BLAST bitscore sums
    blast_votes = (
        block_df
        .dropna(subset=["top_blast_pul", "top_blast_bitscore"])
        .groupby("top_blast_pul")["top_blast_bitscore"].sum()
        .sort_values(ascending=False))
    if not blast_votes.empty:
        best_pul       = blast_votes.index[0]
        best_pul_score = float(blast_votes.iloc[0])
        hit_method     = "blastp"
    else:
        custom_votes = block_df["best_custom_pul"].dropna().astype(str)
        custom_votes = custom_votes[custom_votes != ""]
        if len(custom_votes) > 0:
            best_pul       = custom_votes.value_counts().index[0]
            best_pul_score = float(custom_votes.value_counts().iloc[0])
            hit_method     = "custom_hmm_vote"
        else:
            best_pul       = ""
            best_pul_score = 0.0
            hit_method     = ""

    substrate    = ""
    ref_organism = ""
    if best_pul and best_pul in meta_by_pul:
        row = meta_by_pul[best_pul]
        substrate    = row.get("substrate_final", "") or row.get("Substrate", "")
        ref_organism = row.get("organism_name", "") or row.get("Organism", "")

    # Confidence tiering
    if (counts["susC_count"] >= 1 and counts["susD_count"] >= 1
            and counts["cazyme_count"] >= 1
            and counts["ref_hmm_gene_count"] >= 2):
        confidence = "high"
    elif (counts["susC_count"] >= 1 and counts["susD_count"] >= 1
          and counts["cazyme_count"] >= 1):
        confidence = "medium"
    elif (counts["cazyme_count"] >= 2
          and (counts["tc_count"] + counts["tf_count"]
               + counts["stp_count"]) >= 1):
        confidence = "low"
    elif (counts["cazyme_count"] >= 1
          and counts["ref_hmm_gene_count"] >= 1):
        confidence = "tentative"
    else:
        confidence = "discard"

    return counts, best_pul, best_pul_score, hit_method, substrate, ref_organism, confidence


# ── Output writing ──────────────────────────────────────────────────────────

def write_candidate_outputs(candidate_id, block_df, genome, proteins, cds,
                            outdir):
    cdir = outdir / candidate_id
    cdir.mkdir(parents=True, exist_ok=True)
    contig       = block_df.iloc[0]["contig"]
    region_start = int(block_df["start"].min())
    region_end   = int(block_df["end"].max())
    region_seq   = genome[contig].seq[region_start - 1 : region_end]

    region_rec = genome[contig][:0]
    region_rec.id = f"{candidate_id}|{contig}:{region_start}-{region_end}"
    region_rec.description = ""
    region_rec.seq = region_seq
    SeqIO.write(region_rec, cdir / f"{candidate_id}.region.fna", "fasta")

    prot_recs = [proteins[g] for g in block_df["gene_id"] if g in proteins]
    cds_recs  = [cds[g] for g in block_df["gene_id"] if g in cds]
    if prot_recs:
        SeqIO.write(prot_recs, cdir / f"{candidate_id}.proteins.faa", "fasta")
    if cds_recs:
        SeqIO.write(cds_recs, cdir / f"{candidate_id}.cds.fna", "fasta")

    block_df.to_csv(cdir / f"{candidate_id}.genes.tsv", sep="\t", index=False)

    with open(cdir / f"{candidate_id}.gff3", "w") as fh:
        fh.write("##gff-version 3\n")
        for _, row in block_df.iterrows():
            attrs = f"ID={row['gene_id']};Name={row['gene_id']}"
            fh.write("\t".join([
                str(row["contig"]),
                str(row["source"]),
                str(row["feature"]),
                str(row["start"]),
                str(row["end"]),
                ".",
                str(row["strand"]),
                str(row["phase"]),
                attrs,
            ]) + "\n")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Call candidate PULs from HMM/BLAST annotations "
                    "and export wet-lab-friendly outputs")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--reference-table", required=True)
    ap.add_argument("--cluster-manifest", required=True)
    ap.add_argument("--custom-domtblout", required=True)
    ap.add_argument("--suscd-domtblout", default=None)
    ap.add_argument("--dbcan-domtblout", default=None)
    ap.add_argument("--blast-tsv", default=None)
    ap.add_argument("--custom-ievalue", type=float, default=1e-5)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    candidate_dir = outdir / "candidate_PULs"
    candidate_dir.mkdir(parents=True, exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────
    genes    = read_prodigal_gff(args.gff)
    proteins = read_fasta_dict(args.proteins)
    cds      = read_fasta_dict(args.cds)
    genome   = read_fasta_dict(args.genome)

    # Reconcile gene IDs if needed
    genes = _reconcile_gene_ids(genes, proteins)

    ref_table = pd.read_csv(args.reference_table, sep="\t", dtype=str).fillna("")

    # Pick the most informative row per PUL for substrate/organism lookup
    # (most non-empty fields wins)
    def _row_info_score(row):
        return sum(1 for v in row if str(v).strip())
    meta_by_pul = {}
    for _, r in ref_table.iterrows():
        pid = r["pul_id"]
        if pid not in meta_by_pul or _row_info_score(r) > _row_info_score(meta_by_pul[pid]):
            meta_by_pul[pid] = r

    # ── Parse search results ───────────────────────────────────────────
    custom_hits = parse_hmmer_domtblout(
        args.custom_domtblout, evalue_cutoff=args.custom_ievalue, cov_cutoff=0.40)
    pfam_hits = parse_hmmer_domtblout(
        args.suscd_domtblout, evalue_cutoff=None, cov_cutoff=0.30)
    dbcan_hits = parse_hmmer_domtblout(
        args.dbcan_domtblout, evalue_cutoff=1e-15, cov_cutoff=0.35)
    blast_hits = parse_blast(args.blast_tsv)

    # ── Diagnostics ────────────────────────────────────────────────────
    log(f"Genes from GFF:         {len(genes)}")
    log(f"Proteins from FASTA:    {len(proteins)}")
    log(f"Gene-FASTA overlap:     "
        f"{len(set(genes['gene_id']) & set(proteins.keys()))}")
    log(f"Custom HMM hits (filt): {len(custom_hits)}")
    log(f"Pfam SusCD hits (filt): {len(pfam_hits)}")
    log(f"dbCAN hits (filt):      {len(dbcan_hits)}")
    log(f"BLAST hits:             {len(blast_hits)}")

    # ── Annotate ───────────────────────────────────────────────────────
    genes = join_gene_annotations(
        genes, proteins, custom_hits, args.cluster_manifest,
        pfam_hits, dbcan_hits, blast_hits, ref_table)
    genes.to_csv(outdir / "gene_annotation_master.tsv", sep="\t", index=False)

    n_signal = int(genes["signal_gene"].sum())
    log(f"Signal genes:           {n_signal}")
    log(f"  susC_like:            {int(genes['susC_like'].sum())}")
    log(f"  susD_like:            {int(genes['susD_like'].sum())}")
    log(f"  is_cazyme:            {int(genes['is_cazyme'].sum())}")
    log(f"  has_reference_hmm:    {int(genes['has_reference_hmm_hit'].sum())}")
    log(f"  is_tf:                {int(genes['is_tf'].sum())}")
    log(f"  is_stp:               {int(genes['is_stp'].sum())}")

    if n_signal == 0:
        log("WARNING: No signal genes found – no candidate PULs will be "
            "produced.  Check that gene IDs match between GFF, FASTA, and "
            "HMMER/BLAST outputs.")

    # ── Build and score blocks ─────────────────────────────────────────
    blocks = build_candidate_blocks(genes)
    log(f"Raw candidate blocks:   {len(blocks)}")

    summary_rows = []
    bed_rows = []
    kept = 0
    for contig, left_idx, right_idx in blocks:
        sub = (genes[genes["contig"] == contig]
               .sort_values("gene_order").reset_index(drop=True))
        block_df = sub.loc[left_idx:right_idx].copy().reset_index(drop=True)
        (counts, best_pul, best_pul_score, hit_method,
         substrate, ref_organism, confidence) = score_block(block_df, meta_by_pul)
        if confidence == "discard":
            continue
        kept += 1
        candidate_id = f"PULcand_{kept:04d}"
        region_start = int(block_df["start"].min())
        region_end   = int(block_df["end"].max())
        length_bp    = region_end - region_start + 1

        summary_rows.append({
            "candidate_id":            candidate_id,
            "contig":                  contig,
            "start":                   region_start,
            "end":                     region_end,
            "length_bp":               length_bp,
            **counts,
            "best_reference_pul":      best_pul,
            "best_reference_pul_score": best_pul_score,
            "best_reference_method":   hit_method,
            "predicted_substrate":     substrate,
            "reference_organism":      ref_organism,
            "confidence":              confidence,
        })
        bed_rows.append([contig, region_start - 1, region_end,
                         candidate_id, counts["cazyme_count"], "."])
        write_candidate_outputs(candidate_id, block_df, genome,
                                proteins, cds, candidate_dir)

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        conf_order = pd.Categorical(
            summary["confidence"],
            categories=["high", "medium", "low", "tentative"],
            ordered=True)
        summary = (summary
                   .assign(_conf=conf_order)
                   .sort_values(["_conf", "cazyme_count",
                                 "best_reference_pul_score", "length_bp"],
                                ascending=[True, False, False, False])
                   .drop(columns=["_conf"]))
    summary.to_csv(outdir / "PUL_summary.tsv", sep="\t", index=False)
    pd.DataFrame(
        bed_rows,
        columns=["contig", "start0", "end", "candidate_id", "score", "strand"]
    ).to_csv(outdir / "PUL_regions.bed", sep="\t", header=False, index=False)

    log(f"Candidate PULs kept:    {kept}")
    if not summary.empty:
        for lvl in ["high", "medium", "low", "tentative"]:
            n = int((summary["confidence"] == lvl).sum())
            if n:
                log(f"  {lvl}: {n}")

    with open(outdir / "README_results.txt", "w") as fh:
        fh.write("PUL HMM pipeline results\n")
        fh.write("========================\n")
        fh.write(f"candidate_PUL_count\t{len(summary)}\n")
        if not summary.empty:
            for lvl in ["high", "medium", "low", "tentative"]:
                fh.write(f"{lvl}_confidence\t"
                         f"{(summary['confidence'] == lvl).sum()}\n")
        fh.write("\nMain files:\n")
        fh.write("- PUL_summary.tsv: ranked overview for all candidate loci\n")
        fh.write("- PUL_regions.bed: coordinates of loci\n")
        fh.write("- gene_annotation_master.tsv: per-gene integrated "
                 "annotations\n")
        fh.write("- candidate_PULs/: one folder per locus with region FASTA, "
                 "protein FASTA, CDS FASTA, GFF3, and gene table\n")
        fh.write("\nInterpretation:\n")
        fh.write("high      = susC + susD + CAZyme + >=2 genes with "
                 "reference-HMM support\n")
        fh.write("medium    = susC + susD + CAZyme but weaker reference "
                 "support\n")
        fh.write("low       = CAZyme-rich neighborhood with accessory "
                 "support but incomplete hallmark signature\n")
        fh.write("tentative = at least 1 CAZyme + 1 reference HMM hit, "
                 "but missing canonical SusC/SusD pair\n")


if __name__ == "__main__":
    main()
    