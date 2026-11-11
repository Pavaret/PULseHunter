#!/usr/bin/env python3
"""
call_puls_hmm.py - integrate HMM / BLAST evidence per gene, call candidate
PULs from genomic context and export wet-lab-friendly outputs.

Calling model (v2)
------------------
* SusC-like  = TonB-dependent transporter evidence (Pfam Plug + barrel,
               custom SusC HMM, or dbCAN-PUL 1.B.14 homolog).
* SusD-like  = Pfam SusD family, custom SusD HMM, or 8.A.46 homolog.
* susC/susD pair = a SusC-like and a SusD-like gene on the same strand with
               at most --pair-max-gap genes between them. Only paired
               SusC-like genes are PUL anchors; unpaired TonB-dependent
               transporters are kept as auxiliary transporters.
* Strong signal genes = susC/susD pair members, SusD-like genes and
               degradative CAZymes (GH, PL, CE, CBM, AA). GTs are auxiliary.
* Auxiliary genes = TF, STP, other transporters, GTs.
* Strong genes are chained into a locus when separated by at most
  --max-gap-genes genes, of which at most --max-nonsignal-gap carry no
  signal at all, and by at most --max-gap-bp. Loci are then extended over
  directly adjacent auxiliary genes (regulators, transporters).

Confidence
----------
high   = >=1 susC/susD pair + >=1 degradative CAZyme (canonical PUL)
medium = susC/susD pair without CAZyme, or SusD/TonB transporter + >=2 CAZymes
low    = CAZyme cluster (>=2 CAZymes + transporter / regulator); the
         Gram-positive style "CAZyme gene cluster" without a susC/susD pair
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# Pfam accessions, unversioned so that they match any Pfam release
PFAM_SUSC_EXTRA = {"PF13715", "PF07660"}
PFAM_SUSD = {"PF07980", "PF12741", "PF12771", "PF14322"}

DEGRADATIVE_PREFIXES = ("GH", "PL", "CE", "CBM", "AA")
CONF_LEVELS = ["high", "medium", "low"]

DOMTBL_COLS = [
    "target_name", "target_len", "query_name", "query_acc", "query_len",
    "full_evalue", "full_score", "i_evalue", "domain_score",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "hmm_cov",
]
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def _strip_pfam_version(acc):
    """PF00593.31 -> PF00593"""
    return acc.rsplit(".", 1)[0] if re.match(r"^PF\d+\.\d+$", acc) else acc


def clean_family(name):
    """GH13_10.hmm -> GH13_10"""
    return re.sub(r"\.hmm$", "", str(name))


def cazy_families_from_text(text):
    """Extract CAZy family tokens (GH13_10, CBM20, ...) from free text."""
    return re.findall(r"(?:GH|GT|PL|CE|CBM|AA)\d+(?:_\d+)?", str(text))


def is_degradative(families):
    return any(f.startswith(DEGRADATIVE_PREFIXES) for f in families)


# -- GFF / FASTA reading ----------------------------------------------------

def parse_attrs(attr_str):
    d = {}
    for item in attr_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def read_prodigal_gff(gff_path):
    """Read a Prodigal GFF and build gene ids that match ``prodigal -a``.

    GFF:  NZ_ABC.1  Prodigal  CDS  337  2799  .  +  0  ID=1_1;...
    FAA:  >NZ_ABC.1_1 # 337 # 2799 # 1 # ID=1_1;...

    The GFF ID is seqnum_genenum, the FASTA id is contig_genenum; HMMER and
    BLAST report the FASTA form, so that is the canonical gene_id.
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
            gene_num = raw_id.rsplit("_", 1)[1] if "_" in raw_id else raw_id
            rows.append({
                "contig": parts[0],
                "source": parts[1],
                "feature": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "phase": parts[7],
                "gene_id": f"{parts[0]}_{gene_num}",
                "partial": attrs.get("partial", "00"),
            })
    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit(f"No CDS entries parsed from {gff_path}")
    df = df.sort_values(["contig", "start", "end"]).reset_index(drop=True)
    df["gene_order"] = df.groupby("contig").cumcount() + 1
    df["gene_len_nt"] = df["end"] - df["start"] + 1
    return df


def read_fasta_dict(path):
    return {rec.id: rec for rec in SeqIO.parse(path, "fasta")}


def reconcile_gene_ids(genes, proteins):
    """Fall back to coordinate matching if GFF ids do not match the FASTA."""
    gff_ids = set(genes["gene_id"])
    overlap = gff_ids & set(proteins)
    if len(overlap) >= 0.5 * len(gff_ids):
        return genes
    log(f"WARNING: only {len(overlap)}/{len(gff_ids)} GFF gene ids match the "
        "protein FASTA; reconciling by coordinates")
    coord_to_fasta = {}
    for fid, rec in proteins.items():
        parts = rec.description.split(" # ")
        if len(parts) >= 3:
            try:
                contig = fid.rsplit("_", 1)[0]
                coord_to_fasta[(contig, int(parts[1]), int(parts[2]))] = fid
            except ValueError:
                pass
    genes = genes.copy()
    genes["gene_id"] = [
        coord_to_fasta.get((c, s, e), g)
        for c, s, e, g in zip(genes["contig"], genes["start"],
                              genes["end"], genes["gene_id"])]
    return genes


# -- HMMER / BLAST parsing --------------------------------------------------

def parse_hmmer_domtblout(path, evalue_cutoff=None, cov_cutoff=None):
    """Parse ``hmmsearch --domtblout`` (target = protein, query = HMM)."""
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=DOMTBL_COLS)
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split(maxsplit=22)
            if len(p) < 22:
                continue
            try:
                qlen = int(p[5])
                hmm_from, hmm_to = int(p[15]), int(p[16])
                rows.append((
                    p[0], int(p[2]), p[3], p[4], qlen,
                    float(p[6]), float(p[7]), float(p[12]), float(p[13]),
                    hmm_from, hmm_to, int(p[17]), int(p[18]),
                    (hmm_to - hmm_from + 1) / qlen if qlen else 0.0,
                ))
            except ValueError:
                continue
    df = pd.DataFrame(rows, columns=DOMTBL_COLS)
    if df.empty:
        return df
    if evalue_cutoff is not None:
        df = df[df["i_evalue"] <= evalue_cutoff]
    if cov_cutoff is not None:
        df = df[df["hmm_cov"] >= cov_cutoff]
    return df.copy()


def best_hits(df):
    if df.empty:
        return df
    return (df.sort_values(["target_name", "i_evalue", "full_score", "hmm_cov"],
                           ascending=[True, True, False, False])
              .drop_duplicates("target_name")
              .reset_index(drop=True))


def domain_architecture(df):
    """Non-overlapping domains per protein, N- to C-terminal, joined by '+'.

    Overlapping hits (>50 % of the shorter alignment) are resolved in favour
    of the lower i-Evalue, as dbCAN's hmmscan parser does.
    """
    out = {}
    for gene, sub in df.groupby("target_name", sort=False):
        kept = []
        for r in sub.sort_values("i_evalue").itertuples():
            clash = False
            for k in kept:
                ov = min(r.ali_to, k.ali_to) - max(r.ali_from, k.ali_from) + 1
                shorter = min(r.ali_to - r.ali_from, k.ali_to - k.ali_from) + 1
                if ov > 0.5 * shorter:
                    clash = True
                    break
            if not clash:
                kept.append(r)
        kept.sort(key=lambda r: r.ali_from)
        out[gene] = "+".join(clean_family(r.query_name) for r in kept)
    return out


def classify_suscd(pfam_hits):
    """SusC / SusD calls from Pfam domain architecture (accession based)."""
    out = []
    for gene, sub in pfam_hits.groupby("target_name"):
        domains = {_strip_pfam_version(a) for a in sub["query_acc"].astype(str)}
        # hmmsearch may have been run on a library named by accession
        domains |= {_strip_pfam_version(a) for a in sub["query_name"].astype(str)}
        susc = (("PF07715" in domains
                 and bool(domains & ({"PF00593"} | PFAM_SUSC_EXTRA)))
                or {"PF00593", "PF13715"} <= domains)
        susd = bool(domains & PFAM_SUSD)
        out.append({
            "gene_id": gene,
            "susC_like_pfam": susc,
            "susD_like_pfam": susd,
            "suscd_domains": ";".join(sorted(
                d for d in domains if d.startswith("PF"))),
        })
    return pd.DataFrame(
        out, columns=["gene_id", "susC_like_pfam", "susD_like_pfam",
                      "suscd_domains"])


def parse_blast(path):
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLS)
    return pd.read_csv(path, sep="\t", names=BLAST_COLS, header=None)


def tc_annotation_is_susc(x):
    return bool(re.search(r"(^|[^0-9])1\.B\.14(\.|$|[^0-9])", str(x)))


def tc_annotation_is_susd(x):
    return bool(re.search(r"(^|[^0-9])8\.A\.46(\.|$|[^0-9])", str(x)))


# -- Annotation integration -------------------------------------------------

def join_gene_annotations(genes, proteins, custom_hits, cluster_manifest,
                          pfam_hits, dbcan_hits, blast_hits, ref_table,
                          blast_min_pident=30.0, blast_min_qcov=50.0):
    genes = genes.copy()
    genes["protein_len_aa"] = genes["gene_id"].map(
        lambda x: len(proteins[x].seq.rstrip("*")) if x in proteins else 0)

    # custom marker HMMs: hmmbuild names each profile
    # "<marker_class>__<marker_group>"
    cm = pd.read_csv(cluster_manifest, sep="\t", dtype=str).fillna("")
    cm["hmm_name"] = cm["marker_class"] + "__" + cm["marker_group"]
    cm = cm[["hmm_name", "marker_class", "marker_group",
             "dominant_pul"]].drop_duplicates("hmm_name")
    custom_cols = ["best_custom_hmm", "best_custom_class", "best_custom_group",
                   "best_custom_pul", "best_custom_score",
                   "best_custom_ievalue", "best_custom_cov", "custom_classes"]
    if not custom_hits.empty:
        ch = custom_hits.merge(cm, left_on="query_name", right_on="hmm_name",
                               how="left")
        n_unmapped = int(ch["marker_class"].isna().sum())
        if n_unmapped:
            log(f"WARNING: {n_unmapped} custom HMM hits could not be mapped "
                "to the marker manifest")
        best = best_hits(ch).rename(columns={
            "target_name": "gene_id", "query_name": "best_custom_hmm",
            "marker_class": "best_custom_class",
            "marker_group": "best_custom_group",
            "dominant_pul": "best_custom_pul",
            "full_score": "best_custom_score",
            "i_evalue": "best_custom_ievalue", "hmm_cov": "best_custom_cov"})
        classes = (ch.dropna(subset=["marker_class"])
                     .groupby("target_name")["marker_class"]
                     .apply(lambda x: ";".join(sorted(set(x))))
                     .rename("custom_classes").rename_axis("gene_id")
                     .reset_index())
        genes = (genes.merge(best[["gene_id"] + custom_cols[:-1]],
                             on="gene_id", how="left")
                      .merge(classes, on="gene_id", how="left"))
    else:
        for c in custom_cols:
            genes[c] = pd.NA
    genes["custom_classes"] = genes["custom_classes"].fillna("")

    def has_custom(cls):
        return genes["custom_classes"].map(lambda s: cls in s.split(";"))

    # Pfam SusC / SusD
    genes = genes.merge(classify_suscd(pfam_hits), on="gene_id", how="left")
    for c in ("susC_like_pfam", "susD_like_pfam"):
        genes[c] = genes[c].fillna(False).astype(bool)

    # dbCAN family HMMs
    if not dbcan_hits.empty:
        best = best_hits(dbcan_hits).rename(columns={
            "target_name": "gene_id", "query_name": "best_dbcan_family",
            "full_score": "best_dbcan_score",
            "i_evalue": "best_dbcan_ievalue", "hmm_cov": "best_dbcan_cov"})
        best["best_dbcan_family"] = best["best_dbcan_family"].map(clean_family)
        genes = genes.merge(
            best[["gene_id", "best_dbcan_family", "best_dbcan_score",
                  "best_dbcan_ievalue", "best_dbcan_cov"]],
            on="gene_id", how="left")
        genes["dbcan_families"] = genes["gene_id"].map(
            domain_architecture(dbcan_hits))
    else:
        for c in ["best_dbcan_family", "best_dbcan_score",
                  "best_dbcan_ievalue", "best_dbcan_cov", "dbcan_families"]:
            genes[c] = pd.NA

    # dbCAN-PUL reference homology (BLAST)
    blast_cols = ["top_blast_subject", "top_blast_pul", "top_blast_class",
                  "top_blast_annotation", "top_blast_gene",
                  "top_blast_pident", "top_blast_bitscore",
                  "top_blast_evalue", "top_blast_qcovs"]
    if not blast_hits.empty:
        top = (blast_hits
               .sort_values(["qseqid", "bitscore", "evalue"],
                            ascending=[True, False, True])
               .drop_duplicates("qseqid")
               .merge(ref_table[["full_id", "pul_id", "ref_class",
                                 "ref_annotation", "gene_symbol"]],
                      left_on="sseqid", right_on="full_id", how="left")
               .rename(columns={
                   "qseqid": "gene_id", "sseqid": "top_blast_subject",
                   "pul_id": "top_blast_pul", "ref_class": "top_blast_class",
                   "ref_annotation": "top_blast_annotation",
                   "gene_symbol": "top_blast_gene",
                   "pident": "top_blast_pident",
                   "bitscore": "top_blast_bitscore",
                   "evalue": "top_blast_evalue", "qcovs": "top_blast_qcovs"}))
        genes = genes.merge(top[["gene_id"] + blast_cols], on="gene_id",
                            how="left")
    else:
        for c in blast_cols:
            genes[c] = pd.NA

    # A reference class is only transferred from a convincing homolog
    blast_ok = ((pd.to_numeric(genes["top_blast_pident"], errors="coerce")
                 >= blast_min_pident)
                & (pd.to_numeric(genes["top_blast_qcovs"], errors="coerce")
                   >= blast_min_qcov))
    genes["blast_class"] = genes["top_blast_class"].where(blast_ok)
    blast_annot = genes["top_blast_annotation"].where(blast_ok, "").fillna("")

    # composite calls
    genes["custom_susC_like"] = has_custom("SusC_like")
    genes["custom_susD_like"] = has_custom("SusD_like")
    genes["blast_susC_like"] = (genes["blast_class"].eq("TC")
                                & blast_annot.map(tc_annotation_is_susc))
    genes["blast_susD_like"] = (genes["blast_class"].eq("TC")
                                & blast_annot.map(tc_annotation_is_susd))
    genes["susD_like"] = genes[["susD_like_pfam", "custom_susD_like",
                                "blast_susD_like"]].any(axis=1)
    genes["susC_like"] = (genes[["susC_like_pfam", "custom_susC_like",
                                 "blast_susC_like"]].any(axis=1)
                          & ~genes["susD_like"])
    for kind in ("susC", "susD"):
        genes[f"{kind}_evidence"] = [
            ";".join(src for src, flag in
                     (("pfam", p), ("custom_hmm", c), ("blast", b)) if flag)
            for p, c, b in zip(genes[f"{kind}_like_pfam"],
                               genes[f"custom_{kind}_like"],
                               genes[f"blast_{kind}_like"])]

    # CAZyme family: dbCAN HMM architecture first, then reference labels
    fams, evid = [], []
    for dbcan, c_cls, c_grp, b_cls, b_ann in zip(
            genes["dbcan_families"], genes["best_custom_class"],
            genes["best_custom_group"], genes["blast_class"], blast_annot):
        e = []
        fam = ""
        if isinstance(dbcan, str) and dbcan:
            e.append("dbcan_hmm")
            fam = dbcan
        if c_cls == "CAZyme":
            e.append("custom_hmm")
            fam = fam or "+".join(cazy_families_from_text(c_grp)) or str(c_grp)
        if b_cls == "CAZyme":
            e.append("blast")
            fam = fam or "+".join(cazy_families_from_text(b_ann)) or str(b_ann)
        fams.append(fam)
        evid.append(";".join(e))
    genes["cazyme_family"] = fams
    genes["cazyme_evidence"] = evid
    genes["is_cazyme"] = (genes["cazyme_evidence"].ne("")
                          & ~genes["susC_like"] & ~genes["susD_like"])
    genes["is_degradative_cazyme"] = genes["is_cazyme"] & genes[
        "cazyme_family"].map(lambda f: is_degradative(
            cazy_families_from_text(f)) or not cazy_families_from_text(f))

    genes["is_tf"] = ((genes["best_custom_class"].eq("TF")
                       | genes["blast_class"].eq("TF")) & ~genes["is_cazyme"])
    genes["is_stp"] = ((genes["best_custom_class"].eq("STP")
                        | genes["blast_class"].eq("STP"))
                       & ~genes["is_cazyme"])
    genes["is_tc"] = (genes["best_custom_class"].isin(
                          ["TC_other", "SusC_like", "SusD_like"])
                      | genes["blast_class"].eq("TC")
                      | genes["susC_like"] | genes["susD_like"])
    genes["has_reference_hmm_hit"] = genes["best_custom_hmm"].notna()
    return genes


def mark_suscd_pairs(genes, pair_max_gap=1):
    """Flag SusC-like / SusD-like genes that sit next to each other.

    Each SusD-like gene is paired with the nearest unused SusC-like gene on
    the same strand within ``pair_max_gap`` intervening genes. Tandem
    susC-susD-susC-susD arrangements therefore give two pairs.
    """
    genes = genes.copy()
    genes["susCD_pair"] = ""
    n = 0
    for _, sub in genes.groupby("contig", sort=False):
        idx = sub.sort_values("gene_order").index.tolist()
        used = set()
        for pos, i in enumerate(idx):
            if not genes.at[i, "susD_like"]:
                continue
            best = None
            for d in range(1, pair_max_gap + 2):
                for q in (pos - d, pos + d):
                    if 0 <= q < len(idx):
                        j = idx[q]
                        if (genes.at[j, "susC_like"] and j not in used
                                and genes.at[j, "strand"]
                                == genes.at[i, "strand"]):
                            best = j
                            break
                if best is not None:
                    break
            if best is not None:
                n += 1
                used.add(best)
                genes.at[i, "susCD_pair"] = f"pair{n:03d}"
                genes.at[best, "susCD_pair"] = f"pair{n:03d}"
    paired = genes["susCD_pair"].ne("")

    genes["strong_signal"] = (paired | genes["susD_like"]
                              | genes["is_degradative_cazyme"])
    genes["aux_signal"] = (~genes["strong_signal"]
                           & (genes["is_tf"] | genes["is_stp"] | genes["is_tc"]
                              | genes["is_cazyme"]))
    genes["signal_gene"] = genes["strong_signal"] | genes["aux_signal"]

    def gene_class(r):
        if r["susC_like"] and r["susCD_pair"]:
            return "SusC"
        if r["susD_like"]:
            return "SusD"
        if r["is_cazyme"]:
            return "CAZyme"
        if r["is_tf"]:
            return "TF"
        if r["is_stp"]:
            return "STP"
        if r["is_tc"]:
            return "Transporter"
        return "Other"
    genes["gene_class"] = genes.apply(gene_class, axis=1)

    def label(r):
        cls = r["gene_class"]
        if cls == "CAZyme":
            return r["cazyme_family"] or "CAZyme"
        if cls in ("SusC", "SusD"):
            return f"sus{cls[-1]}"
        if cls == "Transporter":
            return "TBDT" if r["susC_like"] else "TC"
        return "" if cls == "Other" else cls
    genes["gene_label"] = genes.apply(label, axis=1)
    return genes


# -- Candidate locus building -----------------------------------------------

def build_candidate_blocks(genes, max_gap_genes=4, max_nonsignal_gap=2,
                           max_gap_bp=10000, max_aux_extension=3):
    """Chain strong signal genes into loci; return (contig, left, right)
    with positions relative to the contig's gene order (0-based)."""
    blocks = []
    for contig, sub in genes.groupby("contig", sort=False):
        sub = sub.sort_values("gene_order").reset_index(drop=True)
        strong = sub["strong_signal"].to_numpy()
        signal = sub["signal_gene"].to_numpy()
        aux = sub["aux_signal"].to_numpy()
        starts, ends = sub["start"].to_numpy(), sub["end"].to_numpy()
        strong_idx = [i for i, s in enumerate(strong) if s]
        if not strong_idx:
            continue
        chains = [[strong_idx[0]]]
        for idx in strong_idx[1:]:
            prev = chains[-1][-1]
            between = idx - prev - 1
            nonsignal = int((~signal[prev + 1:idx]).sum())
            bp_gap = max(0, int(starts[idx]) - int(ends[prev]) - 1)
            if (between <= max_gap_genes and nonsignal <= max_nonsignal_gap
                    and bp_gap <= max_gap_bp):
                chains[-1].append(idx)
            else:
                chains.append([idx])
        claimed = set()
        for chain in chains:
            claimed.update(range(chain[0], chain[-1] + 1))
        for chain in chains:
            left, right = chain[0], chain[-1]
            n = 0
            while (left > 0 and n < max_aux_extension and aux[left - 1]
                   and (left - 1) not in claimed):
                left -= 1
                n += 1
            n = 0
            while (right < len(sub) - 1 and n < max_aux_extension
                   and aux[right + 1] and (right + 1) not in claimed):
                right += 1
                n += 1
            claimed.update(range(left, right + 1))
            blocks.append((contig, left, right))
    return blocks


# -- Reference PUL assignment -----------------------------------------------

def reference_pul_votes(block_genes, blast_hits, pul_of, min_pident, min_qcov):
    """Score every reference PUL by the summed best bitscore of each gene.

    Unlike a top-hit vote, all retained BLAST hits are used, so a gene whose
    best hit lies in a different-but-similar PUL still supports the PUL that
    the rest of the locus agrees on.
    """
    if blast_hits.empty:
        return pd.DataFrame()
    h = blast_hits[blast_hits["qseqid"].isin(block_genes)
                   & (blast_hits["pident"] >= min_pident)
                   & (blast_hits["qcovs"] >= min_qcov)].copy()
    if h.empty:
        return pd.DataFrame()
    h["pul_id"] = h["sseqid"].map(pul_of)
    best = (h.sort_values("bitscore", ascending=False)
             .drop_duplicates(["qseqid", "pul_id"]))
    return (best.groupby("pul_id")
                .agg(score=("bitscore", "sum"),
                     n_genes=("qseqid", "nunique"),
                     n_ref_proteins=("sseqid", "nunique"),
                     mean_pident=("pident", "mean"))
                .sort_values(["score", "n_genes"], ascending=False))


def score_block(block_df, blast_hits, ref_info, min_pident, min_qcov):
    paired = block_df["susCD_pair"].ne("")
    counts = {
        "num_genes": len(block_df),
        "susCD_pairs": int(block_df.loc[paired, "susCD_pair"].nunique()),
        "susC_count": int((block_df["susC_like"] & paired).sum()),
        "susD_count": int(block_df["susD_like"].sum()),
        "tbdt_unpaired_count": int((block_df["susC_like"] & ~paired).sum()),
        "cazyme_count": int(block_df["is_cazyme"].sum()),
        "degradative_cazyme_count": int(
            block_df["is_degradative_cazyme"].sum()),
        "tc_count": int(block_df["is_tc"].sum()),
        "tf_count": int(block_df["is_tf"].sum()),
        "stp_count": int(block_df["is_stp"].sum()),
    }
    fams = []
    for f in block_df.loc[block_df["is_cazyme"], "cazyme_family"]:
        fams.extend(cazy_families_from_text(f) or [f])
    counts["cazyme_families"] = ",".join(sorted(set(fams)))

    ref = {"best_reference_pul": "", "best_reference_pul_score": 0.0,
           "ref_genes_matched": 0, "ref_pul_size": 0, "ref_pul_coverage": 0.0,
           "ref_mean_pident": 0.0, "predicted_substrate": "",
           "substrate_support": "", "reference_organism": "",
           "reference_verification": ""}
    votes = reference_pul_votes(set(block_df["gene_id"]), blast_hits,
                                ref_info["pul_of"], min_pident, min_qcov)
    if not votes.empty:
        pul = votes.index[0]
        v = votes.iloc[0]
        size = ref_info["pul_size"].get(pul, 0)
        cov = min(1.0, v["n_ref_proteins"] / size) if size else 0.0
        meta = ref_info["meta"].get(pul, {})
        if v["n_genes"] >= 3 and cov >= 0.5 and v["mean_pident"] >= 50:
            support = "strong"
        elif v["n_genes"] >= 2 and cov >= 0.25:
            support = "moderate"
        else:
            support = "weak"
        ref.update({
            "best_reference_pul": pul,
            "best_reference_pul_score": round(float(v["score"]), 1),
            "ref_genes_matched": int(v["n_genes"]),
            "ref_pul_size": int(size),
            "ref_pul_coverage": round(cov, 2),
            "ref_mean_pident": round(float(v["mean_pident"]), 1),
            "predicted_substrate": meta.get("substrate_final", ""),
            "substrate_support": support,
            "reference_organism": meta.get("organism_name", ""),
            "reference_verification": meta.get("verification_final", ""),
        })

    n_caz = counts["degradative_cazyme_count"]
    n_aux = counts["tc_count"] + counts["tf_count"] + counts["stp_count"]
    if counts["susCD_pairs"] >= 1 and n_caz >= 1:
        confidence, pul_type = "high", "canonical_SusCD_PUL"
    elif counts["susCD_pairs"] >= 1:
        confidence, pul_type = "medium", "SusCD_only"
    elif n_caz >= 2 and (counts["susD_count"]
                         + counts["tbdt_unpaired_count"]) >= 1:
        confidence, pul_type = "medium", "incomplete_SusCD_PUL"
    elif n_caz >= 2 and n_aux >= 1:
        confidence, pul_type = "low", "CAZyme_cluster"
    else:
        confidence, pul_type = "discard", ""
    return counts, ref, confidence, pul_type


# -- Output writing ---------------------------------------------------------

GENE_TABLE_COLS = [
    "candidate_id", "gene_id", "contig", "start", "end", "strand",
    "protein_len_aa", "gene_class", "gene_label", "susCD_pair",
    "cazyme_family", "cazyme_evidence", "susC_evidence", "susD_evidence",
    "best_custom_hmm", "top_blast_pul", "top_blast_gene",
    "top_blast_class", "top_blast_annotation", "top_blast_pident",
    "top_blast_qcovs",
]


def write_candidate_outputs(candidate_id, block_df, genome, proteins, cds,
                            outdir, flank_bp=0):
    cdir = outdir / candidate_id
    cdir.mkdir(parents=True, exist_ok=True)
    contig = block_df.iloc[0]["contig"]
    contig_len = len(genome[contig].seq)
    region_start = max(1, int(block_df["start"].min()) - flank_bp)
    region_end = min(contig_len, int(block_df["end"].max()) + flank_bp)

    region_rec = genome[contig][region_start - 1:region_end]
    region_rec.id = f"{candidate_id}|{contig}:{region_start}-{region_end}"
    region_rec.description = ""
    SeqIO.write(region_rec, cdir / f"{candidate_id}.region.fna", "fasta")

    prot_recs = [proteins[g] for g in block_df["gene_id"] if g in proteins]
    cds_recs = [cds[g] for g in block_df["gene_id"] if g in cds]
    if prot_recs:
        SeqIO.write(prot_recs, cdir / f"{candidate_id}.proteins.faa", "fasta")
    if cds_recs:
        SeqIO.write(cds_recs, cdir / f"{candidate_id}.cds.fna", "fasta")

    cols = [c for c in GENE_TABLE_COLS if c in block_df.columns]
    extra = [c for c in block_df.columns if c not in cols]
    block_df[cols + extra].to_csv(cdir / f"{candidate_id}.genes.tsv",
                                  sep="\t", index=False)

    with open(cdir / f"{candidate_id}.gff3", "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region {contig} 1 {contig_len}\n")
        for r in block_df.itertuples():
            attrs = [f"ID={r.gene_id}", f"Name={r.gene_label or r.gene_id}",
                     f"gene_class={r.gene_class}", f"candidate={candidate_id}"]
            fh.write("\t".join(map(str, [
                r.contig, r.source, r.feature, r.start, r.end, ".",
                r.strand, r.phase, ";".join(attrs)])) + "\n")


# -- Main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Call candidate PULs from HMM/BLAST annotations and "
                    "export wet-lab-friendly outputs")
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
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", default="",
                    help="Sample / genome name written into the summary")
    g = ap.add_argument_group("evidence thresholds")
    g.add_argument("--custom-ievalue", type=float, default=1e-5)
    g.add_argument("--custom-cov", type=float, default=0.40)
    g.add_argument("--dbcan-ievalue", type=float, default=1e-15)
    g.add_argument("--dbcan-cov", type=float, default=0.35)
    g.add_argument("--blast-min-pident", type=float, default=30.0)
    g.add_argument("--blast-min-qcov", type=float, default=50.0)
    g = ap.add_argument_group("locus definition")
    g.add_argument("--pair-max-gap", type=int, default=1,
                   help="Genes allowed between susC and susD of one pair")
    g.add_argument("--max-gap-genes", type=int, default=4,
                   help="Max genes between two chained strong signal genes")
    g.add_argument("--max-nonsignal-gap", type=int, default=2,
                   help="Max of those genes that carry no signal at all")
    g.add_argument("--max-gap-bp", type=int, default=10000)
    g.add_argument("--max-aux-extension", type=int, default=3,
                   help="Adjacent regulator/transporter genes added per edge")
    g.add_argument("--flank-bp", type=int, default=0,
                   help="Extra sequence written around each region FASTA")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    candidate_dir = outdir / "candidate_PULs"
    candidate_dir.mkdir(parents=True, exist_ok=True)

    genes = read_prodigal_gff(args.gff)
    proteins = read_fasta_dict(args.proteins)
    cds = read_fasta_dict(args.cds)
    genome = read_fasta_dict(args.genome)
    genes = reconcile_gene_ids(genes, proteins)

    ref_table = pd.read_csv(args.reference_table, sep="\t",
                            dtype=str).fillna("")
    meta_cols = [c for c in ("substrate_final", "organism_name",
                             "verification_final") if c in ref_table.columns]
    ref_info = {
        "pul_of": dict(zip(ref_table["full_id"], ref_table["pul_id"])),
        "pul_size": ref_table.groupby("pul_id")["full_id"].nunique().to_dict(),
        "meta": (ref_table.drop_duplicates("pul_id")
                          .set_index("pul_id")[meta_cols].to_dict("index")),
    }

    custom_hits = parse_hmmer_domtblout(
        args.custom_domtblout, args.custom_ievalue, args.custom_cov)
    pfam_hits = parse_hmmer_domtblout(args.suscd_domtblout, None, 0.30)
    dbcan_hits = parse_hmmer_domtblout(
        args.dbcan_domtblout, args.dbcan_ievalue, args.dbcan_cov)
    blast_hits = parse_blast(args.blast_tsv)

    log(f"Genes from GFF:         {len(genes)}")
    log(f"Gene-FASTA overlap:     "
        f"{len(set(genes['gene_id']) & set(proteins))}")
    log(f"Custom HMM hits (filt): {len(custom_hits)}")
    log(f"Pfam SusCD hits (filt): {len(pfam_hits)}")
    log(f"dbCAN hits (filt):      {len(dbcan_hits)}")
    log(f"BLAST hits:             {len(blast_hits)}")

    genes = join_gene_annotations(
        genes, proteins, custom_hits, args.cluster_manifest, pfam_hits,
        dbcan_hits, blast_hits, ref_table,
        args.blast_min_pident, args.blast_min_qcov)
    genes = mark_suscd_pairs(genes, args.pair_max_gap)

    for name, col in [("susC-like (TBDT)", "susC_like"),
                      ("susD-like", "susD_like"), ("CAZymes", "is_cazyme"),
                      ("  degradative", "is_degradative_cazyme"),
                      ("TF", "is_tf"), ("STP", "is_stp"),
                      ("strong signal genes", "strong_signal")]:
        log(f"{name + ':':24s}{int(genes[col].sum())}")
    log(f"{'susC/susD pairs:':24s}"
        f"{genes.loc[genes['susCD_pair'].ne(''), 'susCD_pair'].nunique()}")

    blocks = build_candidate_blocks(
        genes, args.max_gap_genes, args.max_nonsignal_gap, args.max_gap_bp,
        args.max_aux_extension)
    log(f"Raw candidate loci:     {len(blocks)}")

    genes["candidate_id"] = ""
    summary_rows, bed_rows = [], []
    by_contig = {c: s.sort_values("gene_order")
                 for c, s in genes.groupby("contig", sort=False)}
    for contig, left, right in sorted(blocks, key=lambda b: (b[0], b[1])):
        block_df = by_contig[contig].iloc[left:right + 1].copy()
        counts, ref, confidence, pul_type = score_block(
            block_df, blast_hits, ref_info,
            args.blast_min_pident, args.blast_min_qcov)
        if confidence == "discard":
            continue
        candidate_id = f"PULcand_{len(summary_rows) + 1:04d}"
        genes.loc[block_df.index, "candidate_id"] = candidate_id
        block_df["candidate_id"] = candidate_id
        start, end = int(block_df["start"].min()), int(block_df["end"].max())
        contig_len = len(genome[contig].seq)
        summary_rows.append({
            "sample": args.sample,
            "candidate_id": candidate_id,
            "contig": contig, "start": start, "end": end,
            "length_bp": end - start + 1,
            "confidence": confidence, "pul_type": pul_type,
            **counts, **ref,
            "at_contig_edge": bool(
                left == 0 or right == len(by_contig[contig]) - 1
                or start < 500 or end > contig_len - 500),
            "gene_order": " ".join(
                (lab if lab else "-") + ("<" if s == "-" else ">")
                for lab, s in zip(block_df["gene_label"], block_df["strand"])),
        })
        bed_rows.append([contig, start - 1, end, candidate_id,
                         {"high": 1000, "medium": 600, "low": 300}[confidence],
                         "."])
        write_candidate_outputs(candidate_id, block_df, genome, proteins,
                                cds, candidate_dir, args.flank_bp)

    genes.to_csv(outdir / "gene_annotation_master.tsv", sep="\t", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "PUL_summary.tsv", sep="\t", index=False)
    pd.DataFrame(bed_rows).to_csv(outdir / "PUL_regions.bed", sep="\t",
                                  header=False, index=False)

    log(f"Candidate PULs kept:    {len(summary)}")
    with open(outdir / "README_results.txt", "w") as fh:
        fh.write("PULseHunter results\n===================\n")
        fh.write(f"candidate_PUL_count\t{len(summary)}\n")
        for lvl in CONF_LEVELS:
            n = int((summary["confidence"] == lvl).sum()) if len(summary) else 0
            log(f"  {lvl}: {n}")
            fh.write(f"{lvl}_confidence\t{n}\n")
        fh.write(
            "\nMain files:\n"
            "- PUL_summary.tsv: one row per candidate locus, in genome order\n"
            "- PUL_regions.bed: locus coordinates (score 1000/600/300 = "
            "high/medium/low)\n"
            "- gene_annotation_master.tsv: per-gene integrated annotation\n"
            "- candidate_PULs/: per-locus region FASTA, proteins, CDS, GFF3, "
            "gene table\n"
            "\nConfidence:\n"
            "high   = susC/susD pair + degradative CAZyme (canonical PUL)\n"
            "medium = susC/susD pair without CAZyme, or SusD / TonB-dependent "
            "transporter with >=2 CAZymes\n"
            "low    = CAZyme cluster (>=2 CAZymes + transporter/regulator) "
            "without susC/susD\n"
            "\nsubstrate_support (transfer from the best dbCAN-PUL match):\n"
            "strong   = >=3 genes, >=50 % of the reference PUL recovered, "
            "mean identity >=50 %\n"
            "moderate = >=2 genes and >=25 % of the reference PUL\n"
            "weak     = anything less; treat the substrate as a hint only\n")


if __name__ == "__main__":
    main()
