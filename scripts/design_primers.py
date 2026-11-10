#!/usr/bin/env python3
"""
design_primers.py - design the best N primer sets for every candidate PUL.

Modes
-----
qpcr  (default) 80-200 bp amplicons inside the locus marker genes
      (susC, then susD, then the CAZymes). Suited to RT-qPCR of PUL
      induction and to strain / locus quantification.
pcr   400-1000 bp amplicons spanning the susC/susD junction (or the first
      marker gene when the locus has no pair). Suited to end-point PCR
      confirmation and Sanger sequencing of the locus.

Every candidate primer is checked by in-silico PCR against the target
genome and, optionally, other genomes (--offtarget-fasta), e.g. the other
members of a synthetic community. A primer is predicted to prime wherever
its 3'-terminal 10 bases match exactly with <= 3 mismatches upstream. Sets are ranked by
(1) specificity, (2) marker-gene priority, (3) Primer3 penalty.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq

MODES = {
    "qpcr": {"product": (80, 200), "opt_tm": 60.0, "min_tm": 58.0,
             "max_tm": 62.0},
    "pcr": {"product": (400, 1000), "opt_tm": 60.0, "min_tm": 57.0,
            "max_tm": 63.0},
}
CLASS_PRIORITY = {"SusC": 0, "SusD": 1, "CAZyme": 2}


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def primer3_pairs(template, mode, n_return, target=None):
    cfg = MODES[mode]
    seq_args = {"SEQUENCE_ID": "t", "SEQUENCE_TEMPLATE": template}
    if target:
        seq_args["SEQUENCE_TARGET"] = [target]
    glob = {
        "PRIMER_TASK": "generic",
        "PRIMER_NUM_RETURN": n_return,
        "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 24,
        "PRIMER_OPT_TM": cfg["opt_tm"], "PRIMER_MIN_TM": cfg["min_tm"],
        "PRIMER_MAX_TM": cfg["max_tm"], "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
        "PRIMER_MIN_GC": 40.0, "PRIMER_OPT_GC_PERCENT": 50.0,
        "PRIMER_MAX_GC": 60.0,
        "PRIMER_GC_CLAMP": 1, "PRIMER_MAX_END_GC": 3,
        "PRIMER_MAX_POLY_X": 4, "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY_TH": 45.0, "PRIMER_MAX_SELF_END_TH": 35.0,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0,
        "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
        "PRIMER_MAX_HAIRPIN_TH": 24.0,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(cfg["product"])],
        "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_SALT_DIVALENT": 1.5,
        "PRIMER_DNTP_CONC": 0.6, "PRIMER_DNA_CONC": 50.0,
        "PRIMER_TM_FORMULA": 1, "PRIMER_SALT_CORRECTIONS": 1,
    }
    try:
        res = primer3.design_primers(seq_args, glob)
    except (OSError, ValueError) as e:
        log(f"  primer3 failed: {e}")
        return []
    out = []
    for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        fpos, flen = res[f"PRIMER_LEFT_{i}"]
        rpos, rlen = res[f"PRIMER_RIGHT_{i}"]
        out.append({
            "fwd": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "rev": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "fwd_tm": round(res[f"PRIMER_LEFT_{i}_TM"], 1),
            "rev_tm": round(res[f"PRIMER_RIGHT_{i}_TM"], 1),
            "fwd_gc": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1),
            "rev_gc": round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 1),
            "product_size": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            "penalty": round(res[f"PRIMER_PAIR_{i}_PENALTY"], 3),
            # 0-based offsets within the template
            "t_start": fpos, "t_end": rpos + 1,
            "fwd_len": flen, "rev_len": rlen,
        })
    return out


_CODE = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _i
    _CODE[ord(_b.lower())] = _i


def kmer_code(seq):
    code = 0
    for ch in seq:
        code = code * 4 + "ACGT".index(ch)
    return code


def rolling_codes(seq_bytes, k):
    """Integer code of every k-mer of a sequence (-1 where it contains N)."""
    base = _CODE[np.frombuffer(seq_bytes, dtype=np.uint8)]
    n = len(base) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.int64)
    codes = np.zeros(n, dtype=np.int64)
    bad = np.zeros(n, dtype=bool)
    for j in range(k):
        window = base[j:j + n]
        bad |= window < 0
        codes = codes * 4 + np.where(window < 0, 0, window)
    codes[bad] = -1
    return codes


def iter_specificity_contigs(genome, offtargets):
    for rec in SeqIO.parse(genome, "fasta"):
        yield f"TARGET|{rec.id}", str(rec.seq).upper()
    for i, path in enumerate(offtargets, 1):
        name = Path(path).stem
        for rec in SeqIO.parse(path, "fasta"):
            yield f"OFF{i}:{name}|{rec.id}", str(rec.seq).upper()


def primer_binding_sites(primers, genome, offtargets, seed_len, max_mismatch):
    """In-silico priming sites of every primer on both strands.

    A site needs an exact match of the primer's 3'-terminal ``seed_len``
    bases and at most ``max_mismatch`` mismatches over the rest of the
    primer. Returns primer -> set of (contig, strand, three_prime_pos) with
    1-based positions.
    """
    plus_seed, minus_seed = defaultdict(list), defaultdict(list)
    for p in primers:
        seed = p[-seed_len:]
        plus_seed[kmer_code(seed)].append(p)
        minus_seed[kmer_code(str(Seq(seed).reverse_complement()))].append(p)
    wanted = np.array(sorted(set(plus_seed) | set(minus_seed)),
                      dtype=np.int64)
    revcomp = {p: str(Seq(p).reverse_complement()) for p in primers}
    sites = defaultdict(set)

    def mismatches(a, b):
        return sum(x != y for x, y in zip(a, b))

    for name, seq in iter_specificity_contigs(genome, offtargets):
        codes = rolling_codes(seq.encode(), seed_len)
        for i in np.flatnonzero(np.isin(codes, wanted)):
            i = int(i)
            code = int(codes[i])
            for p in plus_seed.get(code, ()):
                lo = i + seed_len - len(p)
                if lo >= 0 and mismatches(seq[lo:i + seed_len],
                                          p) <= max_mismatch:
                    sites[p].add((name, "+", i + seed_len))
            for p in minus_seed.get(code, ()):
                hi = i + len(p)
                if hi <= len(seq) and mismatches(seq[i:hi],
                                                 revcomp[p]) <= max_mismatch:
                    sites[p].add((name, "-", i + 1))
    return sites


def count_amplicons(fwd, rev, sites, max_product):
    """In-silico PCR: a priming site on the + strand facing one on the -
    strand within max_product (covers fwd/rev, fwd/fwd and rev/rev)."""
    plus, minus = defaultdict(list), defaultdict(list)
    for p in {fwd, rev}:
        for subj, strand, pos in sites.get(p, ()):
            (plus if strand == "+" else minus)[subj].append(pos)
    on = off = 0
    off_names = set()
    for subj, ppos in plus.items():
        for a in ppos:
            for b in minus.get(subj, ()):
                if 0 < b - a <= max_product:
                    if subj.startswith("TARGET|"):
                        on += 1
                    else:
                        off += 1
                        off_names.add(subj.split("|")[0].split(":", 1)[1])
    return on, off, ";".join(sorted(off_names))


def main():
    ap = argparse.ArgumentParser(
        description="Design the best N primer sets per candidate PUL")
    ap.add_argument("--final-dir", required=True,
                    help="PULseHunter final/ directory")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--mode", choices=sorted(MODES), default="qpcr")
    ap.add_argument("--n-sets", type=int, default=3)
    ap.add_argument("--offtarget-fasta", action="append", default=[],
                    help="Other genome(s) the primers must not amplify "
                         "(repeatable)")
    ap.add_argument("--seed-len", type=int, default=10,
                    help="3'-terminal bases that must match exactly for a "
                         "primer to be predicted to prime (default 10)")
    ap.add_argument("--max-mismatch", type=int, default=3,
                    help="Mismatches tolerated in the rest of the primer "
                         "when predicting priming (default 3)")
    ap.add_argument("--max-offtarget-product", type=int, default=3000)
    ap.add_argument("--pairs-per-gene", type=int, default=15)
    ap.add_argument("--max-genes", type=int, default=4,
                    help="Marker genes tried per locus")
    ap.add_argument("--threads", type=int, default=4,
                    help="Accepted for interface compatibility; unused")
    args = ap.parse_args()

    final = Path(args.final_dir)
    summary = pd.read_csv(final / "PUL_summary.tsv", sep="\t")
    outdir = final / "primers"
    outdir.mkdir(exist_ok=True)
    cols = ["sample", "candidate_id", "rank", "mode", "target_gene",
            "target_class", "target_label", "fwd_primer", "rev_primer",
            "fwd_tm", "rev_tm", "fwd_gc", "rev_gc", "product_size", "contig",
            "amplicon_start", "amplicon_end", "primer3_penalty",
            "on_target_amplicons", "off_target_amplicons",
            "off_target_genomes", "specific"]
    if summary.empty:
        pd.DataFrame(columns=cols).to_csv(outdir / "PUL_primers.tsv",
                                          sep="\t", index=False)
        return
    genes = pd.read_csv(final / "gene_annotation_master.tsv", sep="\t",
                        low_memory=False)
    genes = genes[genes["candidate_id"].notna() & (genes["candidate_id"] != "")]
    genome = {r.id: r.seq for r in SeqIO.parse(args.genome, "fasta")}

    # 1) candidate primer pairs -------------------------------------------
    cands = []
    for cid, block in genes.groupby("candidate_id", sort=True):
        block = block.sort_values("start")
        contig = block.iloc[0]["contig"]
        markers = block[block["gene_class"].isin(CLASS_PRIORITY)].copy()
        markers["prio"] = markers["gene_class"].map(CLASS_PRIORITY)
        markers = markers.sort_values(["prio", "gene_len_nt"],
                                      ascending=[True, False])
        if args.mode == "qpcr":
            for rank, g in enumerate(
                    markers.head(args.max_genes).itertuples()):
                tmpl = str(genome[contig][g.start - 1:g.end]).upper()
                for p in primer3_pairs(tmpl, "qpcr", args.pairs_per_gene):
                    p.update(candidate_id=cid, contig=contig,
                             target_gene=g.gene_id, target_class=g.gene_class,
                             target_label=g.gene_label, gene_rank=rank,
                             amplicon_start=g.start + p["t_start"],
                             amplicon_end=g.start + p["t_end"] - 1)
                    cands.append(p)
        else:
            r_start, r_end = int(block["start"].min()), int(block["end"].max())
            tmpl = str(genome[contig][r_start - 1:r_end]).upper()
            pair = block[block["susCD_pair"].notna()
                         & (block["susCD_pair"] != "")]
            targets = []
            for _, pg in pair.groupby("susCD_pair"):
                pg = pg.sort_values("start")
                junction = (int(pg.iloc[0]["end"]) + int(pg.iloc[-1]["start"])) // 2
                targets.append((junction - r_start - 10, 20,
                                pg.iloc[0]["gene_id"], "SusC/SusD junction"))
            if not targets and len(markers):
                g = markers.iloc[0]
                mid = (int(g["start"]) + int(g["end"])) // 2
                targets.append((mid - r_start - 10, 20, g["gene_id"],
                                g["gene_label"]))
            for rank, (pos, length, gid, lab) in enumerate(
                    targets[:args.max_genes]):
                pos = max(0, min(pos, len(tmpl) - length))
                for p in primer3_pairs(tmpl, "pcr", args.pairs_per_gene,
                                       target=[pos, length]):
                    p.update(candidate_id=cid, contig=contig, target_gene=gid,
                             target_class="junction", target_label=lab,
                             gene_rank=rank,
                             amplicon_start=r_start + p["t_start"],
                             amplicon_end=r_start + p["t_end"] - 1)
                    cands.append(p)
    log(f"Primer3 candidate pairs: {len(cands)} for "
        f"{summary['candidate_id'].nunique()} loci")
    if not cands:
        pd.DataFrame(columns=cols).to_csv(outdir / "PUL_primers.tsv",
                                          sep="\t", index=False)
        return

    # 2) in-silico PCR specificity ----------------------------------------
    primers = sorted({p["fwd"] for p in cands} | {p["rev"] for p in cands})
    sites = primer_binding_sites(primers, args.genome, args.offtarget_fasta,
                                 args.seed_len, args.max_mismatch)
    for p in cands:
        on, off, names = count_amplicons(p["fwd"], p["rev"], sites,
                                         args.max_offtarget_product)
        p["on_target_amplicons"], p["off_target_amplicons"] = on, off
        p["off_target_genomes"] = names
        p["specific"] = (on == 1 and off == 0)

    # 3) pick the best N non-redundant sets per locus ---------------------
    rows = []
    sample = summary["sample"].iloc[0] if "sample" in summary else ""
    by_cand = defaultdict(list)
    for p in cands:
        by_cand[p["candidate_id"]].append(p)
    for cid in summary["candidate_id"]:
        pool = sorted(by_cand.get(cid, []), key=lambda p: (
            not p["specific"], p["off_target_amplicons"],
            abs(p["on_target_amplicons"] - 1), p["gene_rank"], p["penalty"]))
        chosen = []
        for p in pool:
            # skip sets whose amplicon mostly overlaps one already chosen
            redundant = False
            for c in chosen:
                ov = (min(p["amplicon_end"], c["amplicon_end"])
                      - max(p["amplicon_start"], c["amplicon_start"]) + 1)
                if ov > 0.5 * min(p["product_size"], c["product_size"]):
                    redundant = True
                    break
            if not redundant:
                chosen.append(p)
            if len(chosen) == args.n_sets:
                break
        for rank, p in enumerate(chosen, 1):
            rows.append({
                "sample": sample, "candidate_id": cid, "rank": rank,
                "mode": args.mode, "target_gene": p["target_gene"],
                "target_class": p["target_class"],
                "target_label": p["target_label"],
                "fwd_primer": p["fwd"], "rev_primer": p["rev"],
                "fwd_tm": p["fwd_tm"], "rev_tm": p["rev_tm"],
                "fwd_gc": p["fwd_gc"], "rev_gc": p["rev_gc"],
                "product_size": p["product_size"], "contig": p["contig"],
                "amplicon_start": p["amplicon_start"],
                "amplicon_end": p["amplicon_end"],
                "primer3_penalty": p["penalty"],
                "on_target_amplicons": p["on_target_amplicons"],
                "off_target_amplicons": p["off_target_amplicons"],
                "off_target_genomes": p["off_target_genomes"],
                "specific": p["specific"],
            })
    out = pd.DataFrame(rows, columns=cols)
    out.to_csv(outdir / "PUL_primers.tsv", sep="\t", index=False)
    for cid, sub in out.groupby("candidate_id"):
        cdir = final / "candidate_PULs" / cid
        if cdir.is_dir():
            sub.to_csv(cdir / f"{cid}.primers.tsv", sep="\t", index=False)
    with open(outdir / "PUL_primers.fasta", "w") as fh:
        for r in out.itertuples():
            fh.write(f">{r.candidate_id}_set{r.rank}_F\n{r.fwd_primer}\n"
                     f">{r.candidate_id}_set{r.rank}_R\n{r.rev_primer}\n")

    n_loci = summary["candidate_id"].nunique()
    with_any = out["candidate_id"].nunique()
    full = int((out.groupby("candidate_id")["specific"].sum()
                >= args.n_sets).sum())
    log(f"Loci with >=1 primer set: {with_any}/{n_loci}")
    log(f"Loci with {args.n_sets} specific sets: {full}/{n_loci}")
    # sanity: amplicon really is what the primers bracket
    for r in out.head(50).itertuples():
        amp = genome[r.contig][r.amplicon_start - 1:r.amplicon_end]
        assert str(amp).upper().startswith(r.fwd_primer), r
        assert str(Seq(str(amp)).reverse_complement()).upper().startswith(
            r.rev_primer), r


if __name__ == "__main__":
    main()
