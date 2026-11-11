#!/usr/bin/env python3
"""
plot_puls.py - gene-arrow maps for every candidate PUL plus a genome
overview. Writes PNG and SVG per locus into final/figures/.

Each map shows, top to bottom:
  * the genes, drawn to scale as strand-aware arrows, coloured by class and
    labelled with their CAZy family / role (class is never colour-only),
  * the dbCAN3 CGC call(s) overlapping the locus (if dbCAN3 was run),
  * the amplicons of the designed primer sets (if primers were designed).
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch, Polygon, Rectangle  # noqa: E402

# Fixed class -> colour assignment (colour-blind-validated categorical order)
CLASS_COLORS = {
    "SusC": "#2a78d6",
    "SusD": "#eb6834",
    "CAZyme": "#1baf7a",
    "TF": "#eda100",
    "STP": "#e87ba4",
    "Transporter": "#6250d6",
    "Other": "#d3d2cb",
}
CLASS_NAMES = {
    "SusC": "SusC-like",
    "SusD": "SusD-like",
    "CAZyme": "CAZyme",
    "TF": "Transcription factor",
    "STP": "Signal transduction",
    "Transporter": "Other transporter",
    "Other": "Other",
}
INK, INK2, MUTED, AXIS, SURFACE = ("#0b0b0b", "#52514e", "#898781",
                                   "#c3c2b7", "#fcfcfb")
CONF_COLORS = {"high": "#104281", "medium": "#3987e5", "low": "#9ec5f4"}

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 8,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def arrow(x0, x1, strand, y, h, head):
    """Polygon for a gene from x0 to x1 (x0 < x1) pointing along strand."""
    head = min(head, x1 - x0)
    if strand == "-":
        return [(x1, y - h / 2), (x1, y + h / 2), (x0 + head, y + h / 2),
                (x0, y), (x0 + head, y - h / 2)]
    return [(x0, y - h / 2), (x0, y + h / 2), (x1 - head, y + h / 2),
            (x1, y), (x1 - head, y - h / 2)]


def _s(v):
    return "" if pd.isna(v) else str(v)


def plot_locus(row, block, cgcs, primers, out_base):
    start, end = int(row["start"]), int(row["end"])
    span = end - start + 1
    pad = span * 0.02
    has_cgc = cgcs is not None and len(cgcs) > 0
    has_primers = primers is not None and len(primers) > 0
    width = min(16, max(7.5, span / 2600))
    n_sets = int(primers["rank"].max()) if has_primers else 0
    height = 2.5 + 0.3 * has_cgc + 0.14 * n_sets
    fig, ax = plt.subplots(figsize=(width, height))

    y = 0.0
    ax.plot([start, end], [y, y], color=AXIS, lw=1, zorder=1)
    head = span * 0.012
    for g in block.itertuples():
        cls = g.gene_class if g.gene_class in CLASS_COLORS else "Other"
        ax.add_patch(Polygon(arrow(g.start, g.end, g.strand, y, 0.46, head),
                             closed=True, facecolor=CLASS_COLORS[cls],
                             edgecolor=SURFACE, linewidth=0.8, zorder=3))
        lab = _s(g.gene_label)
        if lab:
            if len(lab) > 22:
                lab = lab[:21] + "…"
            ax.text((g.start + g.end) / 2, y + 0.3, lab, rotation=45,
                    ha="left", va="bottom", fontsize=7, color=INK, zorder=4)

    track_y = y - 0.55
    if has_cgc:
        for c in cgcs.itertuples():
            x0, x1 = max(start - pad, c.start), min(end + pad, c.end)
            ax.add_patch(Rectangle((x0, track_y - 0.06), x1 - x0, 0.12,
                                   facecolor=MUTED, edgecolor="none",
                                   zorder=2))
            ax.text(x0, track_y - 0.1, f"dbCAN3 {c.cgc_id}", fontsize=6.5,
                    color=INK2, ha="left", va="top")
        track_y -= 0.42
    if has_primers:
        # one row per primer set so neighbouring amplicons never collide
        for p in primers.itertuples():
            yy = track_y - 0.17 * (p.rank - 1)
            ax.plot([p.amplicon_start, p.amplicon_end], [yy, yy],
                    color=INK, lw=3, solid_capstyle="butt", zorder=3)
            ax.text(p.amplicon_end + span * 0.004, yy,
                    f"set {p.rank} ({p.product_size} bp)", fontsize=6,
                    color=INK2, ha="left", va="center")
        ax.text(start - pad, track_y, "primers ", fontsize=6.5,
                color=INK2, ha="right", va="center")
        track_y -= 0.17 * n_sets + 0.05

    ax.set_xlim(start - pad, end + pad)
    ax.set_ylim(track_y - 0.1, y + 1.05)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(axis="x", colors=MUTED, labelsize=7, length=3)
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v / 1000:,.1f}"))
    ax.set_xlabel(f"{row['contig']} position (kb)", fontsize=7.5, color=INK2)

    sample = _s(row.get("sample"))
    title = f"{sample + '  ' if sample else ''}{row['candidate_id']}"
    sub = [f"{span / 1000:.1f} kb", f"{int(row['num_genes'])} genes",
           f"confidence: {row['confidence']}"]
    substrate = _s(row.get("predicted_substrate"))
    if substrate:
        sub.append(f"best dbCAN-PUL match: {row['best_reference_pul']} "
                   f"({substrate}; {_s(row.get('substrate_support'))} support)")
    fig.text(0.01, 0.97, title, fontsize=10.5, fontweight="bold",
             color=INK, ha="left", va="top")
    fig.text(0.01, 0.97 - 0.27 / height, "  ·  ".join(sub), fontsize=8,
             color=INK2, ha="left", va="top")

    present = [c for c in CLASS_COLORS if c in set(block["gene_class"])]
    fig.legend(handles=[Patch(facecolor=CLASS_COLORS[c], edgecolor="none",
                              label=CLASS_NAMES[c]) for c in present],
              loc="lower right", bbox_to_anchor=(0.995, 0.0),
              ncol=len(present), fontsize=6.5, frameon=False,
              handlelength=1.1, columnspacing=1.0, labelcolor=INK2)
    fig.subplots_adjust(left=0.07 if has_primers else 0.02, right=0.99,
                        top=1 - 0.55 / height, bottom=0.78 / height)
    for ext in ("png", "svg"):
        fig.savefig(f"{out_base}.{ext}", dpi=200)
    plt.close(fig)


def plot_overview(summary, contig_lengths, out_base, sample):
    """Where the loci sit on the (largest) contigs, shaded by confidence."""
    contigs = [c for c in sorted(contig_lengths, key=contig_lengths.get,
                                 reverse=True)
               if c in set(summary["contig"])][:15]
    if not contigs:
        return
    fig, ax = plt.subplots(figsize=(10, 0.6 + 0.42 * len(contigs)))
    for i, c in enumerate(contigs):
        yy = len(contigs) - i - 1
        ax.plot([0, contig_lengths[c]], [yy, yy], color=AXIS, lw=1.2,
                solid_capstyle="round", zorder=1)
        for r in summary[summary["contig"] == c].itertuples():
            ax.add_patch(Rectangle((r.start, yy - 0.22),
                                   max(r.end - r.start, 1), 0.44,
                                   facecolor=CONF_COLORS[r.confidence],
                                   edgecolor=SURFACE, linewidth=0.3,
                                   zorder=2))
    ax.set_yticks(range(len(contigs)))
    ax.set_yticklabels(list(reversed(
        [c if len(c) <= 28 else c[:27] + "…" for c in contigs])),
        fontsize=7, color=INK2)
    ax.set_ylim(-0.6, len(contigs) - 0.4)
    ax.set_xlim(0, max(contig_lengths[c] for c in contigs) * 1.01)
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v / 1e6:.1f}"))
    ax.set_xlabel("position (Mb)", fontsize=8, color=INK2)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(axis="both", colors=MUTED, length=3)
    n_other = summary["contig"].nunique() - len(contigs)
    note = f" (15 largest contigs; {n_other} more not shown)" if n_other > 0 else ""
    ax.set_title(f"{sample + ': ' if sample else ''}{len(summary)} candidate "
                 f"PULs along the genome{note}", fontsize=10, loc="left",
                 color=INK, fontweight="bold")
    counts = summary["confidence"].value_counts()
    ax.legend(handles=[Patch(facecolor=CONF_COLORS[k], edgecolor="none",
                             label=f"{k} confidence ({counts.get(k, 0)})")
                       for k in CONF_COLORS],
              loc="lower right", fontsize=7, frameon=False, ncol=3,
              bbox_to_anchor=(1.0, 1.0), labelcolor=INK2)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_base}.{ext}", dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot candidate PUL gene maps")
    ap.add_argument("--final-dir", required=True)
    ap.add_argument("--genome", default=None,
                    help="Genome FASTA (for contig lengths in the overview)")
    args = ap.parse_args()

    final = Path(args.final_dir)
    figdir = final / "figures"
    figdir.mkdir(exist_ok=True)
    summary = pd.read_csv(final / "PUL_summary.tsv", sep="\t")
    if summary.empty:
        return
    genes = pd.read_csv(final / "gene_annotation_master.tsv", sep="\t",
                        low_memory=False)
    genes["gene_label"] = genes["gene_label"].fillna("")

    cgc = None
    cgc_path = final / "dbcan" / "dbcan_cgc_regions.tsv"
    if cgc_path.exists():
        cgc = pd.read_csv(cgc_path, sep="\t")
    primers = None
    primer_path = final / "primers" / "PUL_primers.tsv"
    if primer_path.exists():
        primers = pd.read_csv(primer_path, sep="\t")

    for _, row in summary.iterrows():
        cid = row["candidate_id"]
        block = genes[genes["candidate_id"] == cid].sort_values("start")
        c = None
        if cgc is not None and not cgc.empty:
            c = cgc[(cgc["contig"] == row["contig"])
                    & (cgc["end"] >= row["start"])
                    & (cgc["start"] <= row["end"])]
        p = primers[primers["candidate_id"] == cid] if primers is not None else None
        plot_locus(row, block, c, p, figdir / cid)

    if args.genome:
        from Bio import SeqIO
        lengths = {r.id: len(r.seq) for r in SeqIO.parse(args.genome, "fasta")}
    else:
        lengths = genes.groupby("contig")["end"].max().to_dict()
    sample = _s(summary["sample"].iloc[0]) if "sample" in summary else ""
    plot_overview(summary, lengths, figdir / "PUL_genome_overview", sample)


if __name__ == "__main__":
    main()
