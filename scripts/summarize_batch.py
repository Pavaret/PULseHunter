#!/usr/bin/env python3
"""
summarize_batch.py - combine the per-genome PULseHunter results of a batch
run into cohort tables, three summary figures and an index.html.
"""

import argparse
import html
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

INK, INK2, MUTED, AXIS, GRID, SURFACE = ("#0b0b0b", "#52514e", "#898781",
                                         "#c3c2b7", "#e1e0d9", "#fcfcfb")
CONF_COLORS = {"high": "#104281", "medium": "#3987e5", "low": "#9ec5f4"}
SUPPORT_COLORS = {"both": "#2a78d6", "PULseHunter_only": "#eb6834",
                  "dbCAN3_only": "#1baf7a"}
SUPPORT_NAMES = {"both": "Found by both", "PULseHunter_only": "PULseHunter only",
                 "dbCAN3_only": "dbCAN3 only"}
BLUES = LinearSegmentedColormap.from_list(
    "blues", ["#eef4fd", "#9ec5f4", "#3987e5", "#104281"])

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 8.5,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
    "xtick.color": MUTED, "ytick.color": MUTED,
})


def pretty(name):
    return str(name).replace("_", " ")


def clean_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=3)


def stacked_barh(table, colors, names, title, xlabel, out_base):
    """Horizontal stacked bars, one row per genome, with a 2px surface gap
    between segments and the segment value written inside when it fits."""
    fig, ax = plt.subplots(figsize=(8, 0.9 + 0.42 * len(table)))
    left = pd.Series(0, index=table.index, dtype=float)
    total_max = table.sum(axis=1).max()
    for col in table.columns:
        ax.barh([pretty(i) for i in table.index], table[col], left=left,
                color=colors[col], edgecolor=SURFACE, linewidth=1.2,
                height=0.62, label=names.get(col, col))
        for y, (v, l) in enumerate(zip(table[col], left)):
            if v >= 0.04 * total_max:
                dark = colors[col] in ("#104281", "#2a78d6", "#3987e5",
                                       "#eb6834")
                ax.text(l + v / 2, y, f"{int(v)}", ha="center", va="center",
                        fontsize=7.5, color="#ffffff" if dark else INK)
        left += table[col]
    for y, tot in enumerate(left):
        ax.text(tot + total_max * 0.012, y, f"{int(tot)}", va="center",
                fontsize=8, color=INK2)
    ax.invert_yaxis()
    ax.set_xlim(0, total_max * 1.09)
    ax.set_xlabel(xlabel, color=INK2)
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    clean_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelcolor=INK)
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold",
                 color=INK, pad=22)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(table.columns),
              frameon=False, fontsize=8, labelcolor=INK2, handlelength=1.1)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_base}.{ext}", dpi=200)
    plt.close(fig)


def heatmap(matrix, title, out_base):
    fig, ax = plt.subplots(figsize=(1.9 + 0.62 * matrix.shape[1],
                                    1.4 + 0.27 * matrix.shape[0]))
    vmax = max(1, matrix.to_numpy().max())
    ax.imshow(matrix.to_numpy(), cmap=BLUES, vmin=0, vmax=vmax,
              aspect="auto")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = matrix.iat[i, j]
            if v:
                ax.text(j, i, int(v), ha="center", va="center", fontsize=7.5,
                        color="#ffffff" if v > 0.55 * vmax else INK)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels([pretty(c).replace("Bacteroides ", "B. ")
                        for c in matrix.columns], rotation=40, ha="right",
                       color=INK, fontstyle="italic")
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, color=INK)
    ax.set_xticks([x - 0.5 for x in range(1, matrix.shape[1])], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, matrix.shape[0])], minor=True)
    ax.grid(which="minor", color=SURFACE, lw=2)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold",
                 color=INK)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_base}.{ext}", dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Summarise a PULseHunter batch")
    ap.add_argument("--batch-dir", required=True)
    args = ap.parse_args()

    batch = Path(args.batch_dir)
    out = batch / "summary"
    out.mkdir(exist_ok=True)
    samples = sorted(p.parent.parent.name
                     for p in batch.glob("*/final/PUL_summary.tsv"))
    if not samples:
        raise SystemExit(f"No finished genomes under {batch}")

    summ, prim, conc, loci_conc, stats = [], [], [], [], []
    for s in samples:
        final = batch / s / "final"
        df = pd.read_csv(final / "PUL_summary.tsv", sep="\t")
        df["sample"] = s
        summ.append(df)
        genes = pd.read_csv(final / "gene_annotation_master.tsv", sep="\t",
                            low_memory=False)
        row = {
            "sample": s, "genes": len(genes),
            "contigs": genes["contig"].nunique(),
            "candidate_PULs": len(df),
            **{f"{k}_confidence": int((df["confidence"] == k).sum())
               for k in CONF_COLORS},
            "susCD_pairs": genes.loc[genes["susCD_pair"].notna()
                                     & (genes["susCD_pair"] != ""),
                                     "susCD_pair"].nunique(),
            "cazymes": int(genes["is_cazyme"].sum()),
            "cazymes_in_loci": int((genes["is_cazyme"] & genes[
                "candidate_id"].notna() & (genes["candidate_id"] != "")).sum()),
            "loci_strong_substrate": int(
                (df["substrate_support"] == "strong").sum()),
            "loci_at_contig_edge": int(df["at_contig_edge"].sum()),
        }
        p = final / "primers" / "PUL_primers.tsv"
        if p.exists():
            pdf = pd.read_csv(p, sep="\t")
            pdf["sample"] = s
            prim.append(pdf)
            per = pdf.groupby("candidate_id")["specific"].sum()
            row["loci_with_specific_primers"] = int((per >= 1).sum())
            row["loci_with_3_specific_sets"] = int((per >= 3).sum())
        c = final / "dbcan" / "concordance_summary.tsv"
        if c.exists():
            cdf = pd.read_csv(c, sep="\t")
            cdf["sample"] = s
            conc.append(cdf)
            row.update({k: v for k, v in zip(cdf["metric"], cdf["value"])
                        if k.startswith(("dbcan3", "pulsehunter_only",
                                         "pulsehunter_loci_supported",
                                         "cazymes_"))})
            lc = pd.read_csv(final / "dbcan" / "locus_concordance.tsv",
                             sep="\t")
            lc["sample"] = s
            loci_conc.append(lc)
        stats.append(row)

    summary = pd.concat(summ, ignore_index=True)
    stats = pd.DataFrame(stats)
    summary.to_csv(out / "all_PUL_summary.tsv", sep="\t", index=False)
    stats.to_csv(out / "genome_stats.tsv", sep="\t", index=False)
    if prim:
        pd.concat(prim, ignore_index=True).to_csv(
            out / "all_PUL_primers.tsv", sep="\t", index=False)
    if loci_conc:
        pd.concat(loci_conc, ignore_index=True).to_csv(
            out / "all_locus_concordance.tsv", sep="\t", index=False)

    figs = []
    conf = (summary.groupby(["sample", "confidence"]).size()
            .unstack(fill_value=0).reindex(index=samples,
                                           columns=list(CONF_COLORS),
                                           fill_value=0))
    stacked_barh(conf, CONF_COLORS,
                 {k: f"{k} confidence" for k in CONF_COLORS},
                 "Candidate PULs per genome", "candidate loci",
                 out / "loci_per_genome")
    figs.append(("loci_per_genome", "Candidate PULs per genome by confidence"))

    if loci_conc:
        lc = pd.concat(loci_conc, ignore_index=True)
        sup = (lc.groupby(["sample", "category"]).size().unstack(fill_value=0)
               .reindex(index=samples, columns=list(SUPPORT_COLORS),
                        fill_value=0))
        stacked_barh(sup, SUPPORT_COLORS, SUPPORT_NAMES,
                     "PULseHunter loci vs dbCAN3 CGCs", "loci",
                     out / "pulsehunter_vs_dbcan3")
        figs.append(("pulsehunter_vs_dbcan3",
                     "Agreement between PULseHunter and dbCAN3"))

    good = summary[summary["substrate_support"].isin(["strong", "moderate"])
                   & summary["predicted_substrate"].notna()
                   & (summary["predicted_substrate"] != "")]
    if len(good):
        mat = (good.groupby(["predicted_substrate", "sample"]).size()
               .unstack(fill_value=0).reindex(columns=samples, fill_value=0))
        mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index]
        mat.to_csv(out / "substrate_matrix.tsv", sep="\t")
        top = mat.head(30)
        heatmap(top, "Loci per predicted substrate"
                     + (" (top 30)" if len(mat) > 30 else ""),
                out / "substrate_heatmap")
        figs.append(("substrate_heatmap",
                     "Predicted substrates across genomes"))

    # index.html --------------------------------------------------------
    show = [c for c in [
        "sample", "contigs", "genes", "candidate_PULs", "high_confidence",
        "medium_confidence", "low_confidence", "susCD_pairs", "cazymes",
        "loci_strong_substrate", "dbcan3_cgcs",
        "pulsehunter_loci_supported_by_dbcan3", "pulsehunter_only_loci",
        "dbcan3_only_cgcs", "loci_with_specific_primers",
        "loci_with_3_specific_sets"] if c in stats.columns]
    head = "".join(f"<th>{html.escape(pretty(c))}</th>" for c in show)
    rows = []
    for _, r in stats.iterrows():
        tds = []
        for c in show:
            if c == "sample":
                link = f"../{r[c]}/{r[c]}_PULseHunter_report.html"
                tds.append(f'<td><a href="{link}"><i>{html.escape(pretty(r[c]))}'
                           "</i></a></td>")
            else:
                v = r[c]
                tds.append('<td class="num">'
                           + ("" if pd.isna(v) else f"{int(v):,}") + "</td>")
        rows.append("<tr>" + "".join(tds) + "</tr>")
    fig_html = "".join(
        f'<div class="card"><img src="{n}.svg" alt="{html.escape(a)}"></div>'
        for n, a in figs)
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PULseHunter batch summary</title><style>
body{{margin:0;background:#f9f9f7;color:#0b0b0b;font:14px/1.45 system-ui,sans-serif}}
main{{max-width:1200px;margin:0 auto;padding:24px 16px 64px}}
h1{{font-size:22px;margin:0 0 4px}}p{{color:#52514e;margin:0 0 16px}}
.card{{background:#fcfcfb;border:1px solid #e1e0d9;border-radius:8px;padding:12px 14px;margin:12px 0;overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:12.5px}}
th,td{{padding:6px 8px;border-bottom:1px solid #e1e0d9;text-align:left}}
th{{color:#52514e;font-weight:600}}td.num{{text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto}}a{{color:#2a78d6;text-decoration:none}}
</style></head><body><main><h1>PULseHunter batch summary</h1>
<p>{len(samples)} genomes, {len(summary):,} candidate PULs. Click a genome for its full report.
Tables are in this folder (<code>all_PUL_summary.tsv</code>, <code>genome_stats.tsv</code>,
<code>all_PUL_primers.tsv</code>, <code>all_locus_concordance.tsv</code>).</p>
<div class="card"><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>
{fig_html}</main></body></html>"""
    (out / "index.html").write_text(doc)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
