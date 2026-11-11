#!/usr/bin/env python3
"""
make_report.py - one self-contained HTML report per genome: headline
numbers, genome overview, a filterable locus table and, for every candidate
PUL, its gene map, gene table and primer sets.
"""

import argparse
import base64
import html
from pathlib import Path

import pandas as pd

CSS = """
:root{--surface:#fcfcfb;--page:#f9f9f7;--ink:#0b0b0b;--ink2:#52514e;
--muted:#898781;--line:#e1e0d9;--accent:#2a78d6}
*{box-sizing:border-box}
body{margin:0;background:var(--page);color:var(--ink);
font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1280px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 10px}
h3{font-size:15px;margin:0}
p.sub{color:var(--ink2);margin:0 0 18px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
gap:10px;margin:14px 0}
.tile{background:var(--surface);border:1px solid var(--line);
border-radius:8px;padding:10px 12px}
.tile b{display:block;font-size:22px;font-weight:600}
.tile span{color:var(--ink2);font-size:12px}
.card{background:var(--surface);border:1px solid var(--line);
border-radius:8px;padding:12px 14px;margin:12px 0}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left;
vertical-align:top;white-space:nowrap}
td.wrap{white-space:normal;min-width:180px}
th{color:var(--ink2);font-weight:600;position:sticky;top:0;
background:var(--surface)}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
img{max-width:100%;height:auto;display:block}
.fig{overflow-x:auto}.fig img{max-width:none}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
.badge{display:inline-block;padding:1px 7px;border-radius:10px;
font-size:11.5px;border:1px solid var(--line);color:var(--ink2)}
.badge.high{background:#104281;color:#fff;border-color:#104281}
.badge.medium{background:#3987e5;color:#fff;border-color:#3987e5}
.badge.low{background:#cde2fb;color:#0b0b0b;border-color:#9ec5f4}
.controls{display:flex;gap:10px;flex-wrap:wrap;align-items:center;
margin:8px 0}
input[type=search],select{font:inherit;padding:5px 8px;
border:1px solid var(--line);border-radius:6px;background:var(--surface)}
details{margin-top:8px}summary{cursor:pointer;color:var(--ink2)}
code{font:12px ui-monospace,Menlo,Consolas,monospace}
.meta{color:var(--ink2);font-size:12.5px;margin:4px 0 8px}
.yes{color:#006300}.no{color:#b3261e}
"""

JS = """
const q=document.getElementById('q'),conf=document.getElementById('conf');
function apply(){const t=q.value.toLowerCase(),c=conf.value;let n=0;
document.querySelectorAll('[data-locus]').forEach(el=>{
const ok=(!c||el.dataset.conf===c)&&(!t||el.dataset.text.includes(t));
el.style.display=ok?'':'none';if(ok&&el.tagName==='TR')n++;});
document.getElementById('count').textContent=n+' loci shown';}
q.addEventListener('input',apply);conf.addEventListener('change',apply);apply();
"""


def esc(v):
    return "" if pd.isna(v) else html.escape(str(v))


def img_tag(path, alt):
    data = base64.b64encode(Path(path).read_bytes()).decode()
    return (f'<img loading="lazy" alt="{html.escape(alt)}" '
            f'src="data:image/svg+xml;base64,{data}">')


def table(df, cols, numeric=(), wrap=()):
    head = "".join(
        f'<th class="{"num" if c in numeric else ""}">{html.escape(lab)}</th>'
        for c, lab in cols)
    body = []
    for _, r in df.iterrows():
        tds = []
        for c, _ in cols:
            cls = "num" if c in numeric else ("wrap" if c in wrap else "")
            tds.append(f'<td class="{cls}">{esc(r.get(c, ""))}</td>')
        body.append("<tr>" + "".join(tds) + "</tr>")
    return (f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>')


def main():
    ap = argparse.ArgumentParser(description="Build the HTML report")
    ap.add_argument("--final-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    final = Path(args.final_dir)
    summary = pd.read_csv(final / "PUL_summary.tsv", sep="\t")
    genes = pd.read_csv(final / "gene_annotation_master.tsv", sep="\t",
                        low_memory=False)
    primer_path = final / "primers" / "PUL_primers.tsv"
    primers = (pd.read_csv(primer_path, sep="\t") if primer_path.exists()
               else pd.DataFrame())
    conc_path = final / "dbcan" / "concordance_summary.tsv"
    conc = (dict(pd.read_csv(conc_path, sep="\t")[["metric", "value"]].values)
            if conc_path.exists() else {})
    sample = esc(summary["sample"].iloc[0]) if len(summary) else ""
    figdir = final / "figures"

    n_pairs = genes.loc[genes["susCD_pair"].notna()
                        & (genes["susCD_pair"] != ""), "susCD_pair"].nunique()
    tiles = [
        (len(summary), "candidate PULs"),
        (int((summary.get("confidence") == "high").sum()), "high confidence"),
        (int((summary.get("confidence") == "medium").sum()),
         "medium confidence"),
        (int((summary.get("confidence") == "low").sum()), "low confidence"),
        (n_pairs, "susC/susD pairs"),
        (int(genes["is_cazyme"].sum()), "CAZymes (PULseHunter)"),
    ]
    if conc:
        tiles += [
            (conc["dbcan3_cgcs"], "dbCAN3 CGCs"),
            (conc["pulsehunter_loci_supported_by_dbcan3"],
             "PULs also found by dbCAN3"),
            (conc["pulsehunter_only_loci"], "PULseHunter-only loci"),
            (conc["dbcan3_only_cgcs"], "dbCAN3-only CGCs"),
        ]
    if len(primers):
        ok = int((primers.groupby("candidate_id")["specific"].sum() >= 1).sum())
        tiles.append((ok, "loci with a specific primer set"))

    parts = [f"<h1>PULseHunter report: {sample}</h1>",
             '<p class="sub">Candidate polysaccharide utilization loci, '
             'their dbCAN3 support, gene maps and primer sets.</p>',
             '<div class="tiles">' + "".join(
                 f'<div class="tile"><b>{v}</b><span>{html.escape(k)}</span>'
                 '</div>' for v, k in tiles) + "</div>"]
    ov = figdir / "PUL_genome_overview.svg"
    if ov.exists():
        parts.append(f'<div class="card">{img_tag(ov, "Genome overview")}</div>')

    # locus table -----------------------------------------------------
    parts.append("<h2>Candidate loci</h2>")
    parts.append(
        '<div class="controls"><input id="q" type="search" '
        'placeholder="Filter: substrate, CAZy family, contig…" size="38">'
        '<select id="conf"><option value="">All confidence levels</option>'
        '<option>high</option><option>medium</option><option>low</option>'
        '</select><span id="count" class="meta"></span></div>')
    has_db = "call_support" in summary.columns
    head = ["Locus", "Contig", "Start", "End", "Size (kb)", "Genes",
            "Confidence", "Type", "susC/D pairs", "CAZymes",
            "CAZy families", "Best dbCAN-PUL match", "Substrate (support)"]
    if has_db:
        head.append("dbCAN3")
    rows = []
    for r in summary.itertuples():
        text = " ".join(str(x) for x in (
            r.candidate_id, r.contig, r.pul_type, r.cazyme_families,
            r.predicted_substrate, r.best_reference_pul)).lower()
        sub = esc(r.predicted_substrate)
        sub = f"{sub} ({esc(r.substrate_support)})" if sub else ""
        tds = [
            f'<td><a href="#{r.candidate_id}">{r.candidate_id}</a></td>',
            f"<td>{esc(r.contig)}</td>",
            f'<td class="num">{r.start:,}</td>',
            f'<td class="num">{r.end:,}</td>',
            f'<td class="num">{r.length_bp / 1000:.1f}</td>',
            f'<td class="num">{r.num_genes}</td>',
            f'<td><span class="badge {r.confidence}">{r.confidence}</span></td>',
            f"<td>{esc(r.pul_type).replace('_', ' ')}</td>",
            f'<td class="num">{r.susCD_pairs}</td>',
            f'<td class="num">{r.cazyme_count}</td>',
            f'<td class="wrap">{esc(r.cazyme_families).replace(",", ", ")}</td>',
            f"<td>{esc(r.best_reference_pul)}</td>",
            f'<td class="wrap">{sub}</td>',
        ]
        if has_db:
            tds.append(f"<td>{esc(r.dbcan3_cgcs).replace(';', ', ') or '–'}</td>")
        rows.append(f'<tr data-locus data-conf="{r.confidence}" '
                    f'data-text="{html.escape(text)}">{"".join(tds)}</tr>')
    parts.append('<div class="card scroll" style="max-height:520px;'
                 'overflow-y:auto"><table><thead><tr>'
                 + "".join(f"<th>{h}</th>" for h in head)
                 + "</tr></thead><tbody>" + "".join(rows)
                 + "</tbody></table></div>")

    # per-locus sections ----------------------------------------------
    parts.append("<h2>Locus maps and primers</h2>")
    gene_cols = [("gene_id", "Gene"), ("start", "Start"), ("end", "End"),
                 ("strand", "Strand"), ("protein_len_aa", "aa"),
                 ("gene_class", "Class"), ("gene_label", "Label"),
                 ("cazyme_evidence", "CAZyme evidence"),
                 ("top_blast_pul", "Top dbCAN-PUL hit"),
                 ("top_blast_gene", "Ref. gene"),
                 ("top_blast_pident", "% id")]
    if "dbcan3_hmmer" in genes.columns:
        gene_cols += [("dbcan3_hmmer", "dbCAN3 HMMER"),
                      ("dbcan3_diamond", "dbCAN3 DIAMOND")]
    primer_cols = [("rank", "Set"), ("target_label", "Target"),
                   ("target_gene", "Gene"), ("fwd_primer", "Forward 5'→3'"),
                   ("rev_primer", "Reverse 5'→3'"), ("fwd_tm", "Tm F"),
                   ("rev_tm", "Tm R"), ("product_size", "Product (bp)"),
                   ("amplicon_start", "Start"), ("amplicon_end", "End"),
                   ("on_target_amplicons", "On-target"),
                   ("off_target_amplicons", "Off-target"),
                   ("off_target_genomes", "Off-target genomes"),
                   ("specific", "Specific")]
    for r in summary.itertuples():
        text = " ".join(str(x) for x in (
            r.candidate_id, r.contig, r.pul_type, r.cazyme_families,
            r.predicted_substrate, r.best_reference_pul)).lower()
        sec = [f'<div class="card" id="{r.candidate_id}" data-locus '
               f'data-conf="{r.confidence}" data-text="{html.escape(text)}">',
               f'<h3>{r.candidate_id} <span class="badge {r.confidence}">'
               f'{r.confidence}</span></h3>',
               f'<div class="meta">{esc(r.contig)}:{r.start:,}-{r.end:,} · '
               f'<code>{esc(r.gene_order)}</code></div>']
        svg = figdir / f"{r.candidate_id}.svg"
        if svg.exists():
            sec.append(f'<div class="fig">{img_tag(svg, r.candidate_id)}</div>')
        if len(primers):
            p = primers[primers["candidate_id"] == r.candidate_id]
            if len(p):
                sec.append('<details open><summary>Primer sets</summary>'
                           + table(p, primer_cols,
                                   numeric=("fwd_tm", "rev_tm", "product_size",
                                            "amplicon_start", "amplicon_end",
                                            "on_target_amplicons",
                                            "off_target_amplicons"))
                           + "</details>")
            else:
                sec.append('<div class="meta">No primer set passed the '
                           'design constraints for this locus.</div>')
        block = genes[genes["candidate_id"] == r.candidate_id].sort_values(
            "start")
        sec.append("<details><summary>Genes</summary>"
                   + table(block, gene_cols,
                           numeric=("start", "end", "protein_len_aa",
                                    "top_blast_pident")) + "</details>")
        sec.append("</div>")
        parts.append("".join(sec))

    doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
           f'<meta name="viewport" content="width=device-width,'
           f'initial-scale=1"><title>PULseHunter {sample}</title>'
           f"<style>{CSS}</style></head><body><main>{''.join(parts)}</main>"
           f"<script>{JS}</script></body></html>")
    Path(args.out).write_text(doc)


if __name__ == "__main__":
    main()
