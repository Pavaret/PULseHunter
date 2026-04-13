#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd
from Bio import SeqIO


def clean_name(x: str) -> str:
    x = str(x)
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    x = re.sub(r"_+", "_", x).strip("_")
    return x or "unknown"


def tc_family_code(annot: str) -> str:
    m = re.search(r"\|([0-9]+\.[A-Z]\.[0-9]+)(?:\.|$)", annot)
    if m:
        return m.group(1)
    return "TC_unknown"


def classify_marker(row):
    ref_class = str(row.get("ref_class", "")).strip()
    annot = str(row.get("ref_annotation", "")).strip()
    gene_symbol = str(row.get("gene_symbol", "")).strip()

    if ref_class == "TC" and re.search(r"(^|[^0-9])1\.B\.14(\.|$|[^0-9])", annot):
        return ("SusC_like", "SusC_like")
    if ref_class == "TC" and re.search(r"(^|[^0-9])8\.A\.46(\.|$|[^0-9])", annot):
        return ("SusD_like", "SusD_like")

    if ref_class == "CAZyme":
        fam = annot if annot else gene_symbol if gene_symbol else "CAZyme_misc"
        return ("CAZyme", clean_name(fam))

    if ref_class == "TF":
        fam = annot if annot else gene_symbol if gene_symbol else "TF_misc"
        return ("TF", clean_name(fam))

    if ref_class == "STP":
        fam = annot if annot else gene_symbol if gene_symbol else "STP_misc"
        return ("STP", clean_name(fam))

    if ref_class == "TC":
        fam = tc_family_code(annot)
        return ("TC_other", clean_name(fam))

    return (None, None)


def main():
    ap = argparse.ArgumentParser(description="Build focused marker-group FASTAs and a shell script for HMM building")
    ap.add_argument("--faa", required=True)
    ap.add_argument("--reference-table", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-seed-size", type=int, default=5)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    fasta_dir = outdir / "cluster_fastas"
    aln_dir = outdir / "alignments"
    hmm_dir = outdir / "hmms"
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)
    aln_dir.mkdir(parents=True, exist_ok=True)
    hmm_dir.mkdir(parents=True, exist_ok=True)

    seqs = {rec.id: rec for rec in SeqIO.parse(args.faa, "fasta")}
    ref = pd.read_csv(args.reference_table, sep="\t", dtype=str).fillna("")

    groups = defaultdict(list)
    group_puls = defaultdict(Counter)
    for _, row in ref.iterrows():
        marker_class, marker_group = classify_marker(row)
        if marker_class is None:
            continue
        fid = row["full_id"]
        if fid not in seqs:
            continue
        groups[(marker_class, marker_group)].append(fid)
        group_puls[(marker_class, marker_group)][row.get("pul_id", "")] += 1

    manifest_rows = []
    skipped_rows = []
    idx = 0
    for marker_class, marker_group in sorted(groups):
        members = groups[(marker_class, marker_group)]
        if len(members) < args.min_seed_size:
            skipped_rows.append({
                "marker_class": marker_class,
                "marker_group": marker_group,
                "n_members": len(members),
                "reason": f"below_min_seed_size_{args.min_seed_size}",
            })
            continue
        idx += 1
        cluster_id = f"MKR{idx:05d}"
        fasta_path = fasta_dir / f"{cluster_id}.faa"
        with open(fasta_path, "w") as out:
            for mid in members:
                SeqIO.write(seqs[mid], out, "fasta")

        pul_counts = group_puls[(marker_class, marker_group)]
        dominant_pul = pul_counts.most_common(1)[0][0] if pul_counts else ""
        manifest_rows.append({
            "cluster_id": cluster_id,
            "n_members": len(members),
            "mode": "hmm",
            "majority_class": marker_class,
            "dominant_annotation": marker_group,
            "dominant_pul": dominant_pul,
            "pul_diversity": len(pul_counts),
            "members": ";".join(members),
            "marker_class": marker_class,
            "marker_group": marker_group,
            "fasta": str(fasta_path),
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(outdir / "cluster_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(skipped_rows).to_csv(outdir / "skipped_groups.tsv", sep="\t", index=False)

    build_sh = outdir / "build_hmms.sh"
    with open(build_sh, "w") as fh:
        fh.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        fh.write('BASE="$(cd "$(dirname "$0")" && pwd)"\n')
        fh.write('mkdir -p "$BASE/alignments" "$BASE/hmms"\n')
        fh.write('if [[ ! -s "$BASE/cluster_manifest.tsv" ]]; then echo "No cluster_manifest.tsv entries found" >&2; exit 1; fi\n')
        fh.write('tail -n +2 "$BASE/cluster_manifest.tsv" | while IFS=$\'\\t\' read -r cluster_id n_members mode majority_class dominant_annotation dominant_pul pul_diversity members marker_class marker_group fasta; do\n')
        fh.write('  [[ -z "$cluster_id" ]] && continue\n')
        fh.write('  mafft --auto "$fasta" > "$BASE/alignments/${cluster_id}.aln.faa"\n')
        fh.write('  hmmbuild -n "${marker_class}__${marker_group}" "$BASE/hmms/${cluster_id}.hmm" "$BASE/alignments/${cluster_id}.aln.faa" > "$BASE/hmms/${cluster_id}.hmmbuild.log"\n')
        fh.write('done\n')
        fh.write('shopt -s nullglob\n')
        fh.write('hmms=("$BASE"/hmms/*.hmm)\n')
        fh.write('(( ${#hmms[@]} > 0 )) || { echo "No HMM files produced" >&2; exit 1; }\n')
        fh.write('cat "${hmms[@]}" > "$BASE/custom_pul_profiles.hmm"\n')
        fh.write('hmmpress "$BASE/custom_pul_profiles.hmm" > "$BASE/hmmpress.log"\n')
    build_sh.chmod(0o755)


if __name__ == "__main__":
    main()
