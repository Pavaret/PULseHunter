#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import Counter
import pandas as pd
from Bio import SeqIO


def parse_header(full_id: str):
    parts = full_id.split(":")
    parts += [""] * max(0, 7 - len(parts))
    ref_member_id, pul_id, gene_symbol, locus_tag, protein_acc, ref_class, ref_annotation = parts[:7]
    return {
        "full_id": full_id,
        "ref_member_id": ref_member_id,
        "pul_id": pul_id,
        "gene_symbol": gene_symbol,
        "locus_tag": locus_tag,
        "protein_accession": protein_acc,
        "ref_class": ref_class or "unknown",
        "ref_annotation": ref_annotation,
    }


def main():
    ap = argparse.ArgumentParser(description="Prepare dbCAN-PUL reference tables from FASTA + metadata TSV")
    ap.add_argument("--pul-faa", required=True)
    ap.add_argument("--pul-meta", default=None)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = []
    seq_count = 0
    pul_ids = set()
    classes = Counter()
    with open(outdir / "reference_proteins.faa", "w") as out_faa:
        for rec in SeqIO.parse(args.pul_faa, "fasta"):
            seq_count += 1
            meta = parse_header(rec.id)
            meta["seq_len_aa"] = len(rec.seq)
            records.append(meta)
            pul_ids.add(meta["pul_id"])
            classes[meta["ref_class"]] += 1
            SeqIO.write(rec, out_faa, "fasta")

    ref_df = pd.DataFrame(records)

    meta_df = None
    if args.pul_meta:
        meta_df = pd.read_csv(args.pul_meta, sep="\t", dtype=str).fillna("")
        if "ID" in meta_df.columns and "pul_id" not in meta_df.columns:
            meta_df = meta_df.rename(columns={"ID": "pul_id"})
        merged = ref_df.merge(meta_df, on="pul_id", how="left", suffixes=("", "_meta"))
    else:
        merged = ref_df.copy()

    merged.to_csv(outdir / "reference_proteins.tsv", sep="\t", index=False)

    ref_puls = pd.DataFrame({
        "metric": [
            "reference_sequences",
            "reference_puls_in_faa",
            "reference_classes",
        ],
        "value": [
            seq_count,
            len(pul_ids),
            "; ".join(f"{k}={v}" for k, v in classes.most_common()),
        ]
    })

    if meta_df is not None:
        meta_puls = set(meta_df["pul_id"].astype(str))
        overlap = pul_ids & meta_puls
        ref_puls = pd.concat([
            ref_puls,
            pd.DataFrame({
                "metric": [
                    "metadata_rows",
                    "metadata_puls",
                    "faa_meta_overlap_puls",
                    "faa_only_puls",
                    "meta_only_puls",
                ],
                "value": [
                    len(meta_df),
                    len(meta_puls),
                    len(overlap),
                    len(pul_ids - meta_puls),
                    len(meta_puls - pul_ids),
                ],
            })
        ], ignore_index=True)

        pd.DataFrame(sorted(pul_ids - meta_puls), columns=["pul_id"]).to_csv(
            outdir / "puls_present_only_in_faa.tsv", sep="\t", index=False
        )
        pd.DataFrame(sorted(meta_puls - pul_ids), columns=["pul_id"]).to_csv(
            outdir / "puls_present_only_in_metadata.tsv", sep="\t", index=False
        )

    ref_puls.to_csv(outdir / "reference_qc.tsv", sep="\t", index=False)

    summary_lines = [
        f"reference_sequences\t{seq_count}",
        f"reference_puls_in_faa\t{len(pul_ids)}",
        f"reference_classes\t{' ; '.join(f'{k}={v}' for k, v in classes.most_common())}",
    ]
    if meta_df is not None:
        summary_lines.extend([
            f"metadata_rows\t{len(meta_df)}",
            f"metadata_puls\t{len(meta_puls)}",
            f"faa_meta_overlap_puls\t{len(overlap)}",
            f"faa_only_puls\t{len(pul_ids - meta_puls)}",
            f"meta_only_puls\t{len(meta_puls - pul_ids)}",
        ])

    with open(outdir / "README_reference.txt", "w") as fh:
        fh.write("Prepared PUL reference database\n")
        fh.write("================================\n")
        fh.write("\n".join(summary_lines) + "\n")
        fh.write("\nreference_proteins.tsv is the master lookup table used by the downstream pipeline.\n")
        fh.write("When FASTA and metadata disagree, sequence IDs in the FASTA are treated as authoritative for HMM building.\n")


if __name__ == "__main__":
    main()
