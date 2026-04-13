# PUL HMM pipeline (MMseqs-free, header-driven marker HMM version)

This version removes the MMseqs2 clustering step that failed with `Illegal instruction` and replaces it with a reference-header-driven HMM builder.

## Core idea

1. Predict genes with Prodigal.
2. Parse `PUL_12112023.faa` headers into a structured reference table.
3. Build focused custom HMMs directly from reference labels:
   - `SusC_like` from `TC` records carrying `1.B.14`
   - `SusD_like` from `TC` records carrying `8.A.46`
   - `CAZyme` families from the CAZyme annotations in the FASTA headers
   - optional `TF`, `STP`, and other transporter groups
4. Search the target proteome with those custom HMMs.
5. Optionally add dbCAN-family HMM hits and Pfam SusC/SusD evidence.
6. Use genomic context to call candidate PULs and export wet-lab-friendly outputs.

## Required software

- python3
- prodigal
- mafft
- hmmer (`hmmbuild`, `hmmpress`, `hmmsearch`; `hmmfetch` only if using `--pfam-hmm`)
- BLAST+ (`makeblastdb`, `blastp`)
- python packages: `pandas`, `biopython`

## Main command

```bash
bash run_pul_hmm_pipeline.sh \
  --genome /vol/projects/psivapor/PMIG_project/BioinfoHelper/PUL_BINE/PUL_db/GCF_014131755.1_ASM1413175v1_genomic.fna \
  --pul-faa /vol/projects/psivapor/PMIG_project/BioinfoHelper/PUL_BINE/PUL_db/PUL_12112023.faa \
  --pul-meta /path/to/dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm /path/to/dbCAN.hmm \
  --pfam-hmm /path/to/Pfam-A.hmm \
  --threads 24 \
  --outdir GCF_014131755.1_PUL_HMM
```

## Output highlights

- `final/PUL_summary.tsv`
- `final/gene_annotation_master.tsv`
- `final/PUL_regions.bed`
- `final/candidate_PULs/`
  - one folder per candidate locus with region FASTA, proteins, CDS, GFF3, and gene table

## Notes

- `--dbcan-hmm` is strongly recommended because it improves CAZyme calling.
- `--pfam-hmm` is optional; the custom HMMs and reference homology already support SusC/SusD discovery.
- This version is better aligned to your reference FASTA because the headers already contain structured biological classes.
