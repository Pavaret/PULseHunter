# PULseHunter

Find, annotate, visualise and design primers for polysaccharide utilization
loci (PULs) in bacterial genomes.

PULseHunter runs two complementary arms on the same gene calls and reports
where they agree:

| Arm | What it does |
|---|---|
| **PULseHunter** | Marker HMMs built from dbCAN-PUL, BLASTP against dbCAN-PUL proteins, dbCAN family HMMs and Pfam SusC/SusD domains, combined with genomic context (susC/susD pairs) to call loci and transfer a substrate from the best-matching characterised PUL. |
| **dbCAN3** (optional) | `run_dbcan` (HMMER + DIAMOND) with CGCFinder, giving CAZyme gene clusters (CGCs). |

For every candidate locus you get a gene map, a gene table, sequences, and
the three best primer sets, plus one HTML report per genome.

## Install

```bash
micromamba create -f envs/pulsehunter.yml      # or mamba / conda
micromamba activate pulsehunter
```

The dbCAN3 arm uses a separate `run_dbcan` installation (its own conda
environment is fine; pass it with `--dbcan-env`).

Databases:

- dbCAN-PUL proteins and metadata: `files/resources.tar.xz`
  (`PUL_12112023.faa`, `dbCAN-PUL_Feb-2025.tsv`)
- dbCAN family HMMs (`dbCAN-HMMdb-V*.hmm`) for CAZyme calling
- `Pfam-A.hmm` for SusC/SusD domain evidence
- the `run_dbcan` database directory, only for `--run-dbcan`

## Run one genome

```bash
bash pulsehunter.sh \
  --genome genome.fna \
  --outdir results/my_genome \
  --pul-faa PUL_12112023.faa \
  --pul-meta dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm dbCAN-HMMdb-V14.hmm \
  --pfam-hmm Pfam-A.hmm \
  --run-dbcan --dbcan-db /path/to/dbcan/db --dbcan-env /path/to/envs/run_dbcan \
  --threads 16
```

## Run a folder of genomes

```bash
bash pulsehunter_batch.sh \
  --genome-dir genomes/ \
  --outdir results/ \
  --pul-faa PUL_12112023.faa --pul-meta dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm dbCAN-HMMdb-V14.hmm --pfam-hmm Pfam-A.hmm \
  --run-dbcan --dbcan-db /path/to/dbcan/db --dbcan-env /path/to/envs/run_dbcan \
  --community-specificity \
  --jobs 8 --threads 4
```

The reference (marker HMMs, BLAST database) is built once in
`results/reference/` and shared by all genomes. Finished steps are skipped
on re-run; use `--force` to recompute.

`--community-specificity` checks every primer against all other genomes of
the batch, so the primers are specific within a defined community.

## Outputs

Per genome (`results/<genome>/`):

| File | Content |
|---|---|
| `<genome>_PULseHunter_report.html` | Self-contained report: headline numbers, genome overview, filterable locus table, gene map + primers + gene table per locus |
| `final/PUL_summary.tsv` | One row per candidate locus: coordinates, confidence, type, susC/susD pairs, CAZyme families, best dbCAN-PUL match, predicted substrate and its support, dbCAN3 overlap |
| `final/gene_annotation_master.tsv` | Every gene with all evidence columns |
| `final/PUL_regions.bed` | Locus coordinates |
| `final/candidate_PULs/<locus>/` | Region FASTA, proteins, CDS, GFF3, gene table, primers |
| `final/figures/` | Gene map per locus (PNG + SVG) and a genome overview |
| `final/primers/PUL_primers.tsv`, `.fasta` | Top primer sets per locus with in-silico PCR results |
| `final/dbcan/` | dbCAN3 CGC regions, locus and CAZyme concordance |
| `dbcan3/run/` | Raw `run_dbcan` output |

Batch (`results/summary/`): `index.html`, `all_PUL_summary.tsv`,
`genome_stats.tsv`, `all_PUL_primers.tsv`, `all_locus_concordance.tsv`,
`substrate_matrix.tsv` and summary figures.

## How loci are called

1. **Genes**: Prodigal.
2. **Evidence per gene**
   - SusC-like: TonB-dependent transporter evidence (Pfam Plug + barrel
     domains, the dbCAN-PUL SusC marker HMM, or a 1.B.14 homolog).
   - SusD-like: Pfam SusD families, the SusD marker HMM, or an 8.A.46 homolog.
   - CAZyme: dbCAN family HMM (i-Evalue <= 1e-15, coverage >= 0.35), CAZyme
     marker HMM, or a dbCAN-PUL CAZyme homolog (>= 30 % identity, >= 50 %
     query coverage).
   - TF / STP / other transporters: marker HMMs and dbCAN-PUL homologs.
3. **susC/susD pairs**: a SusC-like and a SusD-like gene on the same strand,
   at most one gene apart. Unpaired TonB-dependent transporters are not
   treated as SusC.
4. **Loci**: pair members, SusD-like genes and degradative CAZymes (GH, PL,
   CE, CBM, AA) are chained when at most 4 genes apart (at most 2 of them
   without any signal) and at most 10 kb apart. Loci are extended over
   directly adjacent regulators and transporters.
5. **Confidence**

   | Level | Rule |
   |---|---|
   | high | susC/susD pair + degradative CAZyme (canonical PUL) |
   | medium | susC/susD pair without CAZyme, or SusD / TonB-dependent transporter with >= 2 CAZymes |
   | low | CAZyme cluster (>= 2 CAZymes + transporter or regulator) without susC/susD |

6. **Substrate**: every dbCAN-PUL entry is scored by the summed best
   bitscore of the locus genes. `substrate_support` is `strong` (>= 3 genes,
   >= 50 % of the reference PUL recovered, mean identity >= 50 %),
   `moderate` (>= 2 genes, >= 25 %) or `weak` (a hint only).

All thresholds are options of `scripts/call_puls_hmm.py`.

## Primers

`--primer-mode qpcr` (default) designs 80-200 bp amplicons in the locus
marker genes (susC, then susD, then CAZymes); `--primer-mode pcr` designs
400-1000 bp amplicons across the susC/susD junction. Primer3 proposes
candidates (Tm 60 °C, 18-24 nt, 40-60 % GC); each is then tested by
in-silico PCR against the target genome and any `--offtarget-fasta`
genomes. A primer is predicted to prime wherever its 3'-terminal 10 bases
match exactly with at most 3 mismatches upstream. The best
`--n-primer-sets` (default 3) non-overlapping sets per locus are reported,
ranked by specificity, marker priority and Primer3 penalty. `specific` means
exactly one predicted product in the target genome and none elsewhere.

In-silico specificity is a prediction; validate primers experimentally.

## Requirements

`python >= 3.10`, `prodigal`, `hmmer`, `mafft`, `blast`, `pandas`,
`biopython`, `matplotlib`, `primer3-py` (all in `envs/pulsehunter.yml`).
