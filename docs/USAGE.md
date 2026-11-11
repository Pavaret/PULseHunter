# Usage

Activate the environment first (`micromamba activate pulsehunter`). See
[INSTALL.md](INSTALL.md) for the databases.

- [One genome](#one-genome)
- [A folder of genomes](#a-folder-of-genomes)
- [Options](#options)
- [Outputs](#outputs)
- [Reading the results](#reading-the-results)
- [Re-running and tuning](#re-running-and-tuning)
- [Troubleshooting](#troubleshooting)

## One genome

```bash
bash pulsehunter.sh \
  --genome genome.fna \
  --outdir results/my_genome \
  --pul-faa resources/PUL_12112023.faa \
  --pul-meta resources/dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm /path/to/dbCAN-HMMdb-V14.hmm \
  --pfam-hmm /path/to/Pfam-A.hmm \
  --threads 16
```

With the dbCAN3 arm:

```bash
bash pulsehunter.sh \
  --genome genome.fna \
  --outdir results/my_genome \
  --pul-faa resources/PUL_12112023.faa \
  --pul-meta resources/dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm /path/to/dbCAN-HMMdb-V14.hmm \
  --pfam-hmm /path/to/Pfam-A.hmm \
  --run-dbcan --dbcan-db /path/to/dbcan/db --dbcan-env /path/to/envs/run_dbcan \
  --threads 16
```

The input is a nucleotide FASTA (complete genome or contigs). Genes are
predicted by Prodigal, so no annotation file is needed.

### Sharing the reference between runs

The marker HMMs and the BLAST database are built from `--pul-faa` on the
first run. To build them once and reuse them:

```bash
# build only
bash pulsehunter.sh --prepare-ref-only --ref-dir refs/pulsehunter_ref \
  --pul-faa resources/PUL_12112023.faa --pul-meta resources/dbCAN-PUL_Feb-2025.tsv \
  --pfam-hmm /path/to/Pfam-A.hmm --threads 32

# later runs: --pul-faa is no longer needed
bash pulsehunter.sh --genome genome.fna --outdir results/my_genome \
  --ref-dir refs/pulsehunter_ref \
  --dbcan-hmm /path/to/dbCAN-HMMdb-V14.hmm --pfam-hmm /path/to/Pfam-A.hmm
```

## A folder of genomes

```bash
bash pulsehunter_batch.sh \
  --genome-dir genomes/ \
  --outdir results/ \
  --pul-faa resources/PUL_12112023.faa --pul-meta resources/dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm /path/to/dbCAN-HMMdb-V14.hmm --pfam-hmm /path/to/Pfam-A.hmm \
  --run-dbcan --dbcan-db /path/to/dbcan/db --dbcan-env /path/to/envs/run_dbcan \
  --community-specificity \
  --jobs 8 --threads 4
```

- Every `*.fna`, `*.fa` or `*.fasta` file in the directory is one genome;
  the file name (without extension) becomes the sample name.
- `--genome-dir` can be given several times.
- The reference is built once in `results/reference/` and shared.
- `--jobs` genomes run in parallel, each with `--threads` threads, so the
  run uses about `jobs x threads` cores.
- `--community-specificity` checks every primer against all other genomes
  of the batch. Use it when the genomes form a defined community and the
  primers must distinguish its members.
- All other options are passed on to `pulsehunter.sh`.

## Options

### `pulsehunter.sh`

| Option | Default | Meaning |
|---|---|---|
| `--genome FILE` | required | Target genome FASTA |
| `--outdir DIR` | required | Output directory |
| `--pul-faa FILE` | required unless `--ref-dir` is built | dbCAN-PUL protein FASTA |
| `--pul-meta FILE` | none | dbCAN-PUL metadata TSV (substrates, organisms) |
| `--dbcan-hmm FILE` | none | dbCAN family HMM database |
| `--pfam-hmm FILE` | none | `Pfam-A.hmm` |
| `--ref-dir DIR` | `OUTDIR/reference` | Shared reference directory |
| `--sample NAME` | genome file name | Sample name in tables and figures |
| `--threads INT` | 8 | CPU threads |
| `--run-dbcan` | off | Also run dbCAN3 and compare |
| `--dbcan-db DIR` | none | run_dbcan database directory |
| `--dbcan-env DIR` | `run_dbcan` from `PATH` | Conda environment prefix containing `bin/run_dbcan` |
| `--primer-mode MODE` | `qpcr` | `qpcr`, `pcr` or `none` |
| `--n-primer-sets INT` | 3 | Primer sets reported per locus |
| `--offtarget-fasta FILE` | none | Genome the primers must not amplify (repeatable) |
| `--no-plots` | off | Skip figures and the HTML report |
| `--prodigal-mode MODE` | `single` | `single` for isolate genomes, `meta` for MAGs or short contigs |
| `--min-seed-size INT` | 5 | Minimum sequences per marker HMM |
| `--custom-ievalue FLOAT` | 1e-5 | Marker HMM i-Evalue cutoff |
| `--max-gap-genes INT` | 4 | Maximum genes between two chained signal genes |
| `--prepare-ref-only` | off | Build the reference and exit |
| `--force` | off | Recompute steps whose outputs exist |

### `pulsehunter_batch.sh`

| Option | Default | Meaning |
|---|---|---|
| `--genome-dir DIR` | required | Directory with genome FASTA files (repeatable) |
| `--outdir DIR` | required | One sub-directory per genome plus `summary/` |
| `--jobs INT` | 2 | Genomes processed in parallel |
| `--threads INT` | 8 | Threads per genome |
| `--community-specificity` | off | Check primers against all other genomes of the batch |

### Primer modes

| Mode | Amplicon | Target | Use |
|---|---|---|---|
| `qpcr` | 80 to 200 bp | inside susC, then susD, then CAZymes | qPCR / RT-qPCR of PUL expression or abundance |
| `pcr` | 400 to 1000 bp | across the susC/susD junction | end-point PCR, presence/absence, sequencing |
| `none` | - | - | skip primer design |

## Outputs

Per genome (`OUTDIR/`):

| Path | Content |
|---|---|
| `<sample>_PULseHunter_report.html` | Self-contained report: headline numbers, genome overview, filterable locus table, and gene map, primers and gene table per locus |
| `final/PUL_summary.tsv` | One row per candidate locus |
| `final/gene_annotation_master.tsv` | Every gene with all evidence columns |
| `final/PUL_regions.bed` | Locus coordinates |
| `final/candidate_PULs/<locus>/` | `region.fna`, `proteins.faa`, `cds.fna`, `gff3`, `genes.tsv`, `primers.tsv` |
| `final/figures/` | Gene map per locus (PNG and SVG) and `PUL_genome_overview` |
| `final/primers/PUL_primers.tsv`, `.fasta` | Primer sets with in-silico PCR results |
| `final/dbcan/` | `dbcan_cgc_regions.tsv`, `locus_concordance.tsv`, `cazyme_concordance.tsv`, `concordance_summary.tsv` |
| `dbcan3/run/` | Raw `run_dbcan` output (`overview.txt`, `cgc_standard.out`) |
| `target/` | Prodigal genes: `proteins.faa`, `cds.fna`, `genes.gff` |
| `search/` | Raw HMMER and BLAST tables |
| `logs/` | One log per step |
| `RUN_FINISHED.txt` | Written when all steps completed |

Batch (`OUTDIR/summary/`):

| File | Content |
|---|---|
| `index.html` | Overview with links to every genome report |
| `genome_stats.tsv` | One row per genome: loci, confidence, pairs, CAZymes, primers, dbCAN3 agreement |
| `all_PUL_summary.tsv` | All loci of all genomes |
| `all_PUL_primers.tsv` | All primer sets |
| `all_locus_concordance.tsv` | PULseHunter loci and dbCAN3 CGCs, with category |
| `substrate_matrix.tsv` | Loci per predicted substrate and genome |
| `loci_per_genome`, `pulsehunter_vs_dbcan3`, `substrate_heatmap` | Summary figures (PNG and SVG) |

## Reading the results

### `PUL_summary.tsv`

| Column | Meaning |
|---|---|
| `candidate_id` | Locus name, numbered along the genome |
| `contig`, `start`, `end`, `length_bp` | Coordinates |
| `confidence` | `high`, `medium` or `low` (rules in [WORKFLOW.md](WORKFLOW.md)) |
| `pul_type` | `canonical_SusCD_PUL` (pair + CAZyme), `SusCD_only` (pair, no CAZyme), `incomplete_SusCD_PUL` (SusD or TonB-dependent transporter + CAZymes, no pair), `CAZyme_cluster` (no SusC/SusD) |
| `susCD_pairs`, `cazyme_count`, `tf_count`, `stp_count`, `tc_count` | Gene content |
| `cazyme_families` | CAZy families in the locus |
| `best_reference_pul` | Best-matching dbCAN-PUL entry |
| `ref_genes_matched`, `ref_pul_size`, `ref_pul_coverage`, `ref_mean_pident` | How much of that reference PUL is recovered and at what identity |
| `predicted_substrate`, `substrate_support` | Substrate of the best reference PUL and how well it is supported (`strong`, `moderate`, `weak`) |
| `at_contig_edge` | The locus touches a contig end and may be truncated |
| `gene_order` | Compact string of the locus, e.g. `GH105 > susC > susD > PL1_2` |
| `dbcan3_cgcs`, `dbcan3_overlap_frac`, `call_support` | Overlapping dbCAN3 CGCs; `call_support` is `both` or `PULseHunter_only` |

### `PUL_primers.tsv`

| Column | Meaning |
|---|---|
| `rank` | 1 is the best set of the locus |
| `target_gene`, `target_class`, `target_label` | Gene the amplicon lies in |
| `fwd_primer`, `rev_primer`, `*_tm`, `*_gc` | Sequences (5' to 3'), melting temperature, GC % |
| `product_size`, `amplicon_start`, `amplicon_end` | Amplicon |
| `primer3_penalty` | Lower is better |
| `on_target_amplicons` | Predicted products in the target genome (1 is wanted) |
| `off_target_amplicons`, `off_target_genomes` | Predicted products in the off-target genomes |
| `specific` | `True` when exactly one product in the target genome and none elsewhere |

### Suggested order

1. Open the HTML report and filter the locus table to `high` confidence.
2. Check `substrate_support`; follow up `strong` and `moderate` loci first.
3. Loci with `call_support = both` are found by both methods.
   `PULseHunter_only` loci are mostly susC/susD pairs with few or no CAZymes;
   dbCAN3-only CGCs (in `locus_concordance.tsv`) are mostly CAZyme clusters
   without a susC/susD pair.
4. Take primers with `specific = True`, rank 1 first.

## Re-running and tuning

- Finished steps are skipped when the same command is run again, so an
  interrupted run can simply be restarted. `--force` recomputes everything.
- To change the calling rules without repeating the searches, run
  `scripts/call_puls_hmm.py` directly on the existing `target/` and
  `search/` files (`--help` lists all thresholds: `--max-gap-genes`,
  `--max-nonsignal-gap`, `--max-gap-bp`, `--pair-max-gap`,
  `--dbcan-ievalue`, `--dbcan-cov`, `--blast-min-pident`,
  `--blast-min-qcov`, ...), then re-run `design_primers.py`,
  `plot_puls.py` and `make_report.py`.
- To redesign primers only:

  ```bash
  python scripts/design_primers.py --final-dir results/my_genome/final \
    --genome genome.fna --mode pcr --n-sets 3 \
    --offtarget-fasta other_genome.fna
  ```

## Troubleshooting

| Problem | What to check |
|---|---|
| A step fails | The message names the log in `OUTDIR/logs/`; batch logs are in `OUTDIR/logs/<sample>.log` |
| `run_dbcan` not found | Pass `--dbcan-env` with the environment prefix (the directory that contains `bin/run_dbcan`) |
| Zero SusC/SusD calls | `--pfam-hmm` missing or not a Pfam-A HMM file |
| Few CAZymes | `--dbcan-hmm` missing |
| Many loci with `at_contig_edge = True` | Fragmented assembly; loci may be split across contigs |
| Fragmented or metagenome-assembled genome | Use `--prodigal-mode meta` |
| No primers for a locus | No pair passed the specificity test; try `--primer-mode pcr` or relax `--max-mismatch` in `design_primers.py` |
| The dbCAN3 arm is slow | Expected with run_dbcan 3; run several genomes in parallel with `--jobs` |
