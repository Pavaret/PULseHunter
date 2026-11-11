# Installation

PULseHunter is a set of Bash and Python scripts. There is nothing to compile:
clone the repository, create the conda environment and point the pipeline at
the databases.

## 1. Get the code

```bash
git clone https://github.com/Pavaret/PULseHunter.git
cd PULseHunter
```

## 2. Create the environment

```bash
micromamba create -f envs/pulsehunter.yml
micromamba activate pulsehunter
```

`mamba env create -f envs/pulsehunter.yml` or
`conda env create -f envs/pulsehunter.yml` work the same way.

The environment contains:

| Tool | Used for |
|---|---|
| python 3.11, pandas, biopython | all scripts |
| prodigal | gene prediction |
| hmmer | marker HMMs, dbCAN family HMMs, Pfam SusC/SusD |
| mafft | alignments for the marker HMMs |
| blast | BLASTP against dbCAN-PUL proteins |
| primer3-py | primer design |
| matplotlib | figures |

## 3. Databases

| Database | Option | Where to get it | Needed |
|---|---|---|---|
| dbCAN-PUL proteins (`PUL_12112023.faa`) | `--pul-faa` | `files/resources.tar.xz` in this repository | yes |
| dbCAN-PUL metadata (`dbCAN-PUL_Feb-2025.tsv`) | `--pul-meta` | `files/resources.tar.xz` | recommended (substrates) |
| dbCAN family HMMs (`dbCAN-HMMdb-V14.hmm`) | `--dbcan-hmm` | dbCAN download site (`dbCAN-HMMdb-V14.txt`, rename to `.hmm`) | recommended (CAZyme calling) |
| Pfam (`Pfam-A.hmm`) | `--pfam-hmm` | EBI Pfam FTP, `current_release/Pfam-A.hmm.gz` | recommended (SusC/SusD domains) |
| run_dbcan database directory | `--dbcan-db` | see the run_dbcan documentation | only for `--run-dbcan` |

Unpack the bundled resources:

```bash
tar -xJf files/resources.tar.xz
ls resources/
# PUL_12112023.faa  dbCAN-PUL_Feb-2025.tsv  dbCAN-PUL_Feb-2025.xlsx
# GCF_014131755.1_ASM1413175v1_genomic.fna   (Bacteroides thetaiotaomicron, example)
# GCF_900537995.1_Roseburia_intestinalis_strain_L1-82_genomic.fna   (example)
```

The HMM databases do not need `hmmpress`; PULseHunter uses `hmmsearch`.
Without `--dbcan-hmm` or `--pfam-hmm` the pipeline still runs, but with less
evidence per gene, so both are strongly recommended.

## 4. dbCAN3 arm (optional)

The dbCAN3 arm calls `run_dbcan` (tested with version 3.0.7). Its
dependencies conflict with the main environment, so install it in its own
environment and pass the environment prefix with `--dbcan-env`:

```bash
micromamba create -n run_dbcan -c conda-forge -c bioconda dbcan=3.0.7
```

Then download and format the run_dbcan database as described in the
run_dbcan documentation, and pass that directory with `--dbcan-db`.

If `run_dbcan` is already on your `PATH`, `--dbcan-env` can be omitted.

## 5. Check the installation

```bash
micromamba activate pulsehunter
bash pulsehunter.sh --help

bash pulsehunter.sh \
  --genome resources/GCF_014131755.1_ASM1413175v1_genomic.fna \
  --outdir test_run/Btheta \
  --pul-faa resources/PUL_12112023.faa \
  --pul-meta resources/dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm /path/to/dbCAN-HMMdb-V14.hmm \
  --pfam-hmm /path/to/Pfam-A.hmm \
  --threads 16
```

The run is finished when `test_run/Btheta/RUN_FINISHED.txt` exists. For this
*Bacteroides thetaiotaomicron* genome expect roughly 90 candidate loci, most
of them high confidence.

## Resources

- The first run builds the reference (about 130 marker HMMs and a BLAST
  database). This takes about 20 minutes on 32 cores and is done only once
  per reference directory.
- One Bacteroides genome takes a few minutes for the PULseHunter arm with 8
  to 16 threads once the reference exists. The Pfam search is the slowest
  step.
- The dbCAN3 arm adds about 20 to 60 minutes per genome, because
  `run_dbcan` 3 scans the HMM database mostly on one core.
- Memory use stays below 8 GB for a typical bacterial genome.
