# PULseHunter

Find, annotate, visualise and design primers for polysaccharide utilization
loci (PULs) in bacterial genomes.

PULseHunter runs two complementary arms on the same gene calls and reports
where they agree:

| Arm | What it does |
|---|---|
| **PULseHunter** | Marker HMMs built from dbCAN-PUL, BLASTP against dbCAN-PUL proteins, dbCAN family HMMs and Pfam SusC/SusD domains, combined with genomic context (susC/susD pairs) to call loci and transfer a substrate from the best-matching characterised PUL. |
| **dbCAN3** (optional) | `run_dbcan` (HMMER + DIAMOND) with CGCFinder, giving CAZyme gene clusters (CGCs). |

For every candidate locus you get a gene map, a gene table, sequences and
the three best primer sets, plus one HTML report per genome and a
cross-genome summary for batches.

## Documentation

| Document | Content |
|---|---|
| [docs/INSTALL.md](docs/INSTALL.md) | Environment, databases, checking the installation |
| [docs/USAGE.md](docs/USAGE.md) | Running one genome or a folder, all options, outputs, troubleshooting |
| [docs/WORKFLOW.md](docs/WORKFLOW.md) | Core idea, workflow diagram, the steps and the calling rules |

## Quick start

```bash
git clone https://github.com/Pavaret/PULseHunter.git
cd PULseHunter
micromamba create -f envs/pulsehunter.yml     # or mamba / conda env create
micromamba activate pulsehunter
tar -xJf files/resources.tar.xz               # dbCAN-PUL proteins + metadata

bash pulsehunter.sh \
  --genome genome.fna \
  --outdir results/my_genome \
  --pul-faa resources/PUL_12112023.faa \
  --pul-meta resources/dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm /path/to/dbCAN-HMMdb-V14.hmm \
  --pfam-hmm /path/to/Pfam-A.hmm \
  --threads 16
```

Open `results/my_genome/my_genome_PULseHunter_report.html`.

Add `--run-dbcan --dbcan-db /path/to/dbcan/db --dbcan-env /path/to/envs/run_dbcan`
to run the dbCAN3 arm as well. For many genomes use `pulsehunter_batch.sh`
(see [docs/USAGE.md](docs/USAGE.md)).

## Repository layout

```
pulsehunter.sh            one genome, all steps
pulsehunter_batch.sh      a folder of genomes + cross-genome summary
envs/pulsehunter.yml      conda environment
files/resources.tar.xz    dbCAN-PUL proteins, metadata, two example genomes
scripts/
  prepare_reference_db.py             parse dbCAN-PUL headers
  build_marker_hmms_from_headers.py   group markers, write the HMM build script
  call_puls_hmm.py                    evidence integration and locus calling
  integrate_dbcan.py                  dbCAN3 concordance
  design_primers.py                   primer design + in-silico PCR
  plot_puls.py                        gene maps and genome overview
  make_report.py                      per-genome HTML report
  summarize_batch.py                  cross-genome tables, figures, index.html
docs/                     installation, usage, workflow
```

## Notes

- Substrates are transferred by homology to characterised PULs; treat `weak`
  support as a hint only.
- Primer specificity is an in-silico prediction; validate primers
  experimentally.

## License

See [LICENSE](LICENSE).
