#!/usr/bin/env bash
# PULseHunter - PUL discovery for one bacterial genome.
#   arm 1: PULseHunter (dbCAN-PUL marker HMMs + reference homology + susC/susD context)
#   arm 2: dbCAN3 run_dbcan CGCFinder (optional, --run-dbcan)
#   then : concordance of the two arms, locus figures, primer design, HTML report
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  pulsehunter.sh --genome genome.fna --outdir results/sample \
      --pul-faa PUL_12112023.faa --pul-meta dbCAN-PUL_Feb-2025.tsv \
      --dbcan-hmm dbCAN-HMMdb-V14.hmm --pfam-hmm Pfam-A.hmm \
      [--run-dbcan --dbcan-db /path/to/dbcan/db --dbcan-env /path/to/conda/env] \
      [--offtarget-fasta other_genome.fna ...] --threads 16

Required:
  --genome FILE          Target genome (FASTA)
  --outdir DIR           Output directory
  --pul-faa FILE         dbCAN-PUL protein FASTA (not needed if --ref-dir is already built)

Recommended:
  --pul-meta FILE        dbCAN-PUL metadata TSV (substrates, organisms)
  --dbcan-hmm FILE       dbCAN family HMM database (CAZyme calling)
  --pfam-hmm FILE        Pfam-A.hmm (SusC/SusD domain evidence)
  --ref-dir DIR          Shared reference directory. Marker HMMs and BLAST db are
                         built once here and reused by every genome (default: OUTDIR/reference)
  --sample NAME          Sample name (default: genome file name)
  --threads INT          CPU threads (default: 8)

dbCAN3 arm:
  --run-dbcan            Also run dbCAN3 (run_dbcan + CGCFinder) and compare
  --dbcan-db DIR         run_dbcan database directory
  --dbcan-env DIR        Conda env prefix that contains bin/run_dbcan
                         (default: run_dbcan from PATH)

Primers and figures:
  --primer-mode MODE     qpcr (80-200 bp in marker genes, default) | pcr (400-1000 bp
                         across the susC/susD junction) | none
  --n-primer-sets INT    Primer sets per locus (default: 3)
  --offtarget-fasta FILE Genome the primers must NOT amplify (repeatable)
  --no-plots             Skip figures and the HTML report

Tuning:
  --prodigal-mode single|meta   (default: single)
  --min-seed-size INT           Min sequences per marker HMM (default: 5)
  --custom-ievalue FLOAT        Marker HMM i-Evalue cutoff (default: 1e-5)
  --max-gap-genes INT           Max genes between chained signal genes (default: 4)
  --prepare-ref-only            Build the reference directory and exit
  --force                       Recompute steps whose outputs already exist
USAGE
}

GENOME=""; PUL_FAA=""; PUL_META=""; OUTDIR=""; PFAM_HMM=""; DBCAN_HMM=""
REF_DIR=""; SAMPLE=""; THREADS=8; PRODIGAL_MODE="single"; MIN_SEED_SIZE=5
CUSTOM_IEVALUE="1e-5"; MAX_GAP_GENES=4
RUN_DBCAN=0; DBCAN_DB=""; DBCAN_ENV=""
PRIMER_MODE="qpcr"; N_PRIMER_SETS=3; OFFTARGETS=(); PLOTS=1
REF_ONLY=0; FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2 ;;
    --pul-faa) PUL_FAA="$2"; shift 2 ;;
    --pul-meta) PUL_META="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --pfam-hmm) PFAM_HMM="$2"; shift 2 ;;
    --dbcan-hmm) DBCAN_HMM="$2"; shift 2 ;;
    --ref-dir) REF_DIR="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --prodigal-mode) PRODIGAL_MODE="$2"; shift 2 ;;
    --min-seed-size) MIN_SEED_SIZE="$2"; shift 2 ;;
    --custom-ievalue) CUSTOM_IEVALUE="$2"; shift 2 ;;
    --max-gap-genes) MAX_GAP_GENES="$2"; shift 2 ;;
    --run-dbcan) RUN_DBCAN=1; shift ;;
    --dbcan-db) DBCAN_DB="$2"; shift 2 ;;
    --dbcan-env) DBCAN_ENV="$2"; shift 2 ;;
    --primer-mode) PRIMER_MODE="$2"; shift 2 ;;
    --n-primer-sets) N_PRIMER_SETS="$2"; shift 2 ;;
    --offtarget-fasta) OFFTARGETS+=("$2"); shift 2 ;;
    --no-plots) PLOTS=0; shift ;;
    --prepare-ref-only) REF_ONLY=1; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

die() { echo "ERROR: $*" >&2; exit 1; }
step() { echo "[$(date +%H:%M:%S)] ${SAMPLE:+$SAMPLE: }$*" >&2; }
# done_or_force FILE -> true when the step still has to run
todo() { [[ "$FORCE" -eq 1 || ! -s "$1" ]]; }

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS="$BASE_DIR/scripts"

if [[ "$REF_ONLY" -eq 1 ]]; then
  [[ -n "$REF_DIR" ]] || die "--prepare-ref-only needs --ref-dir"
else
  [[ -n "$GENOME" && -n "$OUTDIR" ]] || { usage; exit 1; }
  [[ -s "$GENOME" ]] || die "genome not found or empty: $GENOME"
  [[ -n "$REF_DIR" ]] || REF_DIR="$OUTDIR/reference"
  [[ -n "$SAMPLE" ]] || { SAMPLE="$(basename "$GENOME")"; SAMPLE="${SAMPLE%.*}"; }
fi
[[ "$PRODIGAL_MODE" =~ ^(single|meta)$ ]] || die "--prodigal-mode must be single or meta"
[[ "$PRIMER_MODE" =~ ^(qpcr|pcr|none)$ ]] || die "--primer-mode must be qpcr, pcr or none"

for cmd in python3 prodigal mafft hmmbuild hmmpress hmmsearch makeblastdb blastp; do
  command -v "$cmd" >/dev/null 2>&1 || die "required command not found: $cmd"
done

# ---------------------------------------------------------------------------
# Reference: marker HMMs, BLAST db, Pfam SusC/SusD profiles (built once)
# ---------------------------------------------------------------------------
mkdir -p "$REF_DIR/logs"
if todo "$REF_DIR/REFERENCE_READY"; then
  [[ -s "$PUL_FAA" ]] || die "--pul-faa is required to build the reference in $REF_DIR"
  step "building reference in $REF_DIR"
  REF_ARGS=(--pul-faa "$PUL_FAA" --outdir "$REF_DIR")
  [[ -n "$PUL_META" ]] && REF_ARGS+=(--pul-meta "$PUL_META")
  python3 "$SCRIPTS/prepare_reference_db.py" "${REF_ARGS[@]}" > "$REF_DIR/logs/prepare_reference_db.log" 2>&1

  python3 "$SCRIPTS/build_marker_hmms_from_headers.py" \
    --faa "$REF_DIR/reference_proteins.faa" \
    --reference-table "$REF_DIR/reference_proteins.tsv" \
    --outdir "$REF_DIR/custom_hmm" \
    --min-seed-size "$MIN_SEED_SIZE" > "$REF_DIR/logs/build_marker_hmms.log" 2>&1
  [[ -s "$REF_DIR/custom_hmm/cluster_manifest.tsv" ]] || die "no marker groups produced; see $REF_DIR/logs/build_marker_hmms.log"
  THREADS="$THREADS" bash "$REF_DIR/custom_hmm/build_hmms.sh" > "$REF_DIR/logs/custom_hmm_build.log" 2>&1
  [[ -s "$REF_DIR/custom_hmm/custom_pul_profiles.hmm" ]] || die "marker HMM database was not built; see $REF_DIR/logs/custom_hmm_build.log"

  makeblastdb -in "$REF_DIR/reference_proteins.faa" -dbtype prot \
    -out "$REF_DIR/reference_proteins" > "$REF_DIR/logs/makeblastdb.log" 2>&1

  if [[ -n "$PFAM_HMM" ]]; then
    command -v hmmfetch >/dev/null 2>&1 || die "hmmfetch required when --pfam-hmm is given"
    # Resolve the versioned accessions of this Pfam release (PF00593 -> PF00593.31)
    grep -E '^ACC +PF(00593|13715|07715|07660|07980|12741|12771|14322)\.' "$PFAM_HMM" \
      | awk '{print $2}' > "$REF_DIR/suscd_pfam_ids.txt"
    [[ "$(wc -l < "$REF_DIR/suscd_pfam_ids.txt")" -eq 8 ]] \
      || echo "WARNING: only $(wc -l < "$REF_DIR/suscd_pfam_ids.txt")/8 SusC/SusD Pfam models found in $PFAM_HMM" >&2
    hmmfetch -f "$PFAM_HMM" "$REF_DIR/suscd_pfam_ids.txt" > "$REF_DIR/suscd.hmm"
  fi
  date > "$REF_DIR/REFERENCE_READY"
fi
[[ "$REF_ONLY" -eq 1 ]] && { step "reference ready: $REF_DIR"; exit 0; }
[[ -s "$REF_DIR/suscd.hmm" || -z "$PFAM_HMM" ]] || die "$REF_DIR was built without --pfam-hmm; rebuild it with --force"

mkdir -p "$OUTDIR"/{logs,target,search,final}

# ---------------------------------------------------------------------------
# 1) Gene prediction
# ---------------------------------------------------------------------------
if todo "$OUTDIR/target/proteins.faa"; then
  step "gene prediction (prodigal -p $PRODIGAL_MODE)"
  cp "$GENOME" "$OUTDIR/target/genome.fna"
  prodigal -i "$OUTDIR/target/genome.fna" -a "$OUTDIR/target/proteins.faa.tmp" \
    -d "$OUTDIR/target/cds.fna" -o "$OUTDIR/target/genes.gff" \
    -f gff -p "$PRODIGAL_MODE" > "$OUTDIR/logs/prodigal.log" 2>&1
  mv "$OUTDIR/target/proteins.faa.tmp" "$OUTDIR/target/proteins.faa"
fi

# ---------------------------------------------------------------------------
# 2) Searches
# ---------------------------------------------------------------------------
if todo "$OUTDIR/search/custom_pul.domtblout"; then
  step "marker HMM search"
  hmmsearch --cpu "$THREADS" --domtblout "$OUTDIR/search/custom_pul.domtblout.tmp" \
    "$REF_DIR/custom_hmm/custom_pul_profiles.hmm" "$OUTDIR/target/proteins.faa" > /dev/null 2> "$OUTDIR/logs/custom_hmmsearch.log"
  mv "$OUTDIR/search/custom_pul.domtblout.tmp" "$OUTDIR/search/custom_pul.domtblout"
fi

if todo "$OUTDIR/search/reference_vs_target.blastp.tsv"; then
  step "BLASTP against dbCAN-PUL proteins"
  blastp -query "$OUTDIR/target/proteins.faa" -db "$REF_DIR/reference_proteins" \
    -out "$OUTDIR/search/reference_vs_target.blastp.tsv.tmp" \
    -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs' \
    -evalue 1e-10 -max_target_seqs 25 -num_threads "$THREADS" > "$OUTDIR/logs/blastp.log" 2>&1
  mv "$OUTDIR/search/reference_vs_target.blastp.tsv.tmp" "$OUTDIR/search/reference_vs_target.blastp.tsv"
fi

if [[ -n "$DBCAN_HMM" ]]; then
  if todo "$OUTDIR/search/dbcan.domtblout"; then
    step "dbCAN family HMM search"
    hmmsearch --cpu "$THREADS" --domtblout "$OUTDIR/search/dbcan.domtblout.tmp" \
      "$DBCAN_HMM" "$OUTDIR/target/proteins.faa" > /dev/null 2> "$OUTDIR/logs/dbcan_hmmsearch.log"
    mv "$OUTDIR/search/dbcan.domtblout.tmp" "$OUTDIR/search/dbcan.domtblout"
  fi
else
  : > "$OUTDIR/search/dbcan.domtblout"
  echo "WARNING: --dbcan-hmm not supplied; CAZyme calling relies on marker HMMs and BLAST only." >&2
fi

if [[ -s "$REF_DIR/suscd.hmm" ]]; then
  if todo "$OUTDIR/search/suscd.domtblout"; then
    step "Pfam SusC/SusD search"
    hmmsearch --cpu "$THREADS" --cut_ga --domtblout "$OUTDIR/search/suscd.domtblout.tmp" \
      "$REF_DIR/suscd.hmm" "$OUTDIR/target/proteins.faa" > /dev/null 2> "$OUTDIR/logs/suscd_hmmsearch.log"
    mv "$OUTDIR/search/suscd.domtblout.tmp" "$OUTDIR/search/suscd.domtblout"
  fi
else
  : > "$OUTDIR/search/suscd.domtblout"
  echo "WARNING: --pfam-hmm not supplied; SusC/SusD calls rely on marker HMMs and BLAST only." >&2
fi

# ---------------------------------------------------------------------------
# 3) Call candidate PULs
# ---------------------------------------------------------------------------
step "calling candidate PULs"
rm -rf "$OUTDIR/final/candidate_PULs" "$OUTDIR/final/figures" "$OUTDIR/final/primers"
python3 "$SCRIPTS/call_puls_hmm.py" \
  --gff "$OUTDIR/target/genes.gff" \
  --proteins "$OUTDIR/target/proteins.faa" \
  --cds "$OUTDIR/target/cds.fna" \
  --genome "$OUTDIR/target/genome.fna" \
  --reference-table "$REF_DIR/reference_proteins.tsv" \
  --cluster-manifest "$REF_DIR/custom_hmm/cluster_manifest.tsv" \
  --custom-domtblout "$OUTDIR/search/custom_pul.domtblout" \
  --suscd-domtblout "$OUTDIR/search/suscd.domtblout" \
  --dbcan-domtblout "$OUTDIR/search/dbcan.domtblout" \
  --blast-tsv "$OUTDIR/search/reference_vs_target.blastp.tsv" \
  --custom-ievalue "$CUSTOM_IEVALUE" \
  --max-gap-genes "$MAX_GAP_GENES" \
  --sample "$SAMPLE" \
  --outdir "$OUTDIR/final" > "$OUTDIR/logs/call_puls.log" 2>&1

# ---------------------------------------------------------------------------
# 4) dbCAN3 arm (run_dbcan + CGCFinder on the same gene calls) and concordance
# ---------------------------------------------------------------------------
if [[ "$RUN_DBCAN" -eq 1 ]]; then
  [[ -d "$DBCAN_DB" ]] || die "--run-dbcan needs --dbcan-db DIR"
  DBCAN_PATH="$PATH"
  [[ -n "$DBCAN_ENV" ]] && DBCAN_PATH="$DBCAN_ENV/bin:$PATH"
  PATH="$DBCAN_PATH" command -v run_dbcan >/dev/null 2>&1 || die "run_dbcan not found (set --dbcan-env)"
  if todo "$OUTDIR/dbcan3/overview.txt"; then
    step "dbCAN3 (run_dbcan: HMMER + DIAMOND + CGCFinder)"
    rm -rf "$OUTDIR/dbcan3/run"; mkdir -p "$OUTDIR/dbcan3"
    # Same proteins and gene ids as the PULseHunter arm, so the arms are comparable
    sed 's/\*$//' "$OUTDIR/target/proteins.faa" > "$OUTDIR/dbcan3/input_proteins.faa"
    awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && $3=="CDS" {n=$9; sub(/;.*/,"",n); sub(/^ID=[0-9]+_/,"",n); $9="ID="$1"_"n; print}' \
      "$OUTDIR/target/genes.gff" > "$OUTDIR/dbcan3/input_genes.gff"
    ( cd "$OUTDIR/dbcan3" && PATH="$DBCAN_PATH" run_dbcan input_proteins.faa protein \
        -c input_genes.gff --tools hmmer diamond \
        --db_dir "${DBCAN_DB%/}/" --out_dir run \
        --dia_cpu "$THREADS" --hmm_cpu "$THREADS" --tf_cpu "$THREADS" --stp_cpu "$THREADS" \
    ) > "$OUTDIR/logs/run_dbcan.log" 2>&1 || die "run_dbcan failed; see $OUTDIR/logs/run_dbcan.log"
    [[ -s "$OUTDIR/dbcan3/run/overview.txt" ]] || die "run_dbcan produced no overview.txt; see $OUTDIR/logs/run_dbcan.log"
    cp "$OUTDIR/dbcan3/run/overview.txt" "$OUTDIR/dbcan3/overview.txt"
  fi
  step "comparing PULseHunter and dbCAN3 calls"
  python3 "$SCRIPTS/integrate_dbcan.py" --final-dir "$OUTDIR/final" \
    --dbcan-dir "$OUTDIR/dbcan3/run" > "$OUTDIR/logs/integrate_dbcan.log" 2>&1
fi

# ---------------------------------------------------------------------------
# 5) Primers, figures, report
# ---------------------------------------------------------------------------
if [[ "$PRIMER_MODE" != "none" ]]; then
  step "designing primers ($PRIMER_MODE, top $N_PRIMER_SETS per locus)"
  PRIMER_ARGS=(--final-dir "$OUTDIR/final" --genome "$OUTDIR/target/genome.fna"
               --mode "$PRIMER_MODE" --n-sets "$N_PRIMER_SETS" --threads "$THREADS")
  for f in "${OFFTARGETS[@]+"${OFFTARGETS[@]}"}"; do PRIMER_ARGS+=(--offtarget-fasta "$f"); done
  python3 "$SCRIPTS/design_primers.py" "${PRIMER_ARGS[@]}" > "$OUTDIR/logs/design_primers.log" 2>&1
fi

if [[ "$PLOTS" -eq 1 ]]; then
  step "drawing locus maps and report"
  python3 "$SCRIPTS/plot_puls.py" --final-dir "$OUTDIR/final" \
    --genome "$OUTDIR/target/genome.fna" > "$OUTDIR/logs/plot_puls.log" 2>&1
  python3 "$SCRIPTS/make_report.py" --final-dir "$OUTDIR/final" \
    --out "$OUTDIR/${SAMPLE}_PULseHunter_report.html" > "$OUTDIR/logs/make_report.log" 2>&1
fi

cat > "$OUTDIR/RUN_FINISHED.txt" <<DONE
PULseHunter finished: $SAMPLE ($(date))
Report:            $OUTDIR/${SAMPLE}_PULseHunter_report.html
Locus table:       $OUTDIR/final/PUL_summary.tsv
Gene table:        $OUTDIR/final/gene_annotation_master.tsv
Per-locus folders: $OUTDIR/final/candidate_PULs/
Figures:           $OUTDIR/final/figures/
Primers:           $OUTDIR/final/primers/PUL_primers.tsv
dbCAN3 comparison: $OUTDIR/final/dbcan/
DONE
step "done -> $OUTDIR/final/PUL_summary.tsv"
