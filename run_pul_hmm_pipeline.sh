#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_pul_hmm_pipeline.sh \
    --genome GCF_014131755.1_ASM1413175v1_genomic.fna \
    --pul-faa PUL_12112023.faa \
    --pul-meta dbCAN-PUL_Feb-2025.tsv \
    --outdir results/GCF_014131755.1_PUL_HMM \
    --threads 24 \
    --pfam-hmm /path/to/Pfam-A.hmm \
    --dbcan-hmm /path/to/dbCAN.hmm

Required:
  --genome     Target bacterial genome in FASTA
  --pul-faa    Reference PUL protein FASTA (your PUL_12112023.faa)
  --outdir     Output directory

Recommended:
  --pul-meta   dbCAN-PUL metadata TSV
  --pfam-hmm   Full Pfam-A.hmm file; optional orthogonal SusC/SusD support
  --dbcan-hmm  dbCAN family HMM database; recommended for CAZyme calling
  --threads    CPU threads (default: 8)

Optional tuning:
  --prodigal-mode single|meta      Prodigal mode (default: single)
  --min-seed-size INT              Minimum sequences per marker group HMM (default: 5)
  --custom-ievalue FLOAT           Custom HMM domain i-Evalue cutoff used downstream (default in caller: 1e-5)
EOF
}

GENOME=""
PUL_FAA=""
PUL_META=""
OUTDIR=""
PFAM_HMM=""
DBCAN_HMM=""
THREADS=8
PRODIGAL_MODE="single"
MIN_SEED_SIZE=5
CUSTOM_IEVALUE="1e-5"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2 ;;
    --pul-faa) PUL_FAA="$2"; shift 2 ;;
    --pul-meta) PUL_META="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --pfam-hmm) PFAM_HMM="$2"; shift 2 ;;
    --dbcan-hmm) DBCAN_HMM="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --prodigal-mode) PRODIGAL_MODE="$2"; shift 2 ;;
    --min-seed-size) MIN_SEED_SIZE="$2"; shift 2 ;;
    --custom-ievalue) CUSTOM_IEVALUE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "$GENOME" || -z "$PUL_FAA" || -z "$OUTDIR" ]] && { usage; exit 1; }
[[ "$PRODIGAL_MODE" =~ ^(single|meta)$ ]] || { echo "ERROR: --prodigal-mode must be single or meta" >&2; exit 1; }

mkdir -p "$OUTDIR"/{logs,reference,target,search,final,tmp}
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS="$BASE_DIR/scripts"

need_cmds=(python3 prodigal mafft hmmbuild hmmpress hmmsearch makeblastdb blastp)
for cmd in "${need_cmds[@]}"; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: required command not found: $cmd" >&2; exit 1; }
done
if [[ -n "$PFAM_HMM" ]]; then
  command -v hmmfetch >/dev/null 2>&1 || { echo "ERROR: hmmfetch required when --pfam-hmm is provided" >&2; exit 1; }
fi

cp "$GENOME" "$OUTDIR/target/genome.fna"

# 1) Gene prediction
prodigal -i "$OUTDIR/target/genome.fna" \
  -a "$OUTDIR/target/proteins.faa" \
  -d "$OUTDIR/target/cds.fna" \
  -o "$OUTDIR/target/genes.gff" \
  -f gff -p "$PRODIGAL_MODE" > "$OUTDIR/logs/prodigal.log" 2>&1

# 2) Prepare reference DB
REF_ARGS=(--pul-faa "$PUL_FAA" --outdir "$OUTDIR/reference")
if [[ -n "$PUL_META" ]]; then
  REF_ARGS+=(--pul-meta "$PUL_META")
fi
python3 "$SCRIPTS/prepare_reference_db.py" "${REF_ARGS[@]}" > "$OUTDIR/logs/prepare_reference_db.log" 2>&1

# 3) Build focused custom marker HMMs directly from annotated reference headers
python3 "$SCRIPTS/build_marker_hmms_from_headers.py" \
  --faa "$OUTDIR/reference/reference_proteins.faa" \
  --reference-table "$OUTDIR/reference/reference_proteins.tsv" \
  --outdir "$OUTDIR/reference/custom_hmm" \
  --min-seed-size "$MIN_SEED_SIZE" > "$OUTDIR/logs/build_marker_hmms.log" 2>&1

if [[ ! -s "$OUTDIR/reference/custom_hmm/cluster_manifest.tsv" ]]; then
  echo "ERROR: no HMM marker groups were produced. See $OUTDIR/logs/build_marker_hmms.log" >&2
  exit 1
fi

bash "$OUTDIR/reference/custom_hmm/build_hmms.sh" > "$OUTDIR/logs/custom_hmm_build.log" 2>&1

if [[ ! -s "$OUTDIR/reference/custom_hmm/custom_pul_profiles.hmm" ]]; then
  echo "ERROR: custom HMM database was not built successfully." >&2
  exit 1
fi

# 4) Search target proteins against custom HMMs and exact protein reference DB
hmmsearch --cpu "$THREADS" --domtblout "$OUTDIR/search/custom_pul.domtblout" \
  "$OUTDIR/reference/custom_hmm/custom_pul_profiles.hmm" "$OUTDIR/target/proteins.faa" > "$OUTDIR/logs/custom_hmmsearch.log" 2>&1

makeblastdb -in "$OUTDIR/reference/reference_proteins.faa" -dbtype prot -out "$OUTDIR/reference/reference_proteins" > "$OUTDIR/logs/makeblastdb.log" 2>&1
blastp -query "$OUTDIR/target/proteins.faa" -db "$OUTDIR/reference/reference_proteins" \
  -out "$OUTDIR/search/reference_vs_target.blastp.tsv" \
  -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs' \
  -evalue 1e-10 -max_target_seqs 10 -num_threads "$THREADS" > "$OUTDIR/logs/blastp.log" 2>&1

# 5) Optional CAZyme and orthogonal Pfam SusC/SusD annotation
if [[ -n "$DBCAN_HMM" ]]; then
  hmmsearch --cpu "$THREADS" --domtblout "$OUTDIR/search/dbcan.domtblout" \
    "$DBCAN_HMM" "$OUTDIR/target/proteins.faa" > "$OUTDIR/logs/dbcan_hmmsearch.log" 2>&1
else
  : > "$OUTDIR/search/dbcan.domtblout"
  echo "WARNING: --dbcan-hmm not supplied. CAZyme calling will rely more heavily on custom HMMs and BLAST support only." >&2
fi

if [[ -n "$PFAM_HMM" ]]; then
  cat > "$OUTDIR/reference/suscd_pfam_ids.txt" <<'EOF'
PF00593.31
PF13715.13
PF07715.21
PF07660.20
PF07980.18
PF12741.14
PF12771.14
PF14322.13
EOF
  hmmfetch -f "$PFAM_HMM" "$OUTDIR/reference/suscd_pfam_ids.txt" > "$OUTDIR/reference/suscd.hmm"
  hmmpress "$OUTDIR/reference/suscd.hmm" > "$OUTDIR/logs/suscd_hmmpress.log" 2>&1
  hmmsearch --cpu "$THREADS" --cut_ga --domtblout "$OUTDIR/search/suscd.domtblout" \
    "$OUTDIR/reference/suscd.hmm" "$OUTDIR/target/proteins.faa" > "$OUTDIR/logs/suscd_hmmsearch.log" 2>&1
else
  : > "$OUTDIR/search/suscd.domtblout"
  echo "WARNING: --pfam-hmm not supplied. SusC/SusD calls will come from custom HMMs and reference homology only." >&2
fi

# 6) Integrate evidence and call candidate PULs
python3 "$SCRIPTS/call_puls_hmm.py" \
  --gff "$OUTDIR/target/genes.gff" \
  --proteins "$OUTDIR/target/proteins.faa" \
  --cds "$OUTDIR/target/cds.fna" \
  --genome "$OUTDIR/target/genome.fna" \
  --reference-table "$OUTDIR/reference/reference_proteins.tsv" \
  --cluster-manifest "$OUTDIR/reference/custom_hmm/cluster_manifest.tsv" \
  --custom-domtblout "$OUTDIR/search/custom_pul.domtblout" \
  --suscd-domtblout "$OUTDIR/search/suscd.domtblout" \
  --dbcan-domtblout "$OUTDIR/search/dbcan.domtblout" \
  --blast-tsv "$OUTDIR/search/reference_vs_target.blastp.tsv" \
  --custom-ievalue "$CUSTOM_IEVALUE" \
  --outdir "$OUTDIR/final" > "$OUTDIR/logs/call_puls.log" 2>&1

cat > "$OUTDIR/RUN_FINISHED.txt" <<EOF
PUL HMM pipeline finished successfully.
Main result table: $OUTDIR/final/PUL_summary.tsv
Per-locus folders: $OUTDIR/final/candidate_PULs/
Gene master table: $OUTDIR/final/gene_annotation_master.tsv
EOF

echo "Done. See $OUTDIR/final/PUL_summary.tsv"
