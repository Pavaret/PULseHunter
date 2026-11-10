#!/usr/bin/env bash
# PULseHunter batch runner: every genome FASTA in one or more directories.
# The reference (marker HMMs, BLAST db) is built once and shared. With
# --community-specificity each genome's primers are also checked against all
# other genomes of the batch (useful for defined / synthetic communities).
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  pulsehunter_batch.sh --genome-dir DIR [--genome-dir DIR2 ...] --outdir DIR \
      --pul-faa PUL.faa --pul-meta dbCAN-PUL.tsv --dbcan-hmm dbCAN.hmm --pfam-hmm Pfam-A.hmm \
      [--run-dbcan --dbcan-db DIR --dbcan-env DIR] [--community-specificity] \
      --jobs 4 --threads 8

  --genome-dir DIR          Directory with *.fna / *.fa / *.fasta genomes (repeatable)
  --outdir DIR              One sub-directory per genome plus summary/
  --jobs INT                Genomes processed in parallel (default: 2)
  --threads INT             Threads per genome (default: 8)
  --community-specificity   Check primers against all other genomes of the batch
  All other options are passed on to pulsehunter.sh (see pulsehunter.sh --help).
USAGE
}

GENOME_DIRS=(); OUTDIR=""; JOBS=2; THREADS=8; COMMUNITY=0; PASS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome-dir) GENOME_DIRS+=("$2"); shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --community-specificity) COMMUNITY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    --run-dbcan|--no-plots|--force) PASS+=("$1"); shift ;;
    --*) PASS+=("$1" "$2"); shift 2 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done
[[ ${#GENOME_DIRS[@]} -gt 0 && -n "$OUTDIR" ]] || { usage; exit 1; }

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTDIR/logs"
OUTDIR="$(cd "$OUTDIR" && pwd)"

GENOMES=()
for d in "${GENOME_DIRS[@]}"; do
  n=0
  while IFS= read -r f; do GENOMES+=("$f"); n=$((n + 1)); done \
    < <(find -L "$d" -maxdepth 1 -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \) | sort)
  [[ $n -gt 0 ]] || echo "WARNING: no genome FASTA files (*.fna, *.fa, *.fasta) in $d" >&2
done
[[ ${#GENOMES[@]} -gt 0 ]] || { echo "ERROR: no genomes found" >&2; exit 1; }
dups=$(for g in "${GENOMES[@]}"; do b=$(basename "$g"); echo "${b%.*}"; done | sort | uniq -d)
[[ -z "$dups" ]] || { echo "ERROR: duplicate genome names across directories: $dups" >&2; exit 1; }
echo "[$(date +%H:%M:%S)] ${#GENOMES[@]} genomes -> $OUTDIR" >&2

# Build the shared reference once, before any parallel work
bash "$BASE_DIR/pulsehunter.sh" --prepare-ref-only --ref-dir "$OUTDIR/reference" \
  --threads $((JOBS * THREADS)) "${PASS[@]+"${PASS[@]}"}"

run_one() {
  local genome="$1" sample off=()
  sample="$(basename "$genome")"; sample="${sample%.*}"
  if [[ "$COMMUNITY" -eq 1 ]]; then
    for g in "${GENOMES[@]}"; do [[ "$g" == "$genome" ]] || off+=(--offtarget-fasta "$g"); done
  fi
  if bash "$BASE_DIR/pulsehunter.sh" --genome "$genome" --sample "$sample" \
      --outdir "$OUTDIR/$sample" --ref-dir "$OUTDIR/reference" --threads "$THREADS" \
      "${PASS[@]+"${PASS[@]}"}" "${off[@]+"${off[@]}"}" > "$OUTDIR/logs/$sample.log" 2>&1; then
    echo "[$(date +%H:%M:%S)] finished $sample" >&2
  else
    echo "[$(date +%H:%M:%S)] FAILED $sample (see $OUTDIR/logs/$sample.log)" >&2
    return 1
  fi
}

fail=0
for genome in "${GENOMES[@]}"; do
  while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || fail=$((fail + 1)); done
  run_one "$genome" &
done
while (( $(jobs -rp | wc -l) > 0 )); do wait -n || fail=$((fail + 1)); done

python3 "$BASE_DIR/scripts/summarize_batch.py" --batch-dir "$OUTDIR" > "$OUTDIR/logs/summarize_batch.log" 2>&1 \
  || { echo "ERROR: batch summary failed; see $OUTDIR/logs/summarize_batch.log" >&2; exit 1; }
echo "[$(date +%H:%M:%S)] batch summary: $OUTDIR/summary/index.html" >&2
[[ "$fail" -eq 0 ]] || { echo "ERROR: $fail genome(s) failed" >&2; exit 1; }
