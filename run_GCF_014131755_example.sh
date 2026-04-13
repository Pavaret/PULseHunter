#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

bash "$BASE_DIR/run_pul_hmm_pipeline.sh" \
  --genome GCF_014131755.1_ASM1413175v1_genomic.fna \
  --pul-faa PUL_12112023.faa \
  --pul-meta dbCAN-PUL_Feb-2025.tsv \
  --dbcan-hmm dbCAN.hmm \
  --pfam-hmm Pfam-A.hmm \
  --threads 24 \
  --min-seed-size 5 \
  --outdir GCF_014131755.1_PUL_HMM
