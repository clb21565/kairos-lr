#!/bin/bash
# Tier 3 short-read aligner benchmark (docs/testing.md). Runs everything in this one job;
# finished steps are skipped if it is resubmitted.
#
# Usage (from the repo root):
#   CATALOGUE=/projects/MA/HGT/catalogue.fa OUTDIR=/projects/MA/HGT/bench sbatch slurm/aligner_benchmark.sbatch
#   ... sbatch slurm/aligner_benchmark.sbatch --aligners strobealign,bwa-mem2   # extra args pass through
#SBATCH --job-name=kairos_bench
#SBATCH --output=kairos_bench_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

set -euo pipefail
REPO=${SLURM_SUBMIT_DIR:-$(pwd)}
: "${CATALOGUE:?set CATALOGUE to the catalogue FASTA}"
OUTDIR=${OUTDIR:-$REPO/bench}

eval "$(conda shell.bash hook)"
conda activate kairos-bench        # conda env create -f tools/aligner_benchmark.yaml

python "$REPO/tools/aligner_benchmark.py" --outdir "$OUTDIR" --catalogue "$CATALOGUE" \
    --families-per-state 40 --background-mb 500 --threads "${SLURM_CPUS_PER_TASK:-32}" "$@"
