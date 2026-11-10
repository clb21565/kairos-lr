#!/bin/bash
# Simulate the planted reads of every sample in a spike-in build and prepend that sample's
# background reads (SLURM array, one task per row of <build>/samples.tsv).
#
# Usage (from the repo root):
#   n=$(( $(wc -l < BUILD/samples.tsv) - 2 ))
#   BUILD=/path/to/build BACKGROUNDS=/path/to/backgrounds.tsv SIM=iss \
#       sbatch --array=0-$n slurm/spikein_reads.sbatch
# BACKGROUNDS (optional): TSV "sample<TAB>files" (files: R1,R2 for short reads, one file for ONT),
#   one real metagenome per spike-in sample (e.g. one SBR run per timepoint x compartment).
# SIM: iss (default; short_pe) or badread (ont).
#SBATCH --job-name=kairos_spikein
#SBATCH --output=kairos_spikein_%A_%a.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

set -euo pipefail
REPO=${SLURM_SUBMIT_DIR:-$(pwd)}
: "${BUILD:?set BUILD to the spikein.py build directory}"
SIM=${SIM:-iss}

eval "$(conda shell.bash hook)"
conda activate kairos-spikein       # conda env create -f tools/spikein.yaml

sample=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR == i + 2 {print $1}' "$BUILD/samples.tsv")
bg=""
if [ -n "${BACKGROUNDS:-}" ]; then
    bg=$(awk -F'\t' -v s="$sample" '$1 == s {print $2}' "$BACKGROUNDS")
    [ -n "$bg" ] || { echo "no background for $sample in $BACKGROUNDS" >&2; exit 1; }
fi
python "$REPO/tools/spikein.py" reads --build-dir "$BUILD" --sample "$sample" --simulator "$SIM" \
    --background "$bg" --threads "${SLURM_CPUS_PER_TASK:-8}"
