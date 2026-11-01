#!/bin/bash
# Alternative to the driver: run the whole workflow inside one big job (no per-rule
# submission). Simpler, but every rule waits for the largest allocation; the memory must
# cover mmseqs taxonomy against GTDB.
#
# Usage (from the repo root):
#   sbatch slurm/kairos_single_node.sbatch [target ...]
#   WORKDIR=/projects/MA/HGT/kairos_run1 sbatch slurm/kairos_single_node.sbatch
#SBATCH --job-name=kairos
#SBATCH --output=kairos_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=64
#SBATCH --mem=650G
#SBATCH --time=5-00:00:00

set -euo pipefail
REPO=${SLURM_SUBMIT_DIR:-$(pwd)}
WORKDIR=${WORKDIR:-$REPO}

eval "$(conda shell.bash hook)"
conda activate kairos-lr

CONFIG="$WORKDIR/config.yaml"
[[ -f "$CONFIG" ]] || CONFIG="$REPO/config/config.yaml"

snakemake -s "$REPO/Snakefile" -d "$WORKDIR" --configfile "$CONFIG" \
    --use-conda --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=$(( SLURM_MEM_PER_NODE )) \
    --config threads="$SLURM_CPUS_PER_TASK" \
    --rerun-incomplete --printshellcmds "$@"
