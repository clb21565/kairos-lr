#!/bin/bash
# Driver job: runs Snakemake with the SLURM profile. Snakemake submits each rule as its own
# job with the resources in profiles/slurm/config.yaml; this job only waits and resubmits,
# so it needs little CPU or memory, but must outlive the whole run.
#
# Usage (from the repo root):
#   sbatch slurm/kairos_driver.sbatch                      # everything the inputs allow
#   sbatch slurm/kairos_driver.sbatch detect               # one stage: detect | assess | call
#   WORKDIR=/projects/MA/HGT/kairos_run1 sbatch slurm/kairos_driver.sbatch
#
# WORKDIR holds config.yaml (copied from config/config.yaml and edited), the sample sheet
# and results/. Default: the current directory.
#SBATCH --job-name=kairos_driver
#SBATCH --output=kairos_driver_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00

set -euo pipefail
REPO=${SLURM_SUBMIT_DIR:-$(pwd)}
WORKDIR=${WORKDIR:-$REPO}

# activate the launcher env (conda env create -f environment.yml)
eval "$(conda shell.bash hook)"
conda activate kairos-lr

CONFIG="$WORKDIR/config.yaml"
[[ -f "$CONFIG" ]] || CONFIG="$REPO/config/config.yaml"

snakemake -s "$REPO/Snakefile" -d "$WORKDIR" --configfile "$CONFIG" \
    --profile "$REPO/profiles/slurm" "$@"
