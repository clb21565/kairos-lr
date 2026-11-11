#!/bin/bash
# Build a dereplicated catalogue with tools/build_catalogue.sh as a SLURM job, so it keeps
# running after you log out.
#
# Usage (from the repo root):
#   LIST=site1.txt OUT=/projects/MA/HGT/cat_site1 sbatch slurm/build_catalogue.sbatch
#   MIN_LEN=1500 MIN_COV=0.99 LIST=... OUT=... sbatch slurm/build_catalogue.sbatch
#   sbatch --mem=800G ...          # more memory; mmseqs is capped at 70% of it (SPLIT_MEM)
#
# Resubmitting with the same OUT resumes: finished filtered assemblies (OUT/filtered) and
# finished mmseqs steps are kept. Disk: OUT needs ~1x the uncompressed input for the mmseqs DB
# while it runs, plus the gzipped filtered assemblies.
# Check on it with: squeue --me ; tail -f kairos_catalogue_<jobid>.out
#SBATCH --job-name=kairos_catalogue
#SBATCH --output=kairos_catalogue_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=3-00:00:00

set -euo pipefail
REPO=${SLURM_SUBMIT_DIR:-$(pwd)}
: "${LIST:?set LIST to a file with one assembly path per line}"
: "${OUT:?set OUT to the output directory}"

eval "$(conda shell.bash hook)"
conda activate kairos-catalogue     # conda env create -f tools/build_catalogue.yaml

bash "$REPO/tools/build_catalogue.sh" "$LIST" "$OUT" "${MIN_LEN:-1000}" "${SLURM_CPUS_PER_TASK:-32}"
