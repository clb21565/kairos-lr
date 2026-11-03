#!/bin/bash
# Long-read SV calls (Sniffles2) vs kairos features (tools/sv_compare.py).
# Run from the kairos run directory after assess: uses results/assess/<sample>/<sample>.bam for
# every ont/hifi sample in the sample sheet. Finished per-sample Sniffles2 runs are skipped if
# the job is resubmitted.
#
# Usage (from the run directory):
#   REPO=/path/to/kairos-lr SAMPLES=samples.tsv sbatch $REPO/slurm/sv_compare.sbatch
#   ... optional: CATALOGUE_NAME=catalogue (results/detect/<name>.*), OUTDIR=results/sv,
#       extra args pass through to sv_compare.py
#SBATCH --job-name=kairos_sv
#SBATCH --output=kairos_sv_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

set -euo pipefail
: "${REPO:?set REPO to the kairos-lr checkout}"
: "${SAMPLES:?set SAMPLES to the sample sheet}"
CAT=${CATALOGUE_NAME:-catalogue}
OUTDIR=${OUTDIR:-results/sv}
T=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$OUTDIR/sniffles"

eval "$(conda shell.bash hook)"
conda activate kairos-sv           # conda env create -f $REPO/tools/sv_compare.yaml

# ont/hifi samples from the sheet (columns by name)
mapfile -t LONG < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
    /^#/ {next} $c["read_type"]=="ont" || $c["read_type"]=="hifi" {print $c["sample"]}' "$SAMPLES")
[ ${#LONG[@]} -gt 0 ] || { echo "no ont/hifi samples in $SAMPLES"; exit 1; }

: > "$OUTDIR/sniffles/snf_list.tsv"
for s in "${LONG[@]}"; do
    snf="$OUTDIR/sniffles/$s.snf"
    if [ ! -s "$snf" ]; then
        # --all-contigs: Sniffles2 skips contigs < 1 Mb by default
        sniffles --input "results/assess/$s/$s.bam" --snf "$snf.tmp" --vcf "$OUTDIR/sniffles/$s.vcf" \
            --mosaic --mosaic-include-germline --all-contigs --minsvlen 300 \
            --sample-id "$s" --threads "$T" --allow-overwrite
        mv "$snf.tmp" "$snf"
    fi
    printf '%s\t%s\n' "$snf" "$s" >> "$OUTDIR/sniffles/snf_list.tsv"
done

sniffles --input "$OUTDIR/sniffles/snf_list.tsv" --vcf "$OUTDIR/sniffles/merged.vcf" \
    --all-contigs --threads "$T" --allow-overwrite

python "$REPO/tools/sv_compare.py" --vcf "$OUTDIR/sniffles/merged.vcf" \
    --features "results/detect/$CAT.features" --catalogue "results/detect/$CAT.fa" \
    --assess-dir results/assess --outdir "$OUTDIR" --threads "$T" "$@"
