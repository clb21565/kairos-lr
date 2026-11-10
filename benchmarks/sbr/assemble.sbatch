#!/bin/bash
# Per-sample metaSPAdes assemblies for the SBR benchmark (SLURM array, one task per run),
# as in the original study but short reads only (its Nanopore reads are not in SRA).
#
# Usage (from benchmarks/sbr, after bash fetch_reads.sh):
#   sbatch --array=0-5 assemble.sbatch
#   # SRR26192775 (117 Gb) may need more memory; rerun just that task with e.g.
#   sbatch --array=3 --mem=900G assemble.sbatch
# Then: ls $PWD/assemblies/*.fasta > assemblies.txt
#       cd ../.. && LIST=benchmarks/sbr/assemblies.txt OUT=$PWD/benchmarks/sbr/catalogue \
#           sbatch slurm/build_catalogue.sbatch
#SBATCH --job-name=sbr_assemble
#SBATCH --output=sbr_assemble_%A_%a.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=4-00:00:00

set -euo pipefail
RUNS=(SRR26192772 SRR26192773 SRR26192774 SRR26192775 SRR26192776 SRR26192777)
run=${RUNS[$SLURM_ARRAY_TASK_ID]}
T=${SLURM_CPUS_PER_TASK:-32}
MEM_GB=$(( ${SLURM_MEM_PER_NODE:-512000} / 1024 ))

eval "$(conda shell.bash hook)"
conda activate kairos-assemble      # conda env create -f ../../tools/assemble.yaml

mkdir -p trimmed assemblies
if [ ! -s trimmed/${run}_1.fq.gz ]; then
    fastp -i reads/${run}_1.fastq.gz -I reads/${run}_2.fastq.gz \
          -o trimmed/${run}_1.fq.gz -O trimmed/${run}_2.fq.gz \
          -w $(( T < 16 ? T : 16 )) -j trimmed/${run}.fastp.json -h trimmed/${run}.fastp.html
fi
if [ -d assemblies/${run} ] && [ ! -s assemblies/${run}/contigs.fasta ]; then
    metaspades.py --restart-from last -o assemblies/${run} -t $T -m $MEM_GB
else
    metaspades.py -1 trimmed/${run}_1.fq.gz -2 trimmed/${run}_2.fq.gz \
        -o assemblies/${run} -t $T -m $MEM_GB
fi
# build_catalogue.sh prefixes contig IDs with the file name, so give each run its own
ln -sfn ${run}/contigs.fasta assemblies/${run}.fasta
echo "done: assemblies/${run}.fasta"
