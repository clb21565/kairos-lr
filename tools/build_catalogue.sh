#!/usr/bin/env bash
# Build a dereplicated contig catalogue for kairos from many assemblies.
#
#   tools/build_catalogue.sh assemblies.txt out_dir [min_len] [threads]
#
# assemblies.txt: one assembly FASTA (.fa/.fasta, gzipped or not) per line.
# Writes out_dir/catalogue.fa (cluster representatives) and out_dir/clusters.tsv.
#
# Contigs are renamed <assembly file name>__<contig id>, because assemblers reuse names
# (NODE_1_..., contig_1) across samples.
#
# Clustering keeps near-identical sisters apart: two contigs only merge if they are >= 99%
# identical over >= 95% of BOTH sequences (--cov-mode 0). So a contig with an insert and its
# sister without it stay separate, which kairos needs to find the event.
#
# Staged use for many assemblies: run once per site (or month), then list the per-site
# catalogue.fa files in a new assemblies.txt and run again to merge.
#
# Needs seqkit and mmseqs on PATH (e.g. the workflow/envs/mmseqs.yaml + orfs.yaml envs).
set -euo pipefail

if [[ $# -lt 2 ]]; then
    sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

list=$1
out=$2
min_len=${3:-1000}
threads=${4:-8}
min_id=${MIN_ID:-0.99}
min_cov=${MIN_COV:-0.95}

mkdir -p "$out"

# 1. length filter + rename, all assemblies into one file
: > "$out/contigs.fa"
while read -r fa; do
    [[ -z "$fa" || "$fa" == \#* ]] && continue
    name=$(basename "$fa")
    name=${name%.gz}
    name=${name%.*}
    echo "adding $fa as $name" >&2
    seqkit seq -m "$min_len" "$fa" | seqkit replace -p '^(\S+).*' -r "${name}__\$1" >> "$out/contigs.fa"
done < "$list"
seqkit stats -T "$out/contigs.fa" >&2

# 2. cluster
mmseqs easy-linclust "$out/contigs.fa" "$out/clu" "$out/tmp" \
    --min-seq-id "$min_id" -c "$min_cov" --cov-mode 0 \
    --threads "$threads" -v 1

seqkit seq -i "$out/clu_rep_seq.fasta" > "$out/catalogue.fa"
mv "$out/clu_cluster.tsv" "$out/clusters.tsv"
rm -rf "$out/tmp" "$out/clu_all_seqs.fasta" "$out/clu_rep_seq.fasta" "$out/contigs.fa"
seqkit stats -T "$out/catalogue.fa" >&2
