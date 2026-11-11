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
# Two contigs merge if they are >= 99% identical over >= 95% of BOTH sequences (--cov-mode 0).
# A contig with an insert and its sister without it are only merged when the insert sits near a
# contig end and is small (insert <~ 5% of the contig); an insert of >= 100 bp inside a contig
# is never aligned across, so the sisters stay separate at any setting
# (tools/derep_sensitivity.py). kairos does not need the empty sister: build_features makes the
# empty site from the carrier, and a sister sharing its flanks lowers flank MAPQ.
#
# Resumable: every step writes to a temporary name and is skipped on rerun when its output
# exists. Filtered assemblies go to out_dir/filtered/<name>.fa.gz (one per assembly, kept).
# Rerunning with the same out_dir picks up where the last run died.
#
# Memory: mmseqs sizes its k-mer table from the node's RAM, not the SLURM limit, so on big
# inputs it gets OOM-killed. SPLIT_MEM (default 70% of SLURM_MEM_PER_NODE, else 0 = all RAM)
# makes it process the k-mers in several passes instead. KMER_SCALE (default 0.2, mmseqs nucl
# default) sets k-mers per contig = 21 + KMER_SCALE x length; lower it (e.g. 0.05) to cut
# memory and time further.
#
# Env overrides: MIN_ID (0.99), MIN_COV (0.95), COV_MODE (0), SPLIT_MEM, KMER_SCALE.
# Needs seqkit and mmseqs on PATH (conda env tools/build_catalogue.yaml).
set -euo pipefail

if [[ $# -lt 2 ]]; then
    sed -n '2,29p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

list=$1
out=$2
min_len=${3:-1000}
threads=${4:-8}
min_id=${MIN_ID:-0.99}
min_cov=${MIN_COV:-0.95}
cov_mode=${COV_MODE:-0}
kmer_scale=${KMER_SCALE:-0.2}
if [[ -z "${SPLIT_MEM:-}" ]]; then
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        SPLIT_MEM="$(( SLURM_MEM_PER_NODE * 7 / 10 ))M"
    elif [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
        SPLIT_MEM="$(( SLURM_MEM_PER_CPU * ${SLURM_CPUS_PER_TASK:-1} * 7 / 10 ))M"
    else
        SPLIT_MEM=0
    fi
fi

mkdir -p "$out/filtered" "$out/db"
log() { echo "[$(date '+%F %T')] $*" >&2; }

# 1. length filter + rename, one file per assembly; finished files are never rewritten
filtered=()
while read -r fa; do
    [[ -z "$fa" || "$fa" == \#* ]] && continue
    name=$(basename "$fa")
    name=${name%.gz}
    name=${name%.*}
    dest="$out/filtered/$name.fa.gz"
    filtered+=("$dest")
    if [[ -s "$dest" ]]; then
        continue
    fi
    log "filtering $fa as $name"
    seqkit seq -m "$min_len" "$fa" \
        | seqkit replace -p '^(\S+).*' -r "${name}__\$1" -o "$dest.tmp.gz"
    mv "$dest.tmp.gz" "$dest"
done < "$list"
log "${#filtered[@]} filtered assemblies in $out/filtered"

# 2. mmseqs database straight from the per-assembly files (no concatenated copy)
if [[ ! -e "$out/db/seq.done" ]]; then
    log "createdb"
    rm -f "$out"/db/seq "$out"/db/seq.* "$out"/db/seq_h*
    mmseqs createdb "${filtered[@]}" "$out/db/seq" --dbtype 2 -v 2
    touch "$out/db/seq.done"
fi

# 3. cluster; mmseqs resumes from its own tmp dir if killed part-way
if [[ ! -e "$out/db/clu.done" ]]; then
    log "linclust (split memory $SPLIT_MEM, k-mer scale $kmer_scale)"
    rm -f "$out"/db/clu "$out"/db/clu.*
    mmseqs linclust "$out/db/seq" "$out/db/clu" "$out/db/tmp" \
        --min-seq-id "$min_id" -c "$min_cov" --cov-mode "$cov_mode" \
        --split-memory-limit "$SPLIT_MEM" --kmer-per-seq-scale "$kmer_scale" \
        --threads "$threads" -v 2
    touch "$out/db/clu.done"
fi

# 4. representatives + membership
if [[ ! -s "$out/catalogue.fa" ]]; then
    log "writing catalogue"
    rm -f "$out"/db/rep "$out"/db/rep.*
    mmseqs createsubdb "$out/db/clu" "$out/db/seq" "$out/db/rep" -v 2
    mmseqs convert2fasta "$out/db/rep" "$out/db/rep.fasta" -v 2
    seqkit seq -i "$out/db/rep.fasta" -o "$out/catalogue.fa.tmp"
    mmseqs createtsv "$out/db/seq" "$out/db/seq" "$out/db/clu" "$out/clusters.tsv.tmp" \
        --threads "$threads" -v 2
    mv "$out/clusters.tsv.tmp" "$out/clusters.tsv"
    mv "$out/catalogue.fa.tmp" "$out/catalogue.fa"
fi

# mmseqs DBs are as large as the input; drop them once the catalogue is written
rm -rf "$out/db"
seqkit stats -T "$out/catalogue.fa" >&2
log "done: $out/catalogue.fa, $out/clusters.tsv (filtered inputs kept in $out/filtered)"
