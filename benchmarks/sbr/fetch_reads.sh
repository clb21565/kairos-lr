#!/bin/bash
# Download the six deep SBR runs (~300 Gb of bases, ~116 GB gzipped) from ENA.
# Run from benchmarks/sbr: bash fetch_reads.sh
set -euo pipefail
mkdir -p reads
for run in SRR26192772 SRR26192773 SRR26192774 SRR26192775 SRR26192776 SRR26192777; do
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" \
        | tail -n +2 | while IFS=$'\t' read -r _ ftps md5s; do
        IFS=';' read -ra F <<< "$ftps"; IFS=';' read -ra M <<< "$md5s"
        for i in "${!F[@]}"; do
            f=reads/$(basename "${F[$i]}")
            [ -s "$f" ] && echo "${M[$i]}  $f" | md5sum -c --quiet && continue
            wget -q -c -O "$f" "https://${F[$i]}"
            echo "${M[$i]}  $f" | md5sum -c
        done
    done
done
