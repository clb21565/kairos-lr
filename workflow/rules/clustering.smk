# ─────────────────────────────────────────────────────────────────────────────
# clustering.smk
# Cluster predicted ORFs (nucleotide) with mmseqs2 as the basis for HGT region detection.
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, DET, MINSEQID, MINQUERYCOVER
# ─────────────────────────────────────────────────────────────────────────────

rule cluster_orfs:
    input:
        ffn=f"{DET}.ffn"
    output:
        f"{DET}.ffn.clusts_cluster.tsv"
    params:
        prefix=f"{DET}.ffn.clusts",
        tmp=f"{DET}.mmseqs_cluster_tmp"
    threads: THREADS
    conda:
        "../envs/mmseqs.yaml"
    shell:
        """
        mmseqs easy-cluster {input.ffn} {params.prefix} {params.tmp} --min-seq-id {MINSEQID} -c {MINQUERYCOVER} --cov-mode 1 --threads {threads}
        rm -rf {params.tmp}
        """
