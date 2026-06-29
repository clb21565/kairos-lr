# ─────────────────────────────────────────────────────────────────────────────
# clustering.smk
# Cluster predicted proteins with mmseqs2 as the basis for HGT region detection.
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, MINSEQID, MINQUERYCOVER
# ─────────────────────────────────────────────────────────────────────────────

rule cluster_proteins:
    input:
        ffn="result/{sample}.ffn"
    output:
        "result/{sample}.ffn.clusts_cluster.tsv"
    threads: THREADS
    shell:
        """
        mmseqs easy-cluster {input.ffn} {input}.clusts clusts.tmp --min-seq-id {MINSEQID} -c {MINQUERYCOVER} --cov-mode 1 --threads {threads}
        """

