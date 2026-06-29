# ─────────────────────────────────────────────────────────────────────────────
# annotation.smk
# Functional annotation of predicted proteins against reference databases:
#   - RIPDB     (replication initiation proteins, used for plasmid typing)
#   - mobileOG  (mobile genetic element ontology)
#   - CARD      (antibiotic resistance genes)
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, RIPDB, MOBILEOG, CARD, RIPid, RIPqc, RIPsc, MGEid, MGEqc, MGEsc,
#   ARGid, ARGqc, ARGsc
# ─────────────────────────────────────────────────────────────────────────────

rule diamond_ripdb:
    input:
        "result/{sample}.faa"
    output:
        "result/{sample}-ripdb.tsv"
    params:
        db = RIPDB,
        identity = RIPid,
        qc = RIPqc,
        sc = RIPsc
    threads: THREADS
    shell:
        """
        diamond blastp \
            -q {input} \
            -d {params.db} \
            -p {threads} \
            -k 1 \
            --id {params.identity} \
            --query-cover {params.qc} \
            --subject-cover {params.sc} \
            --outfmt 6 qtitle stitle pident bitscore evalue \
            -o {output}
        """


rule diamond_mobileog:
    input:
        "result/{sample}.faa"
    output:
        "result/{sample}-mobileog.tsv"
    params:
        db = MOBILEOG,
        identity = MGEid,
        qc = MGEqc,
        sc = MGEsc
    threads: THREADS
    shell:
        """
        diamond blastp \
            -q {input} \
            -d {params.db} \
            -p {threads} \
            -k 1 \
            --id {params.identity} \
            --query-cover {params.qc} \
            --subject-cover {params.sc} \
            --outfmt 6 qtitle stitle pident bitscore evalue \
            -o {output}
        """


rule diamond_card:
    input:
        "result/{sample}.faa"
    output:
        "result/{sample}-card.tsv"
    params:
        db = CARD,
        identity = ARGid,
        qc = ARGqc,
        sc = ARGsc
    threads: THREADS
    shell:
        """
        diamond blastp \
            -q {input} \
            -d {params.db} \
            -p {threads} \
            -k 1 \
            --id {params.identity} \
            --query-cover {params.qc} \
            --subject-cover {params.sc} \
            --outfmt 6 qtitle stitle pident bitscore evalue \
            -o {output}
        """

