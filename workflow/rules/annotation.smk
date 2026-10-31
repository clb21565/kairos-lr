# ─────────────────────────────────────────────────────────────────────────────
# annotation.smk
# Functional annotation of predicted proteins against reference databases:
#   - ripdb     (replication initiation proteins, used for plasmid typing)
#   - mobileog  (mobile genetic element ontology)
#   - card      (antibiotic resistance genes)
# Only the databases set in the config are run (ANNOTATION_DBS).
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, DET, ANNOTATION_DBS
# ─────────────────────────────────────────────────────────────────────────────

rule diamond_annotate:
    input:
        f"{DET}.faa"
    output:
        f"{DET}-{{db}}.tsv"
    wildcard_constraints:
        db = "|".join(ANNOTATION_DBS) or "none"
    params:
        db = lambda wc: ANNOTATION_DBS[wc.db]["db"],
        identity = lambda wc: ANNOTATION_DBS[wc.db]["id"],
        qc = lambda wc: ANNOTATION_DBS[wc.db]["qc"],
        sc = lambda wc: ANNOTATION_DBS[wc.db]["sc"]
    threads: THREADS
    conda:
        "../envs/diamond.yaml"
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
