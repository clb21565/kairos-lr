# ─────────────────────────────────────────────────────────────────────────────
# taxonomy.smk
# Assign taxonomy to whole contigs and to individual candidate HGT regions
# using mmseqs taxonomy against a GTDB-formatted database.
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, GTDB_DB   (path to the mmseqs-formatted GTDB database)
# ─────────────────────────────────────────────────────────────────────────────

rule taxonomic_annotation:
    input:
        fasta="{sample}.fasta",
        region_fasta="result/{sample}.potential_hgt_regions.fasta"
    output:
        taxonomy_result="result/{sample}.assignRes2.tsv",
        region_taxonomy_result="result/{sample}.potential_hgt_regions.taxonomy.tsv"
    threads: THREADS
    params:
        database=GTDB_DB,
        outname="result/{sample}.mmseqs_taxonomy"
    shell:
        """
        mmseqs createdb {input.fasta} {params.outname}.contigs
        mmseqs taxonomy {params.outname}.contigs {params.database} {params.outname}.assignments tmpFolder --tax-lineage 1 --majority 0.5 --vote-mode 1 --lca-mode 3 --orf-filter 1 --threads {threads}
        mmseqs createtsv {params.outname}.contigs {params.outname}.assignments {output.taxonomy_result}
        mmseqs createdb {input.region_fasta} {params.outname}.regions
        mmseqs taxonomy {params.outname}.regions {params.database} {params.outname}.regions.assignments tmpFolder --tax-lineage 1 --majority 0.5 --vote-mode 1 --lca-mode 3 --orf-filter 1 --threads {threads}
        mmseqs createtsv {params.outname}.regions {params.outname}.regions.assignments {output.region_taxonomy_result}
        """

