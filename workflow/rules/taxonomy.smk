# ─────────────────────────────────────────────────────────────────────────────
# taxonomy.smk
# Assign taxonomy to whole contigs and to individual candidate HGT regions
# using mmseqs taxonomy against a GTDB-formatted database.
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, DET, GTDB_DB   (path to the mmseqs-formatted GTDB database)
# ─────────────────────────────────────────────────────────────────────────────

MMSEQS_TAXONOMY_OPTS = "--tax-lineage 1 --majority 0.5 --vote-mode 1 --lca-mode 3 --orf-filter 1"

rule taxonomy_contigs:
    input:
        f"{DET}.fa"
    output:
        f"{DET}.assignRes2.tsv"
    threads: THREADS
    params:
        database=GTDB_DB,
        outname=f"{DET}.mmseqs_taxonomy",
        tmp=f"{DET}.mmseqs_taxonomy_contigs_tmp"
    conda:
        "../envs/mmseqs.yaml"
    shell:
        """
        # clear databases left by an earlier run; mmseqs will not overwrite them
        rm -rf {params.tmp} {params.outname}.contigs {params.outname}.contigs.* {params.outname}.contigs_h*
        rm -f {params.outname}.assignments {params.outname}.assignments.*
        mmseqs createdb {input} {params.outname}.contigs
        mmseqs taxonomy {params.outname}.contigs {params.database} {params.outname}.assignments {params.tmp} {MMSEQS_TAXONOMY_OPTS} --threads {threads}
        mmseqs createtsv {params.outname}.contigs {params.outname}.assignments {output}
        rm -rf {params.tmp}
        """


rule taxonomy_regions:
    input:
        f"{DET}.potential_hgt_regions.fasta"
    output:
        f"{DET}.potential_hgt_regions.taxonomy.tsv"
    threads: THREADS
    params:
        database=GTDB_DB,
        outname=f"{DET}.mmseqs_taxonomy",
        tmp=f"{DET}.mmseqs_taxonomy_regions_tmp"
    conda:
        "../envs/mmseqs.yaml"
    shell:
        """
        # clear databases left by an earlier run; mmseqs will not overwrite them
        rm -rf {params.tmp} {params.outname}.regions {params.outname}.regions.* {params.outname}.regions_h*
        rm -f {params.outname}.regions.assignments {params.outname}.regions.assignments.*
        mmseqs createdb {input} {params.outname}.regions
        mmseqs taxonomy {params.outname}.regions {params.database} {params.outname}.regions.assignments {params.tmp} {MMSEQS_TAXONOMY_OPTS} --threads {threads}
        mmseqs createtsv {params.outname}.regions {params.outname}.regions.assignments {output}
        rm -rf {params.tmp}
        """
