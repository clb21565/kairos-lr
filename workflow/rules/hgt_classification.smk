# ─────────────────────────────────────────────────────────────────────────────
# hgt_classification.smk
# Classify HGT events (criteria 1 + criteria 2), run union-find grouping,
# and generate donor-recipient annotations and group reports.
#
# Expects the following to already be defined by the top-level Snakefile:
#   get_HGT
#
# NOTE: get_hgt-v2.py reads its inputs by joining --path + sample-prefixed
# filenames (e.g. {path}/{sample}.assignRes2.tsv). Since all of these files
# live flat under result/, --path must be the result/ directory itself —
# NOT result/{sample}/. See input: block below for the exact files needed.
# get_hgt-v2.py writes ~25 output TSVs into --out-dir; we declare that whole
# directory as the rule's output rather than listing every file individually.
# ─────────────────────────────────────────────────────────────────────────────

rule get_hgt:
    input:
        region_overlaps  = "result/{sample}.kairos_region_overlaps.tsv",
        region_metadata  = "result/{sample}.kairos_region_metadata.tsv",
        contig_taxonomy  = "result/{sample}.assignRes2.tsv",
        region_taxonomy  = "result/{sample}.potential_hgt_regions.taxonomy.tsv"
    output:
        outdir = directory("result/{sample}.hgt/")
    params:
        path = "result/"
    shell:
        "python {get_HGT} --path {params.path} --sample {wildcards.sample} --out-dir {output.outdir}"
