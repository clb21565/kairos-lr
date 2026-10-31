# ─────────────────────────────────────────────────────────────────────────────
# hgt_classification.smk
# Classify HGT events (criteria 1 + criteria 2), run union-find grouping,
# and generate donor-recipient annotations and group reports.
#
# Expects the following to already be defined by the top-level Snakefile:
#   DET, CAT, get_HGT
#
# NOTE: get_hgt-v2.py reads its inputs by joining --path + prefixed filenames
# (e.g. {path}/{CAT}.assignRes2.tsv), so --path is results/detect/ and --sample
# is the catalogue name. It writes ~25 output TSVs into --out-dir; we declare
# that whole directory as the rule's output rather than listing every file.
# ─────────────────────────────────────────────────────────────────────────────

rule get_hgt:
    input:
        region_overlaps  = f"{DET}.kairos_region_overlaps.tsv",
        region_metadata  = f"{DET}.kairos_region_metadata.tsv",
        contig_taxonomy  = f"{DET}.assignRes2.tsv",
        region_taxonomy  = f"{DET}.potential_hgt_regions.taxonomy.tsv"
    output:
        outdir = directory(f"{DET}.hgt/")
    params:
        path = "results/detect/"
    conda:
        "../envs/python.yaml"
    shell:
        "python {get_HGT} --path {params.path} --sample {CAT} --out-dir {output.outdir}"
