# ─────────────────────────────────────────────────────────────────────────────
# features.smk
# Junction, flank and empty-site features for every pair of shared regions
# (build_features.py). These are what assess counts reads against.
#
# Expects the following to already be defined by the top-level Snakefile:
#   DET, FEATURES, build_features
# ─────────────────────────────────────────────────────────────────────────────

rule build_features:
    input:
        fasta=f"{DET}.fa",
        fai=f"{DET}.fa.fai",
        region_metadata=f"{DET}.kairos_region_metadata.tsv",
        region_overlaps=f"{DET}.kairos_region_overlaps.tsv"
    output:
        events=f"{DET}.features/events.tsv",
        breakpoints=f"{DET}.features/breakpoints.tsv",
        flanks=f"{DET}.features/flanks.bed",
        empty_sites=f"{DET}.features/empty_sites.tsv",
        junctions=f"{DET}.features/junctions.fa"
    params:
        outdir=f"{DET}.features"
    conda:
        "../envs/python.yaml"
    shell:
        "python {build_features}"
        " --catalogue {input.fasta}"
        " --region-metadata {input.region_metadata}"
        " --region-overlaps {input.region_overlaps}"
        " --outdir {params.outdir}"
        " --window {FEATURES[window]}"
        " --min-identity {FEATURES[min_identity]}"
        " --min-flank {FEATURES[min_flank]}"
        " --flank-len {FEATURES[flank_len]}"
        " --junction-flank {FEATURES[junction_flank]}"
        " --snap {FEATURES[snap]}"
