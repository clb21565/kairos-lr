# ─────────────────────────────────────────────────────────────────────────────
# hgt_detection.smk
# Detect putative HGT regions from clustered ORFs, merge nearby hits, and
# extract/classify the candidate regions (regions-v2.py).
#
# Expects the following to already be defined by the top-level Snakefile:
#   DET, MIN_CONTIG_ORFS, MERGE_DISTANCE, EDGE_PADDING, SLOP, kairos_dd, regions
# ─────────────────────────────────────────────────────────────────────────────

rule run_derep_detect:
    input:
        coords = f"{DET}.ffn_coords.tsv",
        clusts = f"{DET}.ffn.clusts_cluster.tsv",
        fai    = f"{DET}.fa.fai"
    params:
        out_prefix = f"{DET}.kairos"
    output:
        f"{DET}.kairos_potential_hgts.bed",
        f"{DET}.kairos_deduplicated_overlaps.tsv"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {kairos_dd} --input_clust_file {input.clusts} --minimum_orfs {MIN_CONTIG_ORFS} --max_overlap 1 --out_prefix {params.out_prefix} --coords_file {input.coords} --contigMD {input.fai}
        """


rule get_merged_bed:
    input:
        f"{DET}.kairos_potential_hgts.bed"
    output:
        f"{DET}.kairos_merged_potential_hgts.bed"
    conda:
        "../envs/python.yaml"
    shell:
        """
        bedtools sort -i {input} | bedtools merge -c 1,5 -o count,collapse -d {MERGE_DISTANCE} > {output}
        """


rule run_regions:
    input:
        fasta=f"{DET}.fa",
        hgt_bed=f"{DET}.kairos_merged_potential_hgts.bed",
        contig_overlaps=f"{DET}.kairos_deduplicated_overlaps.tsv",
        fai=f"{DET}.fa.fai"
    output:
        region_overlaps=f"{DET}.kairos_region_overlaps.tsv",
        filtered_overlaps=f"{DET}.kairos_filtered_overlaps.tsv",
        region_metadata=f"{DET}.kairos_region_metadata.tsv",
        fasta=f"{DET}.potential_hgt_regions.fasta",
        bed=f"{DET}.potential_hgt_regions.bed"
    conda:
        "../envs/python.yaml"
    shell:
        "python {regions}"
        " -b {input.hgt_bed}"
        " -o {input.contig_overlaps}"
        " -f {input.fai}"
        " -p {EDGE_PADDING}"
        " --region-overlaps-output {output.region_overlaps}"
        " --filtered-overlaps-output {output.filtered_overlaps}"
        " --region-metadata-output {output.region_metadata}"
        " --filter-reciprocals"
        " --input-fasta {input.fasta}"
        " --extract-regions"
        " --extracted-fasta {output.fasta}"
        " --bed-output {output.bed}"
        " --slop {SLOP}"
