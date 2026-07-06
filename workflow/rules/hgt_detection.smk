# ─────────────────────────────────────────────────────────────────────────────
# hgt_detection.smk
# Detect putative HGT regions from clustered ORFs, merge nearby hits, and
# extract/classify the candidate regions (regions2.py).
#
# Expects the following to already be defined by the top-level Snakefile:
#   PADDING, SLOP, kairos_dd, regions
# ─────────────────────────────────────────────────────────────────────────────

rule run_derep_detect:
    input:
        coords = "result/{sample}.ffn_coords.tsv",
        clusts = "result/{sample}.ffn.clusts_cluster.tsv",
        fai    = "result/{sample}.fasta.fai"
    params:
        out_prefix = "result/{sample}.kairos"
    output:
        "result/{sample}.kairos_potential_hgts.bed",
        "result/{sample}.kairos_deduplicated_overlaps.tsv"
    shell:
        """
        python {kairos_dd} --input_clust_file {input.clusts} --minimum_orfs 3 --max_overlap 1 --out_prefix {params.out_prefix} --coords_file {input.coords} --contigMD {input.fai}
        """


rule get_merged_bed:
    input:
        "result/{sample}.kairos_potential_hgts.bed"
    output:
        "result/{sample}.kairos_merged_potential_hgts.bed"
    shell:
        """
        bedtools sort -i {input} | bedtools merge -c 1,5 -o count,collapse -d {PADDING} > {output}
        """


rule run_regions:
    input:
        fasta="{sample}.fasta",
        hgt_bed="result/{sample}.kairos_merged_potential_hgts.bed",
        contig_overlaps="result/{sample}.kairos_deduplicated_overlaps.tsv",
        fai="result/{sample}.fasta.fai"
    params:
        out_prefix = "result/{sample}"
    output:
        region_overlaps="result/{sample}.kairos_region_overlaps.tsv",
        filtered_overlaps="result/{sample}.kairos_filtered_overlaps.tsv",
        region_metadata="result/{sample}.kairos_region_metadata.tsv",
        fasta="result/{sample}.potential_hgt_regions.fasta",
        bed="result/{sample}.potential_hgt_regions.bed"
    shell:
        "python {regions}"
        " -b {input.hgt_bed}"
        " -o {input.contig_overlaps}"
        " -f {input.fai}"
        " -p 1000"
        " --region-overlaps-output {output.region_overlaps}"
        " --filtered-overlaps-output {output.filtered_overlaps}"
        " --region-metadata-output {output.region_metadata}"
        " --filter-reciprocals"
        " --input-fasta {input.fasta}"
        " --extract-regions"
        " --extracted-fasta {output.fasta}"
        " --bed-output {output.bed}"
        " --slop {SLOP}"
