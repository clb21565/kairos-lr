# ─────────────────────────────────────────────────────────────────────────────
# assess.smk
# Map each sample's reads to the whole catalogue (competitive mapping) and count
# reads at every breakpoint, its flank and its empty site (count_features.py).
#
#   short reads (short_pe, short_se): strobealign (default) or bwa-mem2
#   long reads  (ont, hifi):          minimap2, preset per read type
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, DET, SAMPLES, ASSESS, count_features
# ─────────────────────────────────────────────────────────────────────────────

ASM = "results/assess"
SHORT_TYPES = {"short_pe", "short_se"}


def sample_reads(wc):
    return [r.strip() for r in SAMPLES.loc[wc.sample, "reads"].split(",")]

def read_type(wc):
    return SAMPLES.loc[wc.sample, "read_type"]

def aligner_index(wc):
    if read_type(wc) not in SHORT_TYPES:
        return []    # minimap2 indexes the catalogue on the fly
    if ASSESS["short_aligner"] == "bwa-mem2":
        return f"{DET}.fa.bwt.2bit.64"
    return f"{DET}.fa.r{ASSESS['short_read_length']}.sti"


rule index_strobealign:
    input:
        f"{DET}.fa"
    output:
        f"{DET}.fa.r{{rl}}.sti"
    threads: THREADS
    conda:
        "../envs/align.yaml"
    shell:
        "strobealign --create-index -t {threads} -r {wildcards.rl} {input}"


rule index_bwamem2:
    input:
        f"{DET}.fa"
    output:
        f"{DET}.fa.bwt.2bit.64"
    conda:
        "../envs/align.yaml"
    shell:
        "bwa-mem2 index {input}"


def align_command(wc, input, threads):
    rt = read_type(wc)
    reads = " ".join(input.reads)
    if rt in SHORT_TYPES:
        if ASSESS["short_aligner"] == "bwa-mem2":
            return f"bwa-mem2 mem -t {threads} {input.ref} {reads}"
        return (f"strobealign -t {threads} --use-index -r {ASSESS['short_read_length']}"
                f" {input.ref} {reads}")
    preset = ASSESS["long_read_presets"][rt]
    return f"minimap2 -ax {preset} --secondary=yes -N 5 -t {threads} {input.ref} {reads}"


rule align:
    input:
        ref=f"{DET}.fa",
        index=aligner_index,
        reads=sample_reads
    output:
        bam=f"{ASM}/{{sample}}/{{sample}}.bam",
        bai=f"{ASM}/{{sample}}/{{sample}}.bam.bai"
    params:
        cmd=lambda wc, input, threads: align_command(wc, input, threads)
    threads: THREADS
    conda:
        "../envs/align.yaml"
    shell:
        """
        {params.cmd} | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """


rule contig_coverage:
    input:
        bam=f"{ASM}/{{sample}}/{{sample}}.bam",
        bai=f"{ASM}/{{sample}}/{{sample}}.bam.bai"
    output:
        f"{ASM}/{{sample}}/{{sample}}.coverage.tsv"
    params:
        mapq=ASSESS["min_mapq"]
    conda:
        "../envs/align.yaml"
    shell:
        "samtools coverage -q {params.mapq} {input.bam} > {output}"


rule count_features:
    input:
        bam=f"{ASM}/{{sample}}/{{sample}}.bam",
        bai=f"{ASM}/{{sample}}/{{sample}}.bam.bai",
        breakpoints=f"{DET}.features/breakpoints.tsv",
        flanks=f"{DET}.features/flanks.bed",
        empty_sites=f"{DET}.features/empty_sites.tsv"
    output:
        breakpoints=f"{ASM}/{{sample}}/{{sample}}.breakpoint_counts.tsv",
        empty_sites=f"{ASM}/{{sample}}/{{sample}}.empty_site_counts.tsv"
    params:
        read_type=read_type,
        overhang=lambda wc: ASSESS["overhang"]["short" if read_type(wc) in SHORT_TYPES else "long"]
    conda:
        "../envs/python.yaml"
    shell:
        "python {count_features}"
        " --bam {input.bam}"
        " --sample {wildcards.sample}"
        " --read-type {params.read_type}"
        " --breakpoints {input.breakpoints}"
        " --flanks {input.flanks}"
        " --empty-sites {input.empty_sites}"
        " --overhang {params.overhang}"
        " --min-mapq {ASSESS[min_mapq]}"
        " --junction-min-mapq {ASSESS[junction_min_mapq]}"
        " --max-insert {ASSESS[max_insert]}"
        " --empty-tolerance {ASSESS[empty_tolerance]}"
        " --out-breakpoints {output.breakpoints}"
        " --out-empty-sites {output.empty_sites}"
