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


# Memory scales with catalogue size. Measured peaks (40 Mbp test catalogue):
# strobealign index ~4.2 B/bp, minimap2 map-ont index ~8 B/bp; bwa-mem2 index
# build is documented at ~28 B/bp. Requests below are about 2x those peaks, plus
# a fixed allowance for samtools sort (SORT_THREADS x SORT_MEM).
# Real 4.4 Gbp catalogue (2026-11-02): strobealign index 6.2 B/bp; strobealign and
# minimap2 align jobs peaked at their 9 and 16 B/bp requests, so align has more headroom.
SORT_THREADS = 4
SORT_MEM_MB = int(config["sort_mem_mb"])   # per sort thread

def _cat_mb(input):
    import os
    return os.path.getsize(input.ref if hasattr(input, "ref") else input[0]) / 1e6

def _mem(bytes_per_bp, extra_mb=0):
    def f(wc, input, attempt):
        return int(max(16000, bytes_per_bp * _cat_mb(input) + extra_mb) * attempt)
    return f

def align_mem(wc, input, attempt):
    if read_type(wc) not in SHORT_TYPES:
        bpb = 20                         # minimap2 index + mapping buffers
    elif ASSESS["short_aligner"] == "bwa-mem2":
        bpb = 12
    else:
        bpb = 12
    return _mem(bpb, SORT_THREADS * SORT_MEM_MB + 4000)(wc, input, attempt)


rule index_strobealign:
    input:
        f"{DET}.fa"
    output:
        f"{DET}.fa.r{{rl}}.sti"
    threads: THREADS
    resources:
        mem_mb=_mem(9)
    conda:
        "../envs/align.yaml"
    shell:
        "strobealign --create-index -t {threads} -r {wildcards.rl} {input}"


rule index_bwamem2:
    input:
        f"{DET}.fa"
    output:
        f"{DET}.fa.bwt.2bit.64"
    resources:
        mem_mb=_mem(56)
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
        cmd=lambda wc, input, threads: align_command(wc, input, threads),
        sort_threads=lambda wc, threads: min(threads, SORT_THREADS)
    threads: THREADS
    resources:
        mem_mb=align_mem
    conda:
        "../envs/align.yaml"
    shell:
        """
        {params.cmd} | samtools sort -@ {params.sort_threads} -m {SORT_MEM_MB}M -T {output.bam}.tmp -o {output.bam} -
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
        empty_sites=f"{DET}.features/empty_sites.tsv",
        junctions=f"{DET}.features/junctions.fa"
    output:
        breakpoints=f"{ASM}/{{sample}}/{{sample}}.breakpoint_counts.tsv",
        empty_sites=f"{ASM}/{{sample}}/{{sample}}.empty_site_counts.tsv"
    params:
        read_type=read_type,
        overhang=lambda wc: ASSESS["overhang"]["short" if read_type(wc) in SHORT_TYPES else "long"],
        junctions=lambda wc, input: f"--junctions {input.junctions}" if ASSESS["second_pass"] else ""
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
        " {params.junctions}"
        " --out-breakpoints {output.breakpoints}"
        " --out-empty-sites {output.empty_sites}"
