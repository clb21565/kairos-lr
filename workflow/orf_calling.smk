# ─────────────────────────────────────────────────────────────────────────────
# orf_calling.smk
# Split the input assembly for parallel ORF calling, run Prodigal, and merge
# results back into single per-sample .faa / .ffn / .gff files.
#
# Expects the following to already be defined by the top-level Snakefile
# (via `config.get(...)`) before this file is `include`d:
#   THREADS, FASTA, SAMPLE
# ─────────────────────────────────────────────────────────────────────────────

# ── Split input fasta for parallel prodigal ────────────────────────────────────
checkpoint split_fasta:
    input:
        FASTA
    output:
        directory("result/split_{sample}")
    params:
        nparts = THREADS
    shell:
        "seqkit split2 -p {params.nparts} -O {output} {input}"


# ── Index input fasta ──────────────────────────────────────────────────────────
rule index_fasta:
    input:
        FASTA
    output:
        f"result/{SAMPLE}.fasta.fai"
    shell:
        "samtools faidx {input} -o {output}"


# ── ORF calling on each split part ─────────────────────────────────────────────
rule prodigal:
    input:
        "result/split_{sample}/{sample}.part_{splitnum}.fasta"
    output:
        faa="result/split_{sample}/{sample}.part_{splitnum}.faa",
        ffn="result/split_{sample}/{sample}.part_{splitnum}.ffn",
        gff="result/split_{sample}/{sample}.part_{splitnum}.gff"
    shell:
        "prodigal -i {input} -a {output.faa} -d {output.ffn} -f gff -o {output.gff} -p meta"


# ── Helpers to resolve the dynamic split parts from the split_fasta checkpoint ─
def get_split_parts(wildcards, ext):
    split_dir = checkpoints.split_fasta.get(sample=wildcards.sample).output[0]
    splitnums, = glob_wildcards(f"{split_dir}/{wildcards.sample}.part_{{splitnum}}.fasta")  # always glob .fasta
    return expand("result/split_{sample}/{sample}.part_{splitnum}.{ext}",
                  sample=wildcards.sample, splitnum=splitnums, ext=ext)

def get_split_parts_faa(wildcards):
    return get_split_parts(wildcards, "faa")

def get_split_parts_ffn(wildcards):
    return get_split_parts(wildcards, "ffn")

def get_split_parts_gff(wildcards):
    return get_split_parts(wildcards, "gff")


# ── Merge per-part outputs back into single per-sample files ──────────────────
rule merge_results:
    input:
        faa = get_split_parts_faa,
        ffn = get_split_parts_ffn,
        gff = get_split_parts_gff
    output:
        faa="result/{sample}.faa",
        ffn="result/{sample}.ffn",
        gff="result/{sample}.gff"
    shell:
        "cat {input.faa} > {output.faa} && cat {input.ffn} > {output.ffn} && cat {input.gff} > {output.gff}"


# ── Generate ORF coordinate table from Prodigal headers ───────────────────────
rule generate_coords:
    input:
        ffn="result/{sample}.ffn",
        gff="result/{sample}.gff"
    output:
        "result/{sample}.ffn_coords.tsv"
    shell:
        """
        grep ">" {input.ffn} | sed "s,>,,g" | cut -d " " -f 1,3,5,7,9 | sed "s, ,\t,g" > {output}
        """
