# ─────────────────────────────────────────────────────────────────────────────
# orf_calling.smk
# Stage the catalogue, split it for parallel ORF calling, run Prodigal, and merge
# results back into single .faa / .ffn / .gff files.
#
# Expects the following to already be defined by the top-level Snakefile:
#   THREADS, CATALOGUE, CAT, DET
# ─────────────────────────────────────────────────────────────────────────────

# ── Copy the catalogue into results/ so split/ORF names follow catalogue_name ──
rule stage_catalogue:
    input:
        CATALOGUE
    output:
        f"{DET}.fa"
    conda:
        "../envs/orfs.yaml"
    shell:
        # seqkit seq -i keeps only the ID (first word) of each header; duplicate IDs fail
        "seqkit seq -i {input} | seqkit rmdup -n -D {output}.dups > {output}"
        " && if [ -s {output}.dups ]; then echo 'Duplicate contig IDs in catalogue:' >&2;"
        " cat {output}.dups >&2; rm {output}; exit 1; fi; rm -f {output}.dups"


# ── Split catalogue for parallel prodigal ─────────────────────────────────────
checkpoint split_fasta:
    input:
        f"{DET}.fa"
    output:
        directory(f"results/detect/split_{CAT}")
    params:
        nparts = THREADS
    conda:
        "../envs/orfs.yaml"
    shell:
        "seqkit split2 -p {params.nparts} -O {output} {input}"


# ── Index catalogue ────────────────────────────────────────────────────────────
rule index_fasta:
    input:
        f"{DET}.fa"
    output:
        f"{DET}.fa.fai"
    conda:
        "../envs/orfs.yaml"
    shell:
        "samtools faidx {input} -o {output}"


# ── ORF calling on each split part ─────────────────────────────────────────────
rule prodigal:
    input:
        f"results/detect/split_{CAT}/{CAT}.part_{{splitnum}}.fa"
    output:
        faa=f"results/detect/split_{CAT}/{CAT}.part_{{splitnum}}.faa",
        ffn=f"results/detect/split_{CAT}/{CAT}.part_{{splitnum}}.ffn",
        gff=f"results/detect/split_{CAT}/{CAT}.part_{{splitnum}}.gff"
    conda:
        "../envs/orfs.yaml"
    shell:
        "prodigal -i {input} -a {output.faa} -d {output.ffn} -f gff -o {output.gff} -p meta"


# ── Helpers to resolve the dynamic split parts from the split_fasta checkpoint ─
def get_split_parts(wildcards, ext):
    split_dir = checkpoints.split_fasta.get().output[0]
    splitnums, = glob_wildcards(f"{split_dir}/{CAT}.part_{{splitnum}}.fa")  # always glob .fa
    return expand(f"results/detect/split_{CAT}/{CAT}.part_{{splitnum}}.{ext}",
                  splitnum=splitnums)

def get_split_parts_faa(wildcards):
    return get_split_parts(wildcards, "faa")

def get_split_parts_ffn(wildcards):
    return get_split_parts(wildcards, "ffn")

def get_split_parts_gff(wildcards):
    return get_split_parts(wildcards, "gff")


# ── Merge per-part outputs back into single files ─────────────────────────────
rule merge_results:
    input:
        faa = get_split_parts_faa,
        ffn = get_split_parts_ffn,
        gff = get_split_parts_gff
    output:
        faa=f"{DET}.faa",
        ffn=f"{DET}.ffn",
        gff=f"{DET}.gff"
    shell:
        "cat {input.faa} > {output.faa} && cat {input.ffn} > {output.ffn} && cat {input.gff} > {output.gff}"


# ── Generate ORF coordinate table from Prodigal headers ───────────────────────
rule generate_coords:
    input:
        ffn=f"{DET}.ffn"
    output:
        f"{DET}.ffn_coords.tsv"
    shell:
        """
        grep ">" {input.ffn} | sed "s,>,,g" | cut -d " " -f 1,3,5,7,9 | sed "s, ,\t,g" > {output}
        """
