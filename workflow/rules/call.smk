# ─────────────────────────────────────────────────────────────────────────────
# call.smk
# Pool each observation's samples, call every shared segment present / absent /
# undetermined (detect_state.py), then call in situ HGT events over time (insitu.py).
#
# Expects the following to already be defined by the top-level Snakefile:
#   DET, SAMPLES, CALL, detect_state, insitu
# ─────────────────────────────────────────────────────────────────────────────

CALLDIR = "results/call"

rule detect_state:
    input:
        samples=config["samples"],
        events=f"{DET}.features/events.tsv",
        counts=expand("results/assess/{sample}/{sample}.{kind}_counts.tsv",
                      sample=SAMPLES.index if SAMPLES is not None else [],
                      kind=["breakpoint", "empty_site", "contig"])
    output:
        segments=f"{CALLDIR}/segments.tsv",
        states=f"{CALLDIR}/states.tsv",
        event_segments=f"{CALLDIR}/event_segments.tsv"
    conda:
        "../envs/python.yaml"
    shell:
        "python {detect_state}"
        " --samples {input.samples}"
        " --events {input.events}"
        " --counts-dir results/assess"
        " --min-reads {CALL[min_reads]}"
        " --alpha {CALL[alpha]}"
        " --absent-fraction {CALL[absent_fraction]}"
        " --min-boundaries {CALL[min_boundaries]}"
        " --contig-min-breadth {CALL[contig_min_breadth]}"
        " --contig-breadth-ratio {CALL[contig_breadth_ratio]}"
        " --out-segments {output.segments}"
        " --out-states {output.states}"
        " --out-event-segments {output.event_segments}"


rule insitu:
    input:
        states=f"{CALLDIR}/states.tsv",
        segments=f"{CALLDIR}/segments.tsv",
        event_segments=f"{CALLDIR}/event_segments.tsv",
        taxonomy=[f"{DET}.assignRes2.tsv"] if CALL["use_taxonomy"] else []
    output:
        events=f"{CALLDIR}/insitu_events.tsv",
        excluded=f"{CALLDIR}/insitu_excluded.tsv"
    params:
        observed=",".join(CALL["observed_compartments"]),
        prior=",".join(CALL["prior_compartments"]),
        taxonomy=lambda wc, input: f"--contig-taxonomy {input.taxonomy[0]}" if input.taxonomy else ""
    conda:
        "../envs/python.yaml"
    shell:
        "python {insitu}"
        " --states {input.states}"
        " --segments {input.segments}"
        " --event-segments {input.event_segments}"
        " --observed '{params.observed}'"
        " --prior '{params.prior}'"
        " {params.taxonomy}"
        " --lineage-filter {CALL[lineage_filter]}"
        " --out {output.events}"
        " --excluded {output.excluded}"
