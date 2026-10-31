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
                      kind=["breakpoint", "empty_site"])
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
        f"{CALLDIR}/insitu_events.tsv"
    params:
        observed=",".join(CALL["observed_fractions"]),
        prior=",".join(CALL["prior_fractions"]),
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
        " --out {output}"
