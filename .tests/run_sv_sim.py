"""
Local check for tools/sv_compare.py on the tier 1 simulation (run .tests/run_sim.py first).

Adds one long-read sample, t3_long_ins, in which host4 (no cassette in the catalogue) carries
the host1 cassette in 60% of cells. Then Sniffles2 is run on t1_long, t2_long (assess BAMs) and
t3_long_ins, the calls are merged, and sv_compare should report:
  - the host1 empty site as a DEL matching ES (VAF ~1 at t1, ~0.5 at t2), with the same
    read counts as kairos;
  - the host4 insertion as novel_ins_kairos_segment (a known segment at a new site), VAF ~0.6.

Needs sniffles, minimap2, samtools, pysam, pandas on PATH. On this VM:
  MAMBA_ROOT_PREFIX=~/.local/micromamba/root ~/.local/micromamba/bin/micromamba run -n ksv \\
      python .tests/run_sv_sim.py
"""

import argparse
import importlib.util
import os
import subprocess
import sys

import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
INSERT_AT = 4800
HOST4_INS_FRAC = 0.6


def load_sim():
    spec = importlib.util.spec_from_file_location("make_sim_data", os.path.join(HERE, "make_sim_data.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def run(sim, out, threads):
    os.makedirs(out, exist_ok=True)
    cat = os.path.join(sim, "results", "detect", "catalogue.fa")
    feat = os.path.join(sim, "results", "detect", "catalogue.features")
    es = pd.read_csv(os.path.join(feat, "empty_sites.tsv"), sep="\t").set_index("contig")
    fa = pysam.FastaFile(cat)
    h1 = es.loc["asm1_host1"]
    cassette = fa.fetch("asm1_host1", int(h1["shared_start"]), int(h1["shared_end"]))
    host4 = fa.fetch("asm1_host4")
    genomes = {"asm1_host4": host4, "asm2_host2": fa.fetch("asm2_host2"),
               "asm1_host4_ins": host4[:INSERT_AT] + cassette + host4[INSERT_AT:]}
    fq = os.path.join(out, "t3_long_ins.fq.gz")
    load_sim().long_reads(genomes, {"asm1_host4_ins": round(20 * HOST4_INS_FRAC),
                                    "asm1_host4": round(20 * (1 - HOST4_INS_FRAC)), "asm2_host2": 20}, fq)
    bam = os.path.join(out, "t3_long_ins.bam")
    sh(f"minimap2 -ax map-ont --secondary=yes -N 5 -t {threads} {cat} {fq} 2>/dev/null"
       f" | samtools sort -o {bam} - && samtools index {bam}")

    bams = {s: os.path.join(sim, "results", "assess", s, f"{s}.bam") for s in ("t1_long", "t2_long")}
    bams["t3_long_ins"] = bam
    snfs = []
    for s, b in bams.items():
        snf = os.path.join(out, f"{s}.snf")
        sh(f"sniffles --input {b} --snf {snf} --vcf {out}/{s}.vcf --mosaic --mosaic-include-germline"
           f" --all-contigs --minsvlen 300 --sample-id {s} --threads {threads} --allow-overwrite --quiet")
        snfs.append(snf)
    vcf = os.path.join(out, "merged.vcf")
    sh(f"sniffles --input {' '.join(snfs)} --vcf {vcf} --all-contigs --threads {threads} --allow-overwrite --quiet")
    sh(f"{sys.executable} {os.path.join(REPO, 'tools', 'sv_compare.py')} --vcf {vcf} --features {feat}"
       f" --catalogue {cat} --assess-dir {os.path.join(sim, 'results', 'assess')} --outdir {out}"
       f" --bam t3_long_ins={bam} --threads {threads}")


def check(out):
    calls = pd.read_csv(os.path.join(out, "sv_calls.tsv"), sep="\t", dtype=str, keep_default_na=False)
    calls["start"] = calls["start"].astype(int)
    samp = pd.read_csv(os.path.join(out, "sv_samples.tsv"), sep="\t")
    failures = []

    def expect(ok, msg):
        print(("ok   " if ok else "FAIL ") + msg)
        if not ok:
            failures.append(msg)

    d = calls[(calls["svtype"] == "DEL") & (calls["contig"] == "asm1_host1")]
    expect(len(d) == 1 and d.iloc[0]["class"] == "kairos_empty_site",
           f"host1 DEL matches its empty site ({list(d['class'])})")
    if len(d):
        s = samp[samp["sv_id"] == d.iloc[0]["sv_id"]].set_index("sample")
        expect(s.at["t1_long", "vaf"] > 0.9, f"t1_long DEL VAF ~1 ({s.at['t1_long', 'vaf']})")
        expect(0.3 < s.at["t2_long", "vaf"] < 0.7, f"t2_long DEL VAF ~0.5 ({s.at['t2_long', 'vaf']})")
        expect(s.at["t3_long_ins", "dv"] == 0, "no host1 DEL reads at t3 (host1 absent)")
        k = s.dropna(subset=["kairos_empty"])
        expect((k["dv"] - k["kairos_empty"]).abs().max() <= 2,
               f"Sniffles DV matches kairos empty reads ({list(k['dv'])} vs {list(k['kairos_empty'])})")
        expect((k["alt_reads"] - k["kairos_empty"]).abs().max() <= 2 and
               (k["ref_reads"] - k["kairos_span"]).abs().max() <= 2,
               f"recount matches kairos (alt {list(k['alt_reads'])}, ref {list(k['ref_reads'])})")

    i = calls[(calls["svtype"] == "INS") & (calls["contig"] == "asm1_host4")]
    expect(len(i) == 1, f"one INS on host4 ({len(i)})")
    if len(i):
        r = i.iloc[0]
        expect(abs(r["start"] - INSERT_AT) <= 20, f"INS position ~{INSERT_AT} ({r['start']})")
        expect(r["class"] == "novel_ins_kairos_segment" and r["hit_events"] != "",
               f"INS is a known kairos segment at a new site ({r['class']}, {r['hit_events']})")
        s = samp[samp["sv_id"] == r["sv_id"]].set_index("sample")
        expect(abs(s.at["t3_long_ins", "vaf_recount"] - HOST4_INS_FRAC) < 0.15,
               f"t3_long_ins recounted INS VAF ~{HOST4_INS_FRAC} ({s.at['t3_long_ins', 'vaf_recount']};"
               f" Sniffles2 VAF {s.at['t3_long_ins', 'vaf']})")
        expect(s.at["t1_long", "alt_reads"] == 0, "no recounted host4 INS reads at t1")
        expect(s.at["t1_long", "dv"] == 0, "no host4 INS reads at t1")
    other = calls[~calls.index.isin(d.index.union(i.index))]
    expect(len(other) == 0, f"no other SVs ({len(other)})")
    return failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", default=os.path.join(HERE, "sim"), help="run_sim.py output directory")
    ap.add_argument("--outdir", default=os.path.join(HERE, "sim_sv"))
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    run(a.sim, a.outdir, a.threads)
    failures = check(a.outdir)
    print(f"\n{len(failures)} failure(s)" if failures else "\nall checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
