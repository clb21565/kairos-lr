"""
Local integration test: simulate a small study, run detect (up to the features), assess and
call (without taxonomy), and check the counts, states and in situ events against what was
simulated (see make_sim_data.py).

Needs snakemake, mappy, pysam, pandas, prodigal, mmseqs2, strobealign, minimap2, samtools,
wgsim and seqkit on PATH (no --use-conda; on the HPC use run_sim.py --use-conda instead).

Usage: python .tests/run_sim.py [--outdir .tests/sim] [--use-conda] [--keep] [extra snakemake args]
"""

import argparse
import os
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def run(outdir, use_conda, extra):
    subprocess.run([sys.executable, os.path.join(HERE, "make_sim_data.py"), outdir], check=True)
    open(os.path.join(outdir, "gtdb"), "w").close()     # detect taxonomy is not run here
    with open(os.path.join(outdir, "sim_config.yaml"), "w") as f:
        f.write("call:\n  use_taxonomy: false\n")
    cmd = ["snakemake", "assess", "call", "-s", os.path.join(REPO, "Snakefile"), "-d", outdir,
           "--configfile", os.path.join(REPO, "config", "config.yaml"),
           os.path.join(outdir, "sim_config.yaml"),
           "--config", "catalogue=catalogue.fa", "samples=samples.tsv", "gtdb_db=gtdb",
           "ripdb=", "card=", "threads=2", "cluster_method=easy-linclust",
           "-c", "2", "--quiet", "all", *extra]
    if use_conda:
        cmd.append("--use-conda")
    subprocess.run(cmd, check=True)


def check(outdir):
    res = os.path.join(outdir, "results")
    feat = os.path.join(res, "detect", "catalogue.features")
    bps = pd.read_csv(os.path.join(feat, "breakpoints.tsv"), sep="\t")
    es = pd.read_csv(os.path.join(feat, "empty_sites.tsv"), sep="\t")
    samples = pd.read_csv(os.path.join(outdir, "samples.tsv"), sep="\t")

    counts, empties = [], []
    for s in samples["sample"]:
        counts.append(pd.read_csv(os.path.join(res, "assess", s, f"{s}.breakpoint_counts.tsv"), sep="\t"))
        empties.append(pd.read_csv(os.path.join(res, "assess", s, f"{s}.empty_site_counts.tsv"), sep="\t"))
    c = pd.concat(counts).merge(bps[["bp_id", "contig"]], on="bp_id")
    c["junction"] = c["span"] + c["bridge"]
    e = pd.concat(empties).merge(es[["es_id", "contig"]], on="es_id")
    e["empty"] = e["empty_split"] + e["empty_pairs"]

    J = c.groupby(["sample", "contig"])["junction"].min()        # weaker of the two breakpoints
    F = c.groupby(["sample", "contig"])["flank_depth"].mean()
    E = e.set_index(["sample", "contig"])["empty"]

    failures = []
    def expect(ok, msg):
        print(("ok    " if ok else "FAIL  ") + msg)
        if not ok:
            failures.append(msg)

    expect(len(bps) == 6 and len(es) == 3, f"6 breakpoints and 3 empty sites ({len(bps)}, {len(es)})")
    for s in ("t1_short", "t1_long"):
        expect(J[s, "asm1_host1"] == 0, f"{s}: no host1 junction reads before transfer ({J[s, 'asm1_host1']})")
        expect(E[s, "asm1_host1"] >= 5, f"{s}: host1 empty site seen ({E[s, 'asm1_host1']})")
        expect(abs(F[s, "asm1_host1"] - 20) <= 8, f"{s}: host1 flank depth ~20x ({F[s, 'asm1_host1']:.1f})")
    for s in ("t2_short", "t2_long"):
        expect(J[s, "asm1_host1"] >= 3, f"{s}: host1 junction reads after transfer ({J[s, 'asm1_host1']})")
        expect(E[s, "asm1_host1"] >= 3, f"{s}: host1 empty site still seen (50% carriers) ({E[s, 'asm1_host1']})")
    expect(F["t3_short", "asm1_host1"] < 5, f"t3_short: host1 shallow ({F['t3_short', 'asm1_host1']:.1f}x)")
    for s in samples["sample"]:
        expect(J[s, "asm2_host2"] >= 5, f"{s}: host2 junction reads ({J[s, 'asm2_host2']})")
        expect(E[s, "asm2_host2"] == 0, f"{s}: no host2 empty site ({E[s, 'asm2_host2']})")
        expect(J[s, "asm3_host3"] == 0 and F[s, "asm3_host3"] == 0, f"{s}: host3 absent")

    # stage 3: states per observation (timepoint) and the in situ event
    call = os.path.join(res, "call")
    segs = pd.read_csv(os.path.join(call, "segments.tsv"), sep="\t")
    st = pd.read_csv(os.path.join(call, "states.tsv"), sep="\t").merge(
        segs[["segment_id", "contig"]], on="segment_id")
    S = st.set_index(["timepoint", "contig"])["state"]
    expect(S[1, "asm1_host1"] == "absent", f"t1: host1 segment absent ({S[1, 'asm1_host1']})")
    expect(S[2, "asm1_host1"] == "present", f"t2: host1 segment present ({S[2, 'asm1_host1']})")
    expect(all(S[t, "asm2_host2"] == "present" for t in (1, 2, 3)), "host2 segment present at t1-t3")
    expect(all(S[t, "asm3_host3"] == "undetermined" for t in (1, 2, 3)), "host3 (not sequenced) undetermined")
    cf = st.set_index(["timepoint", "contig"])["carrier_frac_depth"][2, "asm1_host1"]
    expect(0.3 <= cf <= 0.75, f"t2: host1 carrier fraction ~0.5 ({cf})")

    ev = pd.read_csv(os.path.join(call, "insitu_events.tsv"), sep="\t")
    expect(len(ev) == 1, f"one in situ event ({len(ev)})")
    if len(ev) == 1:
        e = ev.iloc[0]
        expect(e["recipient_contig"] == "asm1_host1" and e["timepoint"] == 2,
               f"host1 is the recipient at t2 ({e['recipient_contig']}, t{e['timepoint']})")
        expect(e["direction"] == "directed" and e["tier"] == "A",
               f"directed, tier A ({e['direction']}, {e['tier']})")
        expect(e["donor_contig"] in ("asm2_host2",) and e["donor_last_timepoint"] == 1,
               f"donor host2 last seen at t1 ({e['donor_contig']}, t{e['donor_last_timepoint']})")
        expect(e["persistence"] == "yes", f"still present at t3 ({e['persistence']})")
    return failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "sim"))
    ap.add_argument("--use-conda", action="store_true")
    ap.add_argument("--keep", action="store_true", help="reuse an existing outdir (no re-simulation)")
    args, extra = ap.parse_known_args()
    if not args.keep:
        shutil.rmtree(args.outdir, ignore_errors=True)
        os.makedirs(args.outdir)
        run(args.outdir, args.use_conda, extra)
    failures = check(args.outdir)
    print(f"\n{len(failures)} failed" if failures else "\nall checks passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
