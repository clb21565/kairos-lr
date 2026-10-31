"""
Local integration test: simulate a small study, run detect (up to the features) and assess,
and check the counts against what was simulated (see make_sim_data.py).

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
    cmd = ["snakemake", "assess", "-s", os.path.join(REPO, "Snakefile"), "-d", outdir,
           "--configfile", os.path.join(REPO, "config", "config.yaml"),
           "--config", "catalogue=catalogue.fa", "samples=samples.tsv", "gtdb_db=gtdb",
           "ripdb=", "card=", "threads=2", "cluster_method=easy-linclust",
           "-c", "2", "--quiet", "progress", *extra]
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
