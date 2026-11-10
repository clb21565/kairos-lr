"""
Local end-to-end test of the spike-in benchmark: synthetic genomes and elements -> design ->
build -> wgsim reads -> kairos (detect features, assess, call) on the oracle catalogue ->
insitu.py with the planted lineages -> score against truth.tsv.

One event per scenario (two for insitu: 10x and 2x), each with its own element, 4 timepoints, compartments SYS
(observed) and INF (prior only), transfer at the third timepoint. No background reads.

Needs the same tools as run_sim.py on PATH (snakemake, prodigal, mmseqs2, strobealign,
minimap2, samtools, wgsim, seqkit, mappy, pysam, pandas).

Usage: python .tests/run_spikein.py [--outdir .tests/spikein] [--keep]
"""

import argparse
import os
import random
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SPIKE = os.path.join(REPO, "tools", "spikein.py")
rng = random.Random(11)
STOPS = {"TAA", "TAG", "TGA"}
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def gene(n):
    return "ATG" + "".join(rng.choice(CODONS) for _ in range(n)) + "TAA"


def block(n_genes):
    sp = lambda: "".join(rng.choice("ACGT") for _ in range(150))
    return "".join(sp() + gene(rng.randint(150, 350)) for _ in range(n_genes)) + sp()


def fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n{seq}\n")


def py(*args, **kw):
    return subprocess.run([sys.executable, *args], check=True, **kw)


def make_pools(d):
    os.makedirs(os.path.join(d, "pool"), exist_ok=True)
    rows = []
    for i in range(16):                        # 16 genomes, each its own phylum
        lin = f"d__Bacteria;p__P{i};c__C{i};o__O{i};f__F{i};g__G{i};s__G{i} sp"
        rows.append((f"gen{i:02d}", lin))
    same = "d__Bacteria;p__PX;c__CX;o__OX;f__FX;g__GX;s__GX sp"
    rows += [("genX1", same), ("genX2", same)]  # two strains of one species
    g = []
    for gid, lin in rows:
        p = os.path.join(d, "pool", f"{gid}.fa")
        fasta(p, gid, block(45))                # ~40 kb
        g.append({"genome_id": gid, "fasta": p, "lineage": lin})
    pd.DataFrame(g).to_csv(os.path.join(d, "genomes.tsv"), sep="\t", index=False)
    e = []
    els = [(f"ISsim{i}", "is", 8, 3) for i in range(4)] + \
          [(f"Tnsim{i}", "transposon", 5, 5) for i in range(4)] + [("pSim1", "plasmid", 0, 5)]
    for eid, typ, tsd, n in els:
        p = os.path.join(d, "pool", f"{eid}.fa")
        fasta(p, eid, block(n))
        e.append({"element_id": eid, "fasta": p, "type": typ, "tsd": tsd})
    pd.DataFrame(e).to_csv(os.path.join(d, "elements.tsv"), sep="\t", index=False)


def run(d):
    make_pools(d)
    py(SPIKE, "design", "--genomes", f"{d}/genomes.tsv", "--elements", f"{d}/elements.tsv",
       "--out", f"{d}/events.tsv", "--ranks", "p", "--depths", "10,2", "--carrier-fracs", "0.5")
    b = f"{d}/build"
    py(SPIKE, "build", "--events", f"{d}/events.tsv", "--genomes", f"{d}/genomes.tsv",
       "--elements", f"{d}/elements.tsv", "--out-dir", b, "--timepoints", "1,2,3,4",
       "--transfer-index", "2", "--margin", "3000")
    for s in pd.read_csv(f"{b}/samples.tsv", sep="\t")["sample"]:
        py(SPIKE, "reads", "--build-dir", b, "--sample", s, stdout=subprocess.DEVNULL)

    open(f"{b}/gtdb", "w").close()
    with open(f"{b}/sim_config.yaml", "w") as f:
        f.write("call:\n  use_taxonomy: false\n  lineage_filter: none\n"
                "  observed_compartments: [SYS]\n")
    subprocess.run(["snakemake", "assess", "call", "-s", os.path.join(REPO, "Snakefile"), "-d", b,
                    "--configfile", os.path.join(REPO, "config", "config.yaml"), f"{b}/sim_config.yaml",
                    "--config", "catalogue=oracle_catalogue.fa", "samples=samples.tsv",
                    "gtdb_db=gtdb", "ripdb=", "card=", "threads=2", "cluster_method=easy-linclust",
                    "sort_mem_mb=128", "-c", "2", "--quiet", "all"], check=True)
    c = f"{b}/results/call"
    # call again with the planted lineages, as the taxonomy step would provide them
    py(os.path.join(REPO, "scripts", "insitu.py"), "--states", f"{c}/states.tsv",
       "--segments", f"{c}/segments.tsv", "--event-segments", f"{c}/event_segments.tsv",
       "--observed", "SYS", "--contig-taxonomy", f"{b}/oracle_lineages.tsv",
       "--lineage-filter", "discordant", "--out", f"{c}/insitu_events.lineage.tsv",
       "--excluded", f"{c}/insitu_excluded.lineage.tsv")


def check(d):
    b, c = f"{d}/build", f"{d}/build/results/call"
    out = subprocess.run([sys.executable, SPIKE, "score", "--truth", f"{b}/truth.tsv",
                          "--events", f"{c}/insitu_events.lineage.tsv",
                          "--excluded", f"{c}/insitu_excluded.lineage.tsv",
                          "--out", f"{b}/score.tsv"], check=True, capture_output=True, text=True)
    print(out.stdout)
    res = pd.read_csv(f"{b}/score.tsv", sep="\t")
    failures = []

    def expect(ok, msg):
        print(("ok    " if ok else "FAIL  ") + msg)
        if not ok:
            failures.append(msg)
    obs = {(r.scenario, r.depth): r.observed for r in res.itertuples()}
    don = {(r.scenario, r.depth): r.donor_match for r in res.itertuples()}
    # what the current rules must give; rare_bloom is reported, not asserted (see below)
    expect(obs[("insitu", 10)] == "called_A" and don[("insitu", 10)] in ("planted", "same_element"),
           f"insitu 10x called, tier A, carrier of the element as donor ({obs[('insitu', 10)]}, {don[('insitu', 10)]})")
    for sc in ("immigration", "preexisting", "plasmid"):
        expect(obs[(sc, 10)] == "none", f"{sc}: no event ({obs[(sc, 10)]})")
    expect(obs[("donor_late", 10)] == "unresolved", f"donor_late: unresolved ({obs[('donor_late', 10)]})")
    expect(obs[("same_lineage", 10)] == "excluded", f"same_lineage: excluded ({obs[('same_lineage', 10)]})")
    # near the detection limit the outcome depends on the draw of reads: reported, not asserted
    print(f"info  insitu 2x: {obs[('insitu', 2)]}; 50% carriers (5x): {obs[('carrier_frac', 10)]}; "
          f"rare_bloom (carrier at 0.2x before transfer): {obs[('rare_bloom', 10)]}")
    expect("0 reported events with a recipient that was never planted" in out.stdout,
           "no event has a donor or other unplanted contig as its recipient")
    return failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "spikein"))
    ap.add_argument("--keep", action="store_true", help="reuse an existing outdir")
    args = ap.parse_args()
    if not args.keep:
        shutil.rmtree(args.outdir, ignore_errors=True)
        os.makedirs(args.outdir)
        run(args.outdir)
    failures = check(args.outdir)
    print(f"\n{len(failures)} failed" if failures else "\nall checks passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
