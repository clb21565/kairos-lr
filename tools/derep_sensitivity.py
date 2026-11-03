#!/usr/bin/env python3
"""
Does dereplication merge a contig with an insert into its sister without it?

kairos needs both sisters in the catalogue: the contig carrying a mobile segment (insert) and the
contig with the empty site. This runs tools/build_catalogue.sh on synthetic sister pairs and
reports, for each clustering setting, which pairs collapse into one representative.

Pairs: a random backbone of length L (the empty contig) and the same backbone with an insert of
s bp (the carrier). The insert sits in the middle, or near a contig end (`edge`, 100 bp from the
end; assemblies often break at mobile elements). Optionally the carrier's backbone carries a
fraction of point substitutions (`div`, sister ANI = 1 - div). Pairs are independent random
sequence, so nothing clusters across pairs.

  python tools/derep_sensitivity.py --out derep_sens [--reps 3] [--threads 2]
  python tools/derep_sensitivity.py --out derep_sens --backbones real_contigs.fa   # HPC

Writes <out>/pairs.tsv (one row per pair x setting: merged yes/no) and <out>/summary.tsv
(per setting, length, position and divergence: the largest insert that still merges).
Needs seqkit and mmseqs on PATH (same as build_catalogue.sh).
"""

import argparse
import csv
import itertools
import os
import random
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
LENGTHS = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
INSERTS = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
POSITIONS = ["middle", "edge"]
DIVS = [0.0, 0.005]
# (MIN_ID, MIN_COV, COV_MODE): 0 = both sequences, 1 = target only (containment)
SETTINGS = [(0.99, c, 0) for c in (0.8, 0.9, 0.95, 0.99)] + [(0.99, 0.95, 1), (0.95, 0.95, 0)]


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng, seq, div):
    if div <= 0:
        return seq
    s = list(seq)
    for i in rng.sample(range(len(s)), int(round(len(s) * div))):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs.append("".join(buf))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
    if name:
        seqs.append("".join(buf))
    return seqs


def make_pairs(rng, reps, backbones):
    """Yield (pair_id, L, s, position, div, empty_seq, carrier_seq)."""
    pool = sorted(backbones or [], key=len)
    n = 0
    for L, s, pos, div, r in itertools.product(LENGTHS, INSERTS, POSITIONS, DIVS, range(reps)):
        if pool:
            longer = [x for x in pool if len(x) >= L]
            if not longer:
                continue
            src = rng.choice(longer)
            st = rng.randrange(len(src) - L + 1)
            back = src[st:st + L]
        else:
            back = rand_seq(rng, L)
        at = L // 2 if pos == "middle" else L - 100
        carrier = mutate(rng, back, div)
        carrier = carrier[:at] + rand_seq(rng, s) + carrier[at:]
        n += 1
        yield f"p{n:05d}", L, s, pos, div, back, carrier


def run_setting(out, fasta_list, setting, threads):
    min_id, min_cov, mode = setting
    d = os.path.join(out, f"id{min_id}_cov{min_cov}_mode{mode}")
    env = dict(os.environ, MIN_ID=str(min_id), MIN_COV=str(min_cov), COV_MODE=str(mode))
    subprocess.run(["bash", os.path.join(HERE, "build_catalogue.sh"), fasta_list, d, "1", str(threads)],
                   env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    reps = {}
    with open(os.path.join(d, "clusters.tsv")) as fh:
        for rep, member in csv.reader(fh, delimiter="\t"):
            reps[member] = rep
    return reps


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--backbones", help="FASTA of real contigs to cut backbones from (default: random)")
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)
    backbones = read_fasta(args.backbones) if args.backbones else None
    fa = os.path.join(args.out, "pairs.fa")
    meta = {}
    with open(fa, "w") as fh:
        for pid, L, s, pos, div, empty, carrier in make_pairs(rng, args.reps, backbones):
            fh.write(f">{pid}_empty\n{empty}\n>{pid}_carrier\n{carrier}\n")
            meta[pid] = (L, s, pos, div)
    lst = os.path.join(args.out, "list.txt")
    with open(lst, "w") as fh:
        fh.write(fa + "\n")
    print(f"{len(meta)} pairs", file=sys.stderr)

    rows = []
    for setting in SETTINGS:
        reps = run_setting(args.out, lst, setting, args.threads)
        for pid, (L, s, pos, div) in meta.items():
            a, b = reps.get(f"pairs__{pid}_empty"), reps.get(f"pairs__{pid}_carrier")
            rows.append({"min_id": setting[0], "min_cov": setting[1], "cov_mode": setting[2],
                         "pair": pid, "length": L, "insert": s, "position": pos, "div": div,
                         "merged": "yes" if a and a == b else "no"})
        print(f"done {setting}", file=sys.stderr)

    cols = list(rows[0])
    with open(os.path.join(args.out, "pairs.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # per setting x length x position x div: fraction merged per insert size
    frac = defaultdict(lambda: [0, 0])
    for r in rows:
        k = (r["min_id"], r["min_cov"], r["cov_mode"], r["length"], r["position"], r["div"], r["insert"])
        frac[k][0] += r["merged"] == "yes"
        frac[k][1] += 1
    with open(os.path.join(args.out, "summary.tsv"), "w") as fh:
        fh.write("min_id\tmin_cov\tcov_mode\tlength\tposition\tdiv\tinsert\tmerged_frac\n")
        for k in sorted(frac):
            m, n = frac[k]
            fh.write("\t".join(map(str, k)) + f"\t{m / n:.2f}\n")


if __name__ == "__main__":
    main()
