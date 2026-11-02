#!/usr/bin/env python3
"""
Tier 3: short-read aligner benchmark for the assess stage (docs/testing.md).

Which short-read aligner should assess use by default, and which junction MAPQ threshold?
This script builds a catalogue with known transfers, simulates paired reads from strains that
differ from the catalogue by a set ANI, maps them with each aligner, counts features with
scripts/count_features.py (the same code the workflow runs) and scores the counts and the
present/absent calls against the truth.

Families. Each family takes three pieces of sequence (from a real catalogue, or random with
--synthetic): a recipient backbone A, a donor backbone B and a cassette C. The catalogue holds

  F_R  = A[:p] + C + A[p:]       recipient as assembled (cassette in the middle)
  F_D  = B[:q] + C + B[q:]       donor carrying the same cassette
  F_S  = A mutated to --sister-ani, no cassette   (sister families only)

so C is a shared segment with two breakpoints and an empty site on F_R and on F_D. Features
are written straight from these coordinates (build_features.py format), so detection does not
enter the comparison.

Each family has one state, which sets the strains in every sample:

  carrier         R strain (backbone at sample ANI, cassette identical) at d, donor at d
  empty           R strain without the cassette at d, donor at d
  mixed           R strain with and without the cassette at d/2 each, donor at d
  sister_carrier  as carrier, with F_S in the catalogue (competes for the flank reads)
  sister_empty    sister strain (F_S at sample ANI) at d, donor at d; the R lineage is absent

With --background-mb, unused catalogue contigs are added as decoys: they make the mapping
competitive as in a real run and give the index a realistic size for the memory figures.

Samples: one per (ANI, depth). The strain backbones are F_R / F_S / F_D mutated to the sample
ANI (substitutions, ~10% small indels); 150 bp pairs (insert 400 +/- 40) with wgsim.

Aligners: strobealign, bwa-mem2, bowtie2 (end-to-end), bowtie2-local. Index and mapping
time and peak memory are measured per aligner (mapping output goes through samtools view
only, so the peak is the aligner's).

Scoring, per recipient segment (both breakpoints summed, as detect_state.py does):
  n, n_span       junction reads (span + bridge) and spanning reads
  lam_true        expected spanning reads from the true carrier depth: c (L - 2a) / L
  lam_obs         the same from the measured flank depth (what the workflow uses)
  call            present (n >= k), absent (n < k and P(N < k | lam_obs) <= alpha), undetermined
  expected call   present for carrier / mixed / sister_carrier / donor, absent for empty;
                  sister_empty: anything but present
and per empty site the empty-site templates against the non-carrier depth.

Outputs (in --outdir/results): segments.tsv (one row per aligner, MAPQ, sample, family, role),
summary.tsv (grouped), resources.tsv. Finished steps are skipped when rerun.

Local smoke test (tiny):
  python tools/aligner_benchmark.py --outdir .tests/bench --synthetic --families-per-state 2 \
      --backbone-len 20000 --anis 99.5,100 --depths 5,20 --threads 2
HPC (real contigs; see slurm/aligner_benchmark.sbatch):
  python tools/aligner_benchmark.py --outdir bench --catalogue catalogue.fa \
      --families-per-state 40 --background-mb 500 --threads 32
"""

import argparse
import math
import os
import random
import subprocess
import sys
import time

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COUNT = os.path.join(REPO, "scripts", "count_features.py")
STATES = ["carrier", "empty", "mixed", "sister_carrier", "sister_empty"]
CARRIER_SHARE = {"carrier": 1.0, "empty": 0.0, "mixed": 0.5, "sister_carrier": 1.0, "sister_empty": 0.0}
LINEAGE_DEPTH = {"carrier": 1.0, "empty": 1.0, "mixed": 1.0, "sister_carrier": 1.0, "sister_empty": 0.0}
ALIGNERS = ["strobealign", "bwa-mem2", "bowtie2", "bowtie2-local"]
READ_LEN = 150
FLANK_LEN = 1000
JUNCTION_FLANK = 500


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--outdir", required=True)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--catalogue", help="real catalogue FASTA to take backbones and cassettes from")
    src.add_argument("--synthetic", action="store_true", help="random sequence instead")
    p.add_argument("--families-per-state", type=int, default=10)
    p.add_argument("--backbone-len", type=int, default=40000)
    p.add_argument("--cassette-len", type=int, default=5000)
    p.add_argument("--max-n-frac", type=float, default=0.001, help="skip windows with more Ns")
    p.add_argument("--background-mb", type=float, default=0,
                   help="add this many Mbp of unused catalogue contigs as decoys (with --catalogue); "
                        "no reads come from them")
    p.add_argument("--sister-ani", type=float, default=98.5,
                   help="ANI of the sister contig to the recipient backbone (catalogue)")
    p.add_argument("--anis", default="99.0,99.5,99.8,100", help="sample strain ANI to the catalogue")
    p.add_argument("--depths", default="2,5,10,30", help="per-strain depth d")
    p.add_argument("--aligners", default=",".join(ALIGNERS))
    p.add_argument("--junction-mapqs", default="0,1,5,10,20",
                   help="count_features is run once per junction MAPQ threshold")
    p.add_argument("--overhang", type=int, default=25)
    p.add_argument("--min-mapq", type=int, default=20, help="flank depth MAPQ")
    p.add_argument("--min-reads", type=int, default=2, help="k: junction reads to call present")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--no-second-pass", action="store_true")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--keep-bams", action="store_true")
    return p.parse_args(argv)


# ── sequence ─────────────────────────────────────────────────────────────────

def random_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, ani, rng, indel_frac=0.1, max_indel=3):
    """Substitutions and 1-3 bp indels at rate (1 - ani/100) per base."""
    if ani >= 100:
        return seq
    n = round(len(seq) * (1 - ani / 100))
    out, last = [], 0
    for pos in sorted(rng.sample(range(len(seq)), n)):
        if pos < last:
            continue
        out.append(seq[last:pos])
        r = rng.random()
        if r < indel_frac / 2:                               # deletion
            last = pos + rng.randint(1, max_indel)
            continue
        if r < indel_frac:                                   # insertion after pos
            out.append(seq[pos] + random_seq(rng.randint(1, max_indel), rng))
        else:
            out.append(rng.choice([b for b in "ACGT" if b != seq[pos].upper()]))
        last = pos + 1
    out.append(seq[last:])
    return "".join(out)


def write_fasta(path, seqs):
    with open(path, "w") as out:
        for name, seq in seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


def read_fasta(path):
    seqs, name, parts = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name:
        seqs[name] = "".join(parts)
    return seqs


class Source:
    """Windows of real catalogue contigs (each contig used once), or random sequence."""

    def __init__(self, args, rng):
        self.rng, self.args = rng, args
        if args.synthetic:
            self.fa = None
            return
        import pysam
        if not os.path.exists(args.catalogue + ".fai"):
            pysam.faidx(args.catalogue)
        self.fa = pysam.FastaFile(args.catalogue)
        self.pool = [(c, n) for c, n in zip(self.fa.references, self.fa.lengths)
                     if n >= args.cassette_len]
        rng.shuffle(self.pool)
        self.used = set()

    def take(self, length):
        if self.fa is None:
            return random_seq(length, self.rng), "random"
        for contig, n in self.pool:
            if contig in self.used or n < length:
                continue
            self.used.add(contig)
            for _ in range(5):
                s = self.rng.randint(0, n - length)
                seq = self.fa.fetch(contig, s, s + length).upper()
                if seq.count("N") <= self.args.max_n_frac * length:
                    return seq, f"{contig}:{s}-{s + length}"
        sys.exit(f"catalogue has too few contigs >= {length} bp without Ns")

    def background(self, mbp):
        """Whole unused contigs (>= 1 kb) up to mbp Mbp."""
        out, total = {}, 0
        for contig, n in self.pool:
            if total >= mbp * 1e6:
                break
            if contig not in self.used and n >= 1000:
                self.used.add(contig)
                out[f"BG_{contig}"] = self.fa.fetch(contig).upper()
                total += n
        return out


# ── step 1: catalogue, features, truth ───────────────────────────────────────

def build_catalogue(args, d):
    cat_fa = os.path.join(d, "catalogue.fa")
    if os.path.exists(os.path.join(d, "families.tsv")):
        return
    rng = random.Random(args.seed)
    src = Source(args, rng)
    cat, fams, bps, ess, junc = {}, [], [], [], {}

    def add_segment(contig, seq, s0, s1):
        es_id = f"ES{len(ess) + 1:07d}"
        ids = []
        for pos, side in ((s0, "L"), (s1, "R")):
            bp_id = f"BP{len(bps) + 1:07d}"
            f0, f1 = (pos - FLANK_LEN, pos) if side == "L" else (pos, pos + FLANK_LEN)
            w0, w1 = pos - JUNCTION_FLANK, pos + JUNCTION_FLANK
            bps.append(dict(bp_id=bp_id, contig=contig, pos=pos, side=side, window_start=w0,
                            window_end=w1, flank_bp=FLANK_LEN, _f0=f0, _f1=f1))
            junc[bp_id] = seq[w0:w1]
            ids.append(bp_id)
        ess.append(dict(es_id=es_id, contig=contig, shared_start=s0, shared_end=s1,
                        left_start=s0 - JUNCTION_FLANK, right_end=s1 + JUNCTION_FLANK,
                        junction_pos=JUNCTION_FLANK))
        junc[es_id] = seq[s0 - JUNCTION_FLANK:s0] + seq[s1:s1 + JUNCTION_FLANK]
        return ids, es_id

    Lb, Lc = args.backbone_len, args.cassette_len
    i = 0
    for state in STATES:
        for _ in range(args.families_per_state):
            i += 1
            fam = f"F{i:03d}"
            A, a_src = src.take(Lb)
            B, b_src = src.take(Lb)
            C, c_src = src.take(Lc)
            p, q = Lb // 2, Lb // 3
            R = A[:p] + C + A[p:]
            D = B[:q] + C + B[q:]
            cat[f"{fam}_R"], cat[f"{fam}_D"] = R, D
            if state.startswith("sister"):
                S = mutate(A, args.sister_ani, random.Random(f"{args.seed}-{fam}-S"))
                cat[f"{fam}_S"] = S
            for role, contig, seq, s0 in (("R", f"{fam}_R", R, p), ("D", f"{fam}_D", D, q)):
                (bl, br), es = add_segment(contig, seq, s0, s0 + Lc)
                fams.append(dict(family=fam, state=state, role=role, contig=contig, bp_left=bl,
                                 bp_right=br, es_id=es, contig_s0=s0, contig_s1=s0 + Lc,
                                 src_backbone=a_src if role == "R" else b_src,
                                 src_cassette=c_src))
    if args.background_mb and not args.synthetic:
        cat.update(src.background(args.background_mb))
    os.makedirs(os.path.join(d, "features"), exist_ok=True)
    write_fasta(cat_fa, cat)
    bdf = pd.DataFrame(bps)
    bdf.drop(columns=["_f0", "_f1"]).to_csv(os.path.join(d, "features", "breakpoints.tsv"),
                                            sep="\t", index=False)
    bdf[["contig", "_f0", "_f1", "bp_id"]].to_csv(os.path.join(d, "features", "flanks.bed"),
                                                  sep="\t", index=False, header=False)
    pd.DataFrame(ess).to_csv(os.path.join(d, "features", "empty_sites.tsv"), sep="\t", index=False)
    write_fasta(os.path.join(d, "features", "junctions.fa"), junc)
    pd.DataFrame(fams).to_csv(os.path.join(d, "families.tsv"), sep="\t", index=False)
    print(f"catalogue: {len(cat)} contigs, {sum(map(len, cat.values())) / 1e6:.2f} Mbp; "
          f"{len(bps)} breakpoints, {len(ess)} empty sites", file=sys.stderr)


# ── step 2: strains and reads ────────────────────────────────────────────────

def sample_name(ani, depth):
    return f"ani{ani:g}_d{depth:g}"


def strains(args, d, ani):
    """Strain genomes at this ANI, mutated once per family (shared by all depths)."""
    path = os.path.join(d, "strains", f"ani{ani:g}.fa")
    if os.path.exists(path):
        return read_fasta(path)
    cat = read_fasta(os.path.join(d, "catalogue.fa"))
    fams = pd.read_csv(os.path.join(d, "families.tsv"), sep="\t")
    out = {}
    for f in fams.itertuples():
        rng = random.Random(f"{args.seed}-{f.family}-{f.role}-{ani}")
        seq = cat[f.contig]
        s0, s1 = int(f.contig_s0), int(f.contig_s1)
        left, cas, right = mutate(seq[:s0], ani, rng), seq[s0:s1], mutate(seq[s1:], ani, rng)
        if f.role == "D":
            out[f"{f.family}_D"] = left + cas + right
            continue
        out[f"{f.family}_Rcar"] = left + cas + right
        out[f"{f.family}_Remp"] = left + right
        if f.state == "sister_empty":
            out[f"{f.family}_Sx"] = mutate(cat[f"{f.family}_S"], ani, rng)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    write_fasta(path, out)
    return out


def sample_mix(fams, depth):
    """genome -> depth for one sample."""
    mix = {}
    for f in fams.itertuples():
        if f.role == "D":
            mix[f"{f.family}_D"] = depth
        elif f.state in ("carrier", "sister_carrier"):
            mix[f"{f.family}_Rcar"] = depth
        elif f.state == "empty":
            mix[f"{f.family}_Remp"] = depth
        elif f.state == "mixed":
            mix[f"{f.family}_Rcar"] = depth / 2
            mix[f"{f.family}_Remp"] = depth / 2
        elif f.state == "sister_empty":
            mix[f"{f.family}_Sx"] = depth
    return mix


def simulate_reads(args, d, ani, depth, seed):
    name = sample_name(ani, depth)
    rdir = os.path.join(d, "reads")
    r1, r2 = os.path.join(rdir, f"{name}_R1.fq"), os.path.join(rdir, f"{name}_R2.fq")
    if os.path.exists(r2):
        return r1, r2
    os.makedirs(rdir, exist_ok=True)
    genomes = strains(args, d, ani)
    fams = pd.read_csv(os.path.join(d, "families.tsv"), sep="\t")
    by_depth = {}
    for g, dep in sample_mix(fams, depth).items():
        by_depth.setdefault(dep, {})[g] = genomes[g]
    for k, (dep, gs) in enumerate(sorted(by_depth.items())):
        # wgsim spreads pairs over the sequences in proportion to their length
        tmp = os.path.join(rdir, f"{name}.{k}")
        write_fasta(tmp + ".fa", gs)
        n_pairs = int(dep * sum(map(len, gs.values())) / (2 * READ_LEN))
        subprocess.run(["wgsim", "-S", str(seed + k), "-N", str(n_pairs), "-1", str(READ_LEN),
                        "-2", str(READ_LEN), "-d", "400", "-s", "40", "-e", "0.002", "-r", "0",
                        "-R", "0", "-A", "1", tmp + ".fa", tmp + "_1.fq", tmp + "_2.fq"],
                       check=True, capture_output=True)
    for mate, out in (("1", r1), ("2", r2)):
        with open(out + ".part", "w") as o:
            for k in range(len(by_depth)):
                part = os.path.join(rdir, f"{name}.{k}_{mate}.fq")
                with open(part) as f:
                    for line in f:
                        o.write(line)
                os.remove(part)
        os.replace(out + ".part", out)
    for k in range(len(by_depth)):
        os.remove(os.path.join(rdir, f"{name}.{k}.fa"))
    return r1, r2


# ── step 3: index, align, count ──────────────────────────────────────────────

def timed(cmd, log):
    """Run a shell command; return (wall s, peak RSS MB of the largest process)."""
    t0 = time.time()
    with open(log, "a") as lf:
        lf.write(f"$ {cmd}\n")
        lf.flush()
        p = subprocess.Popen(cmd, shell=True, stdout=lf, stderr=lf, executable="/bin/bash")
        _, status, ru = os.wait4(p.pid, 0)
    p.returncode = os.waitstatus_to_exitcode(status)
    if p.returncode != 0:
        sys.exit(f"failed ({p.returncode}): {cmd}\nsee {log}")
    return time.time() - t0, ru.ru_maxrss / 1024           # ru_maxrss is KB on Linux


def index(aligner, d, threads, log):
    idir = os.path.join(d, "index", aligner.replace("-local", ""))
    ref = os.path.join(idir, "catalogue.fa")
    done = ref + ".done"
    if not os.path.exists(done):
        os.makedirs(idir, exist_ok=True)
        if not os.path.exists(ref):
            os.symlink(os.path.abspath(os.path.join(d, "catalogue.fa")), ref)
        cmd = {"strobealign": f"strobealign --create-index -t {threads} -r {READ_LEN} {ref}",
               "bwa-mem2": f"bwa-mem2 index {ref}",
               "bowtie2": f"bowtie2-build --threads {threads} {ref} {ref}"}[aligner.replace("-local", "")]
        wall, rss = timed(cmd, log)
        open(done, "w").write(f"{wall:.1f}\t{rss:.0f}\n")
    wall, rss = open(done).read().split()
    return ref, float(wall), float(rss)


def align_cmd(aligner, ref, r1, r2, threads):
    if aligner == "strobealign":
        return f"strobealign -t {threads} --use-index -r {READ_LEN} {ref} {r1} {r2}"
    if aligner == "bwa-mem2":
        return f"bwa-mem2 mem -t {threads} {ref} {r1} {r2}"
    local = " --local" if aligner == "bowtie2-local" else ""
    return f"bowtie2{local} -p {threads} -X 1000 -x {ref} -1 {r1} -2 {r2}"


def align(args, d, aligner, ref, name, r1, r2, log):
    adir = os.path.join(d, "align", aligner)
    bam = os.path.join(adir, f"{name}.bam")
    done = bam + ".done"
    if not os.path.exists(done):
        os.makedirs(adir, exist_ok=True)
        raw = bam + ".unsorted.bam"
        wall, rss = timed(f"{align_cmd(aligner, ref, r1, r2, args.threads)} | samtools view -b -1 -o {raw} -", log)
        subprocess.run(["samtools", "sort", "-@", str(min(4, args.threads)), "-m", "512M",
                        "-o", bam, raw], check=True, capture_output=True)
        subprocess.run(["samtools", "index", bam], check=True)
        os.remove(raw)
        open(done, "w").write(f"{wall:.1f}\t{rss:.0f}\n")
    wall, rss = open(done).read().split()
    return bam, float(wall), float(rss)


def count(args, d, aligner, name, bam, jmq):
    cdir = os.path.join(d, "counts", aligner, f"mapq{jmq}")
    out_bp = os.path.join(cdir, f"{name}.breakpoint_counts.tsv")
    out_es = os.path.join(cdir, f"{name}.empty_site_counts.tsv")
    if not os.path.exists(out_es):
        os.makedirs(cdir, exist_ok=True)
        feat = os.path.join(d, "features")
        cmd = [sys.executable, COUNT, "--bam", bam, "--sample", name, "--read-type", "short_pe",
               "--breakpoints", f"{feat}/breakpoints.tsv", "--flanks", f"{feat}/flanks.bed",
               "--empty-sites", f"{feat}/empty_sites.tsv", "--overhang", str(args.overhang),
               "--min-mapq", str(args.min_mapq), "--junction-min-mapq", str(jmq),
               "--out-breakpoints", out_bp + ".part", "--out-empty-sites", out_es + ".part"]
        if not args.no_second_pass:
            cmd += ["--junctions", f"{feat}/junctions.fa"]
        subprocess.run(cmd, check=True, capture_output=True)
        os.replace(out_bp + ".part", out_bp)
        os.replace(out_es + ".part", out_es)
    return out_bp, out_es


# ── step 4: score ────────────────────────────────────────────────────────────

def p_below(k, lam):
    if lam <= 0:
        return 1.0
    return sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k))


def score_sample(args, fams, bp_path, es_path, aligner, jmq, ani, depth):
    bp = pd.read_csv(bp_path, sep="\t").set_index("bp_id")
    es = pd.read_csv(es_path, sep="\t").set_index("es_id")
    frac = (READ_LEN - 2 * args.overhang) / READ_LEN
    rows = []
    for f in fams.itertuples():
        b = bp.loc[[f.bp_left, f.bp_right]]
        share = 1.0 if f.role == "D" else CARRIER_SHARE[f.state]
        lineage = 1.0 if f.role == "D" else LINEAGE_DEPTH[f.state]
        n, n_span = int((b["span"] + b["bridge"]).sum()), int(b["span"].sum())
        lam_true = 2 * share * depth * frac
        lam_obs = float((b["flank_depth"] * frac).sum())
        if n >= args.min_reads:
            call = "present"
        elif p_below(args.min_reads, lam_obs) <= args.alpha:
            call = "absent"
        else:
            call = "undetermined"
        if f.role == "D" or share > 0:
            expected = "present"
        elif f.state == "empty":
            expected = "absent"
        else:
            expected = "not_present"          # sister_empty: R lineage not in the sample
        if call == expected or (expected == "not_present" and call != "present"):
            outcome = "correct"
        elif call == "undetermined":
            outcome = "undetermined"
        else:
            outcome = "wrong"
        empty_depth = (1 - share) * depth * lineage
        rows.append(dict(
            aligner=aligner, junction_mapq=jmq, ani=ani, depth=depth, family=f.family,
            state=f.state, role=f.role, n=n, n_span=n_span, lam_true=round(lam_true, 2),
            lam_obs=round(lam_obs, 2), flank_depth=round(float(b["flank_depth"].mean()), 2),
            true_flank_depth=lineage * depth, empty=int(es.loc[f.es_id, "empty"]),
            empty_split=int(es.loc[f.es_id, "empty_split"]),
            empty_pairs=int(es.loc[f.es_id, "empty_pairs"]),
            empty_remap=int(es.loc[f.es_id, "empty_remap"]),
            lam_empty=round(empty_depth * frac, 2), call=call, expected=expected, outcome=outcome))
    return rows


def summarise(seg):
    g = seg.groupby(["aligner", "junction_mapq", "ani", "depth", "state", "role"])
    s = g.agg(families=("family", "size"), n=("n", "mean"), n_span=("n_span", "mean"),
              lam_true=("lam_true", "mean"), lam_obs=("lam_obs", "mean"),
              flank_depth=("flank_depth", "mean"), true_flank_depth=("true_flank_depth", "mean"),
              empty=("empty", "mean"), lam_empty=("lam_empty", "mean"))
    s["span_ratio"] = (s["n_span"] / s["lam_true"]).where(s["lam_true"] > 0)
    s["flank_ratio"] = (s["flank_depth"] / s["true_flank_depth"]).where(s["true_flank_depth"] > 0)
    s["empty_ratio"] = (s["empty"] / s["lam_empty"]).where(s["lam_empty"] > 0)
    for o in ("correct", "undetermined", "wrong"):
        s[f"frac_{o}"] = g["outcome"].apply(lambda x, o=o: (x == o).mean())
    s["false_junction_reads"] = g.apply(
        lambda x: x.loc[x["expected"] != "present", "n"].sum(), include_groups=False)
    return s.round(3).reset_index()


def headline(summary, resources, jmq):
    """Short table: per aligner, call outcomes and ratios over all ANI/depths at one MAPQ."""
    s = summary[(summary["junction_mapq"] == jmq) & (summary["role"] == "R")]
    t = s.groupby("aligner").agg(
        correct=("frac_correct", "mean"), undetermined=("frac_undetermined", "mean"),
        wrong=("frac_wrong", "mean"), false_junction_reads=("false_junction_reads", "sum"))
    car = s[s["state"].isin(["carrier", "sister_carrier"])].groupby("aligner")
    t["span_ratio_carrier"] = car["span_ratio"].mean()
    t["flank_ratio_sister"] = s[s["state"] == "sister_carrier"].groupby("aligner")["flank_ratio"].mean()
    t["empty_ratio"] = s[s["state"] == "empty"].groupby("aligner")["empty_ratio"].mean()
    t = t.join(resources.set_index("aligner")[["index_s", "index_mb", "align_s_total", "align_mb_max"]])
    return t.round(3)


def main(argv=None):
    args = parse_args(argv)
    d = args.outdir
    os.makedirs(d, exist_ok=True)
    log = os.path.join(d, "commands.log")
    anis = [float(x) for x in args.anis.split(",")]
    depths = [float(x) for x in args.depths.split(",")]
    jmqs = [int(x) for x in args.junction_mapqs.split(",")]
    aligners = args.aligners.split(",")
    for a in aligners:
        if a not in ALIGNERS:
            sys.exit(f"unknown aligner {a}; choose from {ALIGNERS}")

    build_catalogue(args, d)
    fams = pd.read_csv(os.path.join(d, "families.tsv"), sep="\t")

    samples = []
    for i, ani in enumerate(anis):
        for j, depth in enumerate(depths):
            r1, r2 = simulate_reads(args, d, ani, depth, seed=args.seed + 100 * i + 10 * j)
            samples.append((ani, depth, sample_name(ani, depth), r1, r2))
    print(f"{len(samples)} samples simulated", file=sys.stderr)

    seg_rows, res_rows = [], []
    for aligner in aligners:
        ref, iw, irss = index(aligner, d, args.threads, log)
        a_wall, a_rss = 0.0, 0.0
        for ani, depth, name, r1, r2 in samples:
            bam, w, rss = align(args, d, aligner, ref, name, r1, r2, log)
            a_wall, a_rss = a_wall + w, max(a_rss, rss)
            for jmq in jmqs:
                bp_path, es_path = count(args, d, aligner, name, bam, jmq)
                seg_rows += score_sample(args, fams, bp_path, es_path, aligner, jmq, ani, depth)
            if not args.keep_bams:
                for ext in ("", ".bai"):
                    os.remove(bam + ext)
            print(f"{aligner} {name}: aligned in {w:.0f} s", file=sys.stderr)
        res_rows.append(dict(aligner=aligner, index_s=iw, index_mb=irss, align_s_total=round(a_wall, 1),
                             align_mb_max=a_rss, samples=len(samples)))

    out = os.path.join(d, "results")
    os.makedirs(out, exist_ok=True)
    seg = pd.DataFrame(seg_rows)
    seg.to_csv(os.path.join(out, "segments.tsv"), sep="\t", index=False)
    summary = summarise(seg)
    summary.to_csv(os.path.join(out, "summary.tsv"), sep="\t", index=False)
    resources = pd.DataFrame(res_rows)
    resources.to_csv(os.path.join(out, "resources.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    for jmq in jmqs:
        print(f"\n== recipient segments, junction MAPQ >= {jmq} (mean over ANI x depth) ==")
        print(headline(summary, resources, jmq).to_string())
    print(f"\nwritten to {out}/")


if __name__ == "__main__":
    main()
