"""
Simulated study for local tests: a catalogue, the true genomes behind it, reads per sample
and the expected answer.

Hosts carry their own genes; a ~4 kb gene cassette is shared by host1, host2 and host3
(host3 in reverse orientation). host1 also exists without the cassette ("host1_empty"),
which is the pre-transfer state. The catalogue holds host1-3 with the cassette plus host4.

Samples (group G, compartment X):
  t1_short / t1_long   host1_empty 20x, host2 20x, host4 10x  -> cassette absent from host1
  t2_short / t2_long   host1 10x + host1_empty 10x, host2 20x  -> host1 gains it (50% carriers)
  t3_short             host1 3x, host2 20x                     -> host1 too shallow to call

Usage: python make_sim_data.py <outdir> [--wgsim wgsim]
"""

import argparse
import gzip
import os
import random
import subprocess

rng = random.Random(7)
STOPS = {"TAA", "TAG", "TGA"}
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def gene(n_codons):
    return "ATG" + "".join(rng.choice(CODONS) for _ in range(n_codons)) + "TAA"


def spacer(n=150):
    return "".join(rng.choice("ACGT") for _ in range(n))


def block(n_genes):
    return "".join(spacer() + gene(rng.randint(150, 350)) for _ in range(n_genes)) + spacer()


def write_fasta(path, seqs):
    with open(path, "w") as out:
        for name, seq in seqs.items():
            out.write(f">{name}\n{seq}\n")


def short_reads(genomes, depths, prefix, wgsim, seed):
    """Paired 150 bp reads (insert 400 +/- 40) with wgsim, pooled over genomes at the given depths."""
    r1s, r2s = [], []
    for name, depth in depths.items():
        g = genomes[name]
        n_pairs = int(depth * len(g) / 300)
        tmp = f"{prefix}.{name}"
        write_fasta(tmp + ".fa", {name: g})
        subprocess.run([wgsim, "-S", str(seed), "-N", str(n_pairs), "-1", "150", "-2", "150",
                        "-d", "400", "-s", "40", "-e", "0.002", "-r", "0", "-R", "0",
                        tmp + ".fa", tmp + "_1.fq", tmp + "_2.fq"],
                       check=True, capture_output=True)
        r1s.append(tmp + "_1.fq")
        r2s.append(tmp + "_2.fq")
    for parts, out in ((r1s, prefix + "_R1.fq.gz"), (r2s, prefix + "_R2.fq.gz")):
        with gzip.open(out, "wt") as o:
            for p in parts:
                o.write(open(p).read())
                os.remove(p)
    for name in depths:
        os.remove(f"{prefix}.{name}.fa")
    return [prefix + "_R1.fq.gz", prefix + "_R2.fq.gz"]


def long_reads(genomes, depths, path, err=0.01):
    """Reads of 3-15 kb from either strand with 1% substitutions."""
    with gzip.open(path, "wt") as out:
        i = 0
        for name, depth in depths.items():
            g = genomes[name]
            target, made = depth * len(g), 0
            while made < target:
                n = min(len(g), rng.randint(3000, 15000))
                s = rng.randint(0, len(g) - n)
                seq = list(g[s:s + n])
                for j in range(n):
                    if rng.random() < err:
                        seq[j] = rng.choice([b for b in "ACGT" if b != seq[j]])
                seq = "".join(seq)
                if rng.random() < 0.5:
                    seq = revcomp(seq)
                i += 1
                out.write(f"@{name}_{i}\n{seq}\n+\n{'5' * n}\n")
                made += n
    return [path]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--wgsim", default="wgsim")
    args = ap.parse_args()
    os.makedirs(os.path.join(args.outdir, "reads"), exist_ok=True)

    cassette = block(4)
    h1l, h1r = block(6), block(6)
    genomes = {
        "asm1_host1": h1l + cassette + h1r,
        "asm1_host1_empty": h1l + h1r,
        "asm2_host2": block(5) + cassette + block(7),
        "asm3_host3": block(8) + revcomp(cassette) + block(4),
        "asm1_host4": block(10),
    }
    catalogue = {k: v for k, v in genomes.items() if k != "asm1_host1_empty"}
    write_fasta(os.path.join(args.outdir, "catalogue.fa"), catalogue)
    write_fasta(os.path.join(args.outdir, "truth_genomes.fa"), genomes)

    samples = [
        ("t1_short", "short_pe", 1, {"asm1_host1_empty": 20, "asm2_host2": 20, "asm1_host4": 10}),
        ("t1_long", "ont", 1, {"asm1_host1_empty": 20, "asm2_host2": 20, "asm1_host4": 10}),
        ("t2_short", "short_pe", 2, {"asm1_host1": 10, "asm1_host1_empty": 10, "asm2_host2": 20}),
        ("t2_long", "ont", 2, {"asm1_host1": 10, "asm1_host1_empty": 10, "asm2_host2": 20}),
        ("t3_short", "short_pe", 3, {"asm1_host1": 3, "asm2_host2": 20}),
    ]
    rows = ["sample\treads\tread_type\ttimepoint\tgroup\tcompartment"]
    truth = ["sample\tgenome\tdepth"]
    for k, (name, rt, t, depths) in enumerate(samples):
        prefix = os.path.join(args.outdir, "reads", name)
        if rt == "short_pe":
            files = short_reads(genomes, depths, prefix, args.wgsim, seed=k + 1)
        else:
            files = long_reads(genomes, depths, prefix + ".fq.gz")
        rel = [os.path.relpath(f, args.outdir) for f in files]
        rows.append(f"{name}\t{','.join(rel)}\t{rt}\t{t}\tG\tX")
        truth += [f"{name}\t{g}\t{d}" for g, d in depths.items()]
    open(os.path.join(args.outdir, "samples.tsv"), "w").write("\n".join(rows) + "\n")
    open(os.path.join(args.outdir, "truth.tsv"), "w").write("\n".join(truth) + "\n")
    print(f"cassette {len(cassette)} bp; host1 cassette at {len(h1l)}-{len(h1l) + len(cassette)}")


if __name__ == "__main__":
    main()
