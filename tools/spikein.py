#!/usr/bin/env python3
"""
Planted in situ HGT events for benchmarking kairos (spike-in into real metagenome reads).

Each event has a donor genome carrying a mobile element and a recipient genome that is first
present without it (the empty site) and, from the transfer timepoint on, with it inserted at a
known position. A scenario sets who is present where and when, and the truth table records
what kairos *should* call, including the cases where the right answer is "no event":

  scenario       set-up (S = observed compartment, I = prior-only compartment, k = transfer)  expect
  insitu         donor everywhere from t0; recipient in S, empty before k, carrier from k     called
  carrier_frac   insitu with only a fraction of the recipient population carrying it        called
  immigration    insitu, but the carrier is already in I before k                           none
  preexisting    carrier (and donor) present in S from t0                                   none
  rare_bloom     insitu, but the carrier is in S before k below detection (rare_depth)      not_A
  donor_late     insitu, but the donor only appears at k                                    unresolved
  same_lineage   insitu between two genomes of the same species                             excluded
  plasmid        element is a separate replicon (no junction in the host): undetectable    none

Subcommands
  design   genomes.tsv + elements.tsv -> events.tsv (factorial over scenarios, ranks, depths)
  build    events.tsv -> variants.fa, composition.tsv, truth.tsv, junctions.fa, samples.tsv,
           oracle_catalogue.fa, oracle_lineages.tsv
  reads    one sample of composition.tsv -> simulated reads (wgsim / iss / badread), with the
           background reads of that sample prepended
  score    insitu_events.tsv (+ insitu_excluded.tsv) vs truth.tsv, for runs on the oracle
           catalogue (contig names = planted sequence names)

genomes.tsv:  genome_id  fasta  lineage (GTDB, d__..;p__..;...;s__..)
elements.tsv: element_id  fasta  type (is, transposon, cassette, prophage, ice, plasmid)  tsd (bp)
Only the longest record of each genome FASTA (the chromosome) is used.
"""

import argparse
import gzip
import itertools
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

RANKS = ["d", "p", "c", "o", "f", "g", "s"]
SCENARIOS = ["insitu", "carrier_frac", "immigration", "preexisting", "rare_bloom",
             "donor_late", "same_lineage", "plasmid"]
EXPECT = {"insitu": "called", "carrier_frac": "called", "immigration": "none",
          "preexisting": "none", "rare_bloom": "not_A", "donor_late": "unresolved",
          "same_lineage": "excluded", "plasmid": "none"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def longest(path):
    seqs = read_fasta(path)
    if not seqs:
        sys.exit(f"{path}: no sequences")
    return max(seqs.values(), key=len)


def write_fasta(path, seqs, width=80):
    with open(path, "w") as out:
        for name, s in seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(s), width):
                out.write(s[i:i + width] + "\n")


def parse_lineage(s):
    """GTDB or mmseqs lineage string -> names per rank d..s ('' where missing)."""
    ranked = {}
    for x in ("" if pd.isna(s) else str(s)).split(";"):
        x = x.strip()
        if len(x) > 2 and x[0] in RANKS and x[1] == "_":
            ranked[x[0]] = x[2:].lstrip("_")
    return [ranked.get(r, "") for r in RANKS]


def discordant_rank(a, b):
    """First rank at which both lineages are named and differ; 'same' if equal to species."""
    for r, x, y in zip(RANKS, a, b):
        if not x or not y:
            return ""
        if x != y:
            return r
    return "same"


# ── design ──────────────────────────────────────────────────────────────────────────────

def design(genomes, elements, scenarios, ranks, depths, fracs, base_depth, donor_depth,
           rare_depth, replicates, seed, reuse_elements=False):
    """
    Events table. 'insitu' gets the full rank x depth grid, 'carrier_frac' the carrier
    fractions at base_depth, every other scenario one event per rank at base_depth
    ('same_lineage' uses same-species pairs, 'plasmid' plasmid elements). Each genome is used
    by one event only, so events never share reads, and by default each element too: a
    second carrier of the same element is a genuine earlier donor for kairos, which changes
    the right answer of scenarios such as donor_late and same_lineage (reuse_elements=True
    allows it; the score subcommand then reports such donors as 'same_element'). Returns (events, skipped).
    """
    rng = random.Random(seed)
    lin = {g: parse_lineage(l) for g, l in zip(genomes["genome_id"], genomes["lineage"])}
    pairs = {}
    ids = list(lin)
    for a, b in itertools.permutations(ids, 2):
        pairs.setdefault(discordant_rank(lin[a], lin[b]), []).append((a, b))
    for v in pairs.values():
        rng.shuffle(v)
    pool = {"mobile": elements[elements["type"] != "plasmid"]["element_id"].tolist(),
            "plasmid": elements[elements["type"] == "plasmid"]["element_id"].tolist()}
    for v in pool.values():
        rng.shuffle(v)
    n_used = {"mobile": 0, "plasmid": 0}

    plan = []
    for _ in range(replicates):
        for sc in scenarios:
            if sc == "insitu":
                plan += [(sc, r, d, 1.0) for r in ranks for d in depths]
            elif sc == "carrier_frac":
                plan += [(sc, ranks[i % len(ranks)], base_depth, f) for i, f in enumerate(fracs)]
            elif sc == "same_lineage":
                plan += [(sc, "same", base_depth, 1.0)]
            else:
                plan += [(sc, r, base_depth, 1.0) for r in ranks]

    # allocate the rarest ranks first (e.g. same-species pairs), number events in plan order
    order = sorted(range(len(plan)), key=lambda i: len(pairs.get(plan[i][1], [])))
    used, rows, skipped = set(), {}, []
    for i in order:
        sc, rank, depth, frac = plan[i]
        kind = "plasmid" if sc == "plasmid" else "mobile"
        if not pool[kind] or (not reuse_elements and n_used[kind] >= len(pool[kind])):
            skipped.append((sc, rank, f"no unused {kind} element"))
            continue
        pair = next(((a, b) for a, b in pairs.get(rank, []) if a not in used and b not in used), None)
        if pair is None:
            skipped.append((sc, rank, "no unused genome pair at this rank"))
            continue
        used.update(pair)
        element = pool[kind][n_used[kind] % len(pool[kind])]
        n_used[kind] += 1
        rows[i] = ({"scenario": sc,
                     "element_id": element, "donor": pair[0], "recipient": pair[1],
                     "rank": rank, "depth": depth, "carrier_frac": frac,
                     "donor_depth": donor_depth, "rare_depth": rare_depth if sc == "rare_bloom" else 0})
    ev = pd.DataFrame([rows[i] for i in sorted(rows)])
    ev.insert(0, "event_id", [f"ev{j + 1:04d}" for j in range(len(ev))])
    return ev, skipped


# ── build ───────────────────────────────────────────────────────────────────────────────

def insert(host, element, pos, tsd):
    """Insert element after host[:pos], duplicating the tsd bp target site on both sides.
    Returns (sequence, element start, element end) in 0-based half-open coordinates."""
    s = pos + tsd
    return host[:s] + element + host[pos:], s, s + len(element)


def pick_window(seq, window, rng):
    if not window or len(seq) <= window:
        return seq
    s = rng.randint(0, len(seq) - window)
    return seq[s:s + window]


def presence(sc, t, k, comp, observed, d, f, dd, rare):
    """Depths of (donor, recipient empty, recipient carrier) at timepoint index t in comp."""
    in_s = comp in observed
    donor = 0.0 if (sc == "donor_late" and t < k) else dd
    empty = carrier = 0.0
    if in_s:
        if sc == "preexisting" or t >= k:
            carrier, empty = d * f, d * (1 - f)
        else:
            empty = d
            if sc == "rare_bloom":
                carrier = rare
    elif sc == "immigration":
        carrier = d                       # the carrier flows in from the prior compartment
    return donor, empty, carrier


def build(events, genomes, elements, timepoints, compartments, observed, k, window, margin,
          read_type, junction_flank, oracle_empty, group, seed):
    rng = random.Random(seed)
    gpath = dict(zip(genomes["genome_id"], genomes["fasta"]))
    glin = dict(zip(genomes["genome_id"], genomes["lineage"]))
    el = {r.element_id: r for r in elements.itertuples(index=False)}
    el_seq = {e: longest(r.fasta) for e, r in el.items()}

    variants, truth, comp_rows, junctions, oracle, lineages = {}, [], [], {}, {}, []
    for ev in events.itertuples(index=False):
        e = el[ev.element_id]
        eseq = el_seq[ev.element_id]
        tsd = int(e.tsd) if not pd.isna(e.tsd) else 0
        don = pick_window(longest(gpath[ev.donor]), window, rng)
        rec = pick_window(longest(gpath[ev.recipient]), window, rng)
        for name, s in (("donor", don), ("recipient", rec)):
            if len(s) < 2 * margin + 1:
                sys.exit(f"{ev.event_id}: {name} sequence ({len(s)} bp) shorter than 2 x margin")
        strand = rng.choice("+-")
        ins = eseq if strand == "+" else revcomp(eseq)
        n = ev.event_id
        row = ev._asdict()
        row.update({"element_type": e.type, "element_len": len(eseq), "strand": strand,
                    "tsd": tsd, "expect": EXPECT[ev.scenario], "transfer_t": timepoints[k],
                    "donor_contig": f"{n}_don", "recipient_contig": f"{n}_rec",
                    "empty_contig": f"{n}_rec_empty"})
        if e.type == "plasmid":
            variants[f"{n}_don"] = don
            variants[f"{n}_rec"] = rec
            variants[f"{n}_plasmid"] = eseq
            row.update({"rec_start": None, "rec_end": None, "donor_start": None, "donor_end": None,
                        "plasmid_contig": f"{n}_plasmid"})
            carrier_ids, donor_ids = [f"{n}_rec", f"{n}_plasmid"], [f"{n}_don", f"{n}_plasmid"]
        else:
            dpos = rng.randint(margin, len(don) - margin)
            rpos = rng.randint(margin, len(rec) - margin)
            variants[f"{n}_don"], ds, de = insert(don, ins, dpos, tsd)
            variants[f"{n}_rec"], rs, re_ = insert(rec, ins, rpos, tsd)
            row.update({"rec_start": rs, "rec_end": re_, "donor_start": ds, "donor_end": de,
                        "plasmid_contig": ""})
            c = variants[f"{n}_rec"]
            junctions[f"{n}_left"] = c[max(0, rs - junction_flank):rs + junction_flank]
            junctions[f"{n}_right"] = c[max(0, re_ - junction_flank):re_ + junction_flank]
            carrier_ids, donor_ids = [f"{n}_rec"], [f"{n}_don"]
        variants[f"{n}_rec_empty"] = rec
        truth.append(row)

        for t_i, t in enumerate(timepoints):
            for comp in compartments:
                dd, de_, dc = presence(ev.scenario, t_i, k, comp, observed, ev.depth,
                                       ev.carrier_frac, ev.donor_depth, ev.rare_depth)
                sample = f"{group}_{t}_{comp}"
                if ev.scenario == "plasmid":
                    # plasmid copy number 1 per host cell: its depth follows its hosts
                    roles = [(f"{n}_don", "donor", dd), (f"{n}_rec_empty", "empty", de_),
                             (f"{n}_rec", "carrier", dc), (f"{n}_plasmid", "plasmid", dd + dc)]
                else:
                    roles = [(f"{n}_don", "donor", dd), (f"{n}_rec_empty", "empty", de_),
                             (f"{n}_rec", "carrier", dc)]
                for seq_id, role, depth in roles:
                    if depth > 0:
                        comp_rows.append({"sample": sample, "timepoint": t, "compartment": comp,
                                          "event_id": n, "seq_id": seq_id, "role": role,
                                          "depth": round(depth, 4)})
        for s in carrier_ids + donor_ids + ([f"{n}_rec_empty"] if oracle_empty else []):
            oracle[s] = variants[s]
        lineages += [(f"{n}_don", glin[ev.donor]), (f"{n}_rec", glin[ev.recipient]),
                     (f"{n}_rec_empty", glin[ev.recipient])]
        if e.type == "plasmid":
            lineages.append((f"{n}_plasmid", ""))

    samples = []
    for t in timepoints:
        for comp in compartments:
            s = f"{group}_{t}_{comp}"
            reads = (f"reads/{s}_R1.fq.gz,reads/{s}_R2.fq.gz" if read_type == "short_pe"
                     else f"reads/{s}.fq.gz")
            samples.append({"sample": s, "reads": reads, "read_type": read_type, "timepoint": t,
                            "group": group, "compartment": comp})
    return (variants, pd.DataFrame(truth), pd.DataFrame(comp_rows), junctions, oracle,
            lineages, pd.DataFrame(samples))


def mmseqs_lineage(s):
    """GTDB d__X;p__Y -> mmseqs-style d_X;p_Y (what insitu.py --contig-taxonomy reads)."""
    names = parse_lineage(s)
    return ";".join(f"{r}_{x}" for r, x in zip(RANKS, names) if x)


# ── reads ───────────────────────────────────────────────────────────────────────────────

def concat(parts, out):
    """Concatenate (gzipped or plain) FASTQ parts into one gzip file."""
    with open(out, "wb") as o:
        for p in parts:
            if p.endswith(".gz"):
                with open(p, "rb") as fh:
                    shutil.copyfileobj(fh, o)
            else:
                with open(p, "rb") as fh, gzip.GzipFile(fileobj=o, mode="wb") as gz:
                    shutil.copyfileobj(fh, gz)


def simulate(comp, variants, sample, read_type, simulator, out_dir, background, seed, threads,
             read_len=150, insert_mean=400, insert_sd=40, error=0.002, iss_model="novaseq"):
    rows = comp[comp["sample"] == sample]
    if rows.empty:
        sys.exit(f"{sample}: not in composition.tsv")
    os.makedirs(out_dir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix=f"{sample}.", dir=out_dir)
    sim_parts = {"R1": [], "R2": [], "long": []}
    try:
        if simulator == "iss":
            fa = os.path.join(tmp, "genomes.fa")
            write_fasta(fa, {r.seq_id: variants[r.seq_id] for r in rows.itertuples()})
            cov = os.path.join(tmp, "coverage.tsv")
            rows[["seq_id", "depth"]].to_csv(cov, sep="\t", header=False, index=False)
            prefix = os.path.join(tmp, "iss")
            subprocess.run(["iss", "generate", "--genomes", fa, "--coverage_file", cov,
                            "--model", iss_model, "--output", prefix, "--cpus", str(threads),
                            "--seed", str(seed), "--compress"], check=True)
            sim_parts["R1"].append(prefix + "_R1.fastq.gz")
            sim_parts["R2"].append(prefix + "_R2.fastq.gz")
        else:
            for i, r in enumerate(rows.itertuples()):
                fa = os.path.join(tmp, f"{r.seq_id}.fa")
                write_fasta(fa, {r.seq_id: variants[r.seq_id]})
                L = len(variants[r.seq_id])
                if simulator == "wgsim":
                    n_pairs = int(round(r.depth * L / (2 * read_len)))
                    if n_pairs == 0:
                        continue
                    p1, p2 = os.path.join(tmp, f"{i}_1.fq"), os.path.join(tmp, f"{i}_2.fq")
                    subprocess.run(["wgsim", "-S", str(seed + i), "-N", str(n_pairs),
                                    "-1", str(read_len), "-2", str(read_len), "-d", str(insert_mean),
                                    "-s", str(insert_sd), "-e", str(error), "-r", "0", "-R", "0",
                                    fa, p1, p2], check=True, capture_output=True)
                    sim_parts["R1"].append(p1)
                    sim_parts["R2"].append(p2)
                elif simulator == "badread":
                    out = os.path.join(tmp, f"{i}.fq")
                    with open(out, "w") as fh:
                        subprocess.run(["badread", "simulate", "--reference", fa,
                                        "--quantity", f"{r.depth}x", "--seed", str(seed + i)],
                                       check=True, stdout=fh, stderr=subprocess.DEVNULL)
                    sim_parts["long"].append(out)
        bg = [b for b in (background or "").split(",") if b]
        if read_type == "short_pe":
            if bg and len(bg) != 2:
                sys.exit("short_pe background needs two files: R1,R2")
            outs = []
            for j, mate in enumerate(("R1", "R2")):
                out = os.path.join(out_dir, f"{sample}_{mate}.fq.gz")
                concat(([bg[j]] if bg else []) + sim_parts[mate], out)
                outs.append(out)
        else:
            out = os.path.join(out_dir, f"{sample}.fq.gz")
            concat(bg + sim_parts["long"], out)
            outs = [out]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return outs


# ── score (oracle catalogue: contigs are named after the planted sequences) ─────────────

OK = {"called": {"called_A", "called_B"},
      "none": {"none", "excluded"},
      "not_A": {"none", "excluded", "called_B", "unresolved", "wrong_time", "wrong_direction"},
      "unresolved": {"unresolved"},
      "excluded": {"excluded", "none"}}


def observe(t, events, excluded, element_of):
    """What kairos reported for one planted event (a truth.tsv row), judged by its recipient:
    (outcome, donor) where donor is 'planted' (the event's own donor), 'same_element' (another
    contig carrying the same element: also a valid donor as far as the sequence can tell),
    'other' or ''."""
    rec, don = t["recipient_contig"], t["donor_contig"]
    recs = {rec} | ({t["plasmid_contig"]} if t["plasmid_contig"] else set())

    def donor_kind(d):
        if d == don:
            return "planted"
        return "same_element" if element_of.get(d) == t["element_id"] else "other"
    if events is not None and not events.empty:
        mine = events[events["recipient_contig"].isin(recs)]
        right = mine[mine["direction"] == "directed"]
        on_time = right[right["timepoint"].astype(str) == str(t["transfer_t"])]
        if not on_time.empty:
            best = on_time.assign(_p=[donor_kind(d) != "planted" for d in on_time["donor_contig"]]
                                  ).sort_values(["_p", "tier"]).iloc[0]
            return f"called_{best['tier']}", donor_kind(best["donor_contig"])
        if not right.empty:
            return "wrong_time", donor_kind(right.iloc[0]["donor_contig"])
        if not mine.empty:
            return "unresolved", donor_kind(mine.iloc[0]["donor_contig"])
        back = events[(events["recipient_contig"] == don) & (events["donor_contig"].isin(recs))]
        if not back.empty:
            return "wrong_direction", "planted"
    if excluded is not None and not excluded.empty and \
            excluded["recipient_contig"].isin(recs | {don}).any():
        return "excluded", ""
    return "none", ""


def score(truth, events, excluded):
    """Per planted event: observed outcome, donor kind and whether the outcome matches the
    expectation. Also the reported events whose recipient is not a planted recipient (false
    positives: e.g. a donor or an empty-site contig called as a recipient)."""
    truth = truth.copy()
    truth["plasmid_contig"] = truth.get("plasmid_contig", pd.Series("", index=truth.index)).fillna("")
    element_of = {}
    for t in truth.itertuples():
        for c in (t.recipient_contig, t.donor_contig):
            element_of[c] = t.element_id
    rows = []
    for t in truth.to_dict("records"):
        obs, donor = observe(t, events, excluded, element_of)
        rows.append({**t, "observed": obs, "donor_match": donor,
                     "correct": obs in OK[t["expect"]]})
    res = pd.DataFrame(rows)
    recipients = set(truth["recipient_contig"]) | set(truth["plasmid_contig"]) - {""}
    fp = events[~events["recipient_contig"].isin(recipients)] if events is not None else None
    return res, fp

# ── CLI ─────────────────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("design", help="events.tsv from genome and element pools")
    d.add_argument("--genomes", required=True)
    d.add_argument("--elements", required=True)
    d.add_argument("--out", required=True)
    d.add_argument("--scenarios", default=",".join(SCENARIOS))
    d.add_argument("--ranks", default="p,c,o,f,g", help="donor-recipient discordant ranks")
    d.add_argument("--depths", default="0.5,1,2,5,10,30", help="recipient depths (insitu grid)")
    d.add_argument("--carrier-fracs", default="0.05,0.2,0.5,1")
    d.add_argument("--base-depth", type=float, default=10, help="recipient depth, other scenarios")
    d.add_argument("--donor-depth", type=float, default=10)
    d.add_argument("--rare-depth", type=float, default=0.2, help="carrier depth before k (rare_bloom)")
    d.add_argument("--replicates", type=int, default=1)
    d.add_argument("--seed", type=int, default=1)
    d.add_argument("--reuse-elements", action="store_true",
                   help="let several events plant the same element (see design() for the caveat)")

    b = sub.add_parser("build", help="sequences, composition and truth from events.tsv")
    b.add_argument("--events", required=True)
    b.add_argument("--genomes", required=True)
    b.add_argument("--elements", required=True)
    b.add_argument("--out-dir", required=True)
    b.add_argument("--timepoints", default="1,2,3,4,5")
    b.add_argument("--transfer-index", type=int, default=2, help="index of the transfer timepoint")
    b.add_argument("--compartments", default="SYS,INF", help="all compartments sampled")
    b.add_argument("--observed", default="SYS", help="compartments where recipients live")
    b.add_argument("--window", type=int, default=0,
                   help="use a random window of this many bp of each genome (0 = whole chromosome)")
    b.add_argument("--margin", type=int, default=5000, help="min distance of an insertion from an end")
    b.add_argument("--read-type", default="short_pe", choices=["short_pe", "ont"])
    b.add_argument("--junction-flank", type=int, default=250)
    b.add_argument("--oracle-empty", action="store_true",
                   help="also put the empty recipient in the oracle catalogue (as a real catalogue would)")
    b.add_argument("--group", default="SPIKE")
    b.add_argument("--seed", type=int, default=1)

    r = sub.add_parser("reads", help="simulate one sample's reads")
    r.add_argument("--build-dir", required=True)
    r.add_argument("--sample", required=True)
    r.add_argument("--simulator", default="wgsim", choices=["wgsim", "iss", "badread"])
    r.add_argument("--background", default="", help="background reads: R1,R2 (short) or one file (long)")
    r.add_argument("--out-dir", default="", help="default: <build-dir>/reads")
    r.add_argument("--threads", type=int, default=1)
    r.add_argument("--seed", type=int, default=1)

    sc = sub.add_parser("score", help="score insitu_events.tsv against truth.tsv (oracle catalogue)")
    sc.add_argument("--truth", required=True)
    sc.add_argument("--events", required=True)
    sc.add_argument("--excluded", default="")
    sc.add_argument("--out", required=True, help="per-event table; summary goes to stdout")
    args = ap.parse_args(argv)

    if args.cmd == "design":
        genomes = pd.read_csv(args.genomes, sep="\t", dtype=str)
        elements = pd.read_csv(args.elements, sep="\t", dtype={"tsd": float})
        ev, skipped = design(genomes, elements, args.scenarios.split(","), args.ranks.split(","),
                             [float(x) for x in args.depths.split(",")],
                             [float(x) for x in args.carrier_fracs.split(",")],
                             args.base_depth, args.donor_depth, args.rare_depth, args.replicates,
                             args.seed, args.reuse_elements)
        ev.to_csv(args.out, sep="\t", index=False)
        for sc, rank, why in skipped:
            print(f"skipped {sc} at rank {rank}: {why}", file=sys.stderr)
        print(f"{len(ev)} events, {len(skipped)} skipped -> {args.out}", file=sys.stderr)

    elif args.cmd == "build":
        events = pd.read_csv(args.events, sep="\t")
        genomes = pd.read_csv(args.genomes, sep="\t", dtype=str)
        elements = pd.read_csv(args.elements, sep="\t", dtype={"tsd": float})
        tps = args.timepoints.split(",")
        comps = args.compartments.split(",")
        observed = set(args.observed.split(","))
        if not 0 < args.transfer_index < len(tps):
            sys.exit("--transfer-index must leave at least one timepoint before the transfer")
        if observed - set(comps):
            sys.exit("--observed must be a subset of --compartments")
        if "immigration" in set(events["scenario"]) and not set(comps) - observed:
            sys.exit("immigration needs a prior-only compartment (e.g. --compartments SYS,INF)")
        os.makedirs(args.out_dir, exist_ok=True)
        variants, truth, comp, junctions, oracle, lineages, samples = build(
            events, genomes, elements, tps, comps, observed, args.transfer_index, args.window,
            args.margin, args.read_type, args.junction_flank, args.oracle_empty, args.group,
            args.seed)
        o = args.out_dir
        write_fasta(os.path.join(o, "variants.fa"), variants)
        write_fasta(os.path.join(o, "junctions.fa"), junctions)
        write_fasta(os.path.join(o, "oracle_catalogue.fa"), oracle)
        truth.to_csv(os.path.join(o, "truth.tsv"), sep="\t", index=False)
        comp.to_csv(os.path.join(o, "composition.tsv"), sep="\t", index=False)
        samples.to_csv(os.path.join(o, "samples.tsv"), sep="\t", index=False)
        with open(os.path.join(o, "oracle_lineages.tsv"), "w") as fh:
            for contig, lin in lineages:
                fh.write(f"{contig}\t0\tno rank\t-\t{mmseqs_lineage(lin)}\n")
        with open(os.path.join(o, "build.json"), "w") as fh:
            json.dump({k: v for k, v in vars(args).items() if k != "cmd"}, fh, indent=1)
        print(f"{len(truth)} events, {len(samples)} samples, {len(variants)} sequences -> {o}",
              file=sys.stderr)

    elif args.cmd == "reads":
        comp = pd.read_csv(os.path.join(args.build_dir, "composition.tsv"), sep="\t",
                           dtype={"timepoint": str})
        samples = pd.read_csv(os.path.join(args.build_dir, "samples.tsv"), sep="\t", dtype=str)
        rt = samples.set_index("sample").loc[args.sample, "read_type"]
        if (rt == "ont") != (args.simulator == "badread"):
            sys.exit(f"{args.sample} is {rt}: use --simulator badread for ont, wgsim/iss for short_pe")
        variants = read_fasta(os.path.join(args.build_dir, "variants.fa"))
        # a different seed per sample, so samples never share simulated reads
        seed = args.seed * 100003 + samples["sample"].tolist().index(args.sample) * 1009
        outs = simulate(comp, variants, args.sample, rt, args.simulator,
                        args.out_dir or os.path.join(args.build_dir, "reads"),
                        args.background, seed, args.threads)
        print(" ".join(outs))

    elif args.cmd == "score":
        truth = pd.read_csv(args.truth, sep="\t", dtype={"transfer_t": str}, keep_default_na=False,
                            na_values=[""])
        events = pd.read_csv(args.events, sep="\t", dtype={"timepoint": str})
        excluded = (pd.read_csv(args.excluded, sep="\t", dtype={"timepoint": str})
                    if args.excluded and os.path.getsize(args.excluded) else None)
        res, fp = score(truth, events, excluded)
        res.to_csv(args.out, sep="\t", index=False)
        tab = res.groupby(["scenario", "expect", "observed", "donor_match"]).size().rename("n").reset_index()
        print(tab.to_string(index=False))
        print(f"\n{int(res['correct'].sum())}/{len(res)} planted events as expected; "
              f"{len(fp)} reported events with a recipient that was never planted")


if __name__ == "__main__":
    main()
