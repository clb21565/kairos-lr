#!/usr/bin/env python3
"""
kairos detect-state | present / absent / undetermined for every shared segment in every observation

A segment is one contig's copy of a shared region (one side of an event in events.tsv), with up
to two breakpoints. An observation is one (group, timepoint, compartment) in the sample sheet; the
samples in it (e.g. short and long reads of the same sample) are pooled.

For each observation and segment:
  n         junction reads: spanning reads + bridging pairs at the segment's breakpoints
  lambda    expected spanning reads if every host cell carried the segment:
              sum over samples and breakpoints of  flank_depth * (L - 2a) / L
            (L: mean read length, a: overhang). Bridging pairs are counted in n but not in
            lambda, which only makes the absent call more conservative.
  empty     empty-site evidence (split reads + read pairs) for this contig without the segment

  present        n >= k
  absent         n < k and P(N < k | lambda * f) <= alpha   (Poisson; f = --absent-fraction)
  undetermined   otherwise (too shallow to tell), or the segment has no breakpoint

"absent" means the segment is not carried by a fraction f or more of the host population at
that depth. Carrier-fraction estimates are reported as n_span / lambda (spanning reads against
the depth expectation; unbiased, noisy at low depth) and n / (n + empty) (against the empty
site, when it exists; rough, because pairs and reads are mixed).
"""

import argparse
import math
import os
import sys

import pandas as pd

OBS = ["group", "timepoint", "compartment"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Present/absent/undetermined per segment and observation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--samples", required=True, help="sample sheet (TSV)")
    parser.add_argument("--events", required=True, help="events.tsv from build_features.py")
    parser.add_argument("--counts-dir", required=True,
                        help="assess directory with <sample>/<sample>.breakpoint_counts.tsv etc.")
    parser.add_argument("--min-reads", "-k", type=int, default=2,
                        help="junction reads needed to call a segment present")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="largest P(fewer than k reads | carried) at which absent is called")
    parser.add_argument("--absent-fraction", type=float, default=1.0,
                        help="carrier fraction that an absent call rules out")
    parser.add_argument("--out-segments", required=True, help="segments.tsv (segment definitions)")
    parser.add_argument("--out-states", required=True, help="states.tsv (one row per observation x segment)")
    parser.add_argument("--out-event-segments", required=True,
                        help="event_segments.tsv (event_id, segment_1, segment_2)")
    return parser.parse_args(argv)


def poisson_cdf_below(k, lam):
    """P(N < k) for N ~ Poisson(lam)."""
    if lam <= 0:
        return 1.0
    return sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k))


def segments_from_events(events):
    """One row per contig side of an event, deduplicated on (contig, breakpoints, empty site)."""
    sides = []
    for i in (1, 2):
        cols = {f"contig_{i}": "contig", f"shared_start_{i}": "shared_start",
                f"shared_end_{i}": "shared_end", f"bp_left_{i}": "bp_left",
                f"bp_right_{i}": "bp_right", f"empty_site_{i}": "empty_site"}
        s = events[list(cols)].rename(columns=cols)
        s["event_id"], s["side"] = events["event_id"], i
        sides.append(s)
    sides = pd.concat(sides, ignore_index=True).fillna({"bp_left": "", "bp_right": "", "empty_site": ""})
    key = ["contig", "bp_left", "bp_right", "empty_site"]
    seg = (sides.groupby(key, sort=False)
           .agg(shared_start=("shared_start", "min"), shared_end=("shared_end", "max"))
           .reset_index())
    seg.insert(0, "segment_id", [f"SG{i:07d}" for i in range(1, len(seg) + 1)])
    sides = sides.merge(seg[key + ["segment_id"]], on=key)
    return seg, sides[["event_id", "side", "segment_id"]]


def load_counts(samples, counts_dir):
    bp, es = [], []
    for s in samples["sample"]:
        d = os.path.join(counts_dir, s)
        bp.append(pd.read_csv(os.path.join(d, f"{s}.breakpoint_counts.tsv"), sep="\t"))
        es.append(pd.read_csv(os.path.join(d, f"{s}.empty_site_counts.tsv"), sep="\t"))
    bp = pd.concat(bp, ignore_index=True)
    es = pd.concat(es, ignore_index=True)
    meta = samples[["sample"] + OBS]
    return bp.merge(meta, on="sample"), es.merge(meta, on="sample")


def call_states(seg, bp, es, args):
    # expected spanning reads per sample and breakpoint
    L, a = bp["mean_read_len"], bp["overhang"]
    bp["lam"] = (bp["flank_depth"].fillna(0) * ((L - 2 * a) / L).clip(lower=0)).fillna(0)
    bp["n"] = bp["span"] + bp["bridge"]
    per_bp = bp.groupby(OBS + ["bp_id"], as_index=False).agg(
        n=("n", "sum"), n_span=("span", "sum"), lam=("lam", "sum"), flank_depth=("flank_depth", "sum"))
    per_es = es.groupby(OBS + ["es_id"], as_index=False)["empty"].sum()

    # segment -> its breakpoints (long form)
    seg_bp = pd.concat([
        seg[["segment_id", "bp_left"]].rename(columns={"bp_left": "bp_id"}),
        seg[["segment_id", "bp_right"]].rename(columns={"bp_right": "bp_id"}),
    ])
    seg_bp = seg_bp[seg_bp["bp_id"] != ""]

    obs = bp[OBS].drop_duplicates()
    grid = obs.merge(seg[["segment_id", "empty_site"]], how="cross")
    agg = (seg_bp.merge(per_bp, on="bp_id")
           .groupby(OBS + ["segment_id"], as_index=False)
           .agg(n=("n", "sum"), n_span=("n_span", "sum"), lam=("lam", "sum"),
                flank_depth=("flank_depth", "mean"),
                n_breakpoints=("bp_id", "nunique")))
    st = grid.merge(agg, on=OBS + ["segment_id"], how="left")
    st = st.merge(per_es.rename(columns={"es_id": "empty_site"}), on=OBS + ["empty_site"], how="left")
    st[["n", "n_span", "lam", "n_breakpoints"]] = st[["n", "n_span", "lam", "n_breakpoints"]].fillna(0)
    st["empty"] = st["empty"].where(st["empty_site"] != "", float("nan"))

    k, f = args.min_reads, args.absent_fraction
    p_miss = st["lam"].apply(lambda lam: poisson_cdf_below(k, lam * f))
    st["p_miss"] = p_miss.round(4)
    st["state"] = "undetermined"
    st.loc[(st["n_breakpoints"] > 0) & (p_miss <= args.alpha), "state"] = "absent"
    st.loc[st["n"] >= k, "state"] = "present"
    st["carrier_frac_depth"] = (st["n_span"] / st["lam"]).where(st["lam"] > 0).round(3)
    tot = st["n"] + st["empty"]
    st["carrier_frac_empty"] = (st["n"] / tot).where(tot > 0).round(3)
    st["lam"] = st["lam"].round(2)
    st["flank_depth"] = st["flank_depth"].round(2)
    cols = OBS + ["segment_id", "state", "n", "n_span", "lam", "p_miss", "empty", "flank_depth",
                  "carrier_frac_depth", "carrier_frac_empty", "n_breakpoints"]
    return st[cols].astype({"n": int, "n_span": int, "n_breakpoints": int})


def main(argv=None):
    args = parse_args(argv)
    samples = pd.read_csv(args.samples, sep="\t", dtype=str, comment="#").fillna("")
    for c in OBS:
        if c not in samples.columns:
            samples[c] = ""
    events = pd.read_csv(args.events, sep="\t", keep_default_na=False,
                         dtype={"bp_left_1": str, "bp_right_1": str, "empty_site_1": str,
                                "bp_left_2": str, "bp_right_2": str, "empty_site_2": str})
    seg, sides = segments_from_events(events)
    bp, es = load_counts(samples, args.counts_dir)
    states = call_states(seg, bp, es, args)

    seg_out = seg.merge(
        sides.groupby("segment_id")["event_id"].apply(lambda x: ";".join(sorted(set(x)))).rename("events"),
        on="segment_id")
    seg_out.to_csv(args.out_segments, sep="\t", index=False)
    states.to_csv(args.out_states, sep="\t", index=False)
    (sides.pivot(index="event_id", columns="side", values="segment_id")
     .rename(columns={1: "segment_1", 2: "segment_2"}).reset_index()
     .to_csv(args.out_event_segments, sep="\t", index=False))
    print(f"{len(seg)} segments x {len(states[OBS].drop_duplicates())} observations: "
          + ", ".join(f"{k} {v}" for k, v in states["state"].value_counts().items()), file=sys.stderr)


if __name__ == "__main__":
    main()
