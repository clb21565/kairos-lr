#!/usr/bin/env python3
"""
kairos detect-state | present / absent / undetermined for every shared segment in every observation

A segment is one contig's copy of a shared region (one side of an event in events.tsv), with up
to two breakpoints. An observation is one (group, timepoint, compartment) in the sample sheet; the
samples in it (e.g. short and long reads of the same sample) are pooled (--pool observation,
the default). With --pool sample every sample is called on its own: one row per sample x segment,
nothing summed across libraries; insitu.py then combines the per-sample states of an observation
(see there).

For each observation and segment:
  n         junction reads: spanning reads + bridging pairs at the segment's breakpoints
  lambda    expected spanning reads if every host cell carried the segment:
              sum over samples and breakpoints of  flank_depth * (L - 2a) / L
            (L: mean read length, a: overhang). Bridging pairs are counted in n but not in
            lambda, which only makes the absent call more conservative.
  empty     empty-site evidence (split reads + read pairs) for this contig without the segment

A segment counts as present only when both the region and its host contig are seen:

  boundaries     the segment's breakpoints with >= k junction reads each (n_left, n_right)
  contig_state   whether the contig's own sequence (contigs.bed: the contig minus all its
                 regions) is present: pooled breadth >= --contig-min-breadth and
                 breadth >= --contig-breadth-ratio x the breadth expected at that depth
                 (1 - exp(-depth)); the ratio rules out reads piled up on a short stretch.
                 Otherwise "not_detected".

  present        the segment has >= --min-boundaries breakpoints, every one of them has >= k
                 junction reads, and contig_state is present
  partial        n >= k junction reads in total, but not present: only one boundary seen, the
                 segment has fewer than --min-boundaries breakpoints (contig-end region), or the
                 contig itself is not detected. In situ never counts partial as seen before,
                 nor as a recipient or donor; a partial previous state only lowers the tier to B.
  absent         n < k and P(N < k | lambda * f) <= alpha   (Poisson; f = --absent-fraction)
  undetermined   otherwise (too shallow to tell), or the segment has no breakpoint

kairos assess readout (reported, not used in the states above):
  n_window            reads in the original kairos assess windows (window_reads from count_features:
                      >= 100 bp aligned within +-75 bp of a breakpoint), summed over the segment's
                      breakpoints
  n_window_tied       reads that pass the window rule but are tied (MAPQ < the junction MAPQ, i.e.
                      equally good elsewhere), summed the same way; metadata
  contig_window_prop  fraction of the contig's breakpoint windows (all of its segments) with >= 1 read
  contig_window_pass  contig_window_prop > --window-min-prop (kairos assess: 0.9)

"absent" means the segment is not carried by a fraction f or more of the host population at
that depth. Carrier-fraction estimates are reported as n_span / lambda (spanning reads against
the depth expectation; unbiased, noisy at low depth) and n / (n + empty) (against the empty
site, when it exists; rough, because pairs and reads are mixed).
"""

import argparse
import math
import os
import sys

import numpy as np
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
    parser.add_argument("--min-boundaries", type=int, default=2,
                        help="breakpoints a segment needs to be called present (2 = internal regions only)")
    parser.add_argument("--contig-min-breadth", type=float, default=0.5,
                        help="breadth of the contig's own sequence needed to call the contig present")
    parser.add_argument("--contig-breadth-ratio", type=float, default=0.5,
                        help="observed / expected breadth needed to call the contig present")
    parser.add_argument("--window-min-prop", type=float, default=0.9,
                        help="kairos assess readout: contig detected if more than this fraction of its "
                             "breakpoint windows have a read")
    parser.add_argument("--pool", choices=["observation", "sample"], default="observation",
                        help="observation: sum each observation's samples before calling; "
                             "sample: call every sample on its own")
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
        if f"region_type_{i}" in events.columns:
            cols[f"region_type_{i}"] = "region_type"
        s = events[list(cols)].rename(columns=cols)
        s["event_id"], s["side"] = events["event_id"], i
        sides.append(s)
    sides = pd.concat(sides, ignore_index=True).fillna({"bp_left": "", "bp_right": "", "empty_site": ""})
    key = ["contig", "bp_left", "bp_right", "empty_site"]
    seg = (sides.groupby(key, sort=False)
           .agg(shared_start=("shared_start", "min"), shared_end=("shared_end", "max"),
                **({"region_type": ("region_type", lambda x: ",".join(sorted(set(x) - {""})))}
                   if "region_type" in sides.columns else {}))
           .reset_index())
    seg.insert(0, "segment_id", [f"SG{i:07d}" for i in range(1, len(seg) + 1)])
    sides = sides.merge(seg[key + ["segment_id"]], on=key)
    return seg, sides[["event_id", "side", "segment_id"]]


def load_counts(samples, counts_dir):
    bp, es, ct = [], [], []
    for s in samples["sample"]:
        d = os.path.join(counts_dir, s)
        bp.append(pd.read_csv(os.path.join(d, f"{s}.breakpoint_counts.tsv"), sep="\t"))
        es.append(pd.read_csv(os.path.join(d, f"{s}.empty_site_counts.tsv"), sep="\t"))
        ct.append(pd.read_csv(os.path.join(d, f"{s}.contig_counts.tsv"), sep="\t"))
    meta = samples[["sample"] + OBS]
    return [pd.concat(x, ignore_index=True).merge(meta, on="sample") for x in (bp, es, ct)]


def contig_states(ct, args, keys=OBS):
    """Pool the samples in each `keys` unit per contig; breadths combine as independent coverage."""
    ct = ct.assign(miss=1 - ct["breadth"].fillna(0))
    c = ct.groupby(keys + ["contig"], as_index=False).agg(
        contig_reads=("reads", "sum"), contig_depth=("depth", "sum"), miss=("miss", "prod"))
    c["contig_breadth"] = 1 - c["miss"]
    expected = 1 - np.exp(-c["contig_depth"])
    ok = ((c["contig_breadth"] >= args.contig_min_breadth) &
          (c["contig_breadth"] >= args.contig_breadth_ratio * expected) & (c["contig_reads"] > 0))
    c["contig_state"] = np.where(ok, "present", "not_detected")
    c["contig_depth"] = c["contig_depth"].round(2)
    c["contig_breadth"] = c["contig_breadth"].round(3)
    return c.drop(columns="miss")


def call_states(seg, bp, es, ct, args, keys=OBS):
    """States per `keys` unit x segment: OBS pools an observation's samples, OBS + [sample] does not."""
    # expected spanning reads per sample and breakpoint
    L, a = bp["mean_read_len"], bp["overhang"]
    bp["lam"] = (bp["flank_depth"].fillna(0) * ((L - 2 * a) / L).clip(lower=0)).fillna(0)
    bp["n"] = bp["span"] + bp["bridge"]
    for c in ("window_reads", "window_reads_tied"):   # counts from before the window readout
        if c not in bp.columns:
            bp[c] = float("nan")
    per_bp = bp.groupby(keys + ["bp_id"], as_index=False).agg(
        n=("n", "sum"), n_span=("span", "sum"), lam=("lam", "sum"), flank_depth=("flank_depth", "sum"),
        n_window=("window_reads", lambda x: x.sum(min_count=1)),
        n_window_tied=("window_reads_tied", lambda x: x.sum(min_count=1)))
    per_es = es.groupby(keys + ["es_id"], as_index=False)["empty"].sum()

    # segment -> its breakpoints (long form)
    seg_bp = pd.concat([
        seg[["segment_id", "bp_left"]].rename(columns={"bp_left": "bp_id"}).assign(end="left"),
        seg[["segment_id", "bp_right"]].rename(columns={"bp_right": "bp_id"}).assign(end="right"),
    ])
    seg_bp = seg_bp[seg_bp["bp_id"] != ""]

    obs = bp[keys].drop_duplicates()
    grid = obs.merge(seg[["segment_id", "contig", "empty_site"]], how="cross")
    agg = (seg_bp.merge(per_bp, on="bp_id")
           .groupby(keys + ["segment_id"], as_index=False)
           .agg(n=("n", "sum"), n_span=("n_span", "sum"), lam=("lam", "sum"),
                flank_depth=("flank_depth", "mean"),
                n_window=("n_window", lambda x: x.sum(min_count=1)),
                n_window_tied=("n_window_tied", lambda x: x.sum(min_count=1)),
                n_breakpoints=("bp_id", "nunique")))
    k, f = args.min_reads, args.absent_fraction
    ends = seg_bp.merge(per_bp, on="bp_id")
    ends = (ends.pivot_table(index=keys + ["segment_id"], columns="end", values="n", aggfunc="sum")
            .reindex(columns=["left", "right"]).add_prefix("n_").reset_index())
    ends.columns.name = None
    st = grid.merge(agg, on=keys + ["segment_id"], how="left")
    st = st.merge(ends, on=keys + ["segment_id"], how="left")
    st = st.merge(contig_states(ct, args, keys), on=keys + ["contig"], how="left")
    # kairos assess readout: fraction of the contig's breakpoint windows with a read
    cbp = seg_bp.merge(seg[["segment_id", "contig"]], on="segment_id")[["contig", "bp_id"]].drop_duplicates()
    cw = cbp.merge(per_bp, on="bp_id")
    cw = (cw.assign(hit=(cw["n_window"] > 0).where(cw["n_window"].notna()))
          .groupby(keys + ["contig"], as_index=False)["hit"].mean()
          .rename(columns={"hit": "contig_window_prop"}))
    st = st.merge(cw, on=keys + ["contig"], how="left")
    st["contig_window_pass"] = (st["contig_window_prop"] > args.window_min_prop).where(
        st["contig_window_prop"].notna())
    st["contig_window_prop"] = st["contig_window_prop"].round(3)
    st["contig_state"] = st["contig_state"].fillna("not_detected")
    st = st.merge(per_es.rename(columns={"es_id": "empty_site"}), on=keys + ["empty_site"], how="left")
    st[["n", "n_span", "lam", "n_breakpoints"]] = st[["n", "n_span", "lam", "n_breakpoints"]].fillna(0)
    st["empty"] = st["empty"].where(st["empty_site"] != "", float("nan"))

    p_miss = st["lam"].apply(lambda lam: poisson_cdf_below(k, lam * f))
    st["p_miss"] = p_miss.round(4)
    st["state"] = "undetermined"
    st.loc[(st["n_breakpoints"] > 0) & (p_miss <= args.alpha), "state"] = "absent"
    st["boundaries"] = (st[["n_left", "n_right"]] >= k).sum(axis=1)
    st.loc[st["n"] >= k, "state"] = "partial"
    strict = ((st["n_breakpoints"] >= args.min_boundaries) & (st["boundaries"] == st["n_breakpoints"])
              & (st["contig_state"] == "present"))
    st.loc[(st["n"] >= k) & strict, "state"] = "present"
    st["carrier_frac_depth"] = (st["n_span"] / st["lam"]).where(st["lam"] > 0).round(3)
    tot = st["n"] + st["empty"]
    st["carrier_frac_empty"] = (st["n"] / tot).where(tot > 0).round(3)
    st["lam"] = st["lam"].round(2)
    st["flank_depth"] = st["flank_depth"].round(2)
    cols = keys + ["segment_id", "state", "n", "n_left", "n_right", "boundaries", "n_breakpoints",
                  "contig_state", "contig_depth", "contig_breadth", "n_span", "lam", "p_miss", "empty",
                  "flank_depth", "carrier_frac_depth", "carrier_frac_empty",
                  "n_window", "n_window_tied", "contig_window_prop", "contig_window_pass"]
    return st[cols].astype({"n": int, "n_span": int, "n_breakpoints": int, "boundaries": int,
                            "n_left": "Int64", "n_right": "Int64", "n_window": "Int64", "n_window_tied": "Int64",
                            "contig_window_pass": "boolean"})


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
    bp, es, ct = load_counts(samples, args.counts_dir)
    if args.pool == "sample":
        states = call_states(seg, bp, es, ct, args, OBS + ["sample"]).rename(columns={"sample": "samples"})
    else:
        states = call_states(seg, bp, es, ct, args)
        members = samples.groupby(OBS)["sample"].apply(lambda x: ",".join(x)).rename("samples").reset_index()
        states = states.merge(members, on=OBS, how="left")
    states.insert(len(OBS), "samples", states.pop("samples"))

    seg_out = seg.merge(
        sides.groupby("segment_id")["event_id"].apply(lambda x: ";".join(sorted(set(x)))).rename("events"),
        on="segment_id")
    seg_out.to_csv(args.out_segments, sep="\t", index=False)
    states.to_csv(args.out_states, sep="\t", index=False)
    (sides.pivot(index="event_id", columns="side", values="segment_id")
     .rename(columns={1: "segment_1", 2: "segment_2"}).reset_index()
     .to_csv(args.out_event_segments, sep="\t", index=False))
    print(f"{len(seg)} segments x {len(states[OBS + ['samples']].drop_duplicates())} "
          f"{'samples' if args.pool == 'sample' else 'observations'}: "
          + ", ".join(f"{k} {v}" for k, v in states["state"].value_counts().items()), file=sys.stderr)


if __name__ == "__main__":
    main()
