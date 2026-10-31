#!/usr/bin/env python3
"""
Compare in situ HGT calls from insitu.py with the R get_isHGTs_weekly() results
(240903-<rank>hgt-res.csv: one row per lineage pair and week ti, `edges` = "recipient->donor;...").

Matching is on (recipient contig, donor contig, timepoint). Timepoints in the sample sheet should
be the R ordered_week, or pass --timepoint-map (TSV: timepoint, week). If the catalogue renamed
contigs (e.g. assembly prefixes), pass --contig-map (TSV: old_contig, new_contig).

Each R-only edge gets a reason from the call outputs:
  no_event            the two contigs share no region pair in events.tsv
  recipient_not_new   the recipient's segment was already present earlier (or in another fraction at t)
  recipient_not_seen  the recipient's segment is not present at t (undetermined or absent)
  donor_not_before    the donor's segment was not present before t
  other_timepoint     the pair is an event, but called at another timepoint

Usage:
  python tools/compare_insitu_r.py --r 240903-phgt-res.csv --insitu results/call/insitu_events.tsv \\
      --call-dir results/call [--rank p] [--timepoint-map map.tsv] [--contig-map map.tsv] --out cmp.tsv
"""

import argparse
import sys

import pandas as pd


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Compare insitu.py events with the R in situ results.")
    ap.add_argument("--r", required=True, help="R result CSV (write.csv of get_isHGTs_weekly output)")
    ap.add_argument("--insitu", required=True, help="insitu_events.tsv")
    ap.add_argument("--call-dir", required=True, help="results/call (states, segments, event_segments)")
    ap.add_argument("--rank", default="", help="only new events discordant at this rank (d,p,c,o,f,g,s)")
    ap.add_argument("--timepoint-map", default="", help="TSV: timepoint, week")
    ap.add_argument("--contig-map", default="", help="TSV: old_contig, new_contig")
    ap.add_argument("--out", required=True, help="per-edge comparison TSV")
    return ap.parse_args(argv)


def r_edges(path, contig_map):
    r = pd.read_csv(path)
    rows = []
    for edges, ti in zip(r["edges"], r["ti"]):
        for e in str(edges).split(";"):
            if "->" not in e:
                continue
            rec, don = e.split("->", 1)
            rows.append((contig_map.get(rec, rec), contig_map.get(don, don), str(ti)))
    return pd.DataFrame(rows, columns=["recipient", "donor", "week"]).drop_duplicates()


def main(argv=None):
    a = parse_args(argv)
    cmap = {}
    if a.contig_map:
        m = pd.read_csv(a.contig_map, sep="\t", header=None, dtype=str)
        cmap = dict(zip(m[0], m[1]))
    R = r_edges(a.r, cmap)

    ev = pd.read_csv(a.insitu, sep="\t", dtype=str, keep_default_na=False)
    if a.rank:
        if "discordant_rank" not in ev:
            sys.exit("--rank needs insitu events with taxonomy (call.use_taxonomy: true)")
        ev = ev[ev["discordant_rank"] == a.rank]
    tmap = {}
    if a.timepoint_map:
        m = pd.read_csv(a.timepoint_map, sep="\t", header=None, dtype=str)
        tmap = dict(zip(m[0], m[1]))
    ev["week"] = ev["timepoint"].map(lambda t: tmap.get(t, t))
    new = ev.rename(columns={"recipient_contig": "recipient", "donor_contig": "donor"})
    # unresolved events count in both directions
    unres = new[new["direction"] == "unresolved"].rename(
        columns={"recipient": "donor", "donor": "recipient"})
    new = pd.concat([new, unres], ignore_index=True)
    N = new[["recipient", "donor", "week", "direction", "tier", "persistence"]].drop_duplicates(
        ["recipient", "donor", "week"])

    cmp = R.assign(in_r=True).merge(N.assign(in_new=True), on=["recipient", "donor", "week"], how="outer")
    cmp[["in_r", "in_new"]] = cmp[["in_r", "in_new"]].fillna(False).astype(bool)

    # reasons for R-only edges
    cd = a.call_dir
    seg = pd.read_csv(f"{cd}/segments.tsv", sep="\t", keep_default_na=False)
    es = pd.read_csv(f"{cd}/event_segments.tsv", sep="\t")
    st = pd.read_csv(f"{cd}/states.tsv", sep="\t", dtype={"timepoint": str}, keep_default_na=False)
    st["week"] = st["timepoint"].map(lambda t: tmap.get(t, t))
    seg_contig = dict(zip(seg["segment_id"], seg["contig"]))
    pairs = {}
    for s1, s2 in zip(es["segment_1"], es["segment_2"]):
        c1, c2 = seg_contig[s1], seg_contig[s2]
        pairs.setdefault((c1, c2), []).append((s1, s2))
        pairs.setdefault((c2, c1), []).append((s2, s1))
    present = st[st["state"] == "present"]
    weeks_present = present.groupby("segment_id")["week"].apply(
        lambda w: set(pd.to_numeric(w, errors="coerce").dropna())).to_dict()
    new_pairs = set(zip(N["recipient"], N["donor"]))

    def reason(row):
        if row.in_new or not row.in_r:
            return ""
        segs = pairs.get((row.recipient, row.donor))
        if not segs:
            return "no_event"
        if (row.recipient, row.donor) in new_pairs:
            return "other_timepoint"
        t = float(row.week)
        rec_weeks = set().union(*(weeks_present.get(r, set()) for r, _ in segs))
        don_weeks = set().union(*(weeks_present.get(d, set()) for _, d in segs))
        if t not in rec_weeks:
            return "recipient_not_seen"
        if any(w < t for w in rec_weeks):
            return "recipient_not_new"
        if not any(w < t for w in don_weeks):
            return "donor_not_before"
        return "recipient_not_new"      # present in another fraction at t (same-day rule)

    cmp["reason"] = cmp.apply(reason, axis=1)
    cmp.sort_values(["week", "recipient", "donor"]).to_csv(a.out, sep="\t", index=False)

    both = (cmp.in_r & cmp.in_new).sum()
    print(f"R edges: {cmp.in_r.sum()}   new events: {cmp.in_new.sum()}   in both: {both}")
    new_only = cmp[cmp.in_new & ~cmp.in_r].groupby(["direction", "tier"]).size()
    r_only = cmp[cmp.in_r & ~cmp.in_new]["reason"].value_counts()
    print("\nnew only, by direction and tier:")
    print(new_only.to_string() if len(new_only) else "  none")
    print("\nR only, by reason:")
    print(r_only.to_string() if len(r_only) else "  none")


if __name__ == "__main__":
    main()
