#!/usr/bin/env python3
"""
kairos in situ | HGT events that happened within the sampled system, ordered in time

Port of get_isHGTs_weekly() (CIWARS influent/effluent study) to segment states from
detect_state.py. "Present" is detect_state's strict call: every boundary of the region has
junction reads and the host contig's own sequence is detected. "Seen" is present or partial
(some junction reads, but not the full region + contig evidence). For each event (two contigs sharing a segment), each contig is tried as the
recipient. Within a group, at timepoint t:

  1. Recipient new: its segment is present in an observed compartment f at t, was not seen in any
     prior compartment at an earlier timepoint, and is not seen at t in any other prior compartment
     (<=: rules out same-day import, e.g. influent -> effluent). t must have an earlier
     observation (nothing is "new" at the first timepoint).
  2. Donor first: the other contig's segment was present in a prior compartment at an earlier
     timepoint (<: this is what gives the direction). If instead it is first seen at t,
     the event is kept with direction "unresolved".
  3. Tier: A if, at the previous timepoint, the recipient's segment was confidently absent
     (detect_state "absent") in a prior compartment; B if it was only undetermined (could be
     immigration or below detection) or not observed there.
  4. Persistence: whether the segment is still present in an observed compartment at the next
     timepoint (yes / partial / no / undetermined / NA).
  5. Lineage (--contig-taxonomy): HGT is between lineages, so events whose two contigs have the
     same lineage are never reported. --lineage-filter discordant (default) keeps only events
     whose contigs are classified and differ at some rank; not_same also keeps events whose
     lineages cannot be compared (unclassified, or identical as far as both go); none keeps all.
     Removed events go to --excluded with the reason.

Timepoints are ordered numerically if they are all numbers, else as dates if they all parse,
else as text. Observed and prior compartments come from --observed / --prior (default: all).
"""

import argparse
import sys

import pandas as pd

RANKS = ["d", "p", "c", "o", "f", "g", "s"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="In situ HGT calls from segment states.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--states", required=True, help="states.tsv from detect_state.py")
    parser.add_argument("--segments", required=True, help="segments.tsv from detect_state.py")
    parser.add_argument("--event-segments", required=True, help="event_segments.tsv from detect_state.py")
    parser.add_argument("--observed", default="",
                        help="comma-separated compartments in which recipients are looked for (default: all)")
    parser.add_argument("--prior", default="",
                        help="comma-separated compartments that count as having seen a contig (default: all)")
    parser.add_argument("--contig-taxonomy", default="",
                        help="mmseqs taxonomy TSV for the catalogue contigs (assignRes2.tsv), optional")
    parser.add_argument("--lineage-filter", choices=["discordant", "not_same", "none"], default="discordant",
                        help="which events to keep by lineage relation (needs --contig-taxonomy)")
    parser.add_argument("--out", required=True, help="insitu_events.tsv")
    parser.add_argument("--excluded", default="", help="events removed by the lineage filter, with reason")
    return parser.parse_args(argv)


def time_order(values):
    """Map each timepoint label to its rank."""
    vals = pd.Series(sorted(set(values)))
    num = pd.to_numeric(vals, errors="coerce")
    if num.notna().all():
        key = num
    else:
        dates = pd.to_datetime(vals, errors="coerce")
        key = dates if dates.notna().all() else vals
    order = vals.iloc[key.argsort(kind="stable").values].tolist()
    return {v: i for i, v in enumerate(order)}


def read_lineages(path):
    """contig -> list of names per rank (d..s); mmseqs createtsv with --tax-lineage 1."""
    tax = pd.read_csv(path, sep="\t", header=None)
    out = {}
    for contig, lineage in zip(tax[0], tax[tax.columns[-1]]):
        # keep the ranks d..s only (lineages may start with root/no-rank entries)
        parts = [] if pd.isna(lineage) else str(lineage).split(";")
        ranked = {x[0]: x[2:] for x in parts if len(x) > 2 and x[1] == "_" and x[0] in RANKS}
        out[contig] = [ranked.get(r, "") for r in RANKS]
    return out


def discordant_rank(a, b):
    """First rank (d..s) at which two lineages both have a name and differ; '' if none."""
    return lineage_relation(a, b)[1]


def lineage_relation(a, b):
    """
    (relation, discordant_rank, shared_rank) for two lineages (names per rank d..s):
      discordant   both named and different at some rank (the first such rank is returned)
      same         identical names at every rank down to species
      unresolved   identical as far as both are named, but one stops before species
                   (includes unclassified contigs)
    shared_rank is the deepest rank at which both are named and equal.
    """
    shared = ""
    for r, x, y in zip(RANKS, a, b):
        if not x or not y:
            return "unresolved", "", shared
        if x != y:
            return "discordant", r, shared
        shared = r
    return "same", "", shared


class Group:
    def __init__(self, st, observed, prior):
        self.observed, self.prior = observed, prior
        self.times = sorted(st["t"].unique())
        self.prior_times = sorted(st.loc[st["compartment"].isin(prior), "t"].unique())
        idx = st.set_index(["segment_id", "t", "compartment"])
        self.state = idx["state"].to_dict()
        self.row = idx
        self.present = self._index(st[st["state"] == "present"])
        self.seen = self._index(st[st["state"].isin(["present", "partial"])])

    @staticmethod
    def _index(d):
        return {seg: list(zip(x["t"], x["compartment"])) for seg, x in d.groupby("segment_id")}

    def present_at(self, seg, t=None, before=None, compartments=None, seen=False):
        """(t, compartment) where seg is present (seen=True: present or partial), filtered."""
        hits = (self.seen if seen else self.present).get(seg, [])
        return [(tt, f) for tt, f in hits
                if (t is None or tt == t) and (before is None or tt < before)
                and (compartments is None or f in compartments)]

    def recipient_new(self, seg):
        """(t, compartment) at which seg first appears as a recipient, or None."""
        hits = self.present_at(seg, compartments=self.observed | self.prior)
        if not hits:
            return None
        t = min(tt for tt, _ in hits)
        if not self.prior_times or t <= self.prior_times[0]:
            return None      # no earlier observation to compare with
        if self.present_at(seg, before=t, compartments=self.prior, seen=True):
            return None
        for f in sorted(self.observed):
            if (t, f) in hits:
                others = [ff for _, ff in self.present_at(seg, t=t, compartments=self.prior, seen=True)
                          if ff != f]
                if not others:
                    return t, f
        return None

    def states_at(self, seg, t, compartments):
        return [self.state.get((seg, t, f)) for f in compartments if (seg, t, f) in self.state]


def call_group(g, st, ev, observed, prior, seg_contig):
    grp = Group(st, observed, prior)
    rows = {}
    for e in ev.itertuples(index=False):
        for rec, don in ((e.segment_1, e.segment_2), (e.segment_2, e.segment_1)):
            new = grp.recipient_new(rec)
            if new is None:
                continue
            t, f = new
            donor_before = grp.present_at(don, before=t, compartments=prior)
            if donor_before:
                direction = "directed"
                t_d = max(tt for tt, _ in donor_before)
                donor_last = (t_d, ",".join(sorted({ff for tt, ff in donor_before if tt == t_d})))
            elif grp.present_at(don, t=t) and not grp.present_at(don, before=t, seen=True):
                direction = "unresolved"
                donor_last = (None, "")
            else:
                continue
            key = (e.event_id, t)
            if direction == "unresolved" and key in rows:
                continue     # same event already reported from the other side
            if key in rows and rows[key]["direction"] == "directed":
                continue

            tp = max(tt for tt in grp.prior_times if tt < t)
            prev = grp.states_at(rec, tp, prior)
            prev_state = ("absent" if "absent" in prev else
                          "undetermined" if prev else "not_observed")
            # NA when the segment has no empty site (nothing to count), not 0
            prev_vals = [grp.row.loc[(rec, tp, ff), "empty"] for ff in prior
                         if (rec, tp, ff) in grp.state and
                         pd.notna(grp.row.loc[(rec, tp, ff), "empty"])]
            prev_empty = sum(prev_vals) if prev_vals else None
            tier = "A" if prev_state == "absent" else "B"

            later = [tt for tt in grp.times if tt > t]
            if later:
                nxt = grp.states_at(rec, min(later), observed)
                persistence = ("yes" if "present" in nxt else "partial" if "partial" in nxt else
                               "no" if "absent" in nxt else "undetermined" if nxt else "NA")
            else:
                persistence = "NA"

            r = grp.row.loc[(rec, t, f)]
            rows[key] = {
                "group": g, "t": t, "event_id": e.event_id, "direction": direction, "tier": tier,
                "recipient_contig": seg_contig[rec], "recipient_segment": rec, "recipient_compartment": f,
                "recipient_reads": int(r["n"]),
                "recipient_reads_left": r.get("n_left"), "recipient_reads_right": r.get("n_right"),
                "recipient_contig_breadth": r.get("contig_breadth"),
                "recipient_carrier_frac_depth": r["carrier_frac_depth"],
                "recipient_carrier_frac_empty": r["carrier_frac_empty"],
                "prev_t": tp, "prev_state": prev_state, "prev_empty_reads": prev_empty,
                "donor_contig": seg_contig[don], "donor_segment": don,
                "donor_last_t": donor_last[0], "donor_last_compartment": donor_last[1],
                "persistence": persistence,
            }
    return list(rows.values())


def main(argv=None):
    args = parse_args(argv)
    st = pd.read_csv(args.states, sep="\t", dtype={"group": str, "timepoint": str, "compartment": str},
                     keep_default_na=False, na_values=[""])
    st[["group", "compartment"]] = st[["group", "compartment"]].fillna("")
    seg = pd.read_csv(args.segments, sep="\t", keep_default_na=False)
    seg_contig = dict(zip(seg["segment_id"], seg["contig"]))
    ev = pd.read_csv(args.event_segments, sep="\t")

    compartments = set(st["compartment"])
    observed = set(args.observed.split(",")) - {""} or compartments
    prior = set(args.prior.split(",")) - {""} or compartments
    for name, fs in (("--observed", observed), ("--prior", prior)):
        if fs - compartments:
            sys.exit(f"{name}: compartments not in the sample sheet: {', '.join(sorted(fs - compartments))}")

    rows = []
    for g, gst in st.groupby("group"):
        order = time_order(gst["timepoint"])
        gst = gst.assign(t=gst["timepoint"].map(order))
        label = {v: k for k, v in order.items()}
        for r in call_group(g, gst, ev, observed, prior, seg_contig):
            r["timepoint"] = label[r.pop("t")]
            r["prev_timepoint"] = label.get(r.pop("prev_t"), "")
            r["donor_last_timepoint"] = label.get(r.pop("donor_last_t"), "")
            rows.append(r)

    cols = ["group", "timepoint", "event_id", "direction", "tier",
            "recipient_contig", "recipient_segment", "recipient_compartment", "recipient_reads",
            "recipient_reads_left", "recipient_reads_right", "recipient_contig_breadth",
            "recipient_carrier_frac_depth", "recipient_carrier_frac_empty",
            "prev_timepoint", "prev_state", "prev_empty_reads",
            "donor_contig", "donor_segment", "donor_last_timepoint", "donor_last_compartment",
            "persistence"]
    out = pd.DataFrame(rows, columns=cols)

    if args.contig_taxonomy:
        lin = read_lineages(args.contig_taxonomy)
        empty = [""] * len(RANKS)

        def fmt(c):
            # unclassified contigs get an empty field, not ";;;;;;"
            parts = lin.get(c, empty)
            return ";".join(parts) if any(parts) else ""
        out["recipient_lineage"] = [fmt(c) for c in out["recipient_contig"]]
        out["donor_lineage"] = [fmt(c) for c in out["donor_contig"]]
        rel = [lineage_relation(lin.get(a, empty), lin.get(b, empty))
               for a, b in zip(out["recipient_contig"], out["donor_contig"])]
        out["lineage_relation"] = [x[0] for x in rel]
        out["discordant_rank"] = [x[1] for x in rel]
        out["shared_rank"] = [x[2] for x in rel]
        keep = {"discordant": {"discordant"}, "not_same": {"discordant", "unresolved"},
                "none": {"discordant", "unresolved", "same"}}[args.lineage_filter]
        drop = ~out["lineage_relation"].isin(keep)
        excluded = out[drop].assign(reason="lineage_" + out.loc[drop, "lineage_relation"])
        out = out[~drop]
        print(f"lineage filter ({args.lineage_filter}): removed {int(drop.sum())} events "
              f"({', '.join(f'{k} {v}' for k, v in excluded['lineage_relation'].value_counts().items()) or 'none'})",
              file=sys.stderr)
        if args.excluded:
            excluded.to_csv(args.excluded, sep="\t", index=False)
    elif args.lineage_filter != "none":
        sys.exit("--lineage-filter needs --contig-taxonomy (or use --lineage-filter none)")
    elif args.excluded:
        out.head(0).assign(reason=[]).to_csv(args.excluded, sep="\t", index=False)

    out.to_csv(args.out, sep="\t", index=False)
    print(f"{len(out)} in situ events: " + (", ".join(
        f"{d} {t}: {n}" for (d, t), n in out.groupby(["direction", "tier"]).size().items()) or "none"),
        file=sys.stderr)


if __name__ == "__main__":
    main()
