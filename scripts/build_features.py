#!/usr/bin/env python3
"""
kairos build-features | junction, flank and empty-site features for every pair of shared regions

For each pair of regions that share ORFs (regions-v2.py region overlaps), the two
regions are aligned (minimap2 asm5 via mappy) with some sequence either side, so the
shared segment is placed where the two contigs actually stop matching rather than at
the ORF boundaries. On each contig that gives up to two breakpoints, where shared
sequence meets the contig's own sequence:

  contig:  [ own flank ][======= shared =======][ own flank ]
                        ^ left breakpoint       ^ right breakpoint

Outputs (in --outdir):
  events.tsv       one row per region pair: shared interval and breakpoint/empty-site IDs on each contig
  breakpoints.tsv  one row per breakpoint: contig, position, side, junction window coordinates
  flanks.bed       the contig's own sequence next to each breakpoint (other regions removed);
                   its read depth sets the detection limit
  empty_sites.tsv  one row per shared segment with both breakpoints: the contig with the segment removed
  junctions.fa     junction windows (BP*) and empty-site sequences (ES*), for the second-pass remapping
  contigs.bed      every event contig's own sequence: the whole contig minus all of its regions and
                   shared segments; its depth and breadth say whether the contig itself is present

See docs/features.md for what each file and column means.

All coordinates are 0-based, half-open. A breakpoint at position p lies between bases p-1 and p.
"""

import argparse
import os
import sys

import mappy
import pandas as pd
import pysam


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build junction, flank and empty-site features for shared regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--catalogue", required=True, help="catalogue FASTA (indexed, .fai next to it)")
    parser.add_argument("--region-metadata", required=True, help="kairos_region_metadata.tsv")
    parser.add_argument("--region-overlaps", required=True, help="kairos_region_overlaps.tsv")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window", type=int, default=5000,
                        help="bp either side of each region included in the pairwise alignment")
    parser.add_argument("--min-identity", type=float, default=0.95,
                        help="minimum identity (matches / alignment length) of the shared segment")
    parser.add_argument("--min-flank", type=int, default=200,
                        help="a breakpoint needs at least this much own sequence beyond it; "
                             "otherwise the segment is treated as reaching the contig end")
    parser.add_argument("--flank-len", type=int, default=1000,
                        help="length of own sequence next to each breakpoint used for flank depth")
    parser.add_argument("--junction-flank", type=int, default=500,
                        help="bp either side of the junction in junction windows and empty-site sequences")
    parser.add_argument("--snap", type=int, default=10,
                        help="breakpoints on the same contig and side within this distance are merged "
                             "(alignment ends vary by a few bp between partners)")
    return parser.parse_args(argv)


def load_regions(path):
    regions = pd.read_csv(path, sep="\t")
    # kairos-dd writes Prodigal's 1-based ORF starts into the BED start column
    regions["s0"] = regions["start"] - 1
    regions["e0"] = regions["stop"]
    return regions.set_index("region_idx")


def load_pairs(path):
    overlaps = pd.read_csv(path, sep="\t")
    pairs = overlaps[["region_idx_contig1", "region_idx_contig2"]].drop_duplicates()
    pairs = pairs[pairs["region_idx_contig1"] != pairs["region_idx_contig2"]]
    key = pairs.apply(lambda r: tuple(sorted(r)), axis=1)
    pairs = pd.DataFrame(key.drop_duplicates().tolist(), columns=["r1", "r2"])
    return pairs.sort_values(["r1", "r2"]).reset_index(drop=True)


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def subtract(start, end, intervals):
    """[start, end) minus a list of (s, e) intervals → list of (s, e)."""
    pieces = [(start, end)]
    for s, e in sorted(intervals):
        nxt = []
        for p0, p1 in pieces:
            if e <= p0 or s >= p1:
                nxt.append((p0, p1))
                continue
            if s > p0:
                nxt.append((p0, s))
            if e < p1:
                nxt.append((e, p1))
        pieces = nxt
    return pieces


def best_hit(aligner, query, w1, w2, reg1, reg2, min_identity):
    """Alignment of region 2's window onto region 1's window that covers both regions best."""
    best, best_score = None, 0
    for h in aligner.map(query):
        if h.blen == 0 or h.mlen / h.blen < min_identity:
            continue
        t0, t1 = w1 + h.r_st, w1 + h.r_en
        q0, q1 = w2 + h.q_st, w2 + h.q_en
        o1 = overlap(t0, t1, *reg1)
        o2 = overlap(q0, q1, *reg2)
        if o1 == 0 or o2 == 0:
            continue
        if o1 + o2 > best_score:
            best, best_score = h, o1 + o2
    return best


class FeatureStore:
    """Deduplicates breakpoints and empty sites shared by several events."""

    def __init__(self, fasta, regions, contig_len, args):
        self.fa = fasta
        self.args = args
        self.contig_len = contig_len
        self.region_iv = regions.groupby("contig")[["s0", "e0"]].apply(
            lambda d: list(zip(d["s0"], d["e0"]))).to_dict()
        self.bp = {}       # (contig, pos, side) -> row
        self.bp_by_side = {}   # (contig, side) -> rows, for snapping
        self.es = {}       # (contig, start, end) -> row

    def breakpoint(self, contig, pos, side, shared):
        """
        side L: shared sequence starts at pos (own flank to the left); side R: it ends at pos.
        Returns the breakpoint row; an existing one within --snap bp is reused.
        """
        a = self.args
        near = [r for r in self.bp_by_side.get((contig, side), []) if abs(r["pos"] - pos) <= a.snap]
        if near:
            return min(near, key=lambda r: abs(r["pos"] - pos))
        L = self.contig_len[contig]
        if side == "L":
            f0, f1 = max(0, pos - a.flank_len), pos
        else:
            f0, f1 = pos, min(L, pos + a.flank_len)
        # own flank, minus any other region on this contig (those reads may map elsewhere)
        others = [iv for iv in self.region_iv.get(contig, []) if overlap(*iv, *shared) == 0]
        flanks = subtract(f0, f1, others)
        row = {
            "bp_id": f"BP{len(self.bp) + 1:07d}",
            "contig": contig, "pos": pos, "side": side,
            "window_start": max(0, pos - a.junction_flank),
            "window_end": min(L, pos + a.junction_flank),
            "flank_bp": sum(e - s for s, e in flanks),
            "_flanks": flanks,
        }
        self.bp[(contig, pos, side)] = row
        self.bp_by_side.setdefault((contig, side), []).append(row)
        return row

    def empty_site(self, contig, start, end):
        key = (contig, start, end)
        if key in self.es:
            return self.es[key]["es_id"]
        j = self.args.junction_flank
        L = self.contig_len[contig]
        left = (max(0, start - j), start)
        right = (end, min(L, end + j))
        row = {
            "es_id": f"ES{len(self.es) + 1:07d}",
            "contig": contig, "shared_start": start, "shared_end": end,
            "left_start": left[0], "right_end": right[1],
            "junction_pos": left[1] - left[0],   # position of the junction in the ES sequence
        }
        self.es[key] = row
        return row["es_id"]

    def contig_side(self, contig, s0, s1):
        """Breakpoints and empty site for a shared segment [s0, s1) on one contig."""
        L = self.contig_len[contig]
        mf = self.args.min_flank
        left = self.breakpoint(contig, s0, "L", (s0, s1)) if s0 >= mf else None
        right = self.breakpoint(contig, s1, "R", (s0, s1)) if L - s1 >= mf else None
        empty = self.empty_site(contig, left["pos"], right["pos"]) if left and right else ""
        return (left["bp_id"] if left else ""), (right["bp_id"] if right else ""), empty

    def write(self, outdir):
        bps = pd.DataFrame(self.bp.values(), columns=[
            "bp_id", "contig", "pos", "side", "window_start", "window_end", "flank_bp", "_flanks"])
        bps.drop(columns="_flanks").to_csv(os.path.join(outdir, "breakpoints.tsv"), sep="\t", index=False)
        with open(os.path.join(outdir, "flanks.bed"), "w") as out:
            for r in self.bp.values():
                for s, e in r["_flanks"]:
                    out.write(f"{r['contig']}\t{s}\t{e}\t{r['bp_id']}\n")
        es = pd.DataFrame(self.es.values(), columns=[
            "es_id", "contig", "shared_start", "shared_end", "left_start", "right_end", "junction_pos"])
        es.to_csv(os.path.join(outdir, "empty_sites.tsv"), sep="\t", index=False)
        with open(os.path.join(outdir, "junctions.fa"), "w") as out:
            for r in self.bp.values():
                seq = self.fa.fetch(r["contig"], r["window_start"], r["window_end"])
                out.write(f">{r['bp_id']}\n{seq}\n")
            for r in self.es.values():
                seq = (self.fa.fetch(r["contig"], r["left_start"], r["shared_start"]) +
                       self.fa.fetch(r["contig"], r["shared_end"], r["right_end"]))
                out.write(f">{r['es_id']}\n{seq}\n")


def write_contig_own(events, store, path):
    """BED of each event contig's own sequence: [0, L) minus every region and shared segment on it."""
    remove = {c: list(iv) for c, iv in store.region_iv.items()}
    for i in (1, 2):
        if events.empty:
            break
        for c, s, e in zip(events[f"contig_{i}"], events[f"shared_start_{i}"], events[f"shared_end_{i}"]):
            remove.setdefault(c, []).append((s, e))
    contigs = sorted(set(events[["contig_1", "contig_2"]].values.ravel())) if len(events) else []
    with open(path, "w") as out:
        for c in contigs:
            for s, e in subtract(0, store.contig_len[c], remove.get(c, [])):
                out.write(f"{c}\t{s}\t{e}\n")


def build_events(pairs, regions, fasta, store, args):
    events = []
    aligner, aligner_region, w1 = None, None, 0
    for i, (r1, r2) in enumerate(pairs.itertuples(index=False), start=1):
        g1, g2 = regions.loc[r1], regions.loc[r2]
        c1, c2 = g1["contig"], g2["contig"]
        reg1, reg2 = (g1["s0"], g1["e0"]), (g2["s0"], g2["e0"])
        if aligner_region != r1:   # pairs are sorted by r1, so each target is indexed once
            w1 = max(0, reg1[0] - args.window)
            seq1 = fasta.fetch(c1, w1, min(store.contig_len[c1], reg1[1] + args.window))
            aligner = mappy.Aligner(seq=seq1, preset="asm5")
            aligner_region = r1
        w2 = max(0, reg2[0] - args.window)
        seq2 = fasta.fetch(c2, w2, min(store.contig_len[c2], reg2[1] + args.window))
        h = best_hit(aligner, seq2, w1, w2, reg1, reg2, args.min_identity)

        if h is not None:
            sh1 = (w1 + h.r_st, w1 + h.r_en)
            sh2 = (w2 + h.q_st, w2 + h.q_en)
            strand, identity, refined = h.strand, round(h.mlen / h.blen, 4), True
        else:
            # no alignment passing the filters: keep the ORF-based region coordinates
            sh1, sh2, strand, identity, refined = reg1, reg2, 0, float("nan"), False

        side1 = store.contig_side(c1, *sh1)
        side2 = store.contig_side(c2, *sh2)
        events.append({
            "event_id": f"EV{i:07d}",
            "region_idx_1": r1, "contig_1": c1, "shared_start_1": sh1[0], "shared_end_1": sh1[1],
            "bp_left_1": side1[0], "bp_right_1": side1[1], "empty_site_1": side1[2],
            "region_idx_2": r2, "contig_2": c2, "shared_start_2": sh2[0], "shared_end_2": sh2[1],
            "bp_left_2": side2[0], "bp_right_2": side2[1], "empty_site_2": side2[2],
            "strand": strand, "identity": identity, "shared_len": sh1[1] - sh1[0],
            "refined": refined,
            "region_type_1": g1.get("region_type", ""), "region_type_2": g2.get("region_type", ""),
        })
    return pd.DataFrame(events)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    regions = load_regions(args.region_metadata)
    pairs = load_pairs(args.region_overlaps)
    fasta = pysam.FastaFile(args.catalogue)
    contig_len = dict(zip(fasta.references, fasta.lengths))
    print(f"{len(regions)} regions, {len(pairs)} region pairs", file=sys.stderr)

    store = FeatureStore(fasta, regions, contig_len, args)
    events = build_events(pairs, regions, fasta, store, args)
    if events.empty:
        events = pd.DataFrame(columns=[
            "event_id", "region_idx_1", "contig_1", "shared_start_1", "shared_end_1",
            "bp_left_1", "bp_right_1", "empty_site_1", "region_idx_2", "contig_2",
            "shared_start_2", "shared_end_2", "bp_left_2", "bp_right_2", "empty_site_2",
            "strand", "identity", "shared_len", "refined", "region_type_1", "region_type_2"])
    events.to_csv(os.path.join(args.outdir, "events.tsv"), sep="\t", index=False)
    store.write(args.outdir)
    write_contig_own(events, store, os.path.join(args.outdir, "contigs.bed"))

    n_ref = int(events["refined"].sum()) if len(events) else 0
    print(f"{len(events)} events ({n_ref} refined by alignment), "
          f"{len(store.bp)} breakpoints, {len(store.es)} empty sites", file=sys.stderr)


if __name__ == "__main__":
    main()
