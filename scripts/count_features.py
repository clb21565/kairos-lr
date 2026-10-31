#!/usr/bin/env python3
"""
kairos count-features | read counts at breakpoints, flanks and empty sites for one sample

Input is the sample's reads mapped to the whole catalogue (sorted, indexed BAM) and the
features from build_features.py. For each breakpoint:

  span         reads aligned across the breakpoint with at least --overhang bp on both sides
  bridge       read pairs with one mate in the contig's own flank and the other in the shared
               segment, neither spanning the breakpoint (paired-end only)
  flank_depth  mean depth over the breakpoint's own flank (sets the detection limit)

and for each empty site (contig with its shared segment removed):

  empty_split  reads whose alignment jumps from the left flank straight to the right flank,
               as one large deletion or as a split (supplementary, SA tag) alignment
               (long reads; short reads from aligners that report SA, e.g. bwa-mem2)
  empty_pairs  read pairs with one mate in each flank and an insert that is the segment
               length longer than normal (paired-end only)

Junction evidence (span, bridge, empty_*) uses --junction-min-mapq; flank depth uses --min-mapq.
Junction reads lie partly in shared sequence, so aligners give them modest MAPQ (strobealign:
~15 for a read with 50 bp of own flank) even when the placement is clear; the junction threshold
only has to exclude ties (MAPQ 0). Only primary alignments (and, for empty sites, supplementary
pieces) are used; secondary, duplicate and QC-fail reads are skipped.
Coordinates are 0-based, half-open, as in build_features.py.
"""

import argparse
import re
import sys

import pandas as pd
import pysam

SHORT_TYPES = {"short_pe", "short_se"}
REF_OPS = {0, 2, 3, 7, 8}            # M, D, N, =, X consume the reference
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Count reads at breakpoints, flanks and empty sites for one sample.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--read-type", required=True, choices=["short_pe", "short_se", "ont", "hifi"])
    parser.add_argument("--breakpoints", required=True, help="breakpoints.tsv from build_features.py")
    parser.add_argument("--flanks", required=True, help="flanks.bed from build_features.py")
    parser.add_argument("--empty-sites", required=True, help="empty_sites.tsv from build_features.py")
    parser.add_argument("--overhang", type=int, default=25,
                        help="aligned bp needed on each side of a breakpoint / empty-site junction")
    parser.add_argument("--min-mapq", type=int, default=20, help="MAPQ for flank depth")
    parser.add_argument("--junction-min-mapq", type=int, default=1,
                        help="MAPQ for junction evidence (spanning reads, bridging pairs, empty sites)")
    parser.add_argument("--max-insert", type=int, default=1000,
                        help="largest insert size for a bridging read pair")
    parser.add_argument("--max-gap", type=int, default=10,
                        help="deletions up to this size do not break a spanning alignment")
    parser.add_argument("--empty-tolerance", type=int, default=50,
                        help="how far (bp) an alignment break may be from the empty-site junction")
    parser.add_argument("--out-breakpoints", required=True)
    parser.add_argument("--out-empty-sites", required=True)
    return parser.parse_args(argv)


def usable(read, min_mapq, supplementary=False):
    if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail:
        return False
    if read.is_supplementary and not supplementary:
        return False
    return read.mapping_quality >= min_mapq


def merged_blocks(read, max_gap):
    """Aligned reference blocks, with deletions of up to max_gap bp bridged."""
    out = []
    for s, e in read.get_blocks():
        if out and s - out[-1][1] <= max_gap:
            out[-1][1] = e
        else:
            out.append([s, e])
    return out


def spans(read, pos, overhang, max_gap):
    return any(s <= pos - overhang and e >= pos + overhang for s, e in merged_blocks(read, max_gap))


def ref_len(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")


def mate_end(read):
    """Reference end of the mate: from the MC tag when present, else assume an ungapped mate."""
    if read.has_tag("MC"):
        return read.next_reference_start + ref_len(read.get_tag("MC"))
    return read.next_reference_start + read.query_length


class Counter:
    def __init__(self, bam, args):
        self.bam = bam
        self.a = args
        self.paired = args.read_type == "short_pe"
        self.read_lengths = []

    def note_length(self, read):
        if len(self.read_lengths) < 100000:
            self.read_lengths.append(read.infer_read_length() or read.query_length)

    def breakpoint(self, contig, pos, side):
        a = self.a
        span = 0
        for read in self.bam.fetch(contig, max(0, pos - 1), pos + 1):
            if usable(read, a.junction_min_mapq) and spans(read, pos, a.overhang, a.max_gap):
                span += 1
        bridge = self.bridges(contig, pos, side) if self.paired else 0
        return span, bridge

    def bridges(self, contig, pos, side):
        """Pairs with the flank mate on the contig's own side and the other mate in the shared segment."""
        a = self.a
        n = 0
        if side == "L":    # own flank left of pos, shared segment from pos on
            start, end = max(0, pos - a.max_insert), pos
        else:              # shared segment ends at pos, own flank from pos on
            start, end = pos, pos + a.max_insert
        for read in self.bam.fetch(contig, start, end):
            if not usable(read, a.junction_min_mapq) or not read.is_paired or read.mate_is_unmapped:
                continue
            if read.next_reference_id != read.reference_id or abs(read.template_length) > a.max_insert:
                continue
            if side == "L":
                ok = (not read.is_reverse and read.mate_is_reverse and
                      read.reference_end <= pos and read.next_reference_start >= pos)
            else:
                ok = (read.is_reverse and not read.mate_is_reverse and
                      read.reference_start >= pos and mate_end(read) <= pos)
            n += ok
        return n

    def flank_depth(self, intervals):
        a = self.a
        total = sum(e - s for _, s, e in intervals)
        if total == 0:
            return float("nan")
        bases = 0
        for contig, s, e in intervals:
            for read in self.bam.fetch(contig, s, e):
                if not usable(read, a.min_mapq):
                    continue
                self.note_length(read)
                for b0, b1 in read.get_blocks():
                    bases += max(0, min(b1, e) - max(b0, s))
        return bases / total

    def empty_site(self, contig, s, e):
        """Evidence for the contig without its segment [s, e): (split/deletion reads, read pairs)."""
        pairs = self.empty_pairs(contig, s, e) if self.paired else 0
        return self.empty_split(contig, s, e), pairs

    def empty_pairs(self, contig, s, e):
        """Pairs with the forward mate in the left flank and the reverse mate in the right flank."""
        a = self.a
        tol, seg = a.empty_tolerance, e - s
        n = 0
        for read in self.bam.fetch(contig, max(0, s - a.max_insert), s + tol):
            if (not usable(read, a.junction_min_mapq) or not read.is_paired or read.mate_is_unmapped
                    or read.is_reverse or not read.mate_is_reverse
                    or read.next_reference_id != read.reference_id):
                continue
            if read.reference_end <= s + tol and read.next_reference_start >= e - tol \
                    and 0 < read.template_length - seg <= a.max_insert:
                n += 1
        return n

    def empty_split(self, contig, s, e):
        """Reads that jump from just before s to just after e (deletion or split alignment)."""
        a = self.a
        tol, ov = a.empty_tolerance, a.overhang
        names = set()
        for read in self.bam.fetch(contig, max(0, s - tol - 1), s + tol):
            if not usable(read, a.junction_min_mapq, supplementary=True):
                continue
            # one alignment with a large deletion over the segment
            ref = read.reference_start
            for op, length in read.cigartuples:
                if op in (2, 3) and abs(ref - s) <= tol and abs(ref + length - e) <= tol:
                    if read.reference_start <= ref - ov and read.reference_end >= ref + length + ov:
                        names.add(read.query_name)
                if op in REF_OPS:
                    ref += length
            # split alignment: this piece ends at s, an SA piece starts at e on the same strand
            if (abs(read.reference_end - s) <= tol and read.reference_start <= s - ov
                    and read.has_tag("SA")):
                strand = "-" if read.is_reverse else "+"
                for sa in read.get_tag("SA").rstrip(";").split(";"):
                    rname, spos, sstrand, cigar, mapq, _ = sa.split(",")
                    sa_start = int(spos) - 1
                    if (rname == contig and sstrand == strand and int(mapq) >= a.junction_min_mapq
                            and abs(sa_start - e) <= tol and ref_len(cigar) >= ov):
                        names.add(read.query_name)
        return len(names)


def main(argv=None):
    args = parse_args(argv)
    bps = pd.read_csv(args.breakpoints, sep="\t")
    flanks = pd.read_csv(args.flanks, sep="\t", header=None, names=["contig", "start", "end", "bp_id"])
    flank_iv = {bp: list(zip(d["contig"], d["start"], d["end"])) for bp, d in flanks.groupby("bp_id")}
    es = pd.read_csv(args.empty_sites, sep="\t")

    bam = pysam.AlignmentFile(args.bam)
    counter = Counter(bam, args)

    rows = []
    for r in bps.itertuples(index=False):
        span, bridge = counter.breakpoint(r.contig, r.pos, r.side)
        rows.append({"sample": args.sample, "read_type": args.read_type, "bp_id": r.bp_id,
                     "span": span, "bridge": bridge,
                     "flank_depth": counter.flank_depth(flank_iv.get(r.bp_id, [])),
                     "flank_bp": r.flank_bp})
    out = pd.DataFrame(rows, columns=["sample", "read_type", "bp_id", "span", "bridge",
                                      "flank_depth", "flank_bp"])
    lengths = counter.read_lengths
    out["mean_read_len"] = round(sum(lengths) / len(lengths), 1) if lengths else float("nan")
    out["overhang"] = args.overhang
    out.to_csv(args.out_breakpoints, sep="\t", index=False, float_format="%.3f")

    erows = []
    for r in es.itertuples(index=False):
        split, pairs = counter.empty_site(r.contig, r.shared_start, r.shared_end)
        erows.append({"sample": args.sample, "read_type": args.read_type, "es_id": r.es_id,
                      "empty_split": split, "empty_pairs": pairs})
    pd.DataFrame(erows, columns=["sample", "read_type", "es_id", "empty_split", "empty_pairs"]).to_csv(
        args.out_empty_sites, sep="\t", index=False)

    print(f"{args.sample}: {len(out)} breakpoints, {len(erows)} empty sites, "
          f"mean read length {out['mean_read_len'].iloc[0] if len(out) else 'NA'}", file=sys.stderr)


if __name__ == "__main__":
    main()
