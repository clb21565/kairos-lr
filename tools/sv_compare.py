#!/usr/bin/env python3
"""
Compare long-read SV calls (Sniffles2) with the kairos features: how many insertions and
deletions does an assembly-independent caller see that kairos already models, and how many
are new?

kairos only sees an event when the catalogue holds both a carrier and a related contig.
Reads can show more: a catalogue contig with an extra segment in some samples (INS), or one
whose segment is missing in some reads (DEL). Run Sniffles2 per long-read sample against the
same catalogue BAMs that assess uses, merge across samples, then run this on the merged VCF:

  sniffles --input <s>.bam --snf <s>.snf --vcf <s>.vcf --mosaic --mosaic-include-germline \\
      --all-contigs --minsvlen 300 --sample-id <s>          # per sample
  sniffles --input <s1>.snf <s2>.snf ... --vcf merged.vcf   # multi-sample
  (--all-contigs is essential: Sniffles2 skips contigs < 1 Mb by default, i.e. nearly all of them)

Each SV gets a class:
  DEL  kairos_empty_site       both ends within --tolerance of a kairos empty site
       kairos_segment_overlap  reciprocal overlap >= --min-overlap with a kairos shared segment
       novel_del_elsewhere     deleted sequence also found on another catalogue contig
                               (looks mobile; kairos has no event for it)
       novel_del_unique        deleted sequence found nowhere else in the catalogue
  INS  kairos_breakpoint       within --tolerance of a kairos breakpoint on that contig
       novel_ins_kairos_segment inserted sequence matches a kairos shared segment
                               (a known mobile segment at a new site)
       novel_ins_catalogue     inserted sequence matches other catalogue sequence
       novel_ins_no_hit        inserted sequence not in the catalogue
       novel_ins_no_seq        symbolic <INS>, no sequence to search
  other types (DUP, INV, BND)  other_type
Without --catalogue, novel SVs are novel_del_unsearched / novel_ins_unsearched.

Outputs (in --outdir):
  sv_calls.tsv      one row per SV: position, length, class, matched kairos feature, best hit
  sv_samples.tsv    one row per SV and sample: Sniffles2 GT, DR, DV, VAF; alt/ref reads recounted
                    from the BAM (alt_reads, ref_reads, vaf_recount; use these: Sniffles2 DR/DV
                    can be wrong for long insertions); for kairos_empty_site SVs also kairos
                    empty reads and junction-spanning reads (with --assess-dir)
  es_recovery.tsv   one row per kairos empty site: matched SV, max kairos empty reads per sample
  summary.txt       counts per class

Usage:
  python tools/sv_compare.py --vcf merged.vcf --features results/detect/catalogue.features \\
      --catalogue results/detect/catalogue.fa --assess-dir results/assess --outdir results/sv
  BAMs are taken from --assess-dir/<sample>/<sample>.bam; others via --bam sample=path.
"""

import argparse
import os
import subprocess
import sys

import pandas as pd
import pysam


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Compare Sniffles2 SV calls with kairos features.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--vcf", required=True, help="Sniffles2 VCF (single or multi-sample)")
    ap.add_argument("--features", required=True, help="<catalogue>.features directory from build_features")
    ap.add_argument("--catalogue", default="", help="catalogue FASTA (indexed); enables the donor search")
    ap.add_argument("--assess-dir", default="", help="results/assess, for per-sample kairos counts")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-svlen", type=int, default=300)
    ap.add_argument("--tolerance", type=int, default=50, help="bp between SV and kairos breakpoints")
    ap.add_argument("--min-overlap", type=float, default=0.5, help="reciprocal overlap with a segment")
    ap.add_argument("--all-filters", action="store_true", help="keep SVs that are not FILTER=PASS")
    ap.add_argument("--hit-min-cov", type=float, default=0.8, help="query coverage of a catalogue hit")
    ap.add_argument("--hit-min-id", type=float, default=0.9, help="identity of a catalogue hit")
    ap.add_argument("--bam", action="append", default=[], metavar="SAMPLE=PATH",
                    help="BAM for a VCF sample not under --assess-dir (repeatable)")
    ap.add_argument("--no-recount", action="store_true", help="skip the per-sample read recount")
    ap.add_argument("--recount-min-mapq", type=int, default=20)
    ap.add_argument("--recount-overhang", type=int, default=100,
                    help="ref reads must align this far either side of the SV")
    ap.add_argument("--recount-min-clip", type=int, default=200,
                    help="a clip of at least this length at a breakpoint counts as alt")
    ap.add_argument("--threads", type=int, default=4)
    return ap.parse_args(argv)


# ── input ────────────────────────────────────────────────────────────────────

def read_vcf(path, min_svlen, all_filters):
    """SVs and per-sample genotypes. Coordinates 0-based: a DEL removes [start, end), an INS
    goes in between bases start-1 and start (start == end)."""
    svs, samples = [], []
    with pysam.VariantFile(path) as vf:
        for rec in vf:
            if not all_filters and "PASS" not in rec.filter.keys():
                continue
            svtype = rec.info.get("SVTYPE", "")
            svlen = abs(int(rec.info.get("SVLEN", 0) or 0))
            if svtype in ("DEL", "INS") and svlen < min_svlen:
                continue
            start = rec.pos                       # VCF POS is the base before the event
            end = start + svlen if svtype == "DEL" else start
            alt = rec.alts[0] if rec.alts else ""
            seq = alt.upper() if svtype == "INS" and not alt.startswith("<") else ""
            svs.append({"sv_id": rec.id, "svtype": svtype, "contig": rec.chrom, "start": start,
                        "end": end, "svlen": svlen, "support": rec.info.get("SUPPORT", ""),
                        "ins_seq": seq})
            for name, s in rec.samples.items():
                gt = s.get("GT")
                gt = "/".join("." if a is None else str(a) for a in gt) if gt else "."
                dr, dv = s.get("DR"), s.get("DV")
                tot = (dr or 0) + (dv or 0)
                samples.append({"sv_id": rec.id, "sample": name, "gt": gt, "dr": dr, "dv": dv,
                                "vaf": round(dv / tot, 3) if tot and dv is not None else float("nan")})
    return pd.DataFrame(svs, columns=["sv_id", "svtype", "contig", "start", "end", "svlen",
                                      "support", "ins_seq"]), \
        pd.DataFrame(samples, columns=["sv_id", "sample", "gt", "dr", "dv", "vaf"])


def read_features(fdir):
    bps = pd.read_csv(os.path.join(fdir, "breakpoints.tsv"), sep="\t")
    es = pd.read_csv(os.path.join(fdir, "empty_sites.tsv"), sep="\t")
    ev = pd.read_csv(os.path.join(fdir, "events.tsv"), sep="\t")
    segs = pd.concat([
        ev[["event_id", f"contig_{i}", f"shared_start_{i}", f"shared_end_{i}"]].set_axis(
            ["event_id", "contig", "start", "end"], axis=1) for i in (1, 2)], ignore_index=True)
    return bps, es, segs


# ── matching ─────────────────────────────────────────────────────────────────

def recip_overlap(a0, a1, b0, b1):
    ov = max(0, min(a1, b1) - max(a0, b0))
    return min(ov / max(1, a1 - a0), ov / max(1, b1 - b0))


def match_features(svs, bps, es, segs, tol, min_overlap):
    """Class for SVs that match kairos features; '' for the rest (resolved by the donor search)."""
    es_by = {c: d for c, d in es.groupby("contig")}
    bp_by = {c: d for c, d in bps.groupby("contig")}
    seg_by = {c: d for c, d in segs.groupby("contig")}
    out = []
    for sv in svs.itertuples(index=False):
        cls, feat, diff = "", "", ""
        if sv.svtype == "DEL":
            d = es_by.get(sv.contig)
            if d is not None:
                dd = (d["shared_start"] - sv.start).abs() + (d["shared_end"] - sv.end).abs()
                ok = ((d["shared_start"] - sv.start).abs() <= tol) & ((d["shared_end"] - sv.end).abs() <= tol)
                if ok.any():
                    i = dd[ok].idxmin()
                    cls, feat, diff = "kairos_empty_site", d.at[i, "es_id"], int(dd[i])
            if not cls and sv.contig in seg_by:
                s = seg_by[sv.contig]
                ro = s.apply(lambda r: recip_overlap(sv.start, sv.end, r["start"], r["end"]), axis=1)
                if (ro >= min_overlap).any():
                    cls = "kairos_segment_overlap"
                    feat = ",".join(sorted(s.loc[ro >= min_overlap, "event_id"].unique()))
        elif sv.svtype == "INS":
            d = bp_by.get(sv.contig)
            if d is not None:
                dist = (d["pos"] - sv.start).abs()
                if (dist <= tol).any():
                    i = dist.idxmin()
                    cls, feat, diff = "kairos_breakpoint", d.at[i, "bp_id"], int(dist[i])
        else:
            cls = "other_type"
        out.append((cls, feat, diff))
    svs = svs.copy()
    svs["class"], svs["kairos_feature"], svs["kairos_bp_diff"] = zip(*out) if out else ([], [], [])
    return svs


# ── donor search ─────────────────────────────────────────────────────────────

def search_queries(svs, catalogue):
    """Query sequences for the unmatched INS (inserted sequence) and DEL (deleted sequence)."""
    q = {}
    fa = pysam.FastaFile(catalogue)
    for sv in svs[svs["class"] == ""].itertuples(index=False):
        if sv.svtype == "INS" and sv.ins_seq:
            q[sv.sv_id] = sv.ins_seq
        elif sv.svtype == "DEL":
            q[sv.sv_id] = fa.fetch(sv.contig, sv.start, sv.end).upper()
    fa.close()
    return q


def run_minimap2(queries, catalogue, outdir, threads):
    qfa = os.path.join(outdir, "sv_queries.fa")
    paf = os.path.join(outdir, "sv_queries.paf")
    with open(qfa, "w") as f:
        for k, s in queries.items():
            f.write(f">{k}\n{s}\n")
    cmd = ["minimap2", "-c", "-x", "map-ont", "--secondary=yes", "-N", "20", "-t", str(threads),
           "--split-prefix", os.path.join(outdir, "mm2split"), catalogue, qfa]
    with open(paf, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    return paf


def read_paf(path):
    cols = ["qname", "qlen", "qs", "qe", "strand", "tname", "tlen", "ts", "te", "nmatch", "alen", "mapq"]
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append(p[:12])
    df = pd.DataFrame(rows, columns=cols)
    for c in cols:
        if c not in ("qname", "strand", "tname"):
            df[c] = df[c].astype(int)
    return df


def resolve_novel(svs, hits, segs, min_cov, min_id):
    """Class and best hit for SVs left unmatched by match_features."""
    svs = svs.copy()
    for c in ("hit_contig", "hit_start", "hit_end", "hit_identity", "hit_qcov", "n_hit_contigs", "hit_events"):
        svs[c] = pd.Series("", index=svs.index, dtype=object)   # mixed str/int (pandas 3 infers str)
    if hits is not None and len(hits):
        hits = hits.assign(qcov=(hits["qe"] - hits["qs"]) / hits["qlen"],
                           identity=hits["nmatch"] / hits["alen"])
        hits = hits[(hits["qcov"] >= min_cov) & (hits["identity"] >= min_id)]
    seg_by = {c: d for c, d in segs.groupby("contig")}
    for i, sv in svs.iterrows():
        if sv["class"]:
            continue
        kind = "del" if sv["svtype"] == "DEL" else "ins"
        if hits is None:
            svs.at[i, "class"] = f"novel_{kind}_unsearched"
            continue
        if kind == "ins" and not sv["ins_seq"]:
            svs.at[i, "class"] = "novel_ins_no_seq"
            continue
        h = hits[hits["qname"] == sv["sv_id"]]
        if kind == "del":   # drop the deleted interval itself
            self_hit = (h["tname"] == sv["contig"]) & (h["ts"] < sv["end"]) & (h["te"] > sv["start"])
            h = h[~self_hit]
        if h.empty:
            svs.at[i, "class"] = "novel_del_unique" if kind == "del" else "novel_ins_no_hit"
            continue
        events = set()
        for t in h.itertuples(index=False):
            s = seg_by.get(t.tname)
            if s is None:
                continue
            ov = (s["end"].clip(upper=t.te) - s["start"].clip(lower=t.ts)).clip(lower=0)
            events.update(s.loc[ov >= 0.5 * (t.te - t.ts), "event_id"])
        best = h.loc[h["nmatch"].idxmax()]
        svs.at[i, "hit_contig"] = best["tname"]
        svs.at[i, "hit_start"] = int(best["ts"])
        svs.at[i, "hit_end"] = int(best["te"])
        svs.at[i, "hit_identity"] = round(best["identity"], 4)
        svs.at[i, "hit_qcov"] = round(best["qcov"], 3)
        svs.at[i, "n_hit_contigs"] = h["tname"].nunique()
        svs.at[i, "hit_events"] = ",".join(sorted(events))
        if kind == "del":
            svs.at[i, "class"] = "novel_del_elsewhere"
        else:
            svs.at[i, "class"] = "novel_ins_kairos_segment" if events else "novel_ins_catalogue"
    return svs


# ── per-sample read recount ──────────────────────────────────────────────────
# Sniffles2 FORMAT DR/DV can be off for long insertions (sim: 18 alt / 11 ref reads, written as
# DV=28, DR=0), so alt and ref reads are recounted from each sample's BAM.

CIGAR_INS, CIGAR_DEL, CIGAR_SOFT, CIGAR_HARD = 1, 2, 4, 5
CONSUMES_REF = {0, 2, 3, 7, 8}


def read_events(aln):
    """(op, ref_pos, length) for I/D ops and for clips (op 'C' at the clip's ref position)."""
    ev, pos = [], aln.reference_start
    cig = aln.cigartuples
    for k, (op, n) in enumerate(cig):
        if op in (CIGAR_SOFT, CIGAR_HARD):
            ev.append(("C", pos, n))
        elif op == CIGAR_INS:
            ev.append(("I", pos, n))
        elif op == CIGAR_DEL:
            ev.append(("D", pos, n))
        if op in CONSUMES_REF:
            pos += n
    return ev


def recount(bam_path, svs, tol, min_mapq, overhang, min_clip):
    """Per SV (INS/DEL): alt = reads with a matching I/D op or a clip at a breakpoint;
    ref = other primary reads aligned across the site with `overhang` bp either side."""
    out = []
    with pysam.AlignmentFile(bam_path) as bam:
        for sv in svs.itertuples(index=False):
            if sv.svtype not in ("INS", "DEL") or sv.contig not in bam.references:
                out.append((sv.sv_id, float("nan"), float("nan")))
                continue
            ends = (sv.start,) if sv.svtype == "INS" else (sv.start, sv.end)
            op = "I" if sv.svtype == "INS" else "D"
            alt, spanning = set(), set()
            for a in bam.fetch(sv.contig, max(0, sv.start - tol - overhang), sv.end + tol + overhang):
                if a.is_unmapped or a.is_secondary or a.mapping_quality < min_mapq:
                    continue
                for o, p, n in read_events(a):
                    if (o == op and abs(p - sv.start) <= tol and n >= 0.7 * sv.svlen) or \
                       (o == "C" and n >= min_clip and any(abs(p - e) <= tol for e in ends)):
                        alt.add(a.query_name)
                        break
                if not a.is_supplementary and a.reference_start <= sv.start - overhang and \
                        a.reference_end >= sv.end + overhang:
                    spanning.add(a.query_name)
            out.append((sv.sv_id, len(alt), len(spanning - alt)))
    return pd.DataFrame(out, columns=["sv_id", "alt_reads", "ref_reads"])


def recount_all(samples, svs, bams, tol, min_mapq, overhang, min_clip):
    rows = []
    for name, path in bams.items():
        r = recount(path, svs, tol, min_mapq, overhang, min_clip)
        r["sample"] = name
        rows.append(r)
    if not rows:
        return samples
    r = pd.concat(rows)
    samples = samples.merge(r, on=["sv_id", "sample"], how="left")
    tot = samples["alt_reads"] + samples["ref_reads"]
    samples["vaf_recount"] = (samples["alt_reads"] / tot).round(3)
    return samples


def bam_paths(assess_dir, samples, extra):
    paths = {}
    for s in samples:
        p = os.path.join(assess_dir, s, f"{s}.bam") if assess_dir else ""
        if p and os.path.exists(p):
            paths[s] = p
    for kv in extra:
        name, path = kv.split("=", 1)
        paths[name] = path
    return {s: p for s, p in paths.items() if s in set(samples)}


# ── per-sample comparison with kairos counts ─────────────────────────────────

def kairos_counts(assess_dir, samples):
    es_rows, bp_rows = [], []
    for s in samples:
        e = os.path.join(assess_dir, s, f"{s}.empty_site_counts.tsv")
        b = os.path.join(assess_dir, s, f"{s}.breakpoint_counts.tsv")
        if os.path.exists(e):
            es_rows.append(pd.read_csv(e, sep="\t"))
        if os.path.exists(b):
            bp_rows.append(pd.read_csv(b, sep="\t"))
    es_c = pd.concat(es_rows) if es_rows else pd.DataFrame(columns=["sample", "es_id", "empty"])
    bp_c = pd.concat(bp_rows) if bp_rows else pd.DataFrame(columns=["sample", "bp_id", "span"])
    return es_c, bp_c


def es_breakpoints(es, bps):
    """es_id -> (left bp_id, right bp_id) on the same contig."""
    key = {(r.contig, r.pos, r.side): r.bp_id for r in bps.itertuples(index=False)}
    return {r.es_id: (key.get((r.contig, r.shared_start, "L")), key.get((r.contig, r.shared_end, "R")))
            for r in es.itertuples(index=False)}


def add_kairos_counts(samples, svs, es, bps, es_c, bp_c):
    es_of = svs.loc[svs["class"] == "kairos_empty_site"].set_index("sv_id")["kairos_feature"].to_dict()
    es_bp = es_breakpoints(es, bps)
    empty = es_c.set_index(["sample", "es_id"])["empty"].to_dict()
    span = bp_c.set_index(["sample", "bp_id"])["span"].to_dict()
    rows = []
    for r in samples.itertuples(index=False):
        e = es_of.get(r.sv_id, "")
        ke, ks = float("nan"), float("nan")
        if e:
            ke = empty.get((r.sample, e), float("nan"))
            sp = [span[(r.sample, b)] for b in es_bp.get(e, ()) if b and (r.sample, b) in span]
            ks = sum(sp) / len(sp) if sp else float("nan")
        rows.append((e, ke, ks))
    samples = samples.copy()
    samples["es_id"], samples["kairos_empty"], samples["kairos_span"] = zip(*rows) if rows else ([], [], [])
    tot = samples["kairos_empty"] + samples["kairos_span"]
    samples["kairos_frac_empty"] = (samples["kairos_empty"] / tot).round(3)
    return samples


def es_recovery(es, svs, es_c):
    m = svs[svs["class"] == "kairos_empty_site"].groupby("kairos_feature")["sv_id"].apply(",".join)
    out = es[["es_id", "contig", "shared_start", "shared_end"]].copy()
    out["sv_id"] = out["es_id"].map(m).fillna("")
    if len(es_c):
        mx = es_c.groupby("es_id")["empty"].max()
        out["kairos_max_empty"] = out["es_id"].map(mx)
    return out


def summary(svs, rec, path):
    lines = [f"SVs: {len(svs)}"]
    for (t, c), n in svs.groupby(["svtype", "class"]).size().items():
        lines.append(f"  {t}\t{c}\t{n}")
    lines.append(f"kairos empty sites: {len(rec)}; with a matching DEL: {(rec['sv_id'] != '').sum()}")
    if "kairos_max_empty" in rec:
        seen = rec[rec["kairos_max_empty"] > 0]
        lines.append(f"  with kairos empty reads > 0 in a long-read sample: {len(seen)}; "
                     f"of those with a matching DEL: {(seen['sv_id'] != '').sum()}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main(argv=None):
    a = parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    svs, samples = read_vcf(a.vcf, a.min_svlen, a.all_filters)
    bps, es, segs = read_features(a.features)
    svs = match_features(svs, bps, es, segs, a.tolerance, a.min_overlap)
    hits = None
    if a.catalogue:
        q = search_queries(svs, a.catalogue)
        hits = read_paf(run_minimap2(q, a.catalogue, a.outdir, a.threads)) if q else pd.DataFrame(
            columns=["qname", "qlen", "qs", "qe", "strand", "tname", "tlen", "ts", "te", "nmatch", "alen", "mapq"])
    svs = resolve_novel(svs, hits, segs, a.hit_min_cov, a.hit_min_id)

    es_c, bp_c = kairos_counts(a.assess_dir, samples["sample"].unique()) if a.assess_dir else \
        (pd.DataFrame(columns=["sample", "es_id", "empty"]), pd.DataFrame(columns=["sample", "bp_id", "span"]))
    samples = add_kairos_counts(samples, svs, es, bps, es_c, bp_c)
    if not a.no_recount:
        bams = bam_paths(a.assess_dir, samples["sample"].unique(), a.bam)
        samples = recount_all(samples, svs, bams, a.tolerance, a.recount_min_mapq,
                              a.recount_overhang, a.recount_min_clip)
    rec = es_recovery(es, svs, es_c[es_c["sample"].isin(samples["sample"].unique())])

    svs.drop(columns="ins_seq").to_csv(os.path.join(a.outdir, "sv_calls.tsv"), sep="\t", index=False)
    samples.to_csv(os.path.join(a.outdir, "sv_samples.tsv"), sep="\t", index=False)
    rec.to_csv(os.path.join(a.outdir, "es_recovery.tsv"), sep="\t", index=False)
    summary(svs, rec, os.path.join(a.outdir, "summary.txt"))


if __name__ == "__main__":
    sys.exit(main())
