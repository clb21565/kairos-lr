"""Tests for tools/sv_compare.py (matching Sniffles2 SVs to kairos features)."""

import importlib.util
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("sv_compare", os.path.join(HERE, "..", "tools", "sv_compare.py"))
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)

BPS = pd.DataFrame([
    {"bp_id": "BP1", "contig": "c1", "pos": 1000, "side": "L"},
    {"bp_id": "BP2", "contig": "c1", "pos": 5000, "side": "R"},
    {"bp_id": "BP3", "contig": "c2", "pos": 200, "side": "L"},
    {"bp_id": "BP4", "contig": "c2", "pos": 4200, "side": "R"},
])
ES = pd.DataFrame([{"es_id": "ES1", "contig": "c1", "shared_start": 1000, "shared_end": 5000}])
SEGS = pd.DataFrame([
    {"event_id": "EV1", "contig": "c1", "start": 1000, "end": 5000},
    {"event_id": "EV1", "contig": "c2", "start": 200, "end": 4200},
])


def svs(*rows):
    cols = ["sv_id", "svtype", "contig", "start", "end", "svlen", "support", "ins_seq"]
    return pd.DataFrame([dict(zip(cols, r)) for r in rows], columns=cols)


def classes(df):
    return dict(zip(df["sv_id"], df["class"]))


def test_del_matches_empty_site_within_tolerance():
    d = svs(("a", "DEL", "c1", 1002, 4998, 3996, 10, ""),
            ("b", "DEL", "c1", 1100, 5000, 3900, 10, ""))     # left end 100 bp off
    m = sv.match_features(d, BPS, ES, SEGS, tol=50, min_overlap=0.5)
    assert classes(m) == {"a": "kairos_empty_site", "b": "kairos_segment_overlap"}
    assert m.set_index("sv_id").at["a", "kairos_feature"] == "ES1"
    assert m.set_index("sv_id").at["b", "kairos_feature"] == "EV1"


def test_del_on_partner_segment_without_empty_site_is_segment_overlap():
    # c2 has no empty site in ES (e.g. one breakpoint at a contig end), but the segment matches
    m = sv.match_features(svs(("a", "DEL", "c2", 200, 4200, 4000, 5, "")), BPS, ES, SEGS, 50, 0.5)
    assert classes(m) == {"a": "kairos_segment_overlap"}


def test_ins_at_breakpoint_and_novel_ins():
    d = svs(("a", "INS", "c2", 4230, 4230, 800, 5, "ACGT"),
            ("b", "INS", "c3", 100, 100, 800, 5, "ACGT"),
            ("c", "DUP", "c1", 10, 900, 890, 5, ""))
    m = sv.match_features(d, BPS, ES, SEGS, 50, 0.5)
    assert classes(m) == {"a": "kairos_breakpoint", "b": "", "c": "other_type"}
    assert m.set_index("sv_id").at["a", "kairos_feature"] == "BP4"


def hits(*rows):
    cols = ["qname", "qlen", "qs", "qe", "strand", "tname", "tlen", "ts", "te", "nmatch", "alen", "mapq"]
    return pd.DataFrame([dict(zip(cols, r)) for r in rows], columns=cols)


def test_resolve_novel_classes():
    d = svs(("ins_seg", "INS", "c3", 100, 100, 4000, 5, "A" * 4000),
            ("ins_cat", "INS", "c3", 900, 900, 1000, 5, "A" * 1000),
            ("ins_none", "INS", "c3", 2000, 2000, 1000, 5, "A" * 1000),
            ("ins_sym", "INS", "c3", 3000, 3000, 1000, 5, ""),
            ("del_mob", "DEL", "c4", 100, 1100, 1000, 5, ""),
            ("del_uniq", "DEL", "c4", 3000, 4000, 1000, 5, ""))
    d = sv.match_features(d, BPS, ES, SEGS, 50, 0.5)
    h = hits(("ins_seg", 4000, 0, 3990, "+", "c2", 5000, 210, 4200, 3900, 3990, 60),
             ("ins_cat", 1000, 0, 1000, "+", "c9", 9000, 10, 1010, 990, 1000, 60),
             ("ins_none", 1000, 0, 300, "+", "c9", 9000, 10, 310, 300, 300, 60),    # low coverage
             ("del_mob", 1000, 0, 1000, "+", "c4", 9000, 100, 1100, 1000, 1000, 60),  # itself
             ("del_mob", 1000, 0, 1000, "-", "c7", 9000, 50, 1050, 980, 1000, 0),
             ("del_uniq", 1000, 0, 1000, "+", "c4", 9000, 3000, 4000, 1000, 1000, 60))
    r = sv.resolve_novel(d, h, SEGS, min_cov=0.8, min_id=0.9)
    assert classes(r) == {"ins_seg": "novel_ins_kairos_segment", "ins_cat": "novel_ins_catalogue",
                          "ins_none": "novel_ins_no_hit", "ins_sym": "novel_ins_no_seq",
                          "del_mob": "novel_del_elsewhere", "del_uniq": "novel_del_unique"}
    r = r.set_index("sv_id")
    assert r.at["ins_seg", "hit_events"] == "EV1"
    assert r.at["del_mob", "hit_contig"] == "c7"


def test_unsearched_without_catalogue():
    d = sv.match_features(svs(("a", "INS", "c3", 1, 1, 500, 2, "AC"), ("b", "DEL", "c3", 1, 501, 500, 2, "")),
                          BPS, ES, SEGS, 50, 0.5)
    assert classes(sv.resolve_novel(d, None, SEGS, 0.8, 0.9)) == {
        "a": "novel_ins_unsearched", "b": "novel_del_unsearched"}


def test_kairos_counts_join():
    d = sv.match_features(svs(("a", "DEL", "c1", 1000, 5000, 4000, 9, "")), BPS, ES, SEGS, 50, 0.5)
    samples = pd.DataFrame([{"sv_id": "a", "sample": "s1", "gt": "0/1", "dr": 6, "dv": 4, "vaf": 0.4}])
    es_c = pd.DataFrame([{"sample": "s1", "es_id": "ES1", "empty": 4}])
    bp_c = pd.DataFrame([{"sample": "s1", "bp_id": "BP1", "span": 5}, {"sample": "s1", "bp_id": "BP2", "span": 7}])
    out = sv.add_kairos_counts(samples, d, ES, BPS, es_c, bp_c).iloc[0]
    assert (out["es_id"], out["kairos_empty"], out["kairos_span"], out["kairos_frac_empty"]) == ("ES1", 4, 6.0, 0.4)


def test_read_vcf_coordinates(tmp_path):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=c1,length=10000>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="">\n'
        '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="">\n'
        '##FILTER=<ID=GT,Description="">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n'
        '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="">\n'
        '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "c1\t1000\td1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-4000;SUPPORT=9\tGT:DR:DV\t0/1:6:4\t0/0:9:0\n"
        "c1\t6000\ti1\tN\t" + "A" * 400 + "\t60\tPASS\tSVTYPE=INS;SVLEN=400;SUPPORT=3\tGT:DR:DV\t0/1:6:3\t./.:.:.\n"
        "c1\t7000\ti2\tN\t" + "A" * 100 + "\t60\tPASS\tSVTYPE=INS;SVLEN=100;SUPPORT=3\tGT:DR:DV\t0/1:6:3\t0/0:5:0\n"
        "c1\t8000\td2\tN\t<DEL>\t60\tGT\tSVTYPE=DEL;SVLEN=-900;SUPPORT=2\tGT:DR:DV\t0/1:6:2\t0/0:5:0\n")
    s, g = sv.read_vcf(str(vcf), min_svlen=300, all_filters=False)
    assert list(s["sv_id"]) == ["d1", "i1"]                     # short INS and non-PASS dropped
    d1 = s.iloc[0]
    assert (d1["start"], d1["end"]) == (1000, 5000)             # 0-based [start, end)
    assert s.iloc[1]["ins_seq"] == "A" * 400
    assert g.set_index(["sv_id", "sample"]).at[("d1", "s1"), "vaf"] == 0.4
