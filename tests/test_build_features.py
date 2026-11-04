"""Tests for scripts/build_features.py on small simulated contigs with a known shared segment."""

import importlib.util
import os
import random

import pandas as pd
import pysam
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "build_features", os.path.join(HERE, "..", "scripts", "build_features.py"))
bf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bf)

COMP = str.maketrans("ACGT", "TGCA")


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(s, n, rng):
    s = list(s)
    for i in rng.sample(range(len(s)), n):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


@pytest.fixture
def catalogue(tmp_path):
    """
    host1: 3000 own | 4000 shared | 3000 own
    host2: 2500 own | 4000 shared (2 SNPs, reverse complemented) | 3500 own
    host3: 5000 own | 4000 shared, running to the contig end
    Regions are given as an ORF span inside the shared segment (not its true extent).
    """
    rng = random.Random(1)
    shared = rand_seq(4000, rng)
    f1l, f1r = rand_seq(3000, rng), rand_seq(3000, rng)
    f2l, f2r = rand_seq(2500, rng), rand_seq(3500, rng)
    f3l = rand_seq(5000, rng)
    seqs = {
        "host1": f1l + shared + f1r,
        "host2": f2l + revcomp(mutate(shared, 2, rng)) + f2r,
        "host3": f3l + shared,
    }
    fa = tmp_path / "cat.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    pysam.faidx(str(fa))

    # regions as kairos-dd/regions-v2 write them: 1-based start, inclusive stop
    meta = pd.DataFrame([
        {"contig": "host1", "start": 3000 + 301, "stop": 3000 + 3700, "region_idx": "host1_region_1"},
        {"contig": "host2", "start": 2500 + 301, "stop": 2500 + 3700, "region_idx": "host2_region_1"},
        {"contig": "host3", "start": 5000 + 301, "stop": 5000 + 3700, "region_idx": "host3_region_1"},
    ])
    meta.to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    ov = pd.DataFrame([
        {"region_idx_contig1": "host1_region_1", "region_idx_contig2": "host2_region_1"},
        {"region_idx_contig1": "host2_region_1", "region_idx_contig2": "host1_region_1"},  # reciprocal
        {"region_idx_contig1": "host1_region_1", "region_idx_contig2": "host3_region_1"},
    ])
    ov.to_csv(tmp_path / "ov.tsv", sep="\t", index=False)
    return tmp_path, seqs, (f1l, f1r)


def run(tmp_path, *extra):
    out = tmp_path / "features"
    bf.main(["--catalogue", str(tmp_path / "cat.fa"), "--region-metadata", str(tmp_path / "meta.tsv"),
             "--region-overlaps", str(tmp_path / "ov.tsv"), "--outdir", str(out), *extra])
    read = lambda f: pd.read_csv(out / f, sep="\t", keep_default_na=False)
    return out, read("events.tsv"), read("breakpoints.tsv"), read("empty_sites.tsv")


def test_breakpoints_refined_to_true_segment(catalogue):
    tmp_path, _, _ = catalogue
    _, ev, bps, _ = run(tmp_path)
    assert len(ev) == 2                       # reciprocal pair collapsed
    e = ev.set_index("region_idx_2").loc["host2_region_1"]
    assert e["refined"] == True and e["strand"] == -1
    assert abs(e["shared_start_1"] - 3000) <= 3 and abs(e["shared_end_1"] - 7000) <= 3
    assert abs(e["shared_start_2"] - 2500) <= 3 and abs(e["shared_end_2"] - 6500) <= 3
    assert e["identity"] > 0.99


def test_edge_segment_has_one_breakpoint_and_no_empty_site(catalogue):
    tmp_path, _, _ = catalogue
    _, ev, bps, _ = run(tmp_path)
    e = ev.set_index("region_idx_2").loc["host3_region_1"]
    assert e["bp_left_2"] != "" and e["bp_right_2"] == "" and e["empty_site_2"] == ""
    assert e["empty_site_1"] != ""


def test_breakpoints_and_empty_sites_are_shared_between_events(catalogue):
    tmp_path, _, _ = catalogue
    _, ev, bps, es = run(tmp_path)
    # host1's segment appears in both events with the same breakpoints
    assert ev["bp_left_1"].nunique() == 1 and ev["empty_site_1"].nunique() == 1
    assert not bps.duplicated(["contig", "pos", "side"]).any()
    assert len(es) == 2                       # host1 and host2; host3 reaches its contig end


def test_empty_site_sequence_is_flanks_joined(catalogue):
    tmp_path, _, (f1l, f1r) = catalogue
    out, ev, _, es = run(tmp_path, "--junction-flank", "300")
    es_id = ev["empty_site_1"].iloc[0]
    seqs = {r.name: r.sequence for r in pysam.FastxFile(str(out / "junctions.fa"))}
    row = es.set_index("es_id").loc[es_id]
    s = seqs[es_id]
    j = row["junction_pos"]
    # refined ends may sit a base or two off the simulated boundary
    assert j == 300 and len(s) == 600
    assert s[:j][:280] in f1l[-320:] and s[j:][-280:] in f1r[:320]


def test_flanks_exclude_other_regions(tmp_path):
    rng = random.Random(2)
    shared, other = rand_seq(3000, rng), rand_seq(1500, rng)
    a = rand_seq(500, rng) + other + rand_seq(300, rng) + shared + rand_seq(3000, rng)
    b = rand_seq(3000, rng) + shared + rand_seq(3000, rng)
    (tmp_path / "cat.fa").write_text(f">a\n{a}\n>b\n{b}\n")
    pysam.faidx(str(tmp_path / "cat.fa"))
    s_a = 500 + 1500 + 300
    pd.DataFrame([
        {"contig": "a", "start": 501, "stop": 2000, "region_idx": "a_region_1"},
        {"contig": "a", "start": s_a + 201, "stop": s_a + 2800, "region_idx": "a_region_2"},
        {"contig": "b", "start": 3201, "stop": 5800, "region_idx": "b_region_1"},
    ]).to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    pd.DataFrame([{"region_idx_contig1": "a_region_2", "region_idx_contig2": "b_region_1"}]).to_csv(
        tmp_path / "ov.tsv", sep="\t", index=False)
    out, ev, bps, _ = run(tmp_path, "--flank-len", "1000")
    bed = pd.read_csv(out / "flanks.bed", sep="\t", header=None, names=["contig", "s", "e", "bp"])
    left_a = bps[(bps.contig == "a") & (bps.side == "L")].iloc[0]
    pieces = bed[bed.bp == left_a.bp_id]
    # the 1000 bp to the left of a's segment overlap a_region_1 (ends at 2000); only 2000-2300 is kept
    assert pieces["s"].min() >= 2000 and left_a.flank_bp <= 300


def test_contigs_bed_is_contig_minus_regions(catalogue):
    tmp, seqs, _ = catalogue
    out, ev, _, _ = run(tmp)
    bed = pd.read_csv(out / "contigs.bed", sep="\t", header=None, names=["contig", "s", "e"])
    assert set(bed.contig) == {"host1", "host2", "host3"}
    h1 = bed[bed.contig == "host1"].sort_values("s")
    # host1: own sequence either side of the 4000 bp shared segment (refined to 3000-7000)
    # (alignment ends can extend a few bp where random flank bases happen to match)
    assert len(h1) == 2 and h1.s.iloc[0] == 0 and abs(h1.e.iloc[0] - 3000) <= 5
    assert abs(h1.s.iloc[1] - 7000) <= 5 and h1.e.iloc[1] == 10000
    h3 = bed[bed.contig == "host3"]
    assert len(h3) == 1 and h3.s.iloc[0] == 0 and abs(h3.e.iloc[0] - 5000) <= 5
    assert set(ev.region_type_1) == {""}   # test metadata has no region_type column
