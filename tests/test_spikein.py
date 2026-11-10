"""Unit tests for tools/spikein.py (planted in situ events)."""
import os
import random
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import spikein as sp  # noqa: E402

rng = random.Random(3)


def rand(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


GTDB = {
    "g1": "d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 a",
    "g2": "d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 a",     # same species as g1
    "g3": "d__Bacteria;p__P2;c__C2;o__O2;f__F2;g__G2;s__G2 a",
    "g4": "d__Bacteria;p__P3;c__C3;o__O3;f__F3;g__G3;s__G3 a",
    "g5": "d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G4;s__G4 a",     # genus differs from g1
    "g6": "d__Bacteria;p__P4;c__C4;o__O4;f__F4;g__G5;s__G5 a",
}


@pytest.fixture
def pools(tmp_path):
    rows = []
    for g, lin in GTDB.items():
        p = tmp_path / f"{g}.fa"
        sp.write_fasta(str(p), {f"{g}_chr": rand(30000), f"{g}_pl": rand(3000)})
        rows.append({"genome_id": g, "fasta": str(p), "lineage": lin})
    genomes = pd.DataFrame(rows)
    els = []
    for e, typ, tsd, n in (("IS1", "is", 8, 1500), ("TN1", "transposon", 5, 5000),
                           ("PL1", "plasmid", 0, 4000)):
        p = tmp_path / f"{e}.fa"
        sp.write_fasta(str(p), {e: rand(n)})
        els.append({"element_id": e, "fasta": str(p), "type": typ, "tsd": tsd})
    return genomes, pd.DataFrame(els)


def test_insert_duplicates_target_site():
    host, el = "AAAACCCCGGGGTTTT", "XXX"
    seq, s, e = sp.insert(host, el, 4, 3)
    assert seq == "AAAACCC" + "XXX" + "CCCCGGGGTTTT"
    assert seq[s:e] == "XXX"
    assert seq[s - 3:s] == seq[e:e + 3] == "CCC"           # target-site duplication
    seq0, s0, e0 = sp.insert(host, el, 4, 0)
    assert seq0 == "AAAA" + "XXX" + "CCCCGGGGTTTT" and (s0, e0) == (4, 7)


def test_lineages():
    a, b = sp.parse_lineage(GTDB["g1"]), sp.parse_lineage(GTDB["g5"])
    assert a[0] == "Bacteria" and a[-1] == "G1 a"
    assert sp.discordant_rank(a, b) == "g"
    assert sp.discordant_rank(a, sp.parse_lineage(GTDB["g2"])) == "same"
    assert sp.discordant_rank(a, sp.parse_lineage("d__Bacteria;p__P2")) == "p"
    assert sp.discordant_rank(a, sp.parse_lineage("")) == ""
    assert sp.mmseqs_lineage(GTDB["g3"]).startswith("d_Bacteria;p_P2;c_C2")


def test_presence_scenarios():
    k, S = 2, {"SYS"}
    P = lambda sc, t, comp="SYS", f=1.0: sp.presence(sc, t, k, comp, S, 10, f, 8, 0.2)
    assert P("insitu", 1) == (8, 10, 0) and P("insitu", 2) == (8, 0, 10)
    assert P("insitu", 1, "INF") == (8, 0, 0)                      # donor in the influent only
    assert P("carrier_frac", 3, f=0.2) == (8, 8, 2)
    assert P("immigration", 0, "INF") == (8, 0, 10)                 # carrier flows in early
    assert P("preexisting", 0) == (8, 0, 10)
    assert P("rare_bloom", 0) == (8, 10, 0.2) and P("rare_bloom", 2) == (8, 0, 10)
    assert P("donor_late", 1) == (0, 10, 0) and P("donor_late", 2) == (8, 0, 10)


def test_design_respects_ranks_and_unique_genomes(pools):
    genomes, elements = pools
    ev, skipped = sp.design(genomes, elements, ["insitu", "same_lineage", "plasmid"],
                            ["p", "g"], [5], [1.0], 10, 10, 0.2, 1, seed=1)
    lin = {g: sp.parse_lineage(l) for g, l in GTDB.items()}
    for r in ev.itertuples():
        assert sp.discordant_rank(lin[r.donor], lin[r.recipient]) == r.rank
    used = list(ev["donor"]) + list(ev["recipient"])
    assert len(used) == len(set(used))                               # one event per genome
    assert set(ev.loc[ev["scenario"] == "plasmid", "element_id"]) <= {"PL1"}
    assert set(ev.loc[ev["scenario"] != "plasmid", "element_id"]) <= {"IS1", "TN1"}
    assert len(ev) + len(skipped) == 2 + 1 + 2                      # pool runs out -> skipped, not reused
    assert skipped
    assert ev["element_id"].is_unique                                # one event per element by default
    ev2, _ = sp.design(genomes, elements, ["insitu"], ["p", "g"], [1, 2], [1.0], 10, 10, 0.2, 1,
                       seed=1, reuse_elements=True)
    assert len(ev2) > len(set(ev2["element_id"]))
    assert "same_lineage" in set(ev["scenario"])                    # scarce same-species pair allocated first
    assert list(ev["event_id"]) == [f"ev{i:04d}" for i in range(1, len(ev) + 1)]


def test_build_truth_and_composition(pools, tmp_path):
    genomes, elements = pools
    events = pd.DataFrame([
        {"event_id": "ev0001", "scenario": "insitu", "element_id": "TN1", "donor": "g3",
         "recipient": "g1", "rank": "p", "depth": 10, "carrier_frac": 0.5, "donor_depth": 8,
         "rare_depth": 0},
        {"event_id": "ev0002", "scenario": "plasmid", "element_id": "PL1", "donor": "g4",
         "recipient": "g6", "rank": "p", "depth": 10, "carrier_frac": 1.0, "donor_depth": 8,
         "rare_depth": 0},
    ])
    variants, truth, comp, junc, oracle, lins, samples = sp.build(
        events, genomes, elements, ["1", "2", "3"], ["SYS", "INF"], {"SYS"}, 1, 0, 2000,
        "short_pe", 100, False, "G", seed=2)
    t = truth.set_index("event_id")
    el = sp.read_fasta(elements.set_index("element_id").loc["TN1", "fasta"])["TN1"]
    rec = variants["ev0001_rec"]
    ins = rec[int(t.loc["ev0001", "rec_start"]):int(t.loc["ev0001", "rec_end"])]
    assert ins == (el if t.loc["ev0001", "strand"] == "+" else sp.revcomp(el))
    assert len(rec) == 30000 + 5000 + 5                              # chromosome only + element + TSD
    assert variants["ev0001_rec_empty"] in (rec[:int(t.loc["ev0001", "rec_start"])] +
                                            rec[int(t.loc["ev0001", "rec_end"]) + 5:],)
    assert len(junc["ev0001_left"]) == 200 and "ev0002_left" not in junc
    assert set(oracle) == {"ev0001_rec", "ev0001_don", "ev0002_rec", "ev0002_don", "ev0002_plasmid"}
    assert t.loc["ev0001", "expect"] == "called" and t.loc["ev0002", "expect"] == "none"

    c = comp.set_index(["sample", "seq_id"])["depth"]
    assert ("G_1_SYS", "ev0001_rec") not in c and c["G_1_SYS", "ev0001_rec_empty"] == 10
    assert c["G_2_SYS", "ev0001_rec"] == 5 and c["G_2_SYS", "ev0001_rec_empty"] == 5
    assert c["G_1_INF", "ev0001_don"] == 8 and ("G_1_INF", "ev0001_rec_empty") not in c
    assert c["G_3_SYS", "ev0002_plasmid"] == 8 + 10                  # donor's + carrier's copies
    assert list(samples["sample"]) == ["G_1_SYS", "G_1_INF", "G_2_SYS", "G_2_INF", "G_3_SYS", "G_3_INF"]
    assert dict(lins)["ev0001_rec"] == GTDB["g1"]


def test_score_outcomes():
    truth = pd.DataFrame([
        {"event_id": "e1", "scenario": "insitu", "expect": "called", "recipient_contig": "r1",
         "donor_contig": "d1", "plasmid_contig": "", "transfer_t": "3", "element_id": "IS"},
        {"event_id": "e2", "scenario": "immigration", "expect": "none", "recipient_contig": "r2",
         "donor_contig": "d2", "plasmid_contig": "", "transfer_t": "3", "element_id": "TN"},
        {"event_id": "e3", "scenario": "rare_bloom", "expect": "not_A", "recipient_contig": "r3",
         "donor_contig": "d3", "plasmid_contig": "", "transfer_t": "3", "element_id": "IS"},
        {"event_id": "e4", "scenario": "same_lineage", "expect": "excluded", "recipient_contig": "r4",
         "donor_contig": "d4", "plasmid_contig": "", "transfer_t": "3", "element_id": "PX"},
    ])
    ev = pd.DataFrame([
        {"timepoint": "3", "direction": "directed", "tier": "A", "recipient_contig": "r1", "donor_contig": "d1"},
        {"timepoint": "3", "direction": "directed", "tier": "A", "recipient_contig": "r2", "donor_contig": "d2"},
        {"timepoint": "3", "direction": "directed", "tier": "B", "recipient_contig": "r3", "donor_contig": "d1"},
        {"timepoint": "4", "direction": "directed", "tier": "A", "recipient_contig": "x", "donor_contig": "d1"},
    ])
    ex = pd.DataFrame([{"timepoint": "3", "direction": "directed", "tier": "A",
                        "recipient_contig": "r4", "donor_contig": "d4"}])
    res, fp = sp.score(truth, ev, ex)
    r = res.set_index("event_id")
    assert r.loc["e1", "observed"] == "called_A" and r.loc["e1", "correct"]
    assert r.loc["e2", "observed"] == "called_A" and not r.loc["e2", "correct"]   # immigration called
    assert r.loc["e3", "observed"] == "called_B" and r.loc["e3", "correct"]
    assert r.loc["e1", "donor_match"] == "planted"
    assert r.loc["e3", "donor_match"] == "same_element"              # d1 carries the same IS
    assert r.loc["e4", "observed"] == "excluded" and r.loc["e4", "correct"]
    assert len(fp) == 1 and fp.iloc[0]["recipient_contig"] == "x"
