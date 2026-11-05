"""Tests for scripts/detect_state.py (detection limit) and scripts/insitu.py (in situ rules)."""

import importlib.util
import os

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ds = load("detect_state")
ins = load("insitu")


# ── detection limit ──────────────────────────────────────────────────────────

def test_poisson_threshold_for_k2():
    # absent needs P(N < 2 | lambda) <= 0.05, i.e. lambda >= ~4.74
    assert ds.poisson_cdf_below(2, 4.7) > 0.05
    assert ds.poisson_cdf_below(2, 4.8) <= 0.05
    assert ds.poisson_cdf_below(1, 3.0) == pytest.approx(0.0498, abs=1e-3)


# ── in situ ─────────────────────────────────────────────────────────────────
# Two contigs: R (segment SR) and D (segment SD), one event EV1 between them.

SEGMENTS = pd.DataFrame({"segment_id": ["SR", "SD"], "contig": ["R", "D"]})
EVENTS = pd.DataFrame({"event_id": ["EV1"], "segment_1": ["SR"], "segment_2": ["SD"]})


def states(rows):
    """rows: (timepoint, compartment, segment, state[, empty])"""
    out = []
    for r in rows:
        t, f, seg, state = r[:4]
        out.append({"group": "G", "timepoint": str(t), "compartment": f, "segment_id": seg,
                    "state": state, "n": 5 if state == "present" else 0, "n_span": 0, "lam": 10,
                    "p_miss": 0, "empty": r[4] if len(r) > 4 else float("nan"), "flank_depth": 10,
                    "carrier_frac_depth": 0.5, "carrier_frac_empty": "", "n_breakpoints": 2})
    return pd.DataFrame(out)


def run_insitu(tmp_path, st, observed="", prior="", lineages=None, lineage_filter="none"):
    st.to_csv(tmp_path / "states.tsv", sep="\t", index=False)
    SEGMENTS.to_csv(tmp_path / "segments.tsv", sep="\t", index=False)
    EVENTS.to_csv(tmp_path / "ev.tsv", sep="\t", index=False)
    extra = []
    if lineages:
        # mmseqs createtsv layout: contig, taxid, rank, name, ..., lineage (last column)
        (tmp_path / "tax.tsv").write_text("".join(
            f"{c}\t1\tspecies\tx\t{lin}\n" for c, lin in lineages.items()))
        extra = ["--contig-taxonomy", str(tmp_path / "tax.tsv")]
    ins.main(["--states", str(tmp_path / "states.tsv"), "--segments", str(tmp_path / "segments.tsv"),
              "--event-segments", str(tmp_path / "ev.tsv"), "--observed", observed, "--prior", prior,
              "--lineage-filter", lineage_filter, *extra,
              "--out", str(tmp_path / "out.tsv"), "--excluded", str(tmp_path / "excluded.tsv")])
    return pd.read_csv(tmp_path / "out.tsv", sep="\t", dtype=str, keep_default_na=False)


def test_directed_tier_a(tmp_path):
    out = run_insitu(tmp_path, states([
        (1, "X", "SD", "present"), (1, "X", "SR", "absent", 12),
        (2, "X", "SD", "present"), (2, "X", "SR", "present"),
        (3, "X", "SD", "present"), (3, "X", "SR", "present"),
    ]))
    assert len(out) == 1
    e = out.iloc[0]
    assert (e.recipient_contig, e.donor_contig, e.timepoint) == ("R", "D", "2")
    assert (e.direction, e.tier, e.prev_state, e.persistence) == ("directed", "A", "absent", "yes")
    assert e.donor_last_timepoint == "1" and e.prev_empty_reads == "12.0"


def test_tier_b_when_previous_undetermined(tmp_path):
    out = run_insitu(tmp_path, states([
        (1, "X", "SD", "present"), (1, "X", "SR", "undetermined"),
        (2, "X", "SD", "present"), (2, "X", "SR", "present"),
    ]))
    assert out.iloc[0].tier == "B" and out.iloc[0].persistence == "NA"


def test_unresolved_when_both_appear_together(tmp_path):
    out = run_insitu(tmp_path, states([
        (1, "X", "SD", "absent"), (1, "X", "SR", "absent"),
        (2, "X", "SD", "present"), (2, "X", "SR", "present"),
    ]))
    assert len(out) == 1 and out.iloc[0].direction == "unresolved"


def test_no_event_without_earlier_donor(tmp_path):
    # D first appears after R: R is the donor, D the recipient
    out = run_insitu(tmp_path, states([
        (1, "X", "SR", "present"), (1, "X", "SD", "absent"),
        (2, "X", "SR", "present"), (2, "X", "SD", "present"),
    ]))
    assert len(out) == 1 and out.iloc[0].recipient_contig == "D"


def test_recipient_seen_before_is_not_new(tmp_path):
    out = run_insitu(tmp_path, states([
        (1, "X", "SD", "present"), (1, "X", "SR", "present"),
        (2, "X", "SD", "present"), (2, "X", "SR", "present"),
    ]))
    assert len(out) == 0


def test_same_day_import_excluded(tmp_path):
    # recipient looked for in EFF; at t2 it is also in INF -> could be imported, no event
    rows = [(1, f, "SD", "present") for f in ("INF", "EFF")] + \
           [(1, f, "SR", "absent") for f in ("INF", "EFF")] + \
           [(2, "INF", "SR", "present"), (2, "EFF", "SR", "present"), (2, "INF", "SD", "present")]
    assert len(run_insitu(tmp_path, states(rows), observed="EFF")) == 0
    # if it is only in EFF at t2, it counts
    rows[-2] = (2, "EFF", "SR", "present")
    rows[-3] = (2, "INF", "SR", "absent")
    out = run_insitu(tmp_path, states(rows), observed="EFF")
    assert len(out) == 1 and out.iloc[0].recipient_compartment == "EFF"


def test_donor_seen_only_in_prior_compartment(tmp_path):
    # old CIWARS setup: recipients in INF, prior = INF + EFF; donor seen only in EFF earlier
    out = run_insitu(tmp_path, states([
        (1, "EFF", "SD", "present"), (1, "INF", "SD", "absent"),
        (1, "INF", "SR", "absent"), (1, "EFF", "SR", "absent"),
        (2, "INF", "SR", "present"), (2, "EFF", "SR", "absent"),
    ]), observed="INF", prior="INF,EFF")
    assert len(out) == 1 and out.iloc[0].donor_last_compartment == "EFF"


def test_dates_ordered_as_dates(tmp_path):
    out = run_insitu(tmp_path, states([
        ("2020-09-01", "X", "SD", "present"), ("2020-09-01", "X", "SR", "absent"),
        ("2020-10-15", "X", "SR", "present"), ("2020-10-15", "X", "SD", "present"),
    ]))
    assert out.iloc[0].timepoint == "2020-10-15" and out.iloc[0].donor_last_timepoint == "2020-09-01"
    assert ins.time_order(["10", "9", "100"]) == {"9": 0, "10": 1, "100": 2}


def test_discordant_rank():
    a = ["Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "", "", ""]
    b = ["Bacteria", "Pseudomonadota", "Betaproteobacteria", "Burkholderiales", "", "", ""]
    assert ins.discordant_rank(a, b) == "c"
    assert ins.discordant_rank(a, a) == ""


def test_lineage_relation():
    a = ["Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "", "", ""]
    full = ["Bacteria", "Bacillota", "Bacilli", "Lactobacillales", "Enterococcaceae", "Enterococcus",
            "Enterococcus faecium"]
    assert ins.lineage_relation(full, full) == ("same", "", "s")
    assert ins.lineage_relation(a, a) == ("unresolved", "", "o")
    assert ins.lineage_relation(full, full[:6] + ["Enterococcus faecalis"]) == ("discordant", "s", "g")
    assert ins.lineage_relation(full, [""] * 7) == ("unresolved", "", "")


# ── in situ: lineage filter ──────────────────────────────────────────────────

DIRECTED = [(1, "X", "SD", "present"), (1, "X", "SR", "absent"),
            (2, "X", "SD", "present"), (2, "X", "SR", "present")]
EFAE = "d_Bacteria;p_Bacillota;c_Bacilli;o_Lactobacillales;f_Enterococcaceae;g_Enterococcus;s_Enterococcus faecium"
ECOLI = "d_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria;o_Enterobacterales;f_Enterobacteriaceae;g_Escherichia;s_Escherichia coli"
ECOLI_G = "d_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria;o_Enterobacterales;f_Enterobacteriaceae;g_Escherichia"


def test_same_lineage_never_called(tmp_path):
    for filt in ("discordant", "not_same"):
        out = run_insitu(tmp_path, states(DIRECTED), lineages={"R": ECOLI, "D": ECOLI}, lineage_filter=filt)
        assert len(out) == 0
        exc = pd.read_csv(tmp_path / "excluded.tsv", sep="\t")
        assert list(exc["reason"]) == ["lineage_same"]


def test_discordant_lineage_kept(tmp_path):
    out = run_insitu(tmp_path, states(DIRECTED), lineages={"R": ECOLI, "D": EFAE}, lineage_filter="discordant")
    assert len(out) == 1
    e = out.iloc[0]
    assert (e.lineage_relation, e.discordant_rank, e.shared_rank) == ("discordant", "p", "d")


def test_unresolved_lineage_needs_not_same(tmp_path):
    # same down to genus, one contig not classified to species: cannot rule out the same species
    lin = {"R": ECOLI, "D": ECOLI_G}
    assert len(run_insitu(tmp_path, states(DIRECTED), lineages=lin, lineage_filter="discordant")) == 0
    out = run_insitu(tmp_path, states(DIRECTED), lineages=lin, lineage_filter="not_same")
    assert len(out) == 1 and out.iloc[0].lineage_relation == "unresolved"


def test_filter_without_taxonomy_is_an_error(tmp_path):
    with pytest.raises(SystemExit):
        run_insitu(tmp_path, states(DIRECTED), lineage_filter="discordant")


# ── in situ: partial evidence (one boundary or contig not detected) ──────────

def test_partial_donor_is_not_a_donor(tmp_path):
    out = run_insitu(tmp_path, states([
        (1, "X", "SD", "partial"), (1, "X", "SR", "absent"),
        (2, "X", "SD", "present"), (2, "X", "SR", "present"),
    ]))
    # donor only partially seen at t1: no earlier donor, and it was seen before t2 -> not unresolved either
    assert len(out) == 0


def test_partial_recipient_earlier_is_not_new(tmp_path):
    out = run_insitu(tmp_path, states([
        (1, "X", "SD", "present"), (1, "X", "SR", "partial"),
        (2, "X", "SD", "present"), (2, "X", "SR", "present"),
    ]))
    assert len(out) == 0


def test_partial_recipient_is_not_called(tmp_path):
    out = run_insitu(tmp_path, states([
        (1, "X", "SD", "present"), (1, "X", "SR", "absent"),
        (2, "X", "SD", "present"), (2, "X", "SR", "partial"),
        (3, "X", "SD", "present"), (3, "X", "SR", "present"),
    ]))
    # first full presence at t3, but it was partially seen at t2 -> not new at t3
    assert len(out) == 0


# ── detect_state: region boundaries and contig presence ──────────────────────

def _call(seg, bp_counts, contig_rows, **kw):
    args = type("A", (), dict(min_reads=2, alpha=0.05, absent_fraction=1.0, min_boundaries=2,
                              contig_min_breadth=0.5, contig_breadth_ratio=0.5, **kw))
    obs = {"group": "G", "timepoint": "1", "compartment": "X"}
    bp = pd.DataFrame([{**obs, "sample": "s1", "bp_id": b, "span": n, "bridge": 0, "flank_depth": 10,
                        "mean_read_len": 150, "overhang": 25} for b, n in bp_counts.items()])
    es = pd.DataFrame(columns=list(obs) + ["sample", "es_id", "empty"])
    ct = pd.DataFrame([{**obs, "sample": "s1", "contig": c, "reads": r, "depth": d, "breadth": b}
                       for c, r, d, b in contig_rows])
    return ds.call_states(seg, bp, es, ct, args).set_index("segment_id")


SEG = pd.DataFrame([
    {"segment_id": "S1", "contig": "C1", "bp_left": "B1", "bp_right": "B2", "empty_site": ""},
    {"segment_id": "S2", "contig": "C2", "bp_left": "B3", "bp_right": "", "empty_site": ""},
])


def test_present_needs_both_boundaries_and_contig():
    st = _call(SEG, {"B1": 5, "B2": 5, "B3": 9}, [("C1", 500, 10, 0.99), ("C2", 500, 10, 0.99)])
    assert st.loc["S1", "state"] == "present" and st.loc["S1", "boundaries"] == 2
    # contig-end segment (one breakpoint) is never more than partial
    assert st.loc["S2", "state"] == "partial"


def test_one_boundary_only_is_partial():
    st = _call(SEG, {"B1": 6, "B2": 0, "B3": 0}, [("C1", 500, 10, 0.99), ("C2", 500, 10, 0.99)])
    assert st.loc["S1", "state"] == "partial"
    assert (st.loc["S1", "n_left"], st.loc["S1", "n_right"]) == (6, 0)


def test_contig_not_detected_is_partial():
    # junction reads at both boundaries but the contig's own sequence is barely covered
    st = _call(SEG, {"B1": 5, "B2": 5, "B3": 0}, [("C1", 40, 2.0, 0.05), ("C2", 0, 0, 0)])
    assert st.loc["S1", "contig_state"] == "not_detected" and st.loc["S1", "state"] == "partial"


def test_contig_breadth_pools_samples():
    args = type("A", (), dict(contig_min_breadth=0.5, contig_breadth_ratio=0.5))
    obs = {"group": "G", "timepoint": "1", "compartment": "X"}
    ct = pd.DataFrame([{**obs, "sample": s, "contig": "C", "reads": 10, "depth": 0.4, "breadth": 0.33}
                       for s in ("a", "b")])
    c = ds.contig_states(ct, args).iloc[0]
    # 1 - 0.67^2 = 0.551 >= 0.5, expected 1 - exp(-0.8) = 0.55
    assert c.contig_breadth == pytest.approx(0.551, abs=1e-3) and c.contig_state == "present"


# ── call.pool: sample ────────────────────────────────────────────────────────

def test_pool_sample_keeps_samples_apart():
    # two libraries of one observation, each with one junction read per boundary:
    # pooled they reach k = 2 at both boundaries; per sample, neither boundary does
    args = type("A", (), dict(min_reads=2, alpha=0.05, absent_fraction=1.0, min_boundaries=2,
                              contig_min_breadth=0.5, contig_breadth_ratio=0.5))
    obs = {"group": "G", "timepoint": "1", "compartment": "X"}
    bp = pd.DataFrame([{**obs, "sample": s, "bp_id": b, "span": 1, "bridge": 0, "flank_depth": 1,
                        "mean_read_len": 150, "overhang": 25} for s in ("a", "b") for b in ("B1", "B2")])
    es = pd.DataFrame(columns=list(obs) + ["sample", "es_id", "empty"])
    ct = pd.DataFrame([{**obs, "sample": s, "contig": "C1", "reads": 500, "depth": 10, "breadth": 0.99}
                       for s in ("a", "b")])
    seg = SEG.iloc[[0]]
    pooled = ds.call_states(seg, bp, es, ct, args)
    assert len(pooled) == 1 and pooled.iloc[0].state == "present"
    per = ds.call_states(seg, bp, es, ct, args, ds.OBS + ["sample"])
    assert len(per) == 2 and set(per["state"]) == {"partial"}   # 1 read per boundary < k


def test_combine_samples_strongest_state():
    st = pd.DataFrame([
        {"group": "G", "timepoint": "1", "compartment": "X", "segment_id": "S", "samples": s,
         "state": state, "n": n}
        for s, state, n in [("a", "absent", 0), ("b", "undetermined", 1), ("c", "partial", 3),
                            ("d", "present", 4), ("e", "present", 9)]])
    one = lambda rows: ins.combine_samples(st[st["samples"].isin(list(rows))]).iloc[0]
    assert one("ab").state == "absent"                     # a deep library rules it out
    assert one("abc").state == "partial"
    best = one("abcde")
    assert (best.state, best.samples, best.samples_present) == ("present", "e", 2)


def test_insitu_with_per_sample_states(tmp_path):
    rows = states(DIRECTED)
    # add a second, shallow library to every observation: it never changes the call
    extra = rows.assign(state="undetermined", n=0, samples="lib2")
    st = pd.concat([rows.assign(samples="lib1"), extra], ignore_index=True)
    out = run_insitu(tmp_path, st)
    assert len(out) == 1 and out.iloc[0].recipient_samples == "lib1"
    assert out.iloc[0].samples_present == "1"
