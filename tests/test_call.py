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


def run_insitu(tmp_path, st, observed="", prior=""):
    st.to_csv(tmp_path / "states.tsv", sep="\t", index=False)
    SEGMENTS.to_csv(tmp_path / "segments.tsv", sep="\t", index=False)
    EVENTS.to_csv(tmp_path / "ev.tsv", sep="\t", index=False)
    ins.main(["--states", str(tmp_path / "states.tsv"), "--segments", str(tmp_path / "segments.tsv"),
              "--event-segments", str(tmp_path / "ev.tsv"), "--observed", observed, "--prior", prior,
              "--out", str(tmp_path / "out.tsv")])
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
