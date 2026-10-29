"""Self-checks for domains.py: paper numbers, parser, overlap rule. Run: pytest tests/"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import domains as D  # noqa: E402

# hmmsearch domtblout: target=protein, query=domain; cols 12/13 = i-Evalue/score, 19/20 = env
DOMTBL = """\
# comment
P1 - 300 dA PF1 50 1e-20 80 0 1 2 1e-21 1e-20 70.0 0 1 50 10 60 10 60 0.9 desc
P1 - 300 dB PF2 50 1e-10 40 0 1 1 1e-11 1e-10 35.0 0 1 50 40 90 40 90 0.9 desc
P1 - 300 dA PF1 50 1e-20 80 0 2 2 1e-3 1e-3 5.0 0 1 50 200 250 200 250 0.9 desc
P2 - 300 dB PF2 50 1e-10 40 0 1 1 1e-11 1e-10 35.0 0 1 50 5 55 5 55 0.9 desc
"""


def test_paper_expectations():
    assert D.paper_expectation(5090, 274, 118894)["expected"] == pytest.approx(11.73, abs=0.01)
    assert D.paper_expectation(5090, 274, 118894)["p_zero"] < 1e-5
    assert D.cooccurrence_test(0, 1324, 414, 30712)["expected"] == pytest.approx(17.85, abs=0.01)


def test_domtblout_evalue_and_overlap():
    hits = D.read_domtblout(DOMTBL)
    # P1: dB (40-90) overlaps better dA (10-60) -> dropped; dA at 200-250 fails E-value
    assert list(zip(hits.protein, hits.domain)) == [("P1", "dA"), ("P2", "dB")]
    scan = D.read_domtblout(DOMTBL, program="hmmscan")
    assert set(scan.protein) == {"dA", "dB"}


def test_cooccurrence_tails():
    r = D.cooccurrence_test(10, 10, 10, 1000)
    assert r["p_attract"] < 1e-20 and r["p_avoid"] == pytest.approx(1.0)
    hits = D.read_pairs(b"p1\tA\np1\tB\np2\tA\np3\tB\n")
    m = D.cooccurrence_matrix(hits, ["A", "B"], 100)
    assert m.observed.tolist() == [1] and "q_bh" in m


def test_compare_counts_identity():
    import pandas as pd
    s = pd.Series({"a": 1, "b": 10, "c": 100, "d": 1000})
    _, res = D.compare_counts(s, s)
    assert res["slope"] == pytest.approx(1) and res["pearson_r"] == pytest.approx(1)


def test_security_limits():
    import gzip
    bomb = gzip.compress(b"\0" * 10_000)
    with pytest.raises(ValueError):
        D._text(bomb, limit=1_000)
    two = gzip.compress(b"a\n") + gzip.compress(b"b\n")  # multi-member, like pfam_demo.hmm.gz
    assert D._text(two) == "a\nb\n"
    hits = D.read_pairs(b"=HYPERLINK(1)\t@SUM(A1)\np2\tA\n")
    assert hits.protein.tolist() == ["'=HYPERLINK(1)", "p2"] and hits.domain[0] == "'@SUM(A1)"
    with pytest.raises(ValueError):
        D.run_hmmsearch(b">a\nMKV\n", b"", cpus=1)  # empty HMM file must fail cleanly
