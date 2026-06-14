import math

import pytest

from crispr_app.benchmark import evaluate, pearson, spearman
from crispr_app.scoring import score_breakdown


def test_spearman_perfect_monotonic():
    x = [1, 2, 3, 4, 5]
    y = [10, 20, 30, 40, 50]
    assert spearman(x, y) == pytest.approx(1.0)
    assert spearman(x, [50, 40, 30, 20, 10]) == pytest.approx(-1.0)


def test_pearson_handles_constant():
    assert math.isnan(pearson([1, 1, 1], [1, 2, 3]))


def test_spearman_with_ties():
    # Ties must not raise and stay within [-1, 1].
    r = spearman([1, 1, 2, 3], [2, 2, 1, 4])
    assert -1.0 <= r <= 1.0


def test_evaluate_returns_metrics():
    records = [
        {"guide": "GACGATCAGTCAGGATCACC", "measured": 0.8},
        {"guide": "TTTTTTTTTTTTTTTTTTTT", "measured": 0.05},
        {"guide": "GCGCGCGCGCGCGCGCGCGC", "measured": 0.4},
        {"guide": "GACCATGCATGCATCAGGAG", "measured": 0.6},
    ]
    out = evaluate(records)
    assert out["n"] == 4
    assert -1.0 <= out["spearman"] <= 1.0


def test_score_breakdown_sums_are_present():
    bd = score_breakdown("GACGATCAGTCAGGATCACC", "A")
    assert "contributions" in bd and "intercept" in bd["contributions"]
    assert 0.0 <= bd["score"] <= 1.0
    assert "gc_content" in bd["contributions"]
