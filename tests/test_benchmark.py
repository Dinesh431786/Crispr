import math

import pytest

from crispr_app.benchmark import auc, evaluate, load_crispor_context, pearson, separation_auc, spearman
from crispr_app.scoring import score_breakdown


def test_auc_perfect_and_random():
    # Perfectly separating scores -> AUC 1.0; reversed -> 0.0.
    assert auc([0.1, 0.2, 0.8, 0.9], [0, 0, 1, 1]) == pytest.approx(1.0)
    assert auc([0.9, 0.8, 0.2, 0.1], [0, 0, 1, 1]) == pytest.approx(0.0)


def test_auc_single_class_is_nan():
    assert math.isnan(auc([0.1, 0.2, 0.3], [1, 1, 1]))


def test_separation_auc_matches_ranking():
    # Predictions that rank-align with truth separate top/bottom perfectly.
    preds = list(range(100))
    measured = list(range(100))
    assert separation_auc(preds, measured, 33, 67) == pytest.approx(1.0)


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


def test_load_crispor_context(tmp_path):
    p = tmp_path / "mini.context.tab"
    p.write_text(
        "guide\tseq\tdb\tpos\tmodFreq\tlongSeq\n"
        "g1\tGGCTGCTTTACCCGCTGTGG\thg19\tchr1:1-2:+\t2.9\t"
        + "AAAAAA" + "GGCTGCTTTACCCGCTGTGG" + "TGGAAAAAA" + "\n"
        "g2\tTCCGGGTTGGCCTTCCACTG\thg19\tchr1:3-4:+\t1.1\t"
        + "TTTTTT" + "TCCGGGTTGGCCTTCCACTG" + "CGGTTTTTT" + "\n"
    )
    recs = load_crispor_context(str(p))
    assert len(recs) == 2
    assert recs[0]["guide"] == "GGCTGCTTTACCCGCTGTGG"
    assert recs[0]["measured"] == 2.9
    assert len(recs[0]["mer35"]) == 35  # 6 + 20 + 3 + 6


def test_score_breakdown_sums_are_present():
    bd = score_breakdown("GACGATCAGTCAGGATCACC", "A")
    assert "contributions" in bd and "intercept" in bd["contributions"]
    assert 0.0 <= bd["score"] <= 1.0
    assert "gc_content" in bd["contributions"]
