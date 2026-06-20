import numpy as np

from crispr_app import models
from crispr_app.conformal import (
    calibrate,
    conformal_quantile,
    empirical_coverage,
    interval,
)


def test_quantile_guarantees_coverage_on_holdout():
    # Split-conformal must achieve >= target coverage on exchangeable data.
    rng = np.random.default_rng(0)
    n = 4000
    pred = rng.uniform(0, 1, n)
    target = pred + rng.normal(0, 0.15, n)  # noisy but exchangeable
    cal, test = slice(0, n // 2), slice(n // 2, n)
    q90 = conformal_quantile(target[cal] - pred[cal], alpha=0.10)
    cov = empirical_coverage(pred[test], target[test], q90)
    assert cov >= 0.88  # ~0.90 with finite-sample slack


def test_calibrate_returns_levels():
    pred = np.linspace(0, 1, 200)
    tgt = pred + 0.05
    c = calibrate(pred, tgt)
    assert "q80" in c and "q90" in c
    assert c["q90"] >= c["q80"] >= 0  # higher coverage -> wider interval


def test_interval_is_clamped():
    lo, hi = interval(0.95, 0.4, lo=0.0, hi=1.0)
    assert lo >= 0.0 and hi <= 1.0


def test_shipped_default_model_has_calibrated_interval():
    models.reset_caches()
    ci = models.predict_interval("GACGATCAGTCAGGATCACC")
    assert ci is not None  # default.json ships conformal calibration
    assert 0.0 <= ci["low"] <= ci["point"] <= ci["high"] <= 1.0
    assert ci["coverage"] == 0.9
    models.reset_caches()
