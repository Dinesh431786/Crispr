"""Distribution-free uncertainty for on-target scores (split conformal prediction).

Problem this solves
--------------------
Every common CRISPR on-target tool returns a single point score and no estimate
of how reliable it is. Our own benchmarks show that, depending on context, even
strong models correlate only ~0.3-0.8 with measured efficiency — so a bare point
estimate is overconfident. Researchers selecting guides deserve to know the
*uncertainty*, not just the number.

Method
------
Split conformal regression (Vovk; Lei et al. 2018, JASA "Distribution-Free
Predictive Inference for Regression"). Given a model trained on a proper-training
split and a held-out *calibration* split, let r_i = |y_i - ŷ_i| be the
calibration residuals. For miscoverage alpha, the prediction-interval half-width
is the ((1-alpha)(n+1))/n empirical quantile of {r_i}. The resulting interval
[ŷ - q, ŷ + q] has *guaranteed* marginal coverage >= 1 - alpha under
exchangeability — no distributional assumptions, no retraining, pure NumPy.

The quantiles q are computed at train time and stored in the model so prediction
is O(1). Coverage is verified empirically in the test suite and on real data
(see BENCHMARKS.md).
"""

from __future__ import annotations

import numpy as np

# Coverage levels we calibrate by default: 80% and 90% intervals.
DEFAULT_ALPHAS = (0.20, 0.10)


def conformal_quantile(residuals, alpha: float) -> float:
    """Finite-sample-valid split-conformal quantile of absolute residuals.

    Uses the ((1-alpha)(n+1))/n empirical quantile (the conservative correction
    that guarantees marginal coverage >= 1 - alpha).
    """
    r = np.abs(np.asarray(residuals, dtype=float))
    r = r[np.isfinite(r)]
    n = len(r)
    if n == 0:
        return float("nan")
    level = min(1.0, np.ceil((1 - alpha) * (n + 1)) / n)
    return float(np.quantile(r, level, method="higher"))


def calibrate(predictions, targets, alphas=DEFAULT_ALPHAS) -> dict:
    """Return {"q<percent>": halfwidth} for each coverage level, e.g. {"q90": ...}."""
    residuals = np.asarray(targets, dtype=float) - np.asarray(predictions, dtype=float)
    out = {}
    for a in alphas:
        out[f"q{int(round((1 - a) * 100))}"] = round(conformal_quantile(residuals, a), 4)
    return out


def interval(point: float, halfwidth: float, lo: float = 0.0, hi: float = 1.0) -> tuple[float, float]:
    """Clamp the symmetric conformal interval to the valid score range."""
    return (round(max(lo, point - halfwidth), 3), round(min(hi, point + halfwidth), 3))


def empirical_coverage(predictions, targets, halfwidth: float) -> float:
    """Fraction of targets falling within +/- halfwidth of the prediction."""
    p = np.asarray(predictions, dtype=float)
    t = np.asarray(targets, dtype=float)
    return float(np.mean(np.abs(t - p) <= halfwidth))
