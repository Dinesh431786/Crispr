"""What-if single-nucleotide saturation sensitivity of the on-target score."""
from crispr_app import models


def test_scan_shape_and_bounds():
    models.reset_caches()
    r = models.sensitivity_scan("GATGTGGCGGTCCGGATCGA")
    assert len(r["substitutions"]) == 20 * 3          # every position × 3 alt bases
    assert len(r["per_position"]) == 20
    assert 0.0 <= r["baseline"] <= 1.0
    for s in r["substitutions"]:
        assert s["to"] != s["from"]
        assert 0.0 <= s["score"] <= 1.0
    models.reset_caches()


def test_best_and_worst_are_extremes():
    models.reset_caches()
    r = models.sensitivity_scan("GATGTGGCGGTCCGGATCGA")
    assert r["best"]["delta"] == max(s["delta"] for s in r["substitutions"])
    assert r["worst"]["delta"] == min(s["delta"] for s in r["substitutions"])
    models.reset_caches()


def test_deterministic():
    models.reset_caches()
    a = models.sensitivity_scan("GCTCGACATCGGCAAGGTGT")
    b = models.sensitivity_scan("GCTCGACATCGGCAAGGTGT")
    assert a == b   # fully deterministic, no randomness
    models.reset_caches()
