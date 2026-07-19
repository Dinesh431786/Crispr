"""Uncertainty-aware ranking + normalized (per-guide-width) conformal.

The point of normalized conformal is that interval WIDTH varies per guide, which
makes uncertainty a genuine ranking dimension (constant-width intervals would make
every strategy collapse to `balanced`). These tests assert both properties.
"""
import numpy as np

from crispr_app import models
from crispr_app.analysis import find_gRNAs

SEQ = ("ATGGCCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCC"
       "GCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAG"
       "CTGCAAGAACTCTTCCTCACGCGCGTCGGGCTCGACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCC")


def test_interval_has_per_guide_width():
    """Normalized conformal ships variable half-widths (not one constant q)."""
    models.reset_caches()
    df = find_gRNAs(SEQ, pam="NGG")
    assert "CI_halfwidth" in df.columns
    widths = df["CI_halfwidth"].dropna().unique()
    assert len(widths) > 1  # genuinely heteroscedastic, not a single constant width
    models.reset_caches()


def test_conservative_ranking_reorders():
    """Conservative strategy must differ from balanced (else it's a no-op)."""
    models.reset_caches()
    bal = find_gRNAs(SEQ, pam="NGG", ranking_strategy="balanced")
    con = find_gRNAs(SEQ, pam="NGG", ranking_strategy="conservative")
    assert not bal.empty
    assert list(bal["gRNA"]) != list(con["gRNA"])  # the ordering genuinely changes
    models.reset_caches()


def test_conservative_prefers_certainty():
    """Between two similar scores, conservative favours the tighter interval."""
    models.reset_caches()
    con = find_gRNAs(SEQ, pam="NGG", ranking_strategy="conservative")
    # The conservative #1 maximises (score - half-width); verify it is the argmax.
    val = con["ConsensusScore"].astype(float) - con["CI_halfwidth"].astype(float)
    assert np.isclose(val.iloc[0], val.max())
    models.reset_caches()


def test_displayed_score_is_strategy_invariant():
    """Strategy changes ORDER, never the displayed ConsensusScore of a guide."""
    models.reset_caches()
    bal = find_gRNAs(SEQ, pam="NGG", ranking_strategy="balanced").set_index("gRNA")["ConsensusScore"]
    opt = find_gRNAs(SEQ, pam="NGG", ranking_strategy="optimistic").set_index("gRNA")["ConsensusScore"]
    common = bal.index.intersection(opt.index)
    assert (bal.loc[common] == opt.loc[common]).all()
    models.reset_caches()


def test_invalid_strategy_falls_back_to_balanced():
    models.reset_caches()
    bad = find_gRNAs(SEQ, pam="NGG", ranking_strategy="nonsense")
    bal = find_gRNAs(SEQ, pam="NGG", ranking_strategy="balanced")
    assert list(bad["gRNA"]) == list(bal["gRNA"])
    models.reset_caches()
