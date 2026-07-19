"""Goal-aware (knockout mode) routing: the out-of-frame model changes ranking."""
from crispr_app import models
from crispr_app.analysis import find_gRNAs

GUIDE = "GACGATCAGTCAGGATCACC"
SEQ = ("ATGGCCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCC"
       "GCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAG")


def test_oof_model_ships():
    models.reset_caches()
    assert models._load_oof() is not None  # crispr_app/models/default_oof.json present


def test_knockout_prediction_differs_from_general():
    models.reset_caches()
    general = models.predict_on_target(GUIDE, goal="general")
    knockout = models.predict_on_target(GUIDE, goal="knockout")
    assert 0.0 <= knockout <= 1.0
    assert knockout != general  # a distinct out-of-frame model, not an alias
    models.reset_caches()


def test_knockout_interval_calibrated():
    models.reset_caches()
    ci = models.predict_interval(GUIDE, goal="knockout")
    assert ci is not None and 0.0 <= ci["low"] <= ci["point"] <= ci["high"] <= 1.0
    models.reset_caches()


def test_find_grnas_goal_reorders():
    models.reset_caches()
    gen = find_gRNAs(SEQ, pam="NGG", goal="general")
    ko = find_gRNAs(SEQ, pam="NGG", goal="knockout")
    assert not gen.empty and not ko.empty
    # OnTargetScore comes from a different model, so at least some scores differ.
    assert list(gen["OnTargetScore"]) != list(ko["OnTargetScore"])
    models.reset_caches()


def test_knockin_ranks_by_proximity():
    # Knock-in mode adds CutDist and ranks by cutting x proximity to the edit site.
    df = find_gRNAs(SEQ, pam="NGG", goal="knockin", target_pos=60)
    assert not df.empty and "CutDist" in df.columns
    # Top-ranked guides cluster nearer the edit site than the overall average.
    assert df.head(3)["CutDist"].mean() < df["CutDist"].mean()
    # Displayed score stays the interpretable cutting score (0-1), not crushed.
    assert 0.0 <= df.iloc[0]["ConsensusScore"] <= 1.0
    # A different target site re-ranks to a different #1 guide.
    df2 = find_gRNAs(SEQ, pam="NGG", goal="knockin", target_pos=180)
    assert df.iloc[0]["gRNA"] != df2.iloc[0]["gRNA"]


def test_crispri_ranks_near_tss():
    # Use a longer target so there are guides on both sides of the TSS (the short
    # SEQ has too few, sparsely placed, to test proximity meaningfully).
    long_seq = SEQ + ("GCCGCGGTGGCGGTCTGGACCACGCCGGAGAGCGTCGAAGCGGGGGCGGTGTTCGCC"
                      "GAGATCGGCCCGCGCATGGCCGAGTACAAGCCCACGGTGCGCCTCGCC")
    df = find_gRNAs(long_seq, pam="NGG", goal="crispri", target_pos=140)
    assert "TSSDist" in df.columns
    # The proximity weighting pulls near-TSS guides up: the top cluster is clearly
    # closer to the TSS than the overall set, and the #1 guide is in the near half.
    assert df.head(5)["TSSDist"].mean() < df["TSSDist"].mean()
    assert df.iloc[0]["TSSDist"] <= df["TSSDist"].median()


def test_baseedit_prioritizes_editable_guides():
    df = find_gRNAs(SEQ, pam="NGG", goal="baseedit", editor="CBE")
    assert "BE_targets" in df.columns
    # The #1 guide must have a CBE-editable base in the window.
    assert df.iloc[0]["BE_targets"] != ""
    # CBE tags only reference C positions (from -> C).
    assert all(t.startswith("C") for t in df.iloc[0]["BE_targets"].split(","))
