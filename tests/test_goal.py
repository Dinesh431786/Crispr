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
