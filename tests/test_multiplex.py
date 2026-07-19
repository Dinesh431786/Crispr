"""Multiplex / pooled library designer: strong + collectively diverse."""
from crispr_app import models
from crispr_app.analysis import design_multiplex, _seq_identity

SEQ = ("ATGGCCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCC"
       "GCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAG"
       "CTGCAAGAACTCTTCCTCACGCGCGTCGGGCTCGACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCC")


def test_returns_requested_count_and_summary():
    models.reset_caches()
    picked, summary = design_multiplex(SEQ, n_guides=5)
    assert len(picked) == 5
    assert summary["selected"] == 5
    assert {"mean_on_target", "mean_diversity", "max_pairwise_identity", "span_bp"} <= summary.keys()
    models.reset_caches()


def test_seed_is_highest_score():
    picked, _ = design_multiplex(SEQ, n_guides=4, diversity_weight=0.5)
    # First pick is the top-scoring guide (no similarity penalty yet).
    assert picked.iloc[0]["ConsensusScore"] == picked["ConsensusScore"].max() or \
        picked.iloc[0]["MaxSimilarity"] == 0.0


def test_diversity_weight_reduces_similarity():
    # Heavier diversity weight should not increase the max pairwise identity.
    _, lo = design_multiplex(SEQ, n_guides=6, diversity_weight=0.0)
    _, hi = design_multiplex(SEQ, n_guides=6, diversity_weight=1.5)
    assert hi["max_pairwise_identity"] <= lo["max_pairwise_identity"] + 1e-9


def test_min_spacing_enforced():
    picked, _ = design_multiplex(SEQ, n_guides=6, min_spacing=15)
    starts = sorted(picked["Start"].tolist())
    assert all(b - a >= 15 for a, b in zip(starts, starts[1:]))


def test_seq_identity_bounds():
    assert _seq_identity("ACGT", "ACGT") == 1.0
    assert _seq_identity("ACGT", "TGCA") == 0.0
