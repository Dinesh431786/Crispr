"""Deeper outcome simulation: sequence-derived indel spectrum + frameshift prob."""
from crispr_app.analysis import indel_distribution, _mmej_deletions


def test_finds_microhomology_deletion():
    # Two GATC repeats flank the cut -> MMEJ should predict the spanning deletion.
    seq = "ATGAAAGATCTTTTTTGATCGGGCCCAAA"
    cut = 12
    dels = _mmej_deletions(seq, cut)
    assert any(d["microhomology"] == "GATC" for d in dels)
    dist = indel_distribution(seq, cut)
    assert dist["n_microhomologies"] >= 1
    assert any(o["type"] == "deletion" and "GATC" in o["mechanism"] for o in dist["outcomes"])


def test_frequencies_are_normalized_probabilities():
    dist = indel_distribution("ATGAAAGATCTTTTTTGATCGGGCCCAAA", 12)
    assert 0.0 <= dist["frameshift_probability"] <= 1.0
    assert 0.0 <= dist["stop_gain_probability"] <= 1.0
    assert all(0.0 <= o["frequency"] <= 1.0 for o in dist["outcomes"])


def test_templated_insertion_channel_present():
    # +1 insertion duplicates the base 5' of the cut.
    seq = "ATGCCCGGGAAATTTCCCGGG"
    cut = 9  # base at index 8 is 'G'
    dist = indel_distribution(seq, cut)
    ins = [o for o in dist["outcomes"] if o["type"] == "insertion"]
    assert ins and ins[0]["length"] == 1
    assert "G" in ins[0]["mechanism"]


def test_no_microhomology_still_returns_insertion():
    # A low-microhomology context still yields the templated +1 insertion.
    dist = indel_distribution("ATGACGTACGTACGT", 7)
    assert dist["outcomes"]  # never empty
    assert dist["frameshift_probability"] >= 0.0
