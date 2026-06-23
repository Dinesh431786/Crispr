from crispr_app.structure import max_base_pairs, self_complementarity


def test_no_pairs_in_unpairable_spacer():
    # Poly-A cannot base-pair with itself.
    assert self_complementarity("AAAAAAAAAAAAAAAAAAAA") == 0.0


def test_palindrome_self_pairs():
    # A perfect inverted repeat folds into a hairpin -> high self-complementarity.
    hairpin = "GGGGGGGGGGCCCCCCCCCC"
    assert self_complementarity(hairpin) > 0.6


def test_score_in_unit_interval():
    for g in ["GACGATCAGTCAGGATCACC", "GCGCGCGCGCGCGCGCGCGC", "ATATATATATATATATATAT"]:
        s = self_complementarity(g)
        assert 0.0 <= s <= 1.0


def test_min_loop_enforced():
    # Adjacent complementary bases cannot pair (loop too small).
    assert max_base_pairs("GC", min_loop=3) == 0
