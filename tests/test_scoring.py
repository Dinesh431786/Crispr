from crispr_app.scoring import on_target_score, gc_fraction


def test_score_in_unit_interval():
    for guide in ["GACGTACGTACGTACGTACG", "TTTTTTTTTTTTTTTTTTTT", "GCGCGCGCGCGCGCGCGCGC"]:
        s = on_target_score(guide)
        assert 0.0 <= s <= 1.0


def test_poly_t_is_penalised():
    balanced = "GACGATCAGTCAGGATCAGG"
    polyt = "GACGATCAGTTTTTGATCAGG"[:20]
    assert on_target_score(polyt) < on_target_score(balanced)


def test_pam_proximal_g_beats_t():
    # Identical except PAM-proximal nucleotide (position 20); G favoured over T.
    base = " GACGATCAGTCAGGATCAG".replace(" ", "")
    with_g = base + "G"
    with_t = base + "T"
    assert on_target_score(with_g) > on_target_score(with_t)


def test_extreme_gc_disfavoured():
    optimal = "GACGATCAGTCAGGATCACC"
    all_at = "ATATATATATATATATATAT"
    assert on_target_score(optimal) > on_target_score(all_at)


def test_short_guide_returns_zero():
    assert on_target_score("ACGT") == 0.0


def test_gc_fraction():
    assert gc_fraction("GGCC") == 1.0
    assert gc_fraction("ATAT") == 0.0
