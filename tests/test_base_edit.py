from crispr_app.base_edit import assess, editable_targets, position_weight, summarize


def test_cbe_target_in_window():
    # C at position 5 (1-based) -> CBE target inside the 4-8 window.
    guide = "AAAACAAAAAAAAAAAAAAA"
    t = editable_targets(guide)
    assert any(x["editor"] == "CBE" and x["pos"] == 5 for x in t)


def test_abe_target_in_window():
    guide = "TTTTTATTTTTTTTTTTTTT"  # A at position 6
    t = editable_targets(guide)
    assert any(x["editor"] == "ABE" and x["pos"] == 6 for x in t)


def test_bases_outside_window_ignored():
    # C at position 1 and 20 are outside the 4-8 window.
    guide = "CAATTTTTTTTTTTTTTTTC"
    assert editable_targets(guide) == []


def test_summarize_shape():
    s = summarize("AAAACAAAAAAAAAAAAAAA")
    assert s["editable"] is True
    assert s["CBE_positions"] == [5]
    assert s["ABE_positions"]  # A's at positions 4,6,7,8...
    assert s["window"] == "4-8"


def test_window_efficiency_peaks_at_center():
    # The position-weight profile peaks near position 6 and decays to the edges.
    assert position_weight(6) == max(position_weight(p) for p in range(4, 9))
    assert position_weight(4) < position_weight(6)
    assert position_weight(1) == 0.0  # outside the window


def test_purity_penalises_bystanders():
    # One C in the window -> clean single edit (purity 1.0).
    clean = assess("AAAAACAAAAAAAAAAAAAA", "CBE")  # C at position 6 only
    assert clean["editable"] and not clean["bystanders"]
    assert clean["purity"] == 1.0
    # Two C's in the window -> bystander -> purity < 1 and lower composite score.
    byst = assess("AAACACAAAAAAAAAAAAAA", "CBE")  # C at positions 4 and 6
    assert byst["bystanders"]
    assert byst["purity"] < 1.0
    assert byst["be_score"] < clean["be_score"]


def test_assess_not_editable():
    a = assess("TTTTTTTTTTTTTTTTTTTT", "CBE")  # no C anywhere
    assert a["editable"] is False and a["be_score"] == 0.0
