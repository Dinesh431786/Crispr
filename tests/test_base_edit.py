from crispr_app.base_edit import editable_targets, summarize


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
