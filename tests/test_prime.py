from crispr_app.analysis import design_prime_editing_pegRNAs

def test_design_prime_editing_pegRNAs():
    # Sequence with a PAM near a target
    # 012345678901234567890123
    # GCTAGCTAGCTAGCTAGCTA TGG
    seq = "GCTAGCTAGCTAGCTAGCTATGG" + "A" * 50
    # Target index 20 (the 'A' after PAM if we were looking at that, but let's target position 30)
    target_pos = 30
    desired_edit = "C"

    df = design_prime_editing_pegRNAs(seq, target_pos, desired_edit)
    assert not df.empty
    assert "Spacer" in df.columns
    assert "PBS_Seq" in df.columns
    assert "RTT_Seq" in df.columns

    # Check if edit is applied in RTT
    # Nick pos for the PAM at 20 is 17.
    # RTT starts at 17.
    # Target pos 30 is at index 13 in RTT.
    for _, row in df.iterrows():
        if row["NickPos"] == 17:
            assert row["RTT_Seq"][13] == "C"
