from crispr_app.analysis import find_gRNAs, check_pam


def test_check_pam_ngg():
    assert check_pam("AGG", "NGG")
    assert not check_pam("AAA", "NGG")


def test_find_guides_returns_scores():
    seq = "ATGAGTCTGCTCTTCGCGTTGGAGTGAAATCTGAGATGATGGGTTGAAATCGCAGTTCGACCTGAACTTTTATCTGCTCTTCGCGTTGAGCGGACCGTGGGAAGTTTCGCGTTGATCAGTTCTTCTGCTCTTCGCGTTTAAGCCTTGCGTTGTTTATCTGCTCTTCGCGTTTATCAGCCTGGGCGTTGATCTTTTATCTGCTCTTCGCGTTAACGGAAGCCGG"
    df = find_gRNAs(seq, pam="NGG")
    assert not df.empty
    assert {"HybridScore", "MLScore", "ConsensusScore"}.issubset(df.columns)
