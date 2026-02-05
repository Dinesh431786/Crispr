import pandas as pd

from crispr_app.analysis import check_pam, find_gRNAs, find_off_targets_detailed


def test_check_pam_ngg():
    assert check_pam("AGG", "NGG")
    assert not check_pam("AAA", "NGG")


def test_find_guides_returns_scores():
    seq = "ATGAGTCTGCTCTTCGCGTTGGAGTGAAATCTGAGATGATGGGTTGAAATCGCAGTTCGACCTGAACTTTTATCTGCTCTTCGCGTTGAGCGGACCGTGGGAAGTTTCGCGTTGATCAGTTCTTCTGCTCTTCGCGTTTAAGCCTTGCGTTGTTTATCTGCTCTTCGCGTTTATCAGCCTGGGCGTTGATCTTTTATCTGCTCTTCGCGTTAACGGAAGCCGG"
    df = find_gRNAs(seq, pam="NGG")
    assert not df.empty
    assert {"HybridScore", "MLScore", "ConsensusScore"}.issubset(df.columns)


def test_off_target_search_returns_expected_schema_and_strands():
    guides = pd.DataFrame({"gRNA": ["GAGTCTGCTCTTCGCGTTGG"]})
    background = "".join(
        [
            "TTTGAGTCTGCTCTTCGCGATGGAT",  # + strand 1 mismatch target
            "AACCGGTTAACCGGTTAA",
            "CCATCACGCGAAGAGCAGACTCAA",  # reverse complement target context
        ]
    )

    res = find_off_targets_detailed(guides, background, max_mismatches=2, pam="NGG", include_reverse_strand=True)
    assert {"gRNA", "OffTargetPos", "Mismatches", "TargetSeq", "PAM", "CFD_Score", "Strand"}.issubset(res.columns)
    assert not res.empty
    assert set(res["Strand"].unique()).issubset({"+", "-"})
