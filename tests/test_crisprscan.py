import pytest

from crispr_app.analysis import find_gRNAs
from crispr_app.crisprscan import PARAMS_CRISPRSCAN, score_35mer, score_from_context


def test_reference_vector_matches_crispor():
    # Verbatim doctest vector from CRISPOR's crisporEffScores.py -> 77.
    assert score_35mer("TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA") == 77


def test_param_table_is_complete():
    assert len(PARAMS_CRISPRSCAN) == 91


def test_requires_35mer():
    with pytest.raises(ValueError):
        score_35mer("ACGT")


def test_score_from_context_bounds_and_none():
    seq = "AAAAAA" + "GACGATCAGTCAGGATCACC" + "TGG" + "AAAAAA"  # 6+20+3+6 = 35
    s = score_from_context(seq, 6, 20)
    assert s is not None and 0.0 <= s <= 1.0
    # Too close to the end -> insufficient context -> None.
    assert score_from_context("ACGT" + "GACGATCAGTCAGGATCACC" + "TGG", 4, 20) is None


def test_find_grnas_includes_crisprscan_column():
    seq = ("ATGAGTCTGCTCTTCGCGTTGGAGTGAAATCTGAGATGATGGGTTGAAATCGCAGTTCGACCTGAACTTTT"
           "ATCTGCTCTTCGCGTTGAGCGGACCGTGGGAAGTTTCGCGTTGATCAGTTCTTCTGCTCTTCGCGTTTAAG"
           "CCTTGCGTTGTTTATCTGCTCTTCGCGTTTATCAGCCTGGGCGTTGATCTTTTATCTGCTCTTCGCGTTAACGGAAGCCGG")
    df = find_gRNAs(seq, pam="NGG")
    assert "CRISPRScanScore" in df.columns
    # At least one interior guide should have a computable CRISPRscan score.
    assert any(isinstance(v, float) for v in df["CRISPRScanScore"])
