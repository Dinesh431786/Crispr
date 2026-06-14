import pandas as pd

from crispr_app.analysis import find_off_targets_detailed, summarize_specificity
from crispr_app.offtarget import (
    aggregate_specificity,
    calculate_cfd_score,
    mit_hit_score,
)


def test_cfd_perfect_match_is_one():
    guide = "GACGATCAGTCAGGATCACC"
    assert calculate_cfd_score(guide, guide, "AGG") == 1.0


def test_cfd_mismatch_reduces_score():
    guide = "GACGATCAGTCAGGATCACG"
    off = "GACGATCAGTCAGGATCACT"  # rG:dT at PAM-proximal position (weight 0.9)
    assert calculate_cfd_score(guide, off, "AGG") < 1.0


def test_mit_perfect_match_is_100():
    guide = "GACGATCAGTCAGGATCACC"
    assert mit_hit_score(guide, guide) == 100.0


def test_mit_more_mismatches_lower_score():
    guide = "GACGATCAGTCAGGATCACC"
    one = "GACGATCAGTCAGGATCACA"
    two = "GACGATCAGTCAGGATCAAA"
    assert mit_hit_score(guide, two) < mit_hit_score(guide, one)


def test_aggregate_specificity_bounds():
    assert aggregate_specificity([]) == 100.0
    assert aggregate_specificity([1.0], scale=100.0) < 100.0


def test_both_strand_detection():
    # An off-target (1 mismatch) lives on the reverse strand: embed the
    # reverse-complement of the mismatched protospacer + NGG PAM in the background.
    from Bio.Seq import Seq
    guide = "GACGATCAGTCAGGATCACC"
    off_proto = "GACGATCAGTCAGGATCACA"  # 1 mismatch vs guide
    rc_site = str(Seq(off_proto + "AGG").reverse_complement())
    bg = "A" * 100 + rc_site + "A" * 100
    df = pd.DataFrame({"gRNA": [guide]})
    res = find_off_targets_detailed(df, bg, max_mismatches=1, pam="NGG")
    assert not res.empty
    assert (res["Strand"] == "-").any()
    assert {"CFD_Score", "MIT_Score", "Strand"}.issubset(res.columns)


def test_summarize_specificity_shape():
    df = pd.DataFrame({"gRNA": ["GACGATCAGTCAGGATCACC"]})
    res = find_off_targets_detailed(df, "A" * 200, max_mismatches=2, pam="NGG")
    summary = summarize_specificity(res, ["GACGATCAGTCAGGATCACC"])
    assert {"gRNA", "OffTargetCount", "CFD_Specificity", "MIT_Specificity"}.issubset(summary.columns)
    assert summary.iloc[0]["CFD_Specificity"] == 100.0
