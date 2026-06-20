import random

from Bio.Seq import Seq

from crispr_app.genome import (
    genome_specificity,
    iter_fasta,
    scan_fasta_text,
    scan_records,
)


def _filler(n, seed=0):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_finds_planted_offtarget_forward():
    guide = "GACGATCAGTCAGGATCACC"
    off = "GACGATCAGTCAGGATCACA"  # 1 mismatch
    genome = _filler(500, 1) + off + "AGG" + _filler(500, 2)
    hits = scan_records([("chr1", genome)], [guide], max_mismatches=2, pam="NGG")
    assert not hits.empty
    row = hits.iloc[0]
    assert row["Chrom"] == "chr1"
    assert row["Mismatches"] == 1
    assert row["Pos"] == 500  # planted exactly after 500 bp of filler
    assert row["Strand"] == "+"


def test_finds_offtarget_on_reverse_strand():
    guide = "GACGATCAGTCAGGATCACC"
    off = "GACGATCAGTCAGGATCACA"  # 1 mismatch
    site = str(Seq(off + "AGG").reverse_complement())  # lives on the - strand
    genome = _filler(300, 3) + site + _filler(300, 4)
    hits = scan_records([("chrX", genome)], [guide], max_mismatches=2, pam="NGG")
    assert (hits["Strand"] == "-").any()
    assert (hits["Mismatches"] == 1).all()


def test_multi_record_and_specificity():
    guide = "GACGATCAGTCAGGATCACC"
    off = "GACGATCAGTCAGGATCACA"
    recs = [("c1", _filler(200, 5) + off + "AGG" + _filler(50, 6)),
            ("c2", _filler(400, 7))]
    hits = scan_records(recs, [guide], max_mismatches=2, pam="NGG")
    assert set(hits["Chrom"]) <= {"c1", "c2"}
    spec = genome_specificity(hits, [guide])
    assert spec.iloc[0]["OffTargetCount"] == len(hits)
    assert 0 <= spec.iloc[0]["CFD_Specificity"] <= 100


def test_chunk_boundary_no_miss_no_double_count():
    # A site straddling a tiny chunk boundary must be found exactly once.
    guide = "GACGATCAGTCAGGATCACC"
    off = "GACGATCAGTCAGGATCACA"
    genome = _filler(95, 8) + off + "AGG" + _filler(200, 9)  # site spans pos 95..117
    big = scan_records([("c", genome)], [guide], max_mismatches=2, pam="NGG", chunk_size=10_000_000)
    small = scan_records([("c", genome)], [guide], max_mismatches=2, pam="NGG", chunk_size=100)
    assert len(big) == len(small) == 1
    assert big.iloc[0]["Pos"] == small.iloc[0]["Pos"] == 95


def test_perfect_match_excluded():
    guide = "GACGATCAGTCAGGATCACC"
    genome = _filler(100, 10) + guide + "AGG" + _filler(100, 11)  # exact on-target
    hits = scan_records([("c", genome)], [guide], max_mismatches=2, pam="NGG")
    assert hits.empty  # 0-mismatch intended target is not reported as off-target


def test_fasta_text_parsing():
    guide = "GACGATCAGTCAGGATCACC"
    off = "GACGATCAGTCAGGATCACA"
    fasta = ">chrTest desc\n" + _filler(60, 12) + off + "AGG" + _filler(60, 13)
    hits = scan_fasta_text(fasta, [guide], max_mismatches=2, pam="NGG")
    assert (hits["Chrom"] == "chrTest").all()
