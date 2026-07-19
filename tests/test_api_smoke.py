"""End-to-end API smoke tests: spin the real app and hit every endpoint, so a
refactor can't silently break one. Asserts status + response shape, not science
(that's covered by the unit tests)."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "crispr_app"))
from fastapi.testclient import TestClient  # noqa: E402
import main  # noqa: E402

client = TestClient(main.app)

SEQ = ("ATGGCCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCC"
       "GCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAG"
       "CTGCAAGAACTCTTCCTCACGCGCGTCGGGCTCGACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCC")


def test_health_and_models():
    assert client.get("/health").json()["status"] == "ok"
    m = client.get("/api/models").json()
    assert "active" in m and "available" in m


@pytest.fixture(scope="module")
def guides():
    r = client.post("/api/design", json={"dna_sequence": SEQ, "pam": "NGG"})
    assert r.status_code == 200
    data = r.json()
    assert data["count"] > 0 and data["top_guides"]
    return [g["gRNA"] for g in data["top_guides"]]


def test_design_ranking_strategies_differ():
    bal = client.post("/api/design", json={"dna_sequence": SEQ, "ranking_strategy": "balanced"}).json()
    con = client.post("/api/design", json={"dna_sequence": SEQ, "ranking_strategy": "conservative"}).json()
    assert con["ranking_strategy"] == "conservative"
    assert [g["gRNA"] for g in bal["top_guides"]] != [g["gRNA"] for g in con["top_guides"]]


def test_explain_has_interval(guides):
    out = client.post("/api/explain", json={"guide": guides[0]}).json()
    assert "interval" in out and 0 <= out["interval"]["low"] <= out["interval"]["high"] <= 1


def test_sensitivity(guides):
    g = guides[0]
    out = client.post("/api/sensitivity", json={"guide": g}).json()
    assert len(out["substitutions"]) == len(g) * 3   # every position × 3 alternative bases
    assert len(out["per_position"]) == len(g)
    assert "best" in out and "worst" in out


def test_simulate_has_spectrum(guides):
    out = client.post("/api/simulate", json={"dna_sequence": SEQ, "guide": guides[0]}).json()
    assert "spectrum" in out
    assert 0.0 <= out["spectrum"]["frameshift_probability"] <= 1.0


def test_base_edit(guides):
    out = client.post("/api/base-edit", json={"guides": guides[:10], "editor": "any"}).json()
    assert "results" in out and isinstance(out["results"], list)


def test_multiplex(guides):
    out = client.post("/api/design-multiplex", json={"dna_sequence": SEQ, "n_guides": 4}).json()
    assert out["summary"]["selected"] >= 1
    assert "max_pairwise_identity" in out["summary"]


def test_prime_design():
    out = client.post("/api/prime-design", json={"dna_sequence": SEQ, "target_pos": 60, "desired_edit": "A"}).json()
    assert "pegRNAs" in out


def test_offtargets(guides):
    out = client.post("/api/offtargets", json={"guides": guides[:5], "background_sequence": SEQ}).json()
    assert "specificity" in out


def test_upload_fasta():
    out = client.post("/api/upload-fasta", json={"contents": ">x\n" + SEQ}).json()
    assert out["length"] == len(SEQ)


def test_bad_sequence_400():
    r = client.post("/api/design", json={"dna_sequence": "ZZZZZZZZZZZZZZZZZZZZZZZZ"})
    assert r.status_code == 400
