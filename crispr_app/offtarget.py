"""Off-target scoring: per-site CFD + MIT/Hsu, and aggregate guide specificity.

Two complementary, literature-standard models are implemented:

* **CFD** (Cutting Frequency Determination; Doench et al. 2016, Nat. Biotechnol.
  34:184). Position- and substitution-specific RNA:DNA mismatch penalties plus a
  PAM penalty. Benchmarks show CFD outperforms the older MIT score.
* **MIT / Hsu score** (Hsu et al. 2013, Nat. Biotechnol. 31:827). The canonical
  position-weighted model (``calcHitScore``) retained for cross-validation and
  because many pipelines still report it.

Per-guide *specificity* is aggregated CRISPOR-style as
``100 / (100 + Σ off_target_scores)`` so that a guide with many strong
off-targets approaches 0 and a unique guide approaches 100.
"""

from __future__ import annotations

# --- CFD mismatch matrix (Doench et al. 2016, Suppl. Table 19) --------------
# Keys: "r<RNA>:d<DNA-complement-of-target>"; value is the retained activity.
CFD_SCORES: dict[str, float] = {
    "rA:dA": 1.0, "rA:dC": 1.0, "rA:dG": 0.857, "rA:dT": 1.0,
    "rC:dA": 1.0, "rC:dC": 0.913, "rC:dG": 1.0, "rC:dT": 1.0,
    "rG:dA": 1.0, "rG:dC": 1.0, "rG:dG": 0.714, "rG:dT": 0.9,
    "rU:dA": 1.0, "rU:dC": 0.957, "rU:dG": 0.857, "rU:dT": 1.0,
}

# PAM activity for the two 3' nucleotides of the NGG-class PAM (Doench 2016).
PAM_SCORES: dict[str, float] = {
    "AA": 0.0, "AC": 0.0, "AG": 0.259, "AT": 0.0,
    "CA": 0.0, "CC": 0.0, "CG": 0.107, "CT": 0.0,
    "GA": 0.069, "GC": 0.022, "GG": 1.0, "GT": 0.016,
    "TA": 0.0, "TC": 0.0, "TG": 0.039, "TT": 0.0,
}

# --- MIT / Hsu 2013 position weights (5'->3', index 19 = PAM-proximal) ------
HSU_WEIGHTS: tuple[float, ...] = (
    0.000, 0.000, 0.014, 0.000, 0.000,
    0.395, 0.317, 0.000, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732,
    0.828, 0.615, 0.804, 0.685, 0.583,
)


def calculate_cfd_score(guide_seq: str, off_target_seq: str, pam: str) -> float:
    """CFD score in [0, 1] for an aligned guide / off-target / PAM triple."""
    score = 1.0
    guide_seq = guide_seq.upper().replace("T", "U")
    off_target_seq = off_target_seq.upper()

    for i in range(min(len(guide_seq), len(off_target_seq))):
        if guide_seq[i] != off_target_seq[i]:
            key = f"r{guide_seq[i]}:d{off_target_seq[i]}"
            score *= CFD_SCORES.get(key, 1.0)

    pam_key = pam[1:3] if len(pam) >= 3 else pam[-2:]
    score *= PAM_SCORES.get(pam_key.upper(), 0.0)
    return round(score, 5)


def mit_hit_score(guide_seq: str, off_target_seq: str) -> float:
    """MIT/Hsu single off-target hit score in [0, 100] (Hsu et al. 2013)."""
    guide_seq = guide_seq.upper()
    off_target_seq = off_target_seq.upper()
    n = min(len(guide_seq), len(off_target_seq), len(HSU_WEIGHTS))

    mm_positions: list[int] = []
    score1 = 1.0
    for pos in range(n):
        if guide_seq[pos] != off_target_seq[pos]:
            score1 *= 1.0 - HSU_WEIGHTS[pos]
            mm_positions.append(pos)

    mm_count = len(mm_positions)
    if mm_count == 0:
        return 100.0

    if mm_count < 2:
        score2 = 1.0
    else:
        dists = [mm_positions[i] - mm_positions[i - 1] for i in range(1, mm_count)]
        avg_dist = sum(dists) / len(dists)
        score2 = 1.0 / (((19 - avg_dist) / 19.0) * 4 + 1)

    score3 = 1.0 / (mm_count ** 2)
    return round(score1 * score2 * score3 * 100.0, 4)


def aggregate_specificity(off_target_scores: list[float], scale: float = 1.0) -> float:
    """CRISPOR-style aggregate specificity in [0, 100].

    ``off_target_scores`` are per-site scores (CFD in [0,1] -> use scale=100;
    MIT hit scores in [0,100] -> use scale=1). Higher = more specific guide;
    a guide with no off-targets returns 100, matching the CRISPOR convention.
    """
    total = sum(s * scale for s in off_target_scores)
    return round(10000.0 / (100.0 + total), 2)
