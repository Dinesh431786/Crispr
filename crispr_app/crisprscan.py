"""CRISPRscan (Moreno-Mateos) on-target efficiency model.

A fully published, peer-reviewed SpCas9 on-target model (Moreno-Mateos et al.,
Nat. Methods 2015, 12:982). Unlike a black-box deep net, its logistic-regression
weights are public, so it can be reproduced *exactly* with no training data.

Coefficients, intercept, 35-mer framing, and scoring logic are transcribed
verbatim from the open-source CRISPOR implementation
(maximilianh/crisporWebsite, crisporEffScores.py) and validated in the test
suite against the reference vector embedded there
(``...GGGACTGTA`` 35-mer -> 77).

Input framing: a 35-nt window = 6 nt 5' flank + 20 nt guide + 3 nt PAM (NGG)
+ 6 nt 3' flank, with 1-based positions.
"""

from __future__ import annotations

INTERCEPT = 0.183930943629

# (motif, 1-based position, weight) — transcribed verbatim from CRISPOR.
PARAMS_CRISPRSCAN: list[tuple[str, int, float]] = [
    ("AA", 18, -0.097377097),
    ("TT", 18, -0.094424075), ("TT", 13, -0.08618771), ("CT", 26, -0.084264893),
    ("GC", 25, -0.073453609), ("T", 21, -0.068730497), ("TG", 23, -0.066388075),
    ("AG", 23, -0.054338456), ("G", 30, -0.046315914), ("A", 4, -0.042153521),
    ("AG", 34, -0.041935908), ("GA", 34, -0.037797707), ("A", 18, -0.033820432),
    ("C", 25, -0.031648353), ("C", 31, -0.030715556), ("G", 1, -0.029693709),
    ("C", 16, -0.021638609), ("A", 14, -0.018487229), ("A", 11, -0.018287292),
    ("T", 34, -0.017647692), ("AA", 10, -0.016905415), ("A", 19, -0.015576499),
    ("G", 34, -0.014167123), ("C", 30, -0.013182733), ("GA", 31, -0.01227989),
    ("T", 24, -0.011996172), ("A", 15, -0.010595296), ("G", 4, -0.005448869),
    ("GG", 9, -0.00157799), ("T", 23, -0.001422243), ("C", 15, -0.000477727),
    ("C", 26, -0.000368973), ("T", 27, -0.000280845), ("A", 31, 0.00158975),
    ("GT", 18, 0.002391744), ("C", 9, 0.002449224), ("GA", 20, 0.009740799),
    ("A", 25, 0.010506405), ("A", 12, 0.011633235), ("A", 32, 0.012435231),
    ("T", 22, 0.013224035), ("C", 20, 0.015089514), ("G", 17, 0.01549378),
    ("G", 18, 0.016457816), ("T", 30, 0.017263162), ("A", 13, 0.017628924),
    ("G", 19, 0.017916844), ("A", 27, 0.019126815), ("G", 11, 0.020929039),
    ("TG", 3, 0.022949996), ("GC", 3, 0.024681785), ("G", 14, 0.025116714),
    ("GG", 10, 0.026802158), ("G", 12, 0.027591138), ("G", 32, 0.03071249),
    ("A", 22, 0.031930909), ("G", 20, 0.033957008), ("C", 21, 0.034262921),
    ("TT", 17, 0.03492881), ("T", 13, 0.035445171), ("G", 26, 0.036146649),
    ("A", 24, 0.037466478), ("C", 22, 0.03763162), ("G", 16, 0.037970942),
    ("GG", 12, 0.041883009), ("TG", 18, 0.045908991), ("TG", 31, 0.048136812),
    ("A", 35, 0.048596259), ("G", 15, 0.051129717), ("C", 24, 0.052972314),
    ("TG", 15, 0.053372822), ("GT", 11, 0.053678436), ("GC", 9, 0.054171402),
    ("CA", 30, 0.057759851), ("GT", 24, 0.060952114), ("G", 13, 0.061360905),
    ("CA", 24, 0.06221937), ("AG", 10, 0.063717093), ("G", 10, 0.067739182),
    ("C", 13, 0.069495944), ("GT", 31, 0.07342535), ("GG", 13, 0.074355848),
    ("C", 27, 0.079933922), ("G", 27, 0.085151052), ("CC", 21, 0.088919601),
    ("CC", 23, 0.095072286), ("G", 22, 0.10114438), ("G", 24, 0.105488325),
    ("GT", 23, 0.106718563), ("GG", 25, 0.111559441), ("G", 9, 0.114600681),
]


def score_35mer(seq: str) -> int:
    """CRISPRscan score (0-100) for a 35-nt window, matching CRISPOR exactly."""
    if len(seq) != 35:
        raise ValueError("CRISPRscan requires a 35-nt window")
    seq = seq.upper()
    score = INTERCEPT
    for motif, pos, weight in PARAMS_CRISPRSCAN:
        if seq[pos - 1:pos - 1 + len(motif)] == motif:
            score += weight
    return int(100 * score)


def score_from_context(full_seq: str, guide_start: int, guide_len: int = 20) -> float | None:
    """Return CRISPRscan score in [0, 1] for a guide given its flanking context.

    ``guide_start`` is the 0-based index of the guide's 5' end within ``full_seq``
    (same strand as the guide). Requires 6 nt upstream and 6 nt downstream of the
    20 nt guide + 3 nt PAM. Returns None when the guide is a non-standard length
    or too close to a sequence end for the 35-mer window.
    """
    if guide_len != 20:
        return None
    start = guide_start - 6
    end = guide_start + 29  # 6 + 20 + 3 + 6
    if start < 0 or end > len(full_seq):
        return None
    window = full_seq[start:end].upper()
    if len(window) != 35 or any(c not in "ACGT" for c in window):
        return None
    return round(score_35mer(window) / 100.0, 3)
