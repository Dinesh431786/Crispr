"""On-target gRNA efficiency scoring engine.

This module implements a deterministic, research-informed surrogate for SpCas9
on-target cleavage efficiency. It is *not* a re-trained black-box model; instead
it combines the dominant, well-characterised sequence determinants reported in
the literature into a transparent, calibrated logistic model:

* Position-specific single-nucleotide preferences along the 20-nt spacer
  (Doench et al. 2014, Nat. Biotechnol.; Doench et al. 2016 "Rule Set 2"/Azimuth,
  Nat. Biotechnol. 34:184; Xu et al. 2015, Genome Res.).
* PAM-proximal / seed effects (G favoured at the PAM-proximal position, U/T
  strongly disfavoured) and the NGGN 3'-context (NGGH > NGGG).
* GC content with a quadratic optimum near 50-65 %.
* Nearest-neighbour melting temperature (Bio.SeqUtils.MeltingTemp) with a
  Gaussian optimum, reflecting the duplex-stability term used by Rule Set 2.
* Homopolymer / poly-U penalties (Pol III termination, reduced loading).

The output is squashed through a logistic function to a calibrated [0, 1]
efficiency probability. Modern gold-standard predictors (Rule Set 3, DeepHF,
DeepSpCas9, CRISPRon) are deep models trained on >10^5 guides; this engine is a
fast, dependency-light approximation intended for ranking, and wet-lab
validation remains essential.
"""

from __future__ import annotations

import math
from functools import lru_cache

from Bio.SeqUtils import MeltingTemp as mt

# --- Calibration constants -------------------------------------------------
# Optimal GC fraction window and nearest-neighbour Tm (deg C) for SpCas9 sgRNAs.
GC_OPTIMUM = 0.575          # midpoint of the favourable 50-65% window
GC_TOLERANCE = 0.18         # width of the quadratic GC penalty
TM_OPTIMUM = 65.0           # favourable spacer Tm (Doench RS2 duplex term)
TM_TOLERANCE = 9.0          # Gaussian width for the Tm term

# Position-specific single-nucleotide log-odds contributions for a 20-nt spacer.
# Index 0 = PAM-distal 5' end, index 19 = PAM-proximal 3' end (adjacent to NGG).
# Values encode the *direction and relative magnitude* of effects repeatedly
# reported across Doench 2014/2016 and Xu 2015; they are intentionally small so
# that GC/Tm structure dominates, matching observed feature importances.
_POS_WEIGHTS: dict[str, list[float]] = {
    "A": [0.04, 0.03, 0.02, 0.02, 0.01, 0.0, 0.0, 0.0, -0.01, -0.01,
          0.0, 0.0, 0.02, 0.03, 0.03, 0.02, 0.04, 0.06, 0.02, -0.02],
    "C": [0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.02, 0.03, 0.03, 0.03,
          0.03, 0.03, 0.02, 0.0, -0.02, -0.03, -0.04, -0.05, -0.06, -0.10],
    "G": [0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05,
          0.05, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.10, 0.12, 0.18],
    "T": [-0.02, -0.02, -0.02, -0.02, -0.02, -0.03, -0.03, -0.03, -0.04, -0.04,
          -0.04, -0.05, -0.06, -0.07, -0.08, -0.09, -0.11, -0.13, -0.16, -0.22],
}

# Dinucleotide log-odds adjustments (sparse; only the strongest reported pairs).
_DINUC_WEIGHTS: dict[str, float] = {
    "GG": 0.03,   # favourable, especially PAM-proximal
    "GC": 0.02,
    "TT": -0.06,  # poly-pyrimidine, disfavoured
    "TA": -0.03,
    "AA": -0.02,
}

# NGGN 3'-context: nucleotide immediately 3' of the NGG PAM (Doench 2016).
_NGGN_WEIGHTS: dict[str, float] = {"A": 0.03, "C": 0.02, "T": 0.02, "G": -0.05}

_INTERCEPT = 0.15


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / max(len(seq), 1)


@lru_cache(maxsize=4096)
def _tm(seq: str) -> float:
    try:
        return float(mt.Tm_NN(seq))
    except Exception:
        # Fallback: Wallace rule approximation if NN fails on odd inputs.
        s = seq.upper()
        return 2 * (s.count("A") + s.count("T")) + 4 * (s.count("G") + s.count("C"))


def on_target_score(guide: str, ngg_context: str | None = None) -> float:
    """Calibrated SpCas9 on-target efficiency in [0, 1].

    ``guide`` is the 20-nt (or len-flexible) spacer, 5'->3'. ``ngg_context`` is
    the optional single nucleotide immediately 3' of the NGG PAM (the "N" of
    NGGN) used for the documented 3'-context effect.
    """
    guide = guide.upper()
    n = len(guide)
    if n < 18:
        return 0.0

    logit = _INTERCEPT

    # Position-specific single-nucleotide terms, anchored at the PAM-proximal 3'
    # end so that guides of length 18-25 share the same seed alignment.
    for offset in range(min(n, 20)):
        pos = 19 - offset           # 19 = PAM-proximal
        base = guide[n - 1 - offset]
        col = _POS_WEIGHTS.get(base)
        if col is not None:
            logit += col[pos]

    # Dinucleotide terms across the spacer.
    for i in range(n - 1):
        logit += _DINUC_WEIGHTS.get(guide[i:i + 2], 0.0)

    # GC content: quadratic penalty away from the optimum window.
    gc = gc_fraction(guide)
    logit -= 2.2 * ((gc - GC_OPTIMUM) / GC_TOLERANCE) ** 2

    # Nearest-neighbour Tm: Gaussian reward around the optimum.
    tm_val = _tm(guide)
    logit += 0.45 * math.exp(-((tm_val - TM_OPTIMUM) / TM_TOLERANCE) ** 2)

    # Homopolymer / poly-U penalties (Pol III termination, synthesis issues).
    for base in "ACGT":
        if base * 4 in guide:
            logit -= 0.25 if base == "T" else 0.15

    # NGGN 3'-context effect.
    if ngg_context:
        logit += _NGGN_WEIGHTS.get(ngg_context[0].upper(), 0.0)

    return round(_sigmoid(logit), 3)


def doench_rule_set2_surrogate(guide: str, ngg_context: str | None = None) -> float:
    """Alias preserved for clarity in API responses / downstream callers."""
    return on_target_score(guide, ngg_context)
