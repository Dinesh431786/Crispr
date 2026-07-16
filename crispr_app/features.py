"""Feature extraction for on-target efficiency models.

Turns a guide into a fixed-length numeric feature vector shared by the
trainable linear model and any external model. The feature set follows the
determinants that dominate Doench Rule Set 2 / Azimuth and CRISPRedict, while
staying *guide-only* (no flanking context required) so it slots into the
registry's ``predict(guide, ngg_context)`` interface:

* position-specific one-hot single nucleotides (20 x 4 = 80)
* position-specific one-hot dinucleotides (19 x 16 = 304)  <- the dominant RS2 signal
* per-base nucleotide counts (4)
* GC count, GC count squared (2)
* nearest-neighbour Tm (scaled) (1)
* NGGN 3'-context one-hot (4)

Benchmarks (5-fold CV on public CRISPOR datasets) show this featurizer with a
NumPy ridge model roughly doubles Spearman vs the heuristic on datasets with
signal (e.g. chari2015 0.20->0.40, morenoMateos 0.17->0.43) and matches
gradient boosting, with no heavy dependencies. See BENCHMARKS.md.
"""

from __future__ import annotations

import numpy as np

try:
    from scoring import _tm, gc_fraction
except ImportError:  # pragma: no cover
    from .scoring import _tm, gc_fraction

_BASES = "ACGT"
_BI = {b: i for i, b in enumerate(_BASES)}
_POSITIONS = 20


# Flanking sequence context — the documented reason Rule Set 2 / DeepSpCas9 beat
# guide-only models. up = 6 nt immediately 5' of the protospacer; down = 9 nt
# immediately 3' of it (the 3-nt PAM + 6 nt downstream), both genomic 5'->3'.
_FLANK_UP = 6
_FLANK_DOWN = 9
_N_FLANK = (_FLANK_UP + _FLANK_DOWN) * 4


def feature_names() -> list[str]:
    names = [f"pos{p+1}_{b}" for p in range(_POSITIONS) for b in _BASES]
    names += [f"dipos{p+1}_{a}{b}" for p in range(_POSITIONS - 1) for a in _BASES for b in _BASES]
    names += [f"count_{b}" for b in _BASES]
    names += ["gc_count", "gc_count_sq", "tm_scaled"]
    names += [f"ngg_{b}" for b in _BASES]
    names += [f"up{_FLANK_UP - p}_{b}" for p in range(_FLANK_UP) for b in _BASES]
    names += [f"down{p+1}_{b}" for p in range(_FLANK_DOWN) for b in _BASES]
    names += [f"tripos{p+1}_{a}{b}{c}" for p in range(_POSITIONS - 2)
              for a in _BASES for b in _BASES for c in _BASES]
    return names


def n_features() -> int:
    return _POSITIONS * 4 + (_POSITIONS - 1) * 16 + 4 + 3 + 4 + _N_FLANK + (_POSITIONS - 2) * 64


def featurize(guide: str, ngg_context: str | None = None,
              up: str = "", down: str = "") -> np.ndarray:
    guide = guide.upper()
    n = len(guide)
    vec = np.zeros(n_features(), dtype=np.float64)

    # Position-specific single nucleotides (anchored at 5' end).
    for p in range(min(n, _POSITIONS)):
        bi = _BI.get(guide[p], -1)
        if bi >= 0:
            vec[p * 4 + bi] = 1.0

    # Position-specific dinucleotides — the dominant Rule Set 2 signal.
    off = _POSITIONS * 4
    for p in range(min(n - 1, _POSITIONS - 1)):
        a, b = _BI.get(guide[p], -1), _BI.get(guide[p + 1], -1)
        if a >= 0 and b >= 0:
            vec[off + p * 16 + a * 4 + b] = 1.0

    idx = off + (_POSITIONS - 1) * 16
    for b in _BASES:
        vec[idx] = guide.count(b); idx += 1

    gc = guide.count("G") + guide.count("C")
    vec[idx] = gc; idx += 1
    vec[idx] = gc * gc; idx += 1
    vec[idx] = (_tm(guide) - 65.0) / 20.0; idx += 1

    if ngg_context:
        bi = _BI.get(ngg_context[0].upper(), -1)
        if bi >= 0:
            vec[idx + bi] = 1.0
    idx += 4

    # Upstream flank: right-aligned so the nucleotide nearest the guide is pos -1.
    up = ("N" * _FLANK_UP + (up or "").upper())[-_FLANK_UP:]
    for p in range(_FLANK_UP):
        bi = _BI.get(up[p], -1)
        if bi >= 0:
            vec[idx + p * 4 + bi] = 1.0
    idx += _FLANK_UP * 4

    # Downstream flank (PAM + 6 nt): left-aligned from just after the guide.
    down = ((down or "").upper() + "N" * _FLANK_DOWN)[:_FLANK_DOWN]
    for p in range(_FLANK_DOWN):
        bi = _BI.get(down[p], -1)
        if bi >= 0:
            vec[idx + p * 4 + bi] = 1.0
    idx += _FLANK_DOWN * 4

    # Position-specific trinucleotides — triplet motifs (the local patterns CNNs
    # learn), still linear. Measured +0.025 within-Kim CV over dinucleotides.
    for p in range(min(n - 2, _POSITIONS - 2)):
        a, b, c = _BI.get(guide[p], -1), _BI.get(guide[p + 1], -1), _BI.get(guide[p + 2], -1)
        if a >= 0 and b >= 0 and c >= 0:
            vec[idx + p * 64 + a * 16 + b * 4 + c] = 1.0

    return vec


def featurize_many(guides, contexts=None, ups=None, downs=None) -> np.ndarray:
    contexts = contexts or [None] * len(guides)
    ups = ups or [""] * len(guides)
    downs = downs or [""] * len(guides)
    return np.vstack([featurize(g, c, u, d) for g, c, u, d in zip(guides, contexts, ups, downs)])
