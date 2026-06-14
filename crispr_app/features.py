"""Feature extraction for on-target efficiency models.

Turns a guide (optionally with its NGGN 3'-context) into a fixed-length numeric
feature vector shared by the interpretable surrogate, the trainable linear
model, and any external model. Features follow the determinants used by Doench
Rule Set 2 / Azimuth and CRISPRedict:

* position-specific one-hot nucleotides (20 positions x {A,C,G,T} = 80)
* GC fraction
* nearest-neighbour Tm (scaled)
* selected dinucleotide counts (GG, TT, GC, AA)
* NGGN 3'-context one-hot (4)

The vector is anchored at the PAM-proximal 3' end so guides of length 18-25
share the same seed alignment.
"""

from __future__ import annotations

import numpy as np

try:
    from scoring import _tm, gc_fraction
except ImportError:  # pragma: no cover
    from .scoring import _tm, gc_fraction

_BASES = "ACGT"
_POSITIONS = 20
_DINUCS = ("GG", "TT", "GC", "AA")


def feature_names() -> list[str]:
    names = [f"pos{p+1}_{b}" for p in range(_POSITIONS) for b in _BASES]
    names += ["gc_fraction", "tm_scaled"]
    names += [f"dinuc_{d}" for d in _DINUCS]
    names += [f"ngg_{b}" for b in _BASES]
    return names


def n_features() -> int:
    return _POSITIONS * 4 + 2 + len(_DINUCS) + 4


def featurize(guide: str, ngg_context: str | None = None) -> np.ndarray:
    guide = guide.upper()
    n = len(guide)
    vec = np.zeros(n_features(), dtype=np.float64)

    # Position one-hot, anchored at the 3' (PAM-proximal) end.
    for offset in range(min(n, _POSITIONS)):
        pos = _POSITIONS - 1 - offset
        base = guide[n - 1 - offset]
        bi = _BASES.find(base)
        if bi >= 0:
            vec[pos * 4 + bi] = 1.0

    idx = _POSITIONS * 4
    vec[idx] = gc_fraction(guide); idx += 1
    vec[idx] = (_tm(guide) - 65.0) / 20.0; idx += 1  # centred/scaled

    for d in _DINUCS:
        vec[idx] = sum(1 for i in range(n - 1) if guide[i:i + 2] == d)
        idx += 1

    if ngg_context:
        bi = _BASES.find(ngg_context[0].upper())
        if bi >= 0:
            vec[idx + bi] = 1.0

    return vec


def featurize_many(guides: list[str], contexts: list[str | None] | None = None) -> np.ndarray:
    contexts = contexts or [None] * len(guides)
    return np.vstack([featurize(g, c) for g, c in zip(guides, contexts)])
