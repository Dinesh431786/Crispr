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


def feature_names() -> list[str]:
    names = [f"pos{p+1}_{b}" for p in range(_POSITIONS) for b in _BASES]
    names += [f"dipos{p+1}_{a}{b}" for p in range(_POSITIONS - 1) for a in _BASES for b in _BASES]
    names += [f"count_{b}" for b in _BASES]
    names += ["gc_count", "gc_count_sq", "tm_scaled"]
    names += [f"ngg_{b}" for b in _BASES]
    return names


def n_features() -> int:
    return _POSITIONS * 4 + (_POSITIONS - 1) * 16 + 4 + 3 + 4


def featurize(guide: str, ngg_context: str | None = None) -> np.ndarray:
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

    return vec


def featurize_many(guides: list[str], contexts: list[str | None] | None = None) -> np.ndarray:
    contexts = contexts or [None] * len(guides)
    return np.vstack([featurize(g, c) for g, c in zip(guides, contexts)])
