"""Feature extraction for on-target efficiency models.

Turns a guide into a fixed-length numeric feature vector shared by the
trainable linear model and any external model. The feature set follows the
determinants that dominate Doench Rule Set 2 / Azimuth and CRISPRedict, then
adds several complementary "angles" on the same sequence. Flanking context is
optional, so it still slots into the registry's ``predict(guide, ngg_context)``
interface:

* position-specific one-hot single nucleotides (20 x 4 = 80)
* position-specific one-hot dinucleotides (19 x 16 = 304)  <- the dominant RS2 signal
* per-base nucleotide counts (4)
* GC count, GC count squared (2)
* nearest-neighbour Tm (scaled) (1)
* NGGN 3'-context one-hot (4)
* flanking sequence context — 6 nt upstream + 9 nt downstream one-hot (60)
* position-specific trinucleotides (18 x 64 = 1152)
* gapped/spaced position-specific dinucleotides, gaps 3-7 (1200)  <- long-range coupling
* RC-canonical tetranucleotide counts (136)
* position-independent trinucleotide counts (64)
* energy summaries: regional GC, homopolymer runs, poly-T flag (8)

Benchmarks (5-fold CV on the CRISPRon/Kim set, held-out CRISPOR datasets for
cross-dataset transfer) show this featurizer with a NumPy ridge model reaches
rho ~= 0.766 within-dataset and leads industry tools on cross-dataset mean, with
no heavy dependencies. See BENCHMARKS.md.
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

# Gapped / spaced dinucleotides — long-range positional coupling (a "different
# angle" on the same sequence): the base pair at positions (p, p+gap) for each
# gap in _GAPS, position-specific one-hot. Standard dinucleotides only see the
# adjacent pair (gap 1); these capture the periodic/longer-range preference an
# adjacent-only model misses. Measured +0.008 within-Kim CV and +0.009 on
# held-out cross-datasets on top of the mono/di/tri set.
_GAPS = (3, 4, 5, 6, 7)

# Reverse-complement-canonical tetranucleotide alphabet (136 distinct 4-mers).
# Collapsing each 4-mer with its reverse complement denoises the count and
# transfers better across datasets than the raw 256-way count.
_COMP = str.maketrans("ACGT", "TGCA")


def _rc(mer: str) -> str:
    return mer.translate(_COMP)[::-1]


_TETRA = sorted({min(a + b + c + d, _rc(a + b + c + d))
                 for a in _BASES for b in _BASES for c in _BASES for d in _BASES})
_TETRA_IDX = {m: i for i, m in enumerate(_TETRA)}

# Position-independent trinucleotide vocabulary (a low-dim complement to the
# position-specific trinucleotide block; helps cross-dataset transfer).
_TRI = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_TRI_IDX = {m: i for i, m in enumerate(_TRI)}


def _max_run(guide: str, base: str) -> int:
    best = cur = 0
    for ch in guide:
        cur = cur + 1 if ch == base else 0
        if cur > best:
            best = cur
    return best


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
    names += [f"gap{gap}_pos{p+1}_{a}{b}" for gap in _GAPS
              for p in range(_POSITIONS - gap) for a in _BASES for b in _BASES]
    names += [f"tetracanon_{m}" for m in _TETRA]
    names += [f"tricount_{m}" for m in _TRI]
    names += ["gc_5p10", "gc_3p10", "gc_seed", "run_A", "run_C", "run_G", "run_T", "tttt"]
    return names


def n_features() -> int:
    return (_POSITIONS * 4 + (_POSITIONS - 1) * 16 + 4 + 3 + 4 + _N_FLANK + (_POSITIONS - 2) * 64
            + sum((_POSITIONS - gap) * 16 for gap in _GAPS) + len(_TETRA) + len(_TRI) + 8)


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
    idx += (_POSITIONS - 2) * 64

    # Gapped/spaced position-specific dinucleotides — long-range positional
    # coupling the adjacent-only blocks miss (see _GAPS note above).
    for gap in _GAPS:
        for p in range(min(n - gap, _POSITIONS - gap)):
            a, b = _BI.get(guide[p], -1), _BI.get(guide[p + gap], -1)
            if a >= 0 and b >= 0:
                vec[idx + p * 16 + a * 4 + b] = 1.0
        idx += (_POSITIONS - gap) * 16

    # RC-canonical tetranucleotide counts (longer local motifs, denoised).
    for p in range(min(n - 3, _POSITIONS - 3)):
        ci = _TETRA_IDX.get(_rc(guide[p:p + 4]) if guide[p:p + 4] > _rc(guide[p:p + 4])
                            else guide[p:p + 4])
        if ci is not None:
            vec[idx + ci] += 1.0
    idx += len(_TETRA)

    # Position-independent trinucleotide counts (low-dim, transfer-friendly).
    for p in range(min(n - 2, _POSITIONS - 2)):
        ti = _TRI_IDX.get(guide[p:p + 3])
        if ti is not None:
            vec[idx + ti] += 1.0
    idx += len(_TRI)

    # Energy / melting summaries: regional GC (5' half, 3' half, seed 13-20),
    # homopolymer run length per base, and a poly-T (TTTT) terminator flag.
    vec[idx] = guide[:10].count("G") + guide[:10].count("C")
    vec[idx + 1] = guide[10:].count("G") + guide[10:].count("C")
    vec[idx + 2] = guide[12:20].count("G") + guide[12:20].count("C")
    vec[idx + 3] = _max_run(guide, "A")
    vec[idx + 4] = _max_run(guide, "C")
    vec[idx + 5] = _max_run(guide, "G")
    vec[idx + 6] = _max_run(guide, "T")
    vec[idx + 7] = 1.0 if "TTTT" in guide else 0.0

    return vec


def featurize_many(guides, contexts=None, ups=None, downs=None) -> np.ndarray:
    contexts = contexts or [None] * len(guides)
    ups = ups or [""] * len(guides)
    downs = downs or [""] * len(guides)
    return np.vstack([featurize(g, c, u, d) for g, c, u, d in zip(guides, contexts, ups, downs)])
