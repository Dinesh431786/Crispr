"""sgRNA self-complementarity / secondary-structure propensity.

Spacers that fold back on themselves sequester the seed region and reduce Cas9
on-target activity (Thyme et al. 2016, eLife; Wong et al. 2015, WU-CRISPR,
Genome Biol. 16:218; Moreno-Mateos 2015). We quantify this with a dependency-free
Nussinov maximum-base-pairing computation over the 20-nt spacer (treated as RNA,
allowing Watson-Crick and G·U wobble pairs), normalised to a [0, 1] propensity.

Higher propensity = more internal secondary structure = expected lower activity.
This is a deterministic, O(n^3) (n=20, trivial) feature with no external RNA
folding dependency.
"""

from __future__ import annotations

from functools import lru_cache

_WC = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G")}


def _pairs(a: str, b: str) -> bool:
    return (a, b) in _WC


@lru_cache(maxsize=4096)
def max_base_pairs(seq: str, min_loop: int = 3) -> int:
    """Nussinov maximum number of nested base pairs in ``seq`` (min hairpin loop)."""
    seq = seq.upper()
    n = len(seq)
    if n < min_loop + 2:
        return 0
    dp = [[0] * n for _ in range(n)]
    for span in range(min_loop + 1, n):
        for i in range(n - span):
            j = i + span
            best = dp[i + 1][j]  # i unpaired
            for t in range(i + min_loop + 1, j + 1):
                if _pairs(seq[i], seq[t]):
                    left = dp[i + 1][t - 1] if t - 1 >= i + 1 else 0
                    right = dp[t + 1][j] if t + 1 <= j else 0
                    if left + 1 + right > best:
                        best = left + 1 + right
            dp[i][j] = best
    return dp[0][n - 1]


def self_complementarity(guide: str) -> float:
    """Self-folding propensity in [0, 1] (fraction of the max possible base pairs)."""
    guide = guide.upper()
    n = len(guide)
    if n < 6:
        return 0.0
    return round(max_base_pairs(guide) / (n // 2), 3)
