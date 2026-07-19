"""Benchmark harness for on-target efficiency predictions.

Reputable CRISPR predictors are judged by rank correlation (Spearman's rho)
against experimentally measured editing efficiencies on held-out guides
(e.g. Doench 2014/2016, Kim 2019/DeepSpCas9, Wang/Xu). This module provides the
metrics and an evaluation entry point so accuracy is *measured*, not asserted.

Usage:
    from benchmark import evaluate
    records = [{"guide": "ACGT...", "measured": 0.73}, ...]
    print(evaluate(records))   # -> {"n": ..., "spearman": ..., "pearson": ...}

Reference Spearman's rho on common test sets (literature):
    Azimuth / Rule Set 2 ~0.56, DeepSpCas9 ~0.73, DeepHF ~0.74-0.87,
    CRISPRon ~0.80, AttCRISPR ~0.87.
No SciPy dependency: rank/linear correlations are computed with NumPy.
"""

from __future__ import annotations

import numpy as np

try:
    from scoring import on_target_score
except ImportError:  # pragma: no cover
    from .scoring import on_target_score


def _rankdata(a: np.ndarray) -> np.ndarray:
    """Average ranks (ties shared), matching scipy.stats.rankdata('average')."""
    a = np.asarray(a, dtype=float)
    order = a.argsort()
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)
    # Resolve ties to the average rank.
    _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    return (sums / counts)[inv]


def pearson(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y) -> float:
    if len(x) < 2:
        return float("nan")
    return pearson(_rankdata(np.asarray(x, dtype=float)), _rankdata(np.asarray(y, dtype=float)))


def auc(scores, labels) -> float:
    """Rank-based ROC-AUC (Mann-Whitney U), NumPy-only, ties handled by average rank.

    ``labels`` is a 0/1 array. Measures how well ``scores`` separate the two
    classes — the honest metric for the tool's real job (telling effective guides
    from ineffective ones), which is far less noise-limited than rank correlation
    over every guide including the ambiguous middle.
    """
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels).astype(int)
    npos = int((labels == 1).sum())
    nneg = int((labels == 0).sum())
    if npos == 0 or nneg == 0:
        return float("nan")
    r = _rankdata(scores)
    return float((r[labels == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def separation_auc(preds, measured, lo_pct: float = 33.0, hi_pct: float = 67.0) -> float:
    """AUC for separating high- from low-efficiency guides (top vs bottom slice)."""
    preds = np.asarray(preds, dtype=float)
    measured = np.asarray(measured, dtype=float)
    ql, qh = np.percentile(measured, [lo_pct, hi_pct])
    mask = (measured <= ql) | (measured >= qh)
    return auc(preds[mask], (measured[mask] >= qh).astype(int))


def evaluate(records: list[dict], predictor=on_target_score) -> dict:
    """Score a labelled dataset and return correlation metrics.

    ``records`` is a list of ``{"guide": str, "measured": float}`` (optionally a
    "ngg_context" key). Returns sample size and Spearman/Pearson correlations.
    """
    guides, measured = [], []
    for r in records:
        g = r.get("guide", "")
        if "measured" in r and len(g) >= 18:
            guides.append(g)
            measured.append(float(r["measured"]))

    if len(guides) < 2:
        return {"n": len(guides), "spearman": float("nan"), "pearson": float("nan")}

    preds = [predictor(g, r.get("ngg_context")) for g, r in zip(guides, records)]
    return {
        "n": len(guides),
        "spearman": round(spearman(preds, measured), 4),
        "pearson": round(pearson(preds, measured), 4),
    }


def load_crispor_context(path: str) -> list[dict]:
    """Load a CRISPOR-format ``*.context.tab`` efficiency dataset.

    Columns: guide, seq (20mer or 23mer), db, pos, modFreq, longSeq. Returns
    records with ``guide`` (20mer), ``measured`` (modFreq) and, when derivable
    from ``longSeq``, a 35-mer ``mer35`` for CRISPRscan. Download such datasets
    from the CRISPOR paper repo (github.com/maximilianh/crisporPaper, effData/).
    """
    import csv

    records: list[dict] = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("seq") or "").strip().upper()[:20]
            if len(g) != 20 or any(c not in "ACGT" for c in g):
                continue
            try:
                measured = float(row["modFreq"])
            except (KeyError, ValueError):
                continue
            ls = (row.get("longSeq") or "").strip().upper()
            idx = ls.find(g)
            mer35 = None
            if idx >= 6 and idx + 29 <= len(ls):
                w = ls[idx - 6:idx + 29]
                if len(w) == 35 and all(c in "ACGT" for c in w):
                    mer35 = w
            records.append({"guide": g, "measured": measured, "mer35": mer35})
    return records


def _cli() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Benchmark on-target models on a CRISPOR dataset")
    ap.add_argument("context_tab", help="path to a CRISPOR *.context.tab file")
    args = ap.parse_args()

    recs = load_crispor_context(args.context_tab)
    y = [r["measured"] for r in recs]
    h = [on_target_score(r["guide"]) for r in recs]
    print(f"N = {len(recs)}")
    print(f"  heuristic surrogate   Spearman = {spearman(h, y):.3f}")
    try:
        from crisprscan import score_35mer
    except ImportError:  # pragma: no cover
        from .crisprscan import score_35mer
    cs = [(score_35mer(r["mer35"]), r["measured"]) for r in recs if r["mer35"]]
    if len(cs) > 2:
        print(f"  CRISPRscan            Spearman = {spearman([a for a, _ in cs], [b for _, b in cs]):.3f}  (N={len(cs)})")


if __name__ == "__main__":
    _cli()
