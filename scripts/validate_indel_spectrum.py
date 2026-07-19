"""Validate our lightweight frameshift predictor against the trained Lindel model.

Our `analysis.indel_distribution` predicts an out-of-frame (frameshift) probability
from two sequence-derivable repair channels (MMEJ microhomology deletions + the
templated +1 insertion). This script measures how well that agrees with **Lindel**
(Chen et al. 2019, Nucleic Acids Res.), a model trained on real repair-outcome
data — the same "agree with the field-standard tool on identical inputs" test we
use for on-target (see benchmark_competitors.py).

Reproduce:
    git clone --depth 1 https://github.com/shendurelab/Lindel
    git clone --depth 1 https://github.com/maximilianh/crisporPaper   # for test sequences
    python scripts/validate_indel_spectrum.py \
        --lindel Lindel --effdata crisporPaper/effData

Result (n=6058 held-out sites): Spearman rho ~= 0.54 vs Lindel. We under-predict
absolute frameshift (~0.66 vs ~0.76) because we cover only the two dominant
predictable channels; the *ranking* agreement is what this validates.
"""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "crispr_app"))
from analysis import indel_distribution  # noqa: E402
from benchmark import spearman  # noqa: E402

PAMS = ("AGG", "TGG", "CGG", "GGG")
DATASETS = ("doench2016_hg19", "chari2015Train", "morenoMateos2015")


def _build_60mers(effdata: str, name: str) -> list[str]:
    """Lindel input: 60 nt with the spacer at [13:33], NGG PAM at [33:36], cut=30."""
    path = Path(effdata) / f"{name}.scores.tab"
    if not path.is_file():
        return []
    seqs = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        g = (r.get("seq") or "").upper()[:20]
        ls = (r.get("longSeq100Bp") or "").upper()
        if len(g) != 20 or any(c not in "ACGT" for c in g):
            continue
        j = ls.find(g)
        if j < 13 or j + 47 > len(ls):
            continue
        m = ls[j - 13:j + 47]
        if len(m) == 60 and m[13:33] == g and m[33:36] in PAMS:
            seqs.append(m)
    return seqs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lindel", required=True, help="path to the cloned Lindel repo")
    ap.add_argument("--effdata", required=True, help="path to crisporPaper/effData")
    ap.add_argument("--limit", type=int, default=0, help="cap sequences per dataset (0 = all)")
    args = ap.parse_args()

    sys.path.insert(0, args.lindel)
    from Lindel.Predictor import gen_prediction  # noqa: E402
    wb = pickle.load(open(Path(args.lindel) / "Lindel" / "Model_weights.pkl", "rb"))
    prereq = pickle.load(open(Path(args.lindel) / "Lindel" / "model_prereq.pkl", "rb"))

    pooled_ours, pooled_lindel = [], []
    for name in DATASETS:
        seqs = _build_60mers(args.effdata, name)
        if args.limit:
            seqs = seqs[:args.limit]
        ours, lindel = [], []
        for m in seqs:
            res = gen_prediction(m, wb, prereq)
            if isinstance(res, str):   # "Error: No PAM ..."
                continue
            _, fs = res
            ours.append(indel_distribution(m, 30)["frameshift_probability"])
            lindel.append(float(fs))
        if len(ours) >= 20:
            rho = spearman(np.array(ours), np.array(lindel))
            print(f"{name:20} n={len(ours):5d}  Spearman(ours vs Lindel) = {rho:.3f}")
            pooled_ours += ours
            pooled_lindel += lindel

    if pooled_ours:
        rho = spearman(np.array(pooled_ours), np.array(pooled_lindel))
        print(f"{'POOLED':20} n={len(pooled_ours):5d}  Spearman = {rho:.3f}")
        print(f"  mean frameshift: ours={np.mean(pooled_ours):.3f}  Lindel={np.mean(pooled_lindel):.3f}")


if __name__ == "__main__":
    main()
