"""Reproduce the shipped default on-target model (crispr_app/models/default.json).

Trains a NumPy ridge model over the standard feature set on a pool of public
human SpCas9 efficiency datasets, with within-dataset rank normalisation so the
targets are comparable. Reports leave-one-dataset-out Spearman (the honest
estimate of performance on an unseen context) and writes the model.

Usage:
    # 1) get the datasets (one-time)
    git clone --depth 1 https://github.com/maximilianh/crisporPaper
    # 2) build the model
    python scripts/build_default_model.py --effdata crisporPaper/effData

Datasets pooled: doench2016_hg19, chari2015Train, doench2014-Hs.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "crispr_app"))
from benchmark import spearman  # noqa: E402
from features import featurize_many  # noqa: E402
from train import fit_ridge  # noqa: E402

DATASETS = ["doench2016_hg19", "chari2015Train", "doench2014-Hs"]


def load(effdata: str):
    guides, y, group = [], [], []
    for si, name in enumerate(DATASETS):
        gg, yy = [], []
        with open(f"{effdata}/{name}.context.tab") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                g = row["seq"].strip().upper()[:20]
                if len(g) != 20 or any(c not in "ACGT" for c in g):
                    continue
                try:
                    yy.append(float(row["modFreq"]))
                except ValueError:
                    continue
                gg.append(g)
        yy = np.array(yy)
        ranks = yy.argsort().argsort() / (len(yy) - 1)  # within-dataset percentile
        guides += gg; y += list(ranks); group += [si] * len(gg)
        print(f"  {name:<18} N={len(gg)}")
    return np.array(guides), np.array(y), np.array(group)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--effdata", required=True, help="path to crisporPaper/effData")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "crispr_app" / "models" / "default.json"))
    ap.add_argument("--alpha", type=float, default=20.0)
    args = ap.parse_args()

    guides, y, group = load(args.effdata)
    X = featurize_many(list(guides))
    print(f"pooled N = {len(guides)}")

    for si, name in enumerate(DATASETS):
        tr, te = group != si, group == si
        m = fit_ridge(X[tr], y[tr], args.alpha)
        pred = np.clip(X[te] @ m.weights + m.intercept, 0, 1)
        print(f"  LODO held-out {name:<18} rho = {spearman(pred, y[te]):.3f}")

    final = fit_ridge(X, y, args.alpha)
    final.meta.update({"trained_on": DATASETS, "n": int(len(guides)),
                       "target": "within-dataset percentile", "alpha": args.alpha})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    final.save(args.out)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
