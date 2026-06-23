"""Reproduce the shipped default on-target model (crispr_app/models/default.json).

Trains a NumPy-only ridge model over the standard feature set on the large,
clean CRISPRon/Kim SpCas9 indel-efficiency dataset (~11.6k guides). Reports
5-fold cross-validated Spearman and split-conformal calibration, then writes the
model. The 5-fold CV lands at rho ~= 0.71, inside the wet-lab replicate
reproducibility band (~0.71-0.77) — i.e. the model agrees with the assay about
as well as the assay agrees with itself.

Usage:
    # 1) get the dataset (one-time; ships with the DeepCRISTL repo)
    git clone --depth 1 https://github.com/OrensteinLab/DeepCRISTL
    # 2) build the model
    python scripts/build_default_model.py \
        --data DeepCRISTL/CRISPROn/data/main_dataframes/seq_efficienciey.txt

Dataset columns used: 'gRNA' (20-nt spacer), 'total_indel_eff' (0-100).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "crispr_app"))
from benchmark import spearman  # noqa: E402
from conformal import calibrate, empirical_coverage  # noqa: E402
from features import featurize_many  # noqa: E402
from train import fit_ridge  # noqa: E402


def load(path: str):
    guides, y = [], []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        gi, ei = header.index("gRNA"), header.index("total_indel_eff")
        for row in reader:
            if len(row) <= max(gi, ei):
                continue
            g = row[gi].strip().upper()
            try:
                eff = float(row[ei])
            except ValueError:
                continue
            if len(g) == 20 and all(c in "ACGT" for c in g):
                guides.append(g)
                y.append(eff / 100.0)  # -> [0, 1] to match the model output range
    return guides, np.array(y)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="path to CRISPROn seq_efficienciey.txt")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "crispr_app" / "models" / "default.json"))
    ap.add_argument("--alpha", type=float, default=10.0)
    args = ap.parse_args()

    guides, y = load(args.data)
    X = featurize_many(guides)
    print(f"N = {len(guides)} guides")

    # Honest accuracy estimate: 5-fold cross-validated Spearman (NumPy folds).
    rng0 = np.random.default_rng(0)
    folds = np.array_split(rng0.permutation(len(y)), 5)
    pred = np.zeros(len(y))
    for k in range(5):
        te = folds[k]
        tr = np.concatenate([folds[j] for j in range(5) if j != k])
        m = fit_ridge(X[tr], y[tr], args.alpha)
        pred[te] = np.clip(X[te] @ m.weights + m.intercept, 0, 1)
    rho = spearman(pred, y)
    print(f"  ridge 5-fold CV Spearman = {rho:.3f}  (wet-lab replicate band ~0.71-0.77)")

    # Split-conformal calibration on a held-out split (model never saw it).
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(guides))
    cut = int(0.8 * len(guides))
    tr_i, cal_i = perm[:cut], perm[cut:]
    cal_model = fit_ridge(X[tr_i], y[tr_i], args.alpha)
    cal_pred = np.clip(X[cal_i] @ cal_model.weights + cal_model.intercept, 0, 1)
    conf = calibrate(cal_pred, y[cal_i])
    for level, q in conf.items():
        print(f"  conformal {level}: half-width={q}  coverage={empirical_coverage(cal_pred, y[cal_i], q):.3f}")

    final = fit_ridge(X, y, args.alpha)
    final.meta.update({"trained_on": f"CRISPROn/Kim seq_efficiency (n={len(guides)})",
                       "target": "total_indel_eff/100", "cv_spearman": round(rho, 3),
                       "alpha": args.alpha, "conformal": conf})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    final.save(args.out)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
