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


TARGET_COLUMN = {"cutting": "total_indel_eff", "oof": "out_of_frame efficiency"}


def _load_context(data_path: str) -> dict:
    """Surrogate-ID -> construct sequence (carries the flanking genomic context)."""
    ctx = {}
    p = Path(data_path).parent / "full_seq.txt"
    if not p.is_file():
        return ctx
    with open(p) as fh:
        r = csv.reader(fh, delimiter="\t"); h = next(r)
        si = h.index("Surrogate ID")
        qi = next(i for i, c in enumerate(h) if "sequ" in c.lower())
        for row in r:
            if len(row) > max(si, qi):
                ctx[row[si]] = row[qi].strip().upper()
    return ctx


def _flanks(guide: str, construct: str) -> tuple[str, str]:
    """6 nt upstream, 9 nt (PAM+6) downstream from the surrogate-target occurrence."""
    start = 0
    while True:
        j = construct.find(guide, start)
        if j < 0:
            return "", ""
        if j >= 6 and j + 20 + 3 + 6 <= len(construct) and construct[j + 21:j + 23] == "GG":
            return construct[j - 6:j], construct[j + 20:j + 29]
        start = j + 1


def load(path: str, column: str = "total_indel_eff"):
    ctx = _load_context(path)
    guides, y, ups, downs = [], [], [], []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        gi, ei, si = header.index("gRNA"), header.index(column), header.index("Surrogate ID")
        for row in reader:
            if len(row) <= max(gi, ei, si):
                continue
            g = row[gi].strip().upper()
            try:
                eff = float(row[ei])
            except ValueError:
                continue
            if len(g) == 20 and all(c in "ACGT" for c in g):
                up, down = _flanks(g, ctx.get(row[si], ""))
                guides.append(g); y.append(eff / 100.0); ups.append(up); downs.append(down)
    return guides, np.array(y), ups, downs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="path to CRISPROn seq_efficienciey.txt")
    ap.add_argument("--target", choices=["cutting", "oof"], default="cutting",
                    help="cutting = total_indel_eff (general/default.json); "
                         "oof = out_of_frame efficiency (knockout mode / default_oof.json)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--alpha", type=float, default=10.0)
    args = ap.parse_args()

    column = TARGET_COLUMN[args.target]
    if args.out is None:
        fname = "default.json" if args.target == "cutting" else "default_oof.json"
        args.out = str(Path(__file__).resolve().parent.parent / "crispr_app" / "models" / fname)

    guides, y, ups, downs = load(args.data, column)
    X = featurize_many(guides, ups=ups, downs=downs)
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
                       "target": f"{column}/100", "goal": args.target,
                       "cv_spearman": round(rho, 3), "alpha": args.alpha, "conformal": conf})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    final.save(args.out)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
