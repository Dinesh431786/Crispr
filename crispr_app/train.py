"""Train a reproducible linear on-target model from real labelled data.

No heavy ML stack: closed-form ridge regression (NumPy only). Given a CSV with
columns ``guide`` and ``measured`` (efficiency in [0,1]; an optional
``ngg_context`` column is used if present), this fits a model over the shared
feature set (:mod:`features`) and writes ``models/linear.json``, which the
registry (:mod:`models`) then loads automatically.

Usage:
    python train.py path/to/dataset.csv [--out models/linear.json] [--alpha 1.0]

Recommended datasets: Doench 2016 (Azimuth), Kim 2019 (DeepSpCas9), Wang/Xu.
The script reports train and held-out Spearman so the gain is measured, honestly.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

try:
    from benchmark import spearman
    from features import featurize_many, n_features
    from models import LinearModel
except ImportError:  # pragma: no cover
    from .benchmark import spearman
    from .features import featurize_many, n_features
    from .models import LinearModel


def load_csv(path: str) -> tuple[list[str], list[str | None], np.ndarray]:
    guides, contexts, y = [], [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            g = (row.get("guide") or row.get("gRNA") or "").strip().upper()
            m = row.get("measured") or row.get("efficiency") or row.get("score")
            if len(g) >= 18 and m not in (None, ""):
                guides.append(g)
                contexts.append((row.get("ngg_context") or "").strip() or None)
                y.append(float(m))
    return guides, contexts, np.asarray(y, dtype=np.float64)


def fit_ridge(X: np.ndarray, y: np.ndarray, alpha: float = 1.0) -> LinearModel:
    """Closed-form ridge: w = (XtX + alpha I)^-1 Xt y, with separate intercept."""
    Xc = X - X.mean(axis=0, keepdims=True)
    yc = y - y.mean()
    d = Xc.shape[1]
    w = np.linalg.solve(Xc.T @ Xc + alpha * np.eye(d), Xc.T @ yc)
    intercept = float(y.mean() - X.mean(axis=0) @ w)
    return LinearModel(w, intercept, meta={"alpha": alpha, "n_features": d})


def train(path: str, out: str, alpha: float = 1.0, seed: int = 0) -> dict:
    guides, contexts, y = load_csv(path)
    if len(guides) < 10:
        raise SystemExit(f"Need >=10 labelled guides, got {len(guides)}.")

    X = featurize_many(guides, contexts)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(guides))
    cut = max(1, int(0.8 * len(guides)))
    tr, te = idx[:cut], idx[cut:]

    model = fit_ridge(X[tr], y[tr], alpha)
    pred_tr = np.clip(X[tr] @ model.weights + model.intercept, 0, 1)
    report = {
        "n_total": len(guides),
        "n_train": len(tr),
        "n_test": len(te),
        "spearman_train": round(spearman(pred_tr, y[tr]), 4),
    }
    if len(te) >= 2:
        pred_te = np.clip(X[te] @ model.weights + model.intercept, 0, 1)
        report["spearman_test"] = round(spearman(pred_te, y[te]), 4)

    # Refit on all data before saving (use the full dataset for deployment).
    full = fit_ridge(X, y, alpha)
    full.meta.update({"trained_on": Path(path).name, "report": report})
    full.save(out)
    report["saved"] = out
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description="Train linear on-target model")
    ap.add_argument("csv")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "models" / "linear.json"))
    ap.add_argument("--alpha", type=float, default=1.0)
    args = ap.parse_args()
    report = train(args.csv, args.out, args.alpha)
    print("Training complete:")
    for k, v in report.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
