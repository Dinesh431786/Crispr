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
from conformal import calibrate, conformal_quantile, empirical_coverage  # noqa: E402
from features import featurize_many  # noqa: E402
from train import fit_ridge  # noqa: E402


TARGET_COLUMN = {"cutting": "total_indel_eff", "oof": "out_of_frame efficiency"}


def _norm_ppf(p: np.ndarray) -> np.ndarray:
    """Inverse standard-normal CDF (Acklam's rational approximation), NumPy-only
    so the build stays dependency-light. Accurate to ~1e-9 on (0,1)."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    p = np.asarray(p, dtype=np.float64)
    x = np.zeros_like(p)
    lo, hi = 0.02425, 1 - 0.02425
    m = p < lo
    q = np.sqrt(-2 * np.log(p[m]))
    x[m] = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    m = (p >= lo) & (p <= hi)
    q = p[m] - 0.5; r = q * q
    x[m] = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    m = p > hi
    q = np.sqrt(-2 * np.log(1 - p[m]))
    x[m] = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    return x


def _gauss_rank(y: np.ndarray) -> np.ndarray:
    """Gaussian-rank (van der Waerden) transform: map values to normal scores by
    rank. Training ridge on this optimises the *ranking* (Spearman) we're judged
    on rather than squared error on raw efficiency."""
    r = np.argsort(np.argsort(y))
    return _norm_ppf((r + 0.5) / len(y))


def _pav(y: np.ndarray) -> np.ndarray:
    """Pool-adjacent-violators: nearest non-decreasing fit to y (NumPy-only)."""
    sv, sw, sn = [], [], []
    for val in y.astype(np.float64):
        cv, cw, cn = float(val), 1.0, 1
        while sv and sv[-1] >= cv:
            pv, pw, pn = sv.pop(), sw.pop(), sn.pop()
            cv = (pv * pw + cv * cw) / (pw + cw); cw += pw; cn += pn
        sv.append(cv); sw.append(cw); sn.append(cn)
    out = np.empty(len(y)); pos = 0
    for v, n in zip(sv, sn):
        out[pos:pos + n] = v; pos += n
    return out


def _isotonic_knots(raw: np.ndarray, y: np.ndarray, n_knots: int = 256):
    """Monotone lookup (xs, ys) mapping a raw model score to the 0-1 efficiency
    scale, fit by isotonic regression. Monotone => preserves ranking exactly."""
    order = np.argsort(raw, kind="mergesort")
    xs = raw[order]; iso = _pav(y[order])
    uniq, inv = np.unique(xs, return_inverse=True)
    fp = np.zeros(len(uniq)); fp[inv] = iso            # last iso value per unique x
    if len(uniq) > n_knots:
        sel = np.linspace(0, len(uniq) - 1, n_knots).round().astype(int)
        uniq, fp = uniq[sel], fp[sel]
    return uniq, np.clip(fp, 0.0, 1.0)


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
    ap.add_argument("--alpha", type=float, default=150.0)
    ap.add_argument("--objective", choices=["rank", "mse"], default="rank",
                    help="rank = train on the Gaussian-rank target (optimises "
                         "Spearman) with a monotone isotonic remap back to the "
                         "0-1 efficiency scale; mse = classic squared-error ridge")
    args = ap.parse_args()

    column = TARGET_COLUMN[args.target]
    if args.out is None:
        fname = "default.json" if args.target == "cutting" else "default_oof.json"
        args.out = str(Path(__file__).resolve().parent.parent / "crispr_app" / "models" / fname)

    guides, y, ups, downs = load(args.data, column)
    X = featurize_many(guides, ups=ups, downs=downs)
    print(f"N = {len(guides)} guides  (objective={args.objective})")

    # Target seen by ridge: Gaussian-rank for the ranking objective, else raw.
    ytr_of = (lambda yy: _gauss_rank(yy)) if args.objective == "rank" else (lambda yy: yy)

    def fit_scored(Xtr, ytr):
        """Fit ridge on the (possibly rank-transformed) target and return a
        function raw->efficiency: identity for MSE, isotonic remap for rank."""
        m = fit_ridge(Xtr, ytr_of(ytr), args.alpha)
        raw_tr = Xtr @ m.weights + m.intercept
        if args.objective == "rank":
            xs, ys = _isotonic_knots(raw_tr, ytr)
            remap = lambda r: np.interp(r, xs, ys)
            calib = (xs, ys)
        else:
            remap = lambda r: np.clip(r, 0, 1)
            calib = None
        return m, remap, calib

    # Honest accuracy estimate: 5-fold cross-validated Spearman (NumPy folds).
    # Spearman is computed on the raw score (the monotone remap can't change it),
    # so this measures the true ranking quality of the shipped objective.
    rng0 = np.random.default_rng(0)
    folds = np.array_split(rng0.permutation(len(y)), 5)
    pred = np.zeros(len(y))
    for k in range(5):
        te = folds[k]
        tr = np.concatenate([folds[j] for j in range(5) if j != k])
        m, _, _ = fit_scored(X[tr], y[tr])
        pred[te] = X[te] @ m.weights + m.intercept
    rho = spearman(pred, y)
    print(f"  5-fold CV Spearman = {rho:.3f}  (wet-lab replicate band ~0.71-0.77)")

    # Split-conformal calibration on a held-out split (model never saw it), on the
    # efficiency scale (after the remap) so intervals stay in efficiency units.
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(guides))
    cut = int(0.8 * len(guides))
    tr_i, cal_i = perm[:cut], perm[cut:]
    cal_model, cal_remap, _ = fit_scored(X[tr_i], y[tr_i])
    cal_pred = np.clip(cal_remap(X[cal_i] @ cal_model.weights + cal_model.intercept), 0, 1)
    conf = calibrate(cal_pred, y[cal_i])
    for level, q in conf.items():
        print(f"  conformal {level}: half-width={q}  coverage={empirical_coverage(cal_pred, y[cal_i], q):.3f}")

    # Normalized (locally-adaptive) conformal: a lightweight sigma-model predicts
    # each guide's absolute residual, giving a per-guide interval WIDTH. Marginal
    # coverage is still guaranteed (sigma is fit on the proper-training split, so
    # the normalized scores on the calibration split stay exchangeable), but now
    # uncertainty varies per guide -- which makes uncertainty-aware ranking real.
    mu_tr = np.clip(cal_remap(X[tr_i] @ cal_model.weights + cal_model.intercept), 0, 1)
    res_tr = np.abs(y[tr_i] - mu_tr)
    sig_model = fit_ridge(X[tr_i], res_tr, max(args.alpha, 300.0))
    floor = float(np.median(res_tr)) * 0.5
    sig_cal = np.maximum(X[cal_i] @ sig_model.weights + sig_model.intercept, floor)
    qn = {lvl: round(conformal_quantile(np.abs(y[cal_i] - cal_pred) / sig_cal, a), 4)
          for a, lvl in [(0.20, "q80"), (0.10, "q90")]}
    normconf = {"weights": sig_model.weights, "intercept": sig_model.intercept,
                "floor": floor, "q": qn}
    for lvl in ("q80", "q90"):
        half = qn[lvl] * sig_cal
        cov = float(np.mean(np.abs(y[cal_i] - cal_pred) <= half))
        print(f"  norm-conformal {lvl}: q={qn[lvl]} coverage={cov:.3f} "
              f"width[min/med/max]={half.min()*2:.3f}/{np.median(half)*2:.3f}/{half.max()*2:.3f}")

    final, _, calib = fit_scored(X, y)
    if calib is not None:
        final.calibration = (np.asarray(calib[0]), np.asarray(calib[1]))
    final.normconformal = {"weights": np.asarray(normconf["weights"], dtype=np.float64),
                           "intercept": float(normconf["intercept"]),
                           "floor": float(normconf["floor"]), "q": normconf["q"]}
    final.meta.update({"trained_on": f"CRISPROn/Kim seq_efficiency (n={len(guides)})",
                       "target": f"{column}/100", "goal": args.target,
                       "objective": args.objective,
                       "cv_spearman": round(rho, 3), "alpha": args.alpha, "conformal": conf})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    final.save(args.out)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
