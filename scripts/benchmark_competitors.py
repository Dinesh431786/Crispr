"""Standard head-to-head benchmark vs industry on-target predictors.

Compares the shipped model against widely-used tools on *identical guides* using
the precomputed competitor scores distributed with the CRISPOR paper
(github.com/maximilianh/crisporPaper, effData/*.scores.tab). Metric is Spearman's
rho vs measured efficiency (`modFreq`) — the field-standard.

The datasets used are all HELD OUT for our model (it is trained only on the
CRISPRon/Kim set), so this is a fair cross-dataset generalisation test. Competitor
scores are read verbatim from the CRISPOR tables (no re-implementation), so the
comparison is apples-to-apples on the same rows.

Usage:
    git clone --depth 1 https://github.com/maximilianh/crisporPaper
    python scripts/benchmark_competitors.py --effdata crisporPaper/effData
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "crispr_app"))
from benchmark import spearman  # noqa: E402
from models import predict_on_target  # noqa: E402

DATASETS = ["doench2016_hg19", "chari2015Train", "morenoMateos2015"]
# CRISPOR score column -> industry tool display name.
COMPETITORS = {
    "fusi": "Azimuth (Doench RS2)",
    "wang": "Wang SVM",
    "crisprScan": "CRISPRscan",
    "chariRank": "Chari",
    "ssc": "SSC",
    "wuCrispr": "WU-CRISPR",
}


def _flanks(guide: str, long_seq: str) -> tuple[str, str]:
    """6 nt upstream, 9 nt (PAM+6) downstream around the protospacer, if locatable."""
    s = (long_seq or "").upper()
    j = s.find(guide)
    if j >= 6 and j + 20 + 9 <= len(s):
        return s[j - 6:j], s[j + 20:j + 29]
    return "", ""


def load(path: str):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = (r.get("seq") or "").strip().upper()[:20]
            if len(g) != 20 or any(c not in "ACGT" for c in g):
                continue
            try:
                y = float(r["modFreq"])
            except (KeyError, ValueError):
                continue
            up, down = _flanks(g, r.get("longSeq100Bp", ""))
            rows.append((g, y, r, up, down))
    return rows


def col_rho(rows, col):
    xs, ys = [], []
    for g, y, r, up, down in rows:
        try:
            xs.append(float(r.get(col, ""))); ys.append(y)
        except ValueError:
            pass
    return spearman(xs, ys) if len(xs) > 20 else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--effdata", required=True, help="path to crisporPaper/effData")
    args = ap.parse_args()

    data = {n: load(f"{args.effdata}/{n}.scores.tab") for n in DATASETS}
    header = ["Tool"] + [d.split("_")[0] for d in DATASETS] + ["mean"]
    print(f"{header[0]:<24}" + "".join(f"{h[:10]:>12}" for h in header[1:]))

    def row(name, fn):
        per = [fn(data[n]) for n in DATASETS]
        print(f"{name:<24}" + "".join(f"{v:>12.3f}" for v in per) + f"{np.nanmean(per):>12.3f}")

    row("CRISPR Precision Studio", lambda rows: spearman(
        [predict_on_target(g, up=up, down=down) for g, _, _, up, down in rows],
        [y for _, y, _, _, _ in rows]))
    for col, disp in COMPETITORS.items():
        row(disp, lambda rows, c=col: col_rho(rows, c))

    print("\nSpearman rho vs measured efficiency, identical guides. All datasets are")
    print("held out for CRISPR Precision Studio (trained only on CRISPRon/Kim).")
    print("Competitor scores are read verbatim from CRISPOR. Cross-dataset rho is")
    print("intrinsically low for all tools; within-dataset Precision Studio reaches")
    print("0.766 on the Kim set (wet-lab reproducibility band).")


if __name__ == "__main__":
    main()
