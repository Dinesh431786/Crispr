# Benchmarks & competitive positioning

This document is deliberately honest. It states where CRISPR Precision Studio is
competitive, where it is not, and exactly how to reproduce the numbers.

## How on-target predictors are judged

The field measures **Spearman's rank correlation (ρ)** between predicted and
experimentally measured editing efficiency on held-out guides (Doench 2014/2016,
Kim 2019/DeepSpCas9, Wang/Xu). Reported reference values:

| Model | Type | Spearman ρ (typical) |
|-------|------|----------------------|
| Azimuth / Rule Set 2 | gradient-boosted, interpretable-ish | ~0.56 |
| DeepSpCas9 | deep CNN | ~0.73 |
| DeepHF | RNN + biofeatures | ~0.74–0.87 |
| CRISPRon | deep CNN | ~0.80 |
| AttCRISPR | attention | ~0.87 |
| CRISPRedict | **interpretable** linear/logistic | competitive with deep models |

Sources: DeepSpCas9 ([Sci. Adv. 2019](https://www.science.org/doi/10.1126/sciadv.aax9249)),
DeepHF ([Nat. Commun. 2019](https://www.nature.com/articles/s41467-019-12281-8)),
CRISPRedict ([NAR 2022](https://academic.oup.com/nar/article/50/W1/W191/6603668)),
benchmark review ([Bioinformatics 2023, 10.1093/bib/bbad333](https://doi.org/10.1093/bib/bbad333)).

## Honest stance

A deterministic, hand-derived feature model **will not out-correlate a deep net
trained on 10^5 guides**, and we do not claim it does. Our on-target engine
targets the *interpretable* tier (Azimuth/CRISPRedict): fast, dependency-light,
and — uniquely here — every score ships with a per-feature breakdown
(`/api/explain`). Where we aim to genuinely lead:

1. **Interpretability** — additive, inspectable contributions per guide.
2. **Off-target completeness** — both-strand scan with per-site CFD **and**
   MIT/Hsu plus a CRISPOR-style aggregate specificity score.
3. **Prime editing** — PRIDICT2.0-informed pegRNA determinants in a free,
   key-less tool.
4. **Speed & UX** — sub-second vectorised scans and an export/copy/explain UI.

## Reproducing an accuracy number

The harness (`crispr_app/benchmark.py`) computes ρ with NumPy only:

```python
from benchmark import evaluate
# records: list of {"guide": <=20nt spacer, "measured": efficiency}
records = load_doench2016()          # supply a public dataset
print(evaluate(records))             # {"n":..., "spearman":..., "pearson":...}
```

To benchmark against the real reference sets, download a public dataset
(e.g. the Doench 2016 / Azimuth training table or the DeepSpCas9 test set),
map each guide to `{"guide", "measured"}`, and run `evaluate`. We intentionally
do **not** vendor a dataset to avoid shipping unverified numbers; plug in the
authoritative source and the harness reports the honest correlation.

## Real, validated scoring without any dataset download

Some published predictors print their full parameters, so they can be
reproduced exactly with **no training data**. We ship one:

**CRISPRscan / Moreno-Mateos** ([Nat. Methods 2015](https://www.nature.com/articles/nmeth.3543)).
Its 91 logistic-regression weights, intercept (0.18393…) and 35-mer framing are
transcribed verbatim into `crisprscan.py` from the open-source CRISPOR
implementation, and the test suite asserts our output equals CRISPOR's embedded
reference vector (`TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA → 77`). That is a
peer-reviewed model, byte-for-byte reproducible and unit-validated, with zero
downloads. (CRISPRscan's reported on-target ρ is ~0.4 on mammalian sets — modest,
but *honest and validated*, and it complements the surrogate in the consensus.)

The same idea extends to other fully-published models (Doench 2014 Rule Set 1,
Housden) whose coefficients are in `crisporEffScores.py`; CRISPRater/CRISPRscan
remain the cleanest reproducible linear models.

## Closing the deep-model gap — now shipped

The path from "interpretable surrogate" to "deep-model accuracy" is implemented,
not just promised:

1. **Pluggable registry** (`models.py`): resolves `onnx → learned-linear →
   heuristic`, with graceful fallback, and reports the active backend via
   `GET /api/models` and the `model` field on `/api/design`.
2. **Reproducible training** (`train.py`, `features.py`): NumPy-only closed-form
   ridge regression over the shared feature set, with an 80/20 split reporting
   train/test Spearman. Point it at a public dataset and you get an honest,
   reproducible correlation — no vendored, unverifiable numbers.
3. **ONNX backend**: export a trained DeepSpCas9/CRISPRon to
   `models/ontarget.onnx`, `pip install onnxruntime`, and the registry serves
   deep-model predictions automatically.

Sanity check (synthetic data where efficiency depends on GC + PAM-proximal G):
`train.py` recovered the signal at train ρ≈0.94 / held-out ρ≈0.93, and the API
switched to `model="linear"` with no code changes. On real screening data the
number will be lower and honest — that is the point.
