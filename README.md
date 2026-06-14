# 🧬 CRISPR Precision Studio

A modern, research-grounded CRISPR guide-design platform with a FastAPI backend
and a responsive single-page web UI. Version 3.0 rebuilds the scientific core
around the determinants reported in current (2023–2026) literature, splitting
the engine into focused, testable modules while keeping the original API-first
architecture intact.

## Architecture

The 3-layer separation is unchanged — only the science layer grew:

```
Browser (templates/index.html + static/app.js)
        │  JSON over fetch()
        ▼
FastAPI (crispr_app/main.py)  ──►  Pydantic validation + utils.validate_sequence
        │
        ▼
Science layer
   ├── scoring.py    on-target efficiency  (Doench RS2 / Azimuth-informed)
   ├── offtarget.py  CFD + MIT/Hsu + aggregate specificity
   ├── prime.py      pegRNA design         (PRIDICT2.0-informed)
   └── analysis.py   pipeline + vectorised both-strand off-target search
        │  pandas DataFrame → JSON
        ▼
Browser renders ranked tables
```

No database; everything is computed per request. No external LLM/API keys are
required for any workflow.

## Train a real model & beat the interpretable baseline

The on-target score is served through a **pluggable registry** (`models.py`)
that resolves `onnx → learned-linear → heuristic` and reports which backend it
used (`GET /api/models`, and a `model` field on `/api/design`).

```bash
# Fit a reproducible linear model on real labelled data (Doench 2016, Kim 2019, …)
cd crispr_app
python train.py path/to/dataset.csv      # CSV columns: guide,measured[,ngg_context]
# -> writes models/linear.json; the API auto-loads it and reports model="linear"
```

`train.py` uses NumPy-only closed-form ridge regression (no heavy ML stack),
holds out 20 % of guides, and prints train/test Spearman so the gain is
**measured, not asserted**. For deep-model accuracy, export a trained
DeepSpCas9/CRISPRon to `models/ontarget.onnx`, install `onnxruntime`, and the
registry uses it automatically. See [BENCHMARKS.md](BENCHMARKS.md).

## Where this tool competes

Honest positioning (see [BENCHMARKS.md](BENCHMARKS.md)): deep models
(DeepHF/DeepSpCas9/CRISPRon, Spearman ρ≈0.73–0.87) lead on raw on-target
correlation and we do not claim to beat them. We target the **interpretable**
tier (Azimuth/CRISPRedict) and aim to lead on dimensions that matter day-to-day:

- **Interpretability** — every on-target score ships with a per-feature
  breakdown (`POST /api/explain`); most tools are black boxes.
- **Off-target completeness** — both-strand scan, per-site CFD *and* MIT/Hsu,
  plus a CRISPOR-style aggregate specificity score.
- **Prime editing** — PRIDICT2.0-informed pegRNA design, free and key-less.
- **Speed & UX** — sub-second vectorised scans; a professional UI with score
  heatmap badges, PAM-highlighted sequences, FASTA upload, CSV export, copy, and
  inline score explanations.

## What's new in 3.0

- **Modular scientific core.** `scoring.py`, `offtarget.py`, and `prime.py`
  isolate each model so it can be reviewed and unit-tested independently;
  `analysis.py` orchestrates them and preserves the public API.
- **Calibrated on-target model.** A transparent logistic model combining
  position-specific single/dinucleotide preferences, a quadratic GC optimum,
  nearest-neighbour Tm, homopolymer/poly-U penalties, and the NGGN 3'-context —
  the dominant features behind Doench *Rule Set 2*/Azimuth.
- **Both-strand off-target scanning.** The vectorised NumPy scanner now sweeps
  the forward **and** reverse-complement strands (the previous engine missed
  reverse-strand off-targets), reporting per-site **CFD** *and* **MIT/Hsu**
  scores plus a CRISPOR-style **aggregate specificity** score per guide.
- **PRIDICT2.0-informed Prime Editing Studio.** pegRNA ranking now optimises PBS
  melting temperature toward ~37 °C, prefers PBS ≈ 13 nt / RTT ≈ 12 nt, requires
  ≥3 nt of 3' homology past the edit, and penalises RTTs that begin with C.
- **Bug fixes.** `numpy` is now declared in `requirements.txt`; the dev container
  launches `uvicorn` (the removed Streamlit references are gone).

## Scientific basis

| Component | Model / source |
|-----------|----------------|
| On-target efficiency | Doench et al. 2014/2016 *Rule Set 2*/Azimuth (Nat. Biotechnol. 34:184); Xu et al. 2015. Modern benchmarks: Rule Set 3, DeepHF, DeepSpCas9, CRISPRon. |
| Off-target (per site) | **CFD** — Doench et al. 2016; **MIT/Hsu** — Hsu et al. 2013 (Nat. Biotechnol. 31:827). |
| Off-target (per guide) | Aggregate specificity `10000 / (100 + Σ off-target scores)` (CRISPOR convention). |
| Prime editing | **PRIDICT2.0** — Mathis et al. 2024 (Nat. Biotechnol., doi:10.1038/s41587-024-02268-2); Anzalone et al. 2019 (Nature 576:149). |

The on-target and pegRNA scores are **deterministic, research-informed
surrogates** designed for fast ranking — not re-trained deep models. They
reproduce the published *direction and relative magnitude* of each determinant
but should not be read as absolute efficiencies. The "10x" framing remains a
development benchmark for speed and feature depth; **wet-lab validation remains
essential for all CRISPR applications.**

## Run locally

```bash
cd crispr_app
pip install -r requirements.txt
uvicorn main:app --reload
```

Open: `http://127.0.0.1:8000`

## API quick reference

- `GET /health`
- `POST /api/design` — ranked gRNAs with `OnTargetScore`, `HybridScore`, `ConsensusScore`
- `POST /api/offtargets` — per-site CFD/MIT hits + per-guide `specificity` summary
- `POST /api/simulate` — protein/indel outcome of an edit
- `POST /api/prime-design` — ranked pegRNAs (Spacer + RTT + PBS) with `Score`
- `POST /api/explain` — interpretable per-feature breakdown of an on-target score
- `POST /api/upload-fasta` — parse pasted FASTA / plain DNA into a clean sequence
- `GET /api/models` — active and available on-target backends

## Tests

```bash
pip install pytest
python -m pytest tests/ -q
```

Covers on-target scoring behaviour, CFD/MIT scoring, aggregate specificity,
both-strand off-target detection, pegRNA design, performance, and dependency
hygiene.

## License

MIT
