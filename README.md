<div align="center">

# 🧬 CRISPR Precision Studio

**Research-grounded CRISPR guide design — interpretable, fast, and honest.**

On-target scoring · both-strand off-target specificity · prime-editing pegRNAs — one clean score per guide, every score explainable.

[![CI](https://github.com/Dinesh431786/Crispr/actions/workflows/ci.yml/badge.svg)](https://github.com/Dinesh431786/Crispr/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.11-3776AB?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-API--first-009688?logo=fastapi&logoColor=white)
![Dependencies](https://img.shields.io/badge/deps-lightweight%20(NumPy)-blue)
![No API keys](https://img.shields.io/badge/API%20keys-none-9ee7d2)
![License](https://img.shields.io/badge/license-MIT-black)

<sub>CI runs the test suite and renders a live UI screenshot (downloadable as a build artifact) on every push.</sub>

</div>

---

## ✨ Highlights

| | Feature | What it means |
|---|---|---|
| 🎯 | **One Efficiency score** | A single 0–100 number ranks each guide — no column soup. |
| 🔍 | **Explainable** | `POST /api/explain` shows the per-feature breakdown behind every score. |
| 🧬 | **Both-strand off-targets** | Vectorised NumPy scan + per-site **CFD** & **MIT/Hsu** + aggregate specificity. |
| 🌟 | **Prime Editing Studio** | PRIDICT2.0-informed pegRNA design (Spacer + RTT + PBS). |
| 📚 | **Peer-reviewed scoring** | **CRISPRscan** weights reproduced verbatim & unit-validated — zero downloads. |
| 🔌 | **Pluggable models** | `onnx → trained-linear → heuristic`, auto-selected and reported. |
| ⚡ | **Lightweight** | No GPU, no LLM keys, no DB — everything computed per request. |

---

## 🚀 Quickstart

```bash
cd crispr_app
pip install -r requirements.txt
uvicorn main:app --reload
```

➡️  Open **http://127.0.0.1:8000**

---

## 🎯 The score (production)

Each guide gets **one Efficiency score, 0–100** (higher = better), color-coded:

| 🟢 High | 🟡 Moderate | 🔴 Low |
|:---:|:---:|:---:|
| ≥ 60 | 40 – 59 | < 40 |

Click **Details** on any guide to see *why* it scored that way (GC, Tm, position-specific features…). Component sub-scores stay in the API/CSV for power users — never on screen.

---

## 📊 Accuracy — measured, not asserted

Held-out Spearman ρ on real public datasets (full table + method in **[BENCHMARKS.md](BENCHMARKS.md)**):

| Configuration | ρ | Notes |
|---|:---:|---|
| **Ships trained** (pooled human SpCas9) | **0.22 – 0.41** | default model, zero setup (leave-one-dataset-out) |
| Trained on your own data | 0.40 – 0.52 | one command — `train.py` |
| Heuristic fallback | ~0.25 | interpretable, used if no model present |
| CRISPRscan (validated) | 0.58 | on its home dataset |
| ONNX backend (deep model) | ~0.85 | bring a DeepSpCas9/CRISPRon export |

> ⚠️ **Honesty note.** No predictor can exceed the ~0.71–0.77 reproducibility ceiling of the wet-lab data itself; published state-of-the-art tops out around ~0.85–0.88. Our scores are deterministic, research-informed surrogates for *ranking*; **wet-lab validation remains essential.**

### Train a stronger model (NumPy-only, no heavy ML stack)

```bash
cd crispr_app
python train.py dataset.csv          # columns: guide,measured[,ngg_context]
# → writes models/linear.json; the API auto-loads it and reports model="linear"

python benchmark.py data.context.tab # measure Spearman on a CRISPOR-format set
```

Position-specific dinucleotide features roughly **double** Spearman on datasets with signal (chari2015 0.20→0.40, morenoMateos 0.17→0.43); gradient boosting matched ridge to ±0.02, so we stay dependency-free.

---

## 🏗️ Architecture

```
Browser (templates/index.html + static/app.js)
        │  JSON over fetch()
        ▼
FastAPI (main.py)  ──►  Pydantic validation + utils.validate_sequence
        │
        ▼
Science layer
   ├── scoring.py     on-target efficiency   (Doench RS2 / Azimuth-informed)
   ├── crisprscan.py  CRISPRscan             (Moreno-Mateos 2015, verbatim)
   ├── offtarget.py   CFD + MIT/Hsu + aggregate specificity
   ├── prime.py       pegRNA design          (PRIDICT2.0-informed)
   ├── features.py / models.py / train.py    pluggable + trainable models
   └── analysis.py    pipeline + vectorised both-strand off-target search
        │  pandas DataFrame → JSON
        ▼
Browser renders one ranked table
```

---

## 🔌 API reference

| Method & route | Purpose |
|---|---|
| `GET /health` | liveness check |
| `POST /api/design` | ranked gRNAs with the `ConsensusScore` (Efficiency) |
| `POST /api/offtargets` | per-site CFD/MIT hits + per-guide specificity summary |
| `POST /api/simulate` | protein / indel outcome of an edit |
| `POST /api/prime-design` | ranked pegRNAs (Spacer + RTT + PBS) |
| `POST /api/explain` | interpretable per-feature score breakdown |
| `POST /api/upload-fasta` | parse pasted FASTA / plain DNA |
| `GET /api/models` | active & available on-target backends |

---

## 🔬 Scientific basis

| Component | Model / source |
|---|---|
| On-target | Doench 2014/2016 *Rule Set 2*/Azimuth (Nat. Biotechnol. 34:184); CRISPRscan (Moreno-Mateos, Nat. Methods 2015) |
| Off-target (site) | **CFD** (Doench 2016) · **MIT/Hsu** (Hsu 2013, Nat. Biotechnol. 31:827) |
| Off-target (guide) | aggregate specificity `10000 / (100 + Σ scores)` (CRISPOR convention) |
| Prime editing | **PRIDICT2.0** (Mathis 2024, doi:10.1038/s41587-024-02268-2); Anzalone 2019 (Nature 576:149) |

---

## ✅ Tests

```bash
pip install pytest
python -m pytest tests/ -q     # 35 passing
```

Covers on-target scoring, CFD/MIT scoring, aggregate specificity, both-strand
off-target detection, pegRNA design, the model registry & trainer, CRISPRscan
reference-vector validation, performance, and dependency hygiene.

---

<div align="center">

**MIT licensed** · No API keys required · Wet-lab validation always essential

</div>
