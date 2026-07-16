<div align="center">

# 🧬 CRISPR Precision Studio

**Design, rank, explain, and validate CRISPR guides — in one lightweight platform.**

Transparent guide *prioritization*: interpretable on-target scoring, both-strand off-target analysis, and prime-editing support — without GPUs, cloud dependencies, or black-box predictions.

[![CI](https://github.com/Dinesh431786/Crispr/actions/workflows/ci.yml/badge.svg)](https://github.com/Dinesh431786/Crispr/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.11-3776AB?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-API--first-009688?logo=fastapi&logoColor=white)
![Dependencies](https://img.shields.io/badge/deps-lightweight%20(NumPy)-blue)
![Accuracy](https://img.shields.io/badge/on--target%20%CF%81-0.77%20(wet--lab%20band)-2e8b57)
![No API keys](https://img.shields.io/badge/API%20keys-none-9ee7d2)
![License](https://img.shields.io/badge/license-MIT-black)

<sub>CI runs the test suite and renders a live UI screenshot (downloadable as a build artifact) on every push.</sub>

</div>

---

## 📸 The interface

![CRISPR Precision Studio UI](docs/ui.png?v=f82c18de)

<sub>Rendered automatically by CI on every push — one **Score** per guide, with a per-feature **Details** breakdown.</sub>

---

## ✨ Highlights

| | Feature | What it means |
|---|---|---|
| 🎯 | **Goal-aware ranking** | Pick your intent — *General* (cutting) · *Knockout* (frameshift model) · *Knock-in/HDR* (cut-proximity) · *CRISPRi/a* (TSS-proximity) · *Base editing* (ABE/CBE target in window) — and guides are ranked by the outcome you actually want, each with a plain-language verdict fusing efficiency + specificity + uncertainty. |
| 🔍 | **Explainable** | Per-feature breakdown shown *by default*; plus an informational self-folding (secondary-structure) QC flag. |
| 📐 | **Calibrated uncertainty** | Distribution-free **conformal** confidence interval per guide — *verified* 90% coverage. |
| 🧬 | **Both-strand off-targets** | Vectorised NumPy scan + per-site **CFD** & **MIT/Hsu** + aggregate specificity. |
| 🌍 | **Genome-wide search** | Stream any (multi-chromosome) FASTA; memory-safe chunked scan, both strands. |
| 🧫 | **Edit-outcome simulation** | Predict the protein consequence of a cut — frameshift, stop loss, indel panel. |
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

## 🔄 Workflow

```
   DNA sequence (paste or FASTA)
            │
            ▼
   Guide discovery  ── both strands, multi-PAM (NGG/NAG/NG/TTTV)
            │
            ▼
   On-target scoring ── built-in model + CRISPRscan
            │
            ▼
   Off-target analysis ── CFD + MIT/Hsu + aggregate specificity
            │
            ▼
   Ranking ── one 0–100 Score per guide
            │
            ▼
   Explanation ── per-feature breakdown (/api/explain)
```

---

## ⚡ Example

A complete, copy-pasteable request (real 288 bp input):

```bash
curl -s -X POST http://127.0.0.1:8000/api/design \
  -H 'Content-Type: application/json' \
  -d '{"pam": "NGG", "dna_sequence": "ATGGCCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCCGCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAGCTGCAAGAACTCTTCCTCACGCGCGTCGGGCTCGACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCCGCGGTGGCGGTCTGGACCACGCCGGAGAGCGTCGAAGCGGGGGCGGTGTTCGCCGAGATCGGCCCGCGCATGGCC"}'
```

Real output → 21 guides found. Top 3 (the API returns `ConsensusScore` in 0–1; shown here ×100 as the UI does):

| # | Guide (5′→3′) | PAM | Strand | GC% | Score | `ConsensusScore` |
|---|---|:---:|:---:|:---:|:---:|:---:|
| 1 | `AAGGTGTGGGTCGCGGACGA` | CGG | + | 65 | **72** | 0.724 |
| 2 | `GATGTGGCGGTCCGGATCGA` | CGG | − | 65 | **72** | 0.718 |
| 3 | `GCTCGACATCGGCAAGGTGT` | GGG | + | 60 | **68** | 0.684 |

`POST /api/explain` then returns the per-feature breakdown (GC, T<sub>m</sub>, position-specific contributions) **and the calibrated confidence interval** behind any guide's score.

---

## 🎯 The score (production)

Each guide gets **one Score, 0–100** (higher = better) — a *relative prioritization* score combining the on-target predictors, **not** a literal % editing rate. Color-coded:

| 🟢 High | 🟡 Moderate | 🔴 Low |
|:---:|:---:|:---:|
| ≥ 60 | 40 – 59 | < 40 |

Click **Details** on any guide to see *why* it scored that way (GC, Tm, position-specific features…). Component sub-scores stay in the API/CSV for power users — never on screen.

---

## 📐 Calibrated uncertainty (what makes this different)

Every common CRISPR tool returns a bare point score. But our benchmarks show
on-target prediction is intrinsically noisy (ρ ≈ 0.25–0.8 depending on context),
so a lone number is *overconfident*. We attach a **distribution-free conformal
prediction interval** (split conformal; Lei et al. 2018) to each guide — a
confidence interval with a **mathematically guaranteed coverage level**, no
distributional assumptions, computed in pure NumPy.

`POST /api/explain` returns, e.g.:

```json
{ "score": 0.69, "interval": { "point": 0.637, "low": 0.207, "high": 1.0,
                               "level": "q90", "coverage": 0.9 } }
```

**Verified, not asserted** — empirical coverage on held-out real data matches the
target almost exactly:

| Interval | Half-width | Target coverage | Measured coverage |
|---|:---:|:---:|:---:|
| 80% | 0.200 | 0.80 | **0.801** |
| 90% | 0.263 | 0.90 | **0.901** |

Wide intervals are a feature, not a bug: they make the model's genuine
uncertainty explicit so a researcher doesn't over-trust a single number. To our
knowledge no other lightweight CRISPR guide tool ships calibrated, coverage-
guaranteed uncertainty.

---

## 📊 Accuracy — measured vs the field (not asserted)

**Head-to-head on identical guides** (Spearman ρ; competitor scores read verbatim
from CRISPOR; all datasets held out for our Kim-trained model):

| Tool | doench2016 | chari2015 | morenoMateos | **mean** |
|---|:---:|:---:|:---:|:---:|
| **OURS (NumPy ridge)** | 0.263 | 0.440 | 0.220 | **0.307** |
| CRISPRscan | 0.108 | 0.123 | 0.579 | 0.270 |
| Azimuth / Rule Set 2 | 0.269 | 0.381 | 0.120 | 0.257 |
| Wang SVM · Chari · SSC · WU-CRISPR | — | — | — | 0.15–0.24 |

On truly held-out data our lightweight model has the **highest mean ρ**, ahead of
CRISPRscan and Azimuth. Full table + reproduce command in **[BENCHMARKS.md](BENCHMARKS.md)**.

**Within a single clean dataset** (CRISPRon/Kim, 11,617 guides, 5-fold CV):
ridge **ρ=0.766** — approaching the deep-CNN CRISPRon (~0.80), past DeepSpCas9
(~0.73). Four pure-NumPy refinements to the *existing* featurizer got there:
**flanking context** (6 nt up + PAM + 6 nt down: 0.707→0.727), **position-
specific trinucleotides** (the triplet motifs a CNN learns, kept linear:
0.727→0.751), and — found by a parallel measurement swarm — **gapped/spaced
dinucleotides** (long-range positional coupling at distances 3–7, a "different
angle" on the same sequence) with **RC-canonical k-mer + energy summaries**
(0.751→0.766). Every step gated on cross-dataset transfer; the last lifted the
held-out mean +0.015 across all three cross-datasets.

> 🎯 **Wet-lab-grade.** Wet-lab replicates of the *same* guide agree only at
> **ρ≈0.71–0.77** — a hard ceiling for *any* predictor. At **ρ=0.766** the model
> agrees with the measurement **about as well as the assay agrees with itself.**
> We don't claim to beat the wet lab; we *match its reproducibility* — NumPy-only,
> zero setup. Cross-dataset ρ is low (~0.15–0.4) for *every* tool, giants included.

Every score is also **explainable and ships with a calibrated confidence
interval**; wet-lab validation remains essential.

### Train a stronger model (NumPy-only, no heavy ML stack)

```bash
cd crispr_app
python train.py dataset.csv          # columns: guide,measured[,ngg_context]
# → writes models/linear.json; the API auto-loads it and reports model="linear"

python benchmark.py data.context.tab # measure Spearman on a CRISPOR-format set
```

Position-specific dinucleotide features roughly **double** Spearman on datasets with signal (chari2015 0.20→0.40, morenoMateos 0.17→0.43); gradient boosting matched ridge to ±0.02, so we stay dependency-free.

---

## 🆚 How it compares

How we compare on the capabilities that actually differ. Our focus is **transparent, explainable prioritization** in a lightweight, API-first package.

| Capability | CRISPR Precision Studio | CRISPOR | CHOPCHOP | Benchling |
|---|:---:|:---:|:---:|:---:|
| Single explainable prioritization score | ✓ | partial¹ | partial¹ | ✗ |
| Per-feature score breakdown (API) | ✓ | ✗ | ✗ | ✗ |
| Both-strand off-target (CFD + MIT) | ✓ | ✓ (reference) | ✓ | ✓ |
| **Genome-wide** off-target search | ✓⁴ | ✓ | ✓ | ✓ |
| Calibrated uncertainty (conformal CI) | ✓ | ✗ | ✗ | ✗ |
| Prime-editing pegRNA design | ✓ | ✗² | partial | ✗ |
| JSON API-first | ✓ | partial | ✗ | ✓ |
| Runs locally, no GPU / no keys | ✓ | ✓³ | ✓³ | ✗ (SaaS) |

<sub>¹ Report several separate scores rather than one explained number. ² CRISPOR targets Cas9/Cas12a guide design; pegRNA design is usually a separate tool (PrimeDesign / pegFinder). ³ Open-source but heavier to self-host. ⁴ Streaming + chunked FASTA scan (any genome); fast on small/medium genomes, minutes at mammalian scale (see BENCHMARKS.md). Marks reflect typical usage and may change as these tools evolve.</sub>

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
| `POST /api/design` | ranked gRNAs; each carries `ConsensusScore` in **0–1** (the UI displays it ×100 as the 0–100 Score) |
| `POST /api/offtargets` | per-site CFD/MIT hits + per-guide specificity (pasted background) |
| `POST /api/offtargets-genome` | genome-wide scan over a (multi-record) FASTA, both strands |
| `POST /api/simulate` | protein / indel outcome of an edit |
| `POST /api/prime-design` | ranked pegRNAs (Spacer + RTT + PBS) |
| `POST /api/explain` | per-feature breakdown + conformal confidence interval |
| `POST /api/upload-fasta` | parse pasted FASTA / plain DNA |
| `GET /api/models` | active & available on-target backends |

---

## 🌍 Genome-wide off-target search

Scan an entire genome FASTA (any number of chromosomes/contigs), both strands,
with the same CFD + MIT/Hsu + aggregate-specificity scoring. Records are streamed
and scanned in overlapping chunks, so peak memory is bounded by the chunk size,
not the genome size.

```bash
# CLI
python crispr_app/genome.py genome.fasta GACGATCAGTCAGGATCACC --max-mismatches 3

# API
curl -s -X POST http://127.0.0.1:8000/api/offtargets-genome \
  -H 'Content-Type: application/json' \
  -d '{"guides": ["GACGATCAGTCAGGATCACC"], "fasta": ">chr1\nACGT...(your FASTA here)", "max_mismatches": 3}'
```

Verified on a real 230 kb genome (both strands) in **0.02 s** (~11 Mb/s/guide);
unit tests confirm planted off-targets are found at the correct coordinates on
both strands, across chunk boundaries, and that the perfect on-target match is
excluded. Mammalian whole-genome scans run in minutes — a seed/FM-index backend
is the next optimisation.

---

## 🌟 Prime editing — how pegRNAs are chosen

For a target base substitution, `prime.py` enumerates and ranks candidate pegRNAs using determinants from PRIDICT2.0 (Mathis 2024) and Anzalone 2019:

1. **Spacer / nick.** Scan NGG PAMs within ~30 nt of the target; place the Cas9 nick 3 bp 5′ of each PAM. Require the edit to fall 0–15 nt downstream of the nick.
2. **PBS (primer-binding site).** Enumerate lengths 8–17 nt; the PBS is the reverse complement of the sequence immediately 5′ of the nick. Its nearest-neighbour **T<sub>m</sub> is optimised toward ~37 °C** (Gaussian reward), with mild length penalties favouring ~13 nt.
3. **RTT (reverse-transcriptase template).** Enumerate lengths 10–20 nt; the RTT encodes the edit and must retain **≥3 nt of 3′ homology past the edit** for flap resolution. Penalties: RTT that **begins with C** (destabilises the edited flap) and RTT GC far from ~55%. Length term favours ~12 nt.
4. **Ranking.** A calibrated logistic score blends the PBS T<sub>m</sub>, PBS/RTT length terms, 3′-homology constraint, RTT-starts-with-C penalty, and GC term into one 0–1 `Score`.

> The pegRNA score is PRIDICT2.0-*informed*, not the trained PRIDICT2.0 network. It reproduces the published determinants for ranking; it has not yet been numerically benchmarked against a PRIDICT test set (on the roadmap). No secondary-structure (e.g. RNAfold) penalty is applied yet.

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
python -m pytest tests/ -q     # 45 passing
```

Covers on-target scoring, CFD/MIT scoring, aggregate specificity, both-strand
off-target detection, pegRNA design, the model registry & trainer, CRISPRscan
reference-vector validation, performance, and dependency hygiene.

---

<div align="center">

**MIT licensed** · No API keys required · Wet-lab validation always essential

</div>
