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

![CRISPR Precision Studio UI](docs/ui.png)

<sub>Rendered automatically by CI on every push — one **Score** per guide, with a per-feature **Details** breakdown.</sub>

---

## 🧭 What it does (in plain English)

Paste a DNA sequence. **CRISPR Precision Studio** finds every possible guide RNA,
**scores and ranks them** so you know which to try first — and, unlike most tools,
**shows you *why*** each guide scored the way it did and **how confident** it is.
It also checks each guide for off-target risk, and designs base-editing and
prime-editing guides.

**Three things make it different:**

1. **One clear score, explained.** Every guide gets a single **0–100** ranking
   score with a plain-language verdict and a per-feature breakdown — no black box.
2. **It tells you when it's unsure.** Each score carries a calibrated confidence
   interval with *mathematically guaranteed* coverage — not a bare number.
3. **Accurate *and* honest.** ρ≈0.77 on-target (wet-lab-grade) and **#1
   cross-dataset** vs six industry tools — every figure measured and reproducible,
   never asserted.

**Who it's for:** molecular biologists choosing guides for a knockout, knock-in,
CRISPRi/a, base-edit or prime-edit experiment — running locally, with no GPU, no
account and no API keys.

---

## ✨ Highlights

| | Feature | What it means |
|---|---|---|
| 🎯 | **Goal-aware ranking** | Pick your intent — *General* (cutting) · *Knockout* (frameshift model) · *Knock-in/HDR* (cut-proximity) · *CRISPRi/a* (TSS-proximity) · *Base editing* (ABE/CBE target in window) — and guides are ranked by the outcome you actually want, each with a plain-language verdict fusing efficiency + specificity + uncertainty. |
| 🔍 | **Explainable** | Per-feature breakdown shown *by default*; plus an informational self-folding (secondary-structure) QC flag. |
| 🔬 | **What-if sensitivity** | One click maps every single-base change in the spacer to its Δscore — see which positions the model cares about and the best/worst mutation (saturation mutagenesis, instant). |
| 📐 | **Calibrated uncertainty** | Distribution-free **conformal** confidence interval per guide — *verified* 90% coverage, with **per-guide adaptive width** (normalized conformal). |
| ⚖️ | **Uncertainty-aware ranking** | Rank by *balanced*, *conservative* (pessimistic CI bound), *robust* (uncertainty-penalised), or *optimistic* — turn the interval into a decision, not just a display. |
| 🧬 | **Both-strand off-targets** | Vectorised NumPy scan + per-site **CFD** & **MIT/Hsu** + aggregate specificity. |
| 🌍 | **Genome-wide search** | Stream any (multi-chromosome) FASTA; memory-safe chunked scan, both strands. |
| 🧫 | **Edit-outcome simulation** | Protein consequence of a cut + a **sequence-derived repair spectrum** (MMEJ microhomology deletions, templated +1 insertion) with an **out-of-frame / knockout probability**. |
| 🧪 | **Base Editing Studio** | Per-guide window efficiency, **bystander** edits, **editing purity**, and a composite BE score — not just "editable yes/no". |
| 🧬 | **Multiplex library designer** | Greedy marginal-gain selection of a strong **and diverse** guide set — avoids the near-duplicates that recombine in pooled libraries. |
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
| 1 | `GCTCGACATCGGCAAGGTGT` | GGG | + | 60 | **71** | 0.712 |
| 2 | `AAGGTGTGGGTCGCGGACGA` | CGG | + | 65 | **70** | 0.697 |
| 3 | `GATGTGGCGGTCCGGATCGA` | CGG | − | 65 | **68** | 0.684 |

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
| 80% | 0.198 | 0.80 | **0.800** |
| 90% | 0.260 | 0.90 | **0.901** |

Wide intervals are a feature, not a bug: they make the model's genuine
uncertainty explicit so a researcher doesn't over-trust a single number. To our
knowledge no other lightweight CRISPR guide tool ships calibrated, coverage-
guaranteed uncertainty.

**Per-guide adaptive width (normalized conformal).** A lightweight sigma-model
gives every guide its *own* interval width — a hard-to-call guide gets a wider
band than an easy one — while marginal coverage stays guaranteed (0.800 / 0.901
verified). This makes the uncertainty **actionable as a ranking dimension**:

| Ranking | Sorts by | Use when |
|---|---|---|
| **Balanced** | the point score | default |
| **Conservative** | score − CI half-width (pessimistic bound) | "even the worst case must be strong" |
| **Robust** | score − ½·half-width | penalise, don't fully discount, uncertainty |
| **Optimistic** | score + CI half-width (best case) | exploratory / willing to gamble |

Because the width varies per guide, a high-but-uncertain guide is genuinely
demoted under *conservative* in favour of a slightly-lower-but-tighter one — the
ordering really changes (it is not a re-label of the same list). Pick it from the
UI **Ranking** dropdown or pass `ranking_strategy` to `POST /api/design`.

---

## 📊 Accuracy — measured vs the field (not asserted)

**Head-to-head on identical guides** (Spearman ρ; competitor scores read verbatim
from CRISPOR; all datasets held out for our Kim-trained model):

| Tool | doench2016 | chari2015 | morenoMateos | **mean** |
|---|:---:|:---:|:---:|:---:|
| **CRISPR Precision Studio** | 0.271 | 0.437 | 0.221 | **0.310** |
| CRISPRscan | 0.108 | 0.123 | 0.579 | 0.270 |
| Azimuth / Rule Set 2 | 0.269 | 0.381 | 0.120 | 0.257 |
| Wang SVM · Chari · SSC · WU-CRISPR | — | — | — | 0.15–0.24 |

On truly held-out data our lightweight model has the **highest mean ρ**, ahead of
CRISPRscan and Azimuth. Full table + reproduce command in **[BENCHMARKS.md](BENCHMARKS.md)**.

**Within a single clean dataset** (CRISPRon/Kim, 11,617 guides, 5-fold CV):
**CRISPR Precision Studio** reaches **ρ=0.767** — approaching the deep-CNN CRISPRon (~0.80), past DeepSpCas9
(~0.73). Pure-NumPy refinements to the *existing* featurizer got there:
**flanking context** (6 nt up + PAM + 6 nt down: 0.707→0.727), **position-
specific trinucleotides** (the triplet motifs a CNN learns, kept linear:
0.727→0.751), **gapped/spaced dinucleotides** (long-range positional coupling at
distances 3–7, a "different angle" on the same sequence) with **RC-canonical
k-mer + energy summaries** (0.751→0.766), and finally a **ranking-target
objective** — training on the Gaussian-rank of efficiency to optimise Spearman
directly (the metric we're judged on), then a monotone isotonic remap back to the
efficiency scale. That last swap needs *no new features* and lifts knockout mode
0.723→0.728 and cross-dataset doench2016 0.263→**0.271, past Azimuth on its own
home set**. Every step gated on cross-dataset transfer.

> 🎯 **Wet-lab-grade.** Wet-lab replicates of the *same* guide agree only at
> **ρ≈0.71–0.77** — a hard ceiling for *any* predictor. At **ρ=0.767** the model
> agrees with the measurement **about as well as the assay agrees with itself.**
> We don't claim to beat the wet lab; we *match its reproducibility* — NumPy-only,
> zero setup. Cross-dataset ρ is low (~0.15–0.4) for *every* tool, giants included.

> **Separating good guides from bad** — the task you act on — runs at **AUC ≈ 0.95**
> (top⅓ vs bottom⅓; 0.97 quartiles, 0.98 deciles), well above GC-only (0.61) or
> random (0.50). Honest caveat: this is *not* a breakthrough past the noise ceiling —
> a replicate that just adds assay-level noise (ρ≈0.78) scores ~0.97 on the same
> task, so ~0.95 is the assay's own separability of the extremes, the same
> "matches its reproducibility" story as the 0.767 Spearman, on a cleaner metric.
> We report it because separation is what you act on, not to claim we cracked 0.9
> ([details](BENCHMARKS.md#good-vs-bad-separation-auc--a-cleaner-metric-honestly-contextualised)).

Every score is also **explainable and ships with a calibrated confidence
interval**; wet-lab validation remains essential.

> **The real frontier (and why we stop here).** Cross-dataset ρ (~0.31) is capped
> by *information*, not modeling: editing efficiency in a living cell depends on
> **chromatin accessibility** — is the target even open? — which a 20-nt sequence
> cannot contain. We measured this directly: added features, kernels, and a trained
> MLP all plateau at the sequence ceiling ([BENCHMARKS.md](BENCHMARKS.md#methods-tested-and-not-shipped-measured-nulls)).
> The honest next lever is epigenomic context (ATAC-seq/DNase), which would break
> the lightweight, sequence-only, zero-setup promise. We deliberately stay on this
> side of that line and say so, rather than chase +0.001s that don't move real use.

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
| `POST /api/design` | ranked gRNAs; each carries `ConsensusScore` in **0–1** (the UI displays it ×100) + a per-guide `CI_low`/`CI_high`/`CI_halfwidth`. Optional `ranking_strategy`: `balanced`\|`conservative`\|`robust`\|`optimistic` |
| `POST /api/offtargets` | per-site CFD/MIT hits + per-guide specificity (pasted background) |
| `POST /api/offtargets-genome` | genome-wide scan over a (multi-record) FASTA, both strands |
| `POST /api/simulate` | protein / indel outcome of an edit |
| `POST /api/base-edit` | per-guide ABE/CBE window efficiency, bystanders, editing purity, composite BE score |
| `POST /api/design-multiplex` | greedy marginal-gain pooled/multiplex guide set (strong + diverse) |
| `POST /api/prime-design` | ranked pegRNAs (Spacer + RTT + PBS) |
| `POST /api/explain` | per-feature breakdown + conformal confidence interval |
| `POST /api/sensitivity` | what-if: Δscore of every single-base spacer mutation (saturation map) |
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

> The pegRNA score is PRIDICT2.0-*informed*, not the trained PRIDICT2.0 network. It reproduces the published determinants for ranking; it has not yet been numerically benchmarked against a PRIDICT test set (on the roadmap). A **lightweight 3′-extension secondary-structure penalty** is applied (dependency-free Nussinov max-base-pairing, shown as `ExtStructure`/**Fold**) — a structure *propensity*, not an RNAfold free energy.

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
python -m pytest tests/ -q     # all passing
```

Covers on-target scoring, CFD/MIT scoring, aggregate specificity, both-strand
off-target detection, pegRNA + structure scoring, base-edit purity, multiplex
selection, uncertainty-aware ranking, the indel spectrum, what-if sensitivity,
the model registry & trainer, CRISPRscan reference-vector validation, end-to-end
API smoke tests, performance, and dependency hygiene.

---

## 🔁 Reproduce every number

Nothing here is asserted — each headline figure has a one-command reproduction.

```bash
# On-target accuracy (rho=0.767, 5-fold CV) + conformal calibration:
git clone --depth 1 https://github.com/OrensteinLab/DeepCRISTL
python scripts/build_default_model.py --data DeepCRISTL/CRISPROn/data/main_dataframes/seq_efficienciey.txt

# Head-to-head vs 6 industry tools on identical guides (#1 cross-dataset mean 0.310):
git clone --depth 1 https://github.com/maximilianh/crisporPaper
python scripts/benchmark_competitors.py --effdata crisporPaper/effData

# Frameshift predictor validated against the trained Lindel model (rho=0.54, n=6058):
git clone --depth 1 https://github.com/shendurelab/Lindel
python scripts/validate_indel_spectrum.py --lindel Lindel --effdata crisporPaper/effData

# Conformal coverage guarantee + all module unit/behaviour tests:
python -m pytest tests/ -q
```

Full method notes and honest scope for every feature are in
**[BENCHMARKS.md](BENCHMARKS.md)** — including the **measured nulls** (things we
tried, that didn't help, and removed).

---

<div align="center">

**MIT licensed** · No API keys required · Wet-lab validation always essential

</div>
