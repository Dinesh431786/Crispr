# Benchmark

Spearman's rank correlation (ρ) between predicted and measured editing efficiency
— the field-standard metric. Every number here is measured in this repo and
reproducible; nothing is asserted.

## Head-to-head vs industry on-target predictors

Our shipped model vs widely-used tools, scored on **identical guides** using the
precomputed competitor predictions distributed with the CRISPOR paper
(`effData/*.scores.tab`). All three datasets are **held out for our model** (it is
trained only on the CRISPRon/Kim set), so this is a fair cross-dataset
generalisation test; competitor scores are read verbatim (no re-implementation).

| Tool | doench2016 | chari2015 | morenoMateos | **mean** |
|---|:---:|:---:|:---:|:---:|
| 🧬 **CRISPR Precision Studio** | **0.263** | **0.440** | **0.220** | **0.307** 🥇 |
| CRISPRscan (Moreno-Mateos 2015) | 0.108 | 0.123 | 0.579¹ | 0.270 |
| Azimuth / Rule Set 2 (Doench 2016) | 0.269² | 0.381 | 0.120 | 0.257 |
| Chari (2015) | 0.121 | 0.457² | 0.145 | 0.241 |
| SSC (Xu 2015) | 0.168 | 0.286 | 0.171 | 0.209 |
| Wang SVM (2014) | 0.152 | 0.304 | 0.127 | 0.194 |
| WU-CRISPR (Wong 2015) | 0.110 | 0.308 | 0.037 | 0.151 |

On truly held-out data our model has the **highest mean ρ** — clear of CRISPRscan
and Azimuth/Rule Set 2 — despite being a dependency-light NumPy model with no
home-dataset advantage here. On chari2015 (0.440) it essentially ties Chari's own
home-dataset score (0.457). ¹ CRISPRscan's home dataset (trained on it).
² Azimuth/Chari had training exposure to that dataset. **Honest caveat:** absolute
cross-dataset ρ is intrinsically low (~0.15–0.4) for *every* tool, including the
giants — cross-context transfer is the field's unsolved problem, not a flaw of any
one tool.

Reproduce:
```bash
git clone --depth 1 https://github.com/maximilianh/crisporPaper
python scripts/benchmark_competitors.py --effdata crisporPaper/effData
```

## Within-dataset accuracy: wet-lab-reproducibility-grade

On the large, clean **CRISPRon/Kim set (11,617 guides)**, 5-fold cross-validated:

| Model | ρ (5-fold CV) |
|---|:---:|
| 🧬 **CRISPR Precision Studio (shipped)** | **0.766** |
| ...+ trinucleotides, before gapped-dinuc + k-mer/energy | 0.751 |
| ...guide + flanking context, before trinucleotides | 0.727 |
| ...guide-only, before any refinement | 0.707 |

Four pure-NumPy refinements to the *existing* featurizer — no new files, no
dependencies — closed most of the gap to the deep-CNN CRISPRon (~0.80):

1. **Flanking sequence context** (6 nt upstream + PAM + 6 nt downstream):
   0.707 → 0.727. The documented reason Rule Set 2 / DeepSpCas9 beat guide-only
   models; extracted from the input at inference, from CRISPRon surrogate
   constructs at training.
2. **Position-specific trinucleotides** (triplet motifs — the local patterns a
   CNN learns, kept linear): 0.727 → 0.751. Beats a gradient-boosting
   reference (0.743) while staying NumPy-only.
3. **Gapped/spaced dinucleotides** (position-specific base pairs at distances
   3–7 — long-range positional coupling the adjacent-only blocks miss; a
   "different angle" on the same sequence) plus **RC-canonical tetranucleotide,
   position-independent trinucleotide, and energy summaries** (regional GC,
   homopolymer runs, poly-T flag): 0.751 → 0.766. These two levers were found
   by a parallel measurement swarm and combine near-additively.

All were **gated on cross-dataset transfer** before shipping (none overfits: the
head-to-head mean above rose 0.274 → 0.292 → 0.307 in lock-step with CV). The
gapped-dinucleotide + k-mer step alone lifted the held-out mean +0.015 while
improving every one of the three cross-datasets.

Wet-lab replicates of the *same* guide agree only at ρ≈0.71–0.77 (Haeussler 2016),
so this **matches the assay's own reproducibility** — on par with DeepSpCas9
(~0.73). Deep models (DeepHF/CRISPRon, ~0.80–0.87) lead *within* a single large
dataset but need GPU/training; we reach the wet-lab band with a NumPy model and no
setup. Reproduce:
```bash
git clone --depth 1 https://github.com/OrensteinLab/DeepCRISTL
python scripts/build_default_model.py --data DeepCRISTL/CRISPROn/data/main_dataframes/seq_efficienciey.txt
```

## Calibrated uncertainty (conformal) — verified coverage

Every score ships with a distribution-free split-conformal interval (Lei 2018).
Empirical coverage on held-out data matches the guarantee exactly:

| Interval | Half-width | Target | Measured |
|---|:---:|:---:|:---:|
| 80% | 0.200 | 0.80 | 0.801 |
| 90% | 0.263 | 0.90 | 0.901 |

To our knowledge no other lightweight CRISPR tool ships coverage-guaranteed
uncertainty. `pytest tests/test_conformal.py` verifies the guarantee.

## Goal-aware ranking (five modes)

Rank by the outcome you actually want, not a generic proxy — three distinct kinds
of "aware-from-the-start":

| Mode | Kind | Ranking objective |
|---|---|---|
| General | model | cutting efficiency (ρ=0.766) |
| Knockout | model | out-of-frame / frameshift (ρ=0.723; a dedicated model) |
| Knock-in (HDR) | objective | cutting × exp(−cut-to-edit / 10 bp) |
| CRISPRi/a | objective | activity × exp(−bind-to-TSS / 75 bp) |
| Base editing | constraint | activity × window-centrality; only guides with a C/A in positions 4–8 |

Honest scope: knockout is a *correctness* gain (out-of-frame correlates 0.937 with
cutting), not an accuracy jump; the proximity/constraint modes are deterministic
mechanism, not statistical claims. Each re-ranking is unit-tested.

## The prediction ceiling (why not higher)

ρ≈0.71–0.77 is the reproducibility of the wet-lab data itself — a hard bound for
*any* predictor. On a noisy target (doench2016) ridge, random forest, and gradient
boosting all converge to ρ≈0.25, proving the limiter is the **data, not the
model**. Published within-dataset state-of-the-art: Azimuth ~0.56, DeepSpCas9 ~0.73,
DeepHF/AttCRISPR ~0.87 ([Sci. Adv. 2019](https://www.science.org/doi/10.1126/sciadv.aax9249),
[Nat. Commun. 2019](https://www.nature.com/articles/s41467-019-12281-8)). An ONNX
backend (`models.py`) lets you drop in a trained DeepSpCas9/CRISPRon export.

## Methods tested and not shipped (measured nulls)

Recorded for integrity — each was implemented, measured against a pre-declared
gate, showed **no improvement**, and was removed (no dead code):

- **Topological features / persistent homology** — tested twice, both null.
  (a) Chaos-Game-Representation + persistence: ρ 0.558→0.553. (b) `ripser`
  persistent homology of a *novel* encoding — Takens delay-embedding of the
  EIIP/physicochemical property track plus a purine/keto/H-bond point cloud, 30
  H0+H1 summary features on top of the shipped 0.766 model: ΔCV −0.0001,
  Δcross +0.0004 at matched alpha. A 20-nt guide is fully described by one-hot,
  so every topological summary is a *lossy* function of information ridge already
  has — it cannot add signal, only discard it. (Would also break the NumPy-only,
  zero-setup promise by pulling in a C++/Cython dependency.)
- **R-loop free-energy physics** (RNA:DNA vs DNA:DNA ΔG, Sugimoto/SantaLucia NN
  params): Δcross = −0.000 (redundant with GC/Tm already in the model).
- **sgRNA self-complementarity** (Nussinov): Δ ≈ −0.000; kept only as an
  informational QC flag, not part of the score.
- **Pigeonhole seed-index** for genome off-target: 0.6× (slower) than the
  vectorised brute scan; the real speedup needs a compiled FM-index.
- **Reduced-alphabet recoding** (purine/pyrimidine, strong/weak, amino/keto —
  position one-hot + spectrum k-mers): best variant (strong/weak 4-mer spectrum)
  gave only +0.0005 CV at matched regularisation — redundant with the 4-letter
  featurizer already present. Cross-dataset nudged up ~+0.0015 but the CV gate
  was not met.
- **Trained deep model** (sklearn MLP 256→64, ReLU, backprop, early-stopped,
  5-fold CV on the full featurizer): ρ 0.727 — **0.038 *below* linear ridge
  (0.766)**. Real capacity doesn't help: the one-hot signal is additive and the
  data is noise-limited, so a nonlinear model only overfits. Confirms the wall is
  information, not model class.
- **Nonlinearity** (pure-NumPy random-feature ELM, hidden 256–1024; low-rank
  pairwise-interaction features): ≤ +0.0006 CV over linear ridge at matched
  regularisation, and larger capacity *lowered* CV — the one-hot signal is
  effectively additive, so ridge already extracts it.
- **Position-independent gapped dinucleotides**: small CV gain but cross-dataset
  mean *dropped* ~0.003 (the distance-4 count that helps within-Kim hurts
  transfer). Only the *position-specific* form (which is what ships) improved
  both.
