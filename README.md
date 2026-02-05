# CRISPR Precision Studio

CRISPR Precision Studio is an API-first guide RNA design and analysis platform engineered for high-confidence candidate prioritization, rapid iteration, and reproducible scientific workflows.

## Executive Summary

This project provides a modern CRISPR design stack with:

- high-throughput gRNA discovery across multiple PAM families,
- consensus scoring for on-target activity prioritization,
- off-target risk profiling using mismatch and CFD-style scoring,
- edit outcome simulation support for downstream experimental planning,
- a responsive web interface backed by a production-ready FastAPI service.

The goal is simple: deliver a rigorous, extensible platform that can be benchmarked transparently and improved continuously against state-of-the-art baselines.

---

## Core Capabilities

- **Multi-PAM guide discovery**
  - Cas9: `NGG`, `NAG`, `NG`
  - Cas12a: `TTTV`
- **Consensus ranking pipeline**
  - Hybrid sequence-quality score
  - ML-inspired feature score
  - Combined consensus score for candidate ordering
- **Off-target assessment**
  - configurable mismatch thresholds
  - CFD-style risk weighting
- **Edit impact simulation**
  - insertion/deletion scenarios
  - translated protein effect previews
- **API + web app delivery**
  - FastAPI backend
  - lightweight browser UI for bench scientists and translational teams

---


## Production Readiness Upgrades

- strict request validation with typed payload models and PAM allow-lists,
- API-level request tracing headers (`X-Request-Id`, `X-Elapsed-Ms`),
- structured error payloads for stable client integration,
- indexed off-target candidate generation for better scalability at larger background sizes,
- bidirectional (forward + reverse) off-target scanning support for improved sensitivity.

## System Architecture

```text
Frontend (HTML/CSS/JS)
        |
        v
FastAPI service (main.py)
        |
        +--> guide discovery + scoring engine (analysis.py)
        +--> sequence utilities/validation (utils.py)
        +--> simulation and off-target workflows
```

This separation keeps the stack clean, testable, and deployment-ready for local, cloud, or internal lab infrastructure.

---

## Benchmarking Strategy (for top-tier performance)

To evaluate and improve against leading tools and commercial-grade standards, use a reproducible benchmark harness with the following protocol:

1. **Dataset design**
   - curated validated guides with measured editing outcomes,
   - stratification by species, locus type, GC bands, and PAM class.
2. **Primary metrics**
   - ranking metrics: Spearman/Pearson, NDCG@k,
   - classification metrics: AUROC, AUPRC for high-activity guide selection,
   - off-target screening metrics: recall/precision at mismatch thresholds.
3. **Latency metrics**
   - p50/p95 response times for `/api/design` and `/api/offtargets`,
   - throughput at concurrent load.
4. **Reproducibility controls**
   - fixed software versions,
   - immutable test sets,
   - versioned benchmark reports committed with each model/scoring revision.

### Performance Improvement Roadmap

- integrate experimentally grounded feature expansions (position-specific signals, thermodynamic terms),
- add optional pluggable model backends with strict benchmark gating,
- profile and optimize hot loops for larger genomic backgrounds,
- introduce CI benchmark checks to prevent regression in ranking quality and latency.

---

## Getting Started

### 1) Install

```bash
cd crispr_app
pip install -r requirements.txt
```

### 2) Run

```bash
uvicorn main:app --reload
```

Open: `http://127.0.0.1:8000`

---

## API Reference

- `GET /health` — service health
- `POST /api/design` — gRNA design and ranking
- `POST /api/offtargets` — off-target search and scoring
- `POST /api/simulate` — edit simulation output

---

## Repository Layout

```text
crispr_app/
  main.py                # FastAPI application
  analysis.py            # scoring, guide discovery, off-target and simulation logic
  utils.py               # sequence validation and parsing helpers
  templates/index.html   # web interface markup
  static/style.css       # UI styling
  static/app.js          # UI interactions with API
tests/
  test_analysis.py
  test_dependencies.py
```

---

## Scientific and Product Positioning

CRISPR Precision Studio is built for teams that require **measurable quality**, **transparent methodology**, and **rapid operational iteration**.

Comparative superiority claims should always be made from formal benchmark evidence on shared datasets and clearly reported metrics.

---

## License

MIT
