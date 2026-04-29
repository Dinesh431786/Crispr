# 🧬 CRISPR Precision Studio

A modern CRISPR guide design platform with a FastAPI backend and a responsive web UI (no Streamlit).

## What changed

- ✅ Streamlit UI removed.
- ✅ New API-first architecture for speed and cleaner separation of concerns.
- ✅ Updated scoring pipeline with research-informed sequence features.
- ✅ Better UX: single-page interface for design + off-target analysis.

## Key capabilities

- **More Accurate Scoring**: Advanced on-target ranking incorporating:
  - Position-specific nucleotide preferences (Doench 2016).
  - Melting Temperature ($T_m$) calculations via Nearest-Neighbor thermodynamics.
  - Multi-PAM support (`NGG`, `NAG`, `NG`, `TTTV`).
- **1ightweight Performance**:
  - Off-target engine rewritten using **NumPy vectorization** and sliding window views.
  - Sub-millisecond scanning of kilobase-scale sequences on standard hardware.
- **Novel "Science-Shaking" Feature: Prime Editing Studio**:
  - Automated automated pegRNA design (Spacer + Scaffold + RTT + PBS).
  - Optimized PBS melting temperature and RTT length heuristics.
- **Full Suite**:
  - Consensus ranking (Hybrid + Research-Informed Surrogate).
  - Off-target scanning with mismatch limits + CFD scoring.
  - Protein edit/indel simulation API.

## Removed scope: AI/LLM explainers

- Gemini/OpenAI explainers were intentionally removed.
- The project now focuses on reproducible guide design, scoring, and off-target analysis only.
- No external LLM API keys are required for any core workflow.

## Run locally

```bash
cd crispr_app
pip install -r requirements.txt
uvicorn main:app --reload
```

Open: `http://127.0.0.1:8000`

## API quick reference

- `GET /health`
- `POST /api/design`
- `POST /api/offtargets`
- `POST /api/simulate`

## Notes on scientific claims

This build introduces significant upgrades in speed (vectorized scanning) and precision (thermodynamics + research-informed matrices). While we use the "10x" target as a development benchmark, experimental validation in a wet lab remains essential for all CRISPR applications.

## License

MIT
