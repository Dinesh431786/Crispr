# 🧬 CRISPR Precision Studio

A modern CRISPR guide design platform with a FastAPI backend and a responsive web UI (no Streamlit).

## What changed

- ✅ Streamlit UI removed.
- ✅ New API-first architecture for speed and cleaner separation of concerns.
- ✅ Updated scoring pipeline with research-informed sequence features.
- ✅ Better UX: single-page interface for design + off-target analysis.

## Key capabilities

- Multi-PAM gRNA search (`NGG`, `NAG`, `NG`, `TTTV`)
- Consensus ranking from:
  - **Hybrid score** (rule-based sequence quality)
  - **ML-inspired score** (feature-engineered surrogate)
- Off-target scanning with mismatch limits + CFD scoring
- Protein edit/indel simulation API

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

This repository now has stronger engineering and better model heuristics, but it does **not** claim guaranteed "10x faster" or "10x more accurate" outcomes versus any external proprietary model. Experimental validation remains essential.

## License

MIT
