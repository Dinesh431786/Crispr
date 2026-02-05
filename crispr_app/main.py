from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field
from starlette.requests import Request

from analysis import find_gRNAs, find_off_targets_detailed, indel_simulations, simulate_protein_edit
from utils import validate_sequence

BASE_DIR = Path(__file__).resolve().parent

app = FastAPI(title="CRISPR Precision Studio", version="2.0.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))


class DesignRequest(BaseModel):
    dna_sequence: str = Field(..., min_length=23)
    pam: str = "NGG"
    guide_length: int = Field(20, ge=18, le=25)
    min_gc: int = Field(40, ge=30, le=80)
    max_gc: int = Field(70, ge=40, le=90)
    add_5prime_g: bool = False


class OffTargetRequest(BaseModel):
    guides: list[str]
    background_sequence: str
    max_mismatches: int = Field(2, ge=0, le=4)
    pam: str = "NGG"


class SimulateRequest(BaseModel):
    dna_sequence: str
    guide: str
    edit_offset: int = 20
    edit_type: str = "del1"


@app.get("/", response_class=HTMLResponse)
def ui(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/api/design")
def design_guides(payload: DesignRequest):
    ok, cleaned = validate_sequence(payload.dna_sequence)
    if not ok:
        raise HTTPException(status_code=400, detail=cleaned)

    guides = find_gRNAs(
        cleaned,
        pam=payload.pam,
        guide_length=payload.guide_length,
        min_gc=payload.min_gc,
        max_gc=payload.max_gc,
        add_5prime_g=payload.add_5prime_g,
    )

    return {
        "count": int(len(guides)),
        "top_guides": guides.head(100).to_dict(orient="records"),
    }


@app.post("/api/offtargets")
def find_offtargets(payload: OffTargetRequest):
    import pandas as pd

    df = pd.DataFrame({"gRNA": payload.guides})
    results = find_off_targets_detailed(df, payload.background_sequence, payload.max_mismatches, pam=payload.pam)
    if results.empty:
        return {"count": 0, "off_targets": []}

    return {
        "count": int(len(results)),
        "off_targets": results.sort_values("CFD_Score", ascending=False).head(500).to_dict(orient="records"),
    }


@app.post("/api/simulate")
def simulate(payload: SimulateRequest):
    idx = payload.dna_sequence.upper().find(payload.guide.upper())
    if idx < 0:
        raise HTTPException(status_code=400, detail="Guide not found in DNA sequence")

    cut = idx + payload.edit_offset
    before, after, frameshift, stop_lost = simulate_protein_edit(payload.dna_sequence, cut, payload.edit_type)
    indels = indel_simulations(payload.dna_sequence, cut)

    return {
        "protein_before": before,
        "protein_after": after,
        "frameshift": frameshift,
        "stop_lost": stop_lost,
        "indel_panel": indels.to_dict(orient="records"),
    }


@app.get("/health")
def health():
    return {"status": "ok"}
