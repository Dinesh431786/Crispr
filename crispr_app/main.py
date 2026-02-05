from __future__ import annotations

from pathlib import Path
import time
import uuid

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field, field_validator

from analysis import find_gRNAs, find_off_targets_detailed, indel_simulations, simulate_protein_edit
from utils import validate_sequence

BASE_DIR = Path(__file__).resolve().parent
ALLOWED_PAMS = {"NGG", "NAG", "NG", "TTTV"}

app = FastAPI(title="CRISPR Precision Studio", version="2.1.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))


@app.middleware("http")
async def add_request_metrics(request: Request, call_next):
    request_id = str(uuid.uuid4())
    start = time.perf_counter()
    response = await call_next(request)
    elapsed_ms = round((time.perf_counter() - start) * 1000, 2)
    response.headers["X-Request-Id"] = request_id
    response.headers["X-Elapsed-Ms"] = str(elapsed_ms)
    return response


@app.exception_handler(HTTPException)
async def http_exception_handler(_: Request, exc: HTTPException):
    return JSONResponse(
        status_code=exc.status_code,
        content={"error": {"code": exc.status_code, "message": str(exc.detail)}},
    )


class DesignRequest(BaseModel):
    dna_sequence: str = Field(..., min_length=23, max_length=2_000_000)
    pam: str = "NGG"
    guide_length: int = Field(20, ge=18, le=25)
    min_gc: int = Field(40, ge=30, le=80)
    max_gc: int = Field(70, ge=40, le=90)
    add_5prime_g: bool = False

    @field_validator("pam")
    @classmethod
    def validate_pam(cls, value: str) -> str:
        if value not in ALLOWED_PAMS:
            raise ValueError(f"Unsupported PAM: {value}")
        return value


class OffTargetRequest(BaseModel):
    guides: list[str] = Field(..., min_length=1, max_length=500)
    background_sequence: str = Field(..., min_length=23, max_length=2_000_000)
    max_mismatches: int = Field(2, ge=0, le=4)
    pam: str = "NGG"
    include_reverse_strand: bool = True

    @field_validator("pam")
    @classmethod
    def validate_pam(cls, value: str) -> str:
        if value not in ALLOWED_PAMS:
            raise ValueError(f"Unsupported PAM: {value}")
        return value


class SimulateRequest(BaseModel):
    dna_sequence: str = Field(..., min_length=23, max_length=2_000_000)
    guide: str = Field(..., min_length=18, max_length=25)
    edit_offset: int = Field(20, ge=0, le=30)
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
        "metadata": {"pam": payload.pam, "guide_length": payload.guide_length},
    }


@app.post("/api/offtargets")
def find_offtargets(payload: OffTargetRequest):
    import pandas as pd

    normalized_guides = [g.upper().strip() for g in payload.guides if g and g.strip()]
    if not normalized_guides:
        raise HTTPException(status_code=400, detail="No valid guides provided")

    df = pd.DataFrame({"gRNA": normalized_guides})
    results = find_off_targets_detailed(
        df,
        payload.background_sequence,
        payload.max_mismatches,
        pam=payload.pam,
        include_reverse_strand=payload.include_reverse_strand,
    )
    if results.empty:
        return {"count": 0, "off_targets": []}

    return {
        "count": int(len(results)),
        "off_targets": results.head(1000).to_dict(orient="records"),
        "metadata": {
            "max_mismatches": payload.max_mismatches,
            "include_reverse_strand": payload.include_reverse_strand,
        },
    }


@app.post("/api/simulate")
def simulate(payload: SimulateRequest):
    clean_seq = payload.dna_sequence.upper().replace("\n", "").replace(" ", "")
    idx = clean_seq.find(payload.guide.upper())
    if idx < 0:
        raise HTTPException(status_code=400, detail="Guide not found in DNA sequence")

    cut = idx + payload.edit_offset
    before, after, frameshift, stop_lost = simulate_protein_edit(clean_seq, cut, payload.edit_type)
    indels = indel_simulations(clean_seq, cut)

    return {
        "protein_before": before,
        "protein_after": after,
        "frameshift": frameshift,
        "stop_lost": stop_lost,
        "indel_panel": indels.to_dict(orient="records"),
    }


@app.get("/health")
def health():
    return {"status": "ok", "version": app.version}
