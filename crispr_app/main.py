from __future__ import annotations

import hashlib
import os
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field
from starlette.requests import Request

from analysis import (
    design_multiplex,
    design_prime_editing_pegRNAs,
    find_gRNAs,
    find_off_targets_detailed,
    indel_simulations,
    simulate_protein_edit,
    summarize_specificity,
)
from base_edit import summarize as base_edit_summarize
from models import active_backend, available_backends, predict_interval
from scoring import score_breakdown
from structure import self_complementarity
from utils import load_fasta_text, validate_sequence

BASE_DIR = Path(__file__).resolve().parent
MAX_SEQ = 200_000  # input guard: max DNA characters accepted by JSON endpoints

app = FastAPI(title="CRISPR Precision Studio", version="3.5.0")

# CORS: configurable for deployment. Default "*" (open) but WITHOUT credentials,
# which is the only spec-valid combination; set CRISPR_CORS_ORIGINS to a
# comma-separated allowlist to lock it down in production.
_origins_env = os.environ.get("CRISPR_CORS_ORIGINS", "*").strip()
_origins = ["*"] if _origins_env == "*" else [o.strip() for o in _origins_env.split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_origins,
    allow_credentials=_origins != ["*"],
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))


def _asset_version() -> str:
    """Short hash of the static assets, for cache-busting the browser on update."""
    h = hashlib.md5()
    for f in ("static/app.js", "static/style.css"):
        try:
            h.update((BASE_DIR / f).read_bytes())
        except FileNotFoundError:
            pass
    return h.hexdigest()[:8]


ASSET_V = _asset_version()


class DesignRequest(BaseModel):
    dna_sequence: str = Field(..., min_length=23, max_length=MAX_SEQ)
    pam: str = "NGG"
    guide_length: int = Field(20, ge=18, le=25)
    min_gc: int = Field(40, ge=30, le=80)
    max_gc: int = Field(70, ge=40, le=90)
    add_5prime_g: bool = False
    goal: str = "general"   # general | knockout | knockin | crispri | baseedit
    target_pos: int | None = Field(None, ge=0)   # edit site (knock-in) or TSS (CRISPRi/a)
    editor: str = "any"     # base-edit: "ABE" | "CBE" | "any"
    # Uncertainty-aware ranking: balanced (point score) | conservative (pessimistic
    # bound) | robust (uncertainty-penalised) | optimistic (best case).
    ranking_strategy: str = "balanced"


class OffTargetRequest(BaseModel):
    guides: list[str] = Field(..., max_length=2000)
    background_sequence: str = Field(..., max_length=MAX_SEQ)
    max_mismatches: int = Field(2, ge=0, le=4)
    pam: str = "NGG"


class SimulateRequest(BaseModel):
    dna_sequence: str = Field(..., max_length=MAX_SEQ)
    guide: str = Field(..., min_length=10, max_length=40)
    edit_offset: int = 20
    edit_type: str = "del1"


class PrimeDesignRequest(BaseModel):
    dna_sequence: str = Field(..., max_length=MAX_SEQ)
    target_pos: int
    desired_edit: str = Field(..., max_length=1)


@app.get("/", response_class=HTMLResponse)
def ui(request: Request):
    resp = templates.TemplateResponse(request, "index.html", {"asset_v": ASSET_V})
    # Never cache the HTML shell, so updated UI (and versioned assets) always load.
    resp.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    return resp


@app.post("/api/design")
def design_guides(payload: DesignRequest):
    ok, cleaned = validate_sequence(payload.dna_sequence)
    if not ok:
        raise HTTPException(status_code=400, detail=cleaned)

    goal = payload.goal if payload.goal in ("general", "knockout", "knockin", "crispri", "baseedit") else "general"
    strategy = payload.ranking_strategy if payload.ranking_strategy in (
        "balanced", "conservative", "robust", "optimistic") else "balanced"
    guides = find_gRNAs(
        cleaned,
        pam=payload.pam,
        guide_length=payload.guide_length,
        min_gc=payload.min_gc,
        max_gc=payload.max_gc,
        add_5prime_g=payload.add_5prime_g,
        goal=goal,
        target_pos=payload.target_pos,
        editor=payload.editor if payload.editor in ("ABE", "CBE", "any") else "any",
        ranking_strategy=strategy,
    )

    # Confidence intervals (flank-aware, goal-appropriate) are computed in
    # find_gRNAs so they match the score's flanking context.
    top = guides.head(100).to_dict(orient="records")
    return {
        "count": int(len(guides)),
        "model": active_backend(),
        "goal": goal,
        "ranking_strategy": strategy,
        "top_guides": top,
    }


@app.post("/api/offtargets")
def find_offtargets(payload: OffTargetRequest):
    import pandas as pd

    df = pd.DataFrame({"gRNA": payload.guides})
    results = find_off_targets_detailed(df, payload.background_sequence, payload.max_mismatches, pam=payload.pam)
    specificity = summarize_specificity(results, list(payload.guides))

    if results.empty:
        return {"count": 0, "off_targets": [], "specificity": specificity.to_dict(orient="records")}

    return {
        "count": int(len(results)),
        "off_targets": results.sort_values("CFD_Score", ascending=False).head(500).to_dict(orient="records"),
        "specificity": specificity.to_dict(orient="records"),
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


@app.post("/api/prime-design")
def prime_design(payload: PrimeDesignRequest):
    ok, cleaned = validate_sequence(payload.dna_sequence)
    if not ok:
        raise HTTPException(status_code=400, detail=cleaned)

    pe_guides = design_prime_editing_pegRNAs(cleaned, payload.target_pos, payload.desired_edit)

    return {
        "count": int(len(pe_guides)),
        "pegRNAs": pe_guides.head(50).to_dict(orient="records"),
    }


class ExplainRequest(BaseModel):
    guide: str = Field(..., min_length=18)
    ngg_context: str | None = None
    goal: str = "general"


class FastaRequest(BaseModel):
    contents: str = Field(..., max_length=MAX_SEQ)


class GenomeOffTargetRequest(BaseModel):
    guides: list[str] = Field(..., max_length=500)
    fasta: str = Field(..., max_length=20_000_000)  # ~20 MB (e.g. a bacterial genome)
    max_mismatches: int = Field(3, ge=0, le=4)
    pam: str = "NGG"


@app.post("/api/explain")
def explain(payload: ExplainRequest):
    """Interpretable breakdown of an on-target score plus a conformal interval."""
    out = score_breakdown(payload.guide, payload.ngg_context)
    goal = payload.goal if payload.goal in ("general", "knockout") else "general"
    ci = predict_interval(payload.guide, payload.ngg_context, goal=goal)
    if ci is not None:
        out["interval"] = ci  # distribution-free CI with guaranteed coverage
    # Informational structural QC (NOT part of the score: in our benchmarks it
    # showed no independent effect on efficiency — flagged for transparency only).
    out["self_complementarity"] = self_complementarity(payload.guide)
    return out


@app.post("/api/upload-fasta")
def upload_fasta(payload: FastaRequest):
    """Parse pasted FASTA (or plain DNA) into a clean sequence."""
    sequence, message = load_fasta_text(payload.contents)
    if not sequence:
        raise HTTPException(status_code=400, detail=message or "Could not parse FASTA.")
    return {"sequence": sequence, "length": len(sequence), "message": message}


@app.post("/api/offtargets-genome")
def offtargets_genome(payload: GenomeOffTargetRequest):
    """Genome-wide off-target scan over a (multi-record) FASTA, both strands."""
    from genome import genome_specificity, scan_fasta_text

    guides = [g.strip().upper() for g in payload.guides if g.strip()]
    if not guides:
        raise HTTPException(status_code=400, detail="No guides provided.")
    hits = scan_fasta_text(payload.fasta, guides, payload.max_mismatches, payload.pam)
    spec = genome_specificity(hits, guides)
    top = hits.sort_values("CFD_Score", ascending=False).head(500) if not hits.empty else hits
    return {
        "count": int(len(hits)),
        "off_targets": top.to_dict(orient="records"),
        "specificity": spec.to_dict(orient="records"),
    }


class BaseEditRequest(BaseModel):
    guides: list[str]
    editor: str = "any"     # "ABE" | "CBE" | "any"


@app.post("/api/base-edit")
def base_edit(payload: BaseEditRequest):
    """Base-editing (ABE/CBE) studio: per guide, window efficiency, bystander
    edits, editing purity, and a composite base-edit score with an explanation."""
    editor = payload.editor if payload.editor in ("ABE", "CBE", "any") else "any"
    guides = [g.strip().upper() for g in payload.guides if g.strip()][:500]
    results = [base_edit_summarize(g, editor) for g in guides]
    results.sort(key=lambda r: r["be_score"], reverse=True)
    return {"editor": editor, "results": results}


class MultiplexRequest(BaseModel):
    dna_sequence: str = Field(..., min_length=23, max_length=MAX_SEQ)
    n_guides: int = Field(6, ge=2, le=30)
    pam: str = "NGG"
    guide_length: int = Field(20, ge=18, le=25)
    min_gc: int = Field(40, ge=30, le=80)
    max_gc: int = Field(70, ge=40, le=90)
    goal: str = "general"
    ranking_strategy: str = "balanced"
    diversity_weight: float = Field(0.5, ge=0.0, le=2.0)
    min_spacing: int = Field(0, ge=0, le=1000)


@app.post("/api/design-multiplex")
def design_multiplex_endpoint(payload: MultiplexRequest):
    """Greedy marginal-gain designer for a multiplex / pooled guide library:
    individually strong yet collectively diverse (low recombination risk)."""
    ok, cleaned = validate_sequence(payload.dna_sequence)
    if not ok:
        raise HTTPException(status_code=400, detail=cleaned)
    goal = payload.goal if payload.goal in ("general", "knockout", "knockin", "crispri", "baseedit") else "general"
    strategy = payload.ranking_strategy if payload.ranking_strategy in (
        "balanced", "conservative", "robust", "optimistic") else "balanced"
    picked, summary = design_multiplex(
        cleaned, n_guides=payload.n_guides, pam=payload.pam,
        guide_length=payload.guide_length, min_gc=payload.min_gc, max_gc=payload.max_gc,
        goal=goal, ranking_strategy=strategy,
        diversity_weight=payload.diversity_weight, min_spacing=payload.min_spacing,
    )
    return {"summary": summary, "guides": picked.to_dict(orient="records")}


@app.get("/api/models")
def models_info():
    """Report which on-target backends are available and which is active."""
    return {"active": active_backend(), "available": available_backends()}


@app.get("/health")
def health():
    return {"status": "ok"}
