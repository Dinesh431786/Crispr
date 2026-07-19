"""Base-editing (ABE / CBE) target analysis — a first-class editing module.

Base editors install point mutations without a double-strand break:

* **CBE** (cytosine base editor): C->T on the protospacer (spacer) strand.
* **ABE** (adenine base editor): A->G on the protospacer strand.

Editing is efficient only within a narrow *activity window* on the protospacer,
classically positions ~4-8 (1-based, counting the PAM-distal 5' end as 1, PAM at
21-23) for BE3/ABE7.10-class editors (Komor et al. 2016; Gaudelli et al. 2017),
with efficiency peaking near position ~6 and falling off toward the window edges.

Beyond *whether* a base is editable, this module scores *how well*:

* **Window efficiency** — a position-weighted relative profile (peak = 1.0),
  approximating the published activity window. Relative, not an absolute rate.
* **Bystander editing** — when several bases of the target type sit in the window
  the editor edits them too, yielding a mixture of products. We report the
  bystanders and the **editing purity** = intended-site weight / all-editable
  weight (the fraction of edits expected to be the intended one).
* **Composite base-edit score** = window-efficiency x purity, so a guide is
  rewarded for placing the target base at a high-efficiency position *and* for
  being a clean single edit. Multiply by specificity when off-target data exists.

Everything is a transparent, position-weighted calculation — no black box.
"""

from __future__ import annotations

import math

WINDOW_START = 4   # 1-based, inclusive
WINDOW_END = 8

_EDITS = {"CBE": ("C", "T"), "ABE": ("A", "G")}

# Relative editing efficiency across the window, peaked near position ~6 (a
# Gaussian approximation to the published BE3 / ABE7.10-class activity window;
# peak normalised to 1.0). These are RELATIVE weights for ranking, not absolute
# editing rates, which vary by editor variant, locus, and delivery.
_PEAK = 6.0
_SPREAD = 1.6


def position_weight(pos: int, window: tuple[int, int] = (WINDOW_START, WINDOW_END)) -> float:
    """Relative editing efficiency (0-1) at a 1-based protospacer position."""
    lo, hi = window
    if not (lo <= pos <= hi):
        return 0.0
    return round(math.exp(-((pos - _PEAK) / _SPREAD) ** 2), 3)


def editable_targets(guide: str, window: tuple[int, int] = (WINDOW_START, WINDOW_END)) -> list[dict]:
    """Return base-editable positions for a guide within the activity window.

    Each entry: editor (CBE/ABE), 1-based position in the spacer, the original
    base, the edited base, and the position's relative window weight.
    """
    guide = guide.upper()
    lo, hi = window
    targets: list[dict] = []
    for editor, (frm, to) in _EDITS.items():
        for pos in range(lo, hi + 1):
            if pos <= len(guide) and guide[pos - 1] == frm:
                targets.append({"editor": editor, "pos": pos, "from": frm, "to": to,
                                "weight": position_weight(pos, window)})
    return targets


def assess(guide: str, editor: str, target_pos: int | None = None,
           window: tuple[int, int] = (WINDOW_START, WINDOW_END)) -> dict:
    """Full base-edit assessment of a guide for one editor.

    ``target_pos`` (1-based spacer position) names the base you intend to edit;
    if omitted or not editable, the highest-efficiency editable base is used as
    the intended target. Bystanders are the other same-type bases in the window.
    """
    guide = guide.upper()
    frm, to = _EDITS[editor]
    lo, hi = window
    hits = [{"pos": p, "weight": position_weight(p, window)}
            for p in range(lo, hi + 1) if p <= len(guide) and guide[p - 1] == frm]
    if not hits:
        return {"editor": editor, "editable": False, "target_pos": None,
                "efficiency": 0.0, "bystanders": [], "purity": 0.0, "be_score": 0.0,
                "explanation": f"No editable {frm} in positions {lo}-{hi} — not {editor}-editable."}

    if target_pos is not None and any(h["pos"] == target_pos for h in hits):
        intended = next(h for h in hits if h["pos"] == target_pos)
    else:
        intended = max(hits, key=lambda h: h["weight"])
    bystanders = [h for h in hits if h["pos"] != intended["pos"]]

    w_int = intended["weight"]
    w_tot = w_int + sum(h["weight"] for h in bystanders)
    purity = round(w_int / w_tot, 3) if w_tot > 0 else 0.0
    be_score = round(w_int * purity, 3)

    if bystanders:
        by = ", ".join(f"{frm}@{h['pos']} ({h['weight']:.2f})" for h in bystanders)
        expl = (f"{editor}: intended {frm}->{to} at position {intended['pos']} "
                f"(window efficiency {w_int:.2f}); {len(bystanders)} bystander "
                f"{frm} in window [{by}] -> editing purity {purity:.0%}.")
    else:
        expl = (f"{editor}: intended {frm}->{to} at position {intended['pos']} "
                f"(window efficiency {w_int:.2f}); no bystander {frm} in window "
                f"-> clean single edit (purity 100%).")

    return {"editor": editor, "editable": True, "target_pos": intended["pos"],
            "to_base": to, "efficiency": w_int,
            "bystanders": [{"pos": h["pos"], "weight": h["weight"]} for h in bystanders],
            "purity": purity, "be_score": be_score, "explanation": expl}


def best_score(guide: str, editor: str = "any",
               window: tuple[int, int] = (WINDOW_START, WINDOW_END)) -> float:
    """Best composite base-edit score for a guide over the allowed editor(s)."""
    editors = ("CBE", "ABE") if editor == "any" else (editor,)
    return max((assess(guide, e, window=window)["be_score"] for e in editors), default=0.0)


def summarize(guide: str, editor: str = "any",
              window: tuple[int, int] = (WINDOW_START, WINDOW_END)) -> dict:
    """Compact per-guide summary for the API/UI, with efficiency + purity."""
    guide = guide.upper()
    t = editable_targets(guide, window)
    cbe = [x["pos"] for x in t if x["editor"] == "CBE"]
    abe = [x["pos"] for x in t if x["editor"] == "ABE"]
    editors = ("CBE", "ABE") if editor == "any" else (editor,)
    assessments = {e: assess(guide, e, window=window) for e in editors}
    best = max(assessments.values(), key=lambda a: a["be_score"], default=None)
    return {
        "gRNA": guide,
        "window": f"{window[0]}-{window[1]}",
        "CBE_positions": cbe,
        "ABE_positions": abe,
        "editable": bool(cbe or abe),
        "assessments": assessments,
        "best_editor": best["editor"] if best and best["editable"] else None,
        "be_score": best["be_score"] if best else 0.0,
        "explanation": best["explanation"] if best else "",
    }
