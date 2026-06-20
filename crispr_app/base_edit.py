"""Base-editing (ABE / CBE) target analysis.

Base editors install point mutations without a double-strand break:

* **CBE** (cytosine base editor): C->T on the protospacer (spacer) strand.
* **ABE** (adenine base editor): A->G on the protospacer strand.

Editing is efficient only within a narrow *activity window* on the protospacer,
classically positions ~4-8 (1-based, counting the PAM-distal 5' end as 1, PAM at
21-23) for BE3/ABE7.10-class editors (Komor et al. 2016; Gaudelli et al. 2017).
This module reports, for a guide, which bases fall in that window and the edit
each base editor would make.
"""

from __future__ import annotations

WINDOW_START = 4   # 1-based, inclusive
WINDOW_END = 8

_EDITS = {"CBE": ("C", "T"), "ABE": ("A", "G")}


def editable_targets(guide: str, window: tuple[int, int] = (WINDOW_START, WINDOW_END)) -> list[dict]:
    """Return base-editable positions for a guide within the activity window.

    Each entry: editor (CBE/ABE), 1-based position in the spacer, the original
    base and the edited base.
    """
    guide = guide.upper()
    lo, hi = window
    targets: list[dict] = []
    for editor, (frm, to) in _EDITS.items():
        for pos in range(lo, hi + 1):
            if pos <= len(guide) and guide[pos - 1] == frm:
                targets.append({"editor": editor, "pos": pos, "from": frm, "to": to})
    return targets


def summarize(guide: str, window: tuple[int, int] = (WINDOW_START, WINDOW_END)) -> dict:
    """Compact per-guide summary for the API/UI."""
    t = editable_targets(guide, window)
    cbe = [x["pos"] for x in t if x["editor"] == "CBE"]
    abe = [x["pos"] for x in t if x["editor"] == "ABE"]
    return {
        "gRNA": guide,
        "window": f"{window[0]}-{window[1]}",
        "CBE_positions": cbe,
        "ABE_positions": abe,
        "editable": bool(cbe or abe),
    }
