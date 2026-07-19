"""Pluggable on-target predictor registry.

Resolution order (``backend="auto"``):

1. **onnx**     - if ``onnxruntime`` is installed and a model file is available
                  (``$CRISPR_ONNX_MODEL`` or ``models/ontarget.onnx``). Lets users
                  drop in a real DeepSpCas9/CRISPRon export for deep-model accuracy.
2. **linear**   - if learned coefficients exist (``$CRISPR_LINEAR_MODEL`` or
                  ``models/linear.json``), produced by ``train.py`` on real data.
3. **heuristic** - the interpretable surrogate in :mod:`scoring` (always available).

Every prediction reports which backend produced it, so scores are never an
unlabelled black box.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

try:
    from conformal import interval as _conf_interval
    from features import featurize
    from scoring import on_target_score
except ImportError:  # pragma: no cover
    from .conformal import interval as _conf_interval
    from .features import featurize
    from .scoring import on_target_score

_BASE_DIR = Path(__file__).resolve().parent
_LINEAR_PATH = os.environ.get("CRISPR_LINEAR_MODEL", str(_BASE_DIR / "models" / "linear.json"))
# Shipped, version-controlled default model (pooled human SpCas9). A user-trained
# linear.json takes precedence over it.
_DEFAULT_PATH = str(_BASE_DIR / "models" / "default.json")
# Goal-aware "knockout mode": a model trained on out-of-frame (frameshift)
# efficiency instead of total cutting — ranks guides by functional knockout.
_OOF_PATH = str(_BASE_DIR / "models" / "default_oof.json")
_ONNX_PATH = os.environ.get("CRISPR_ONNX_MODEL", str(_BASE_DIR / "models" / "ontarget.onnx"))

_linear_cache: "LinearModel | None" = None
_oof_cache: "LinearModel | None | bool" = False   # False = not yet tried
_onnx_cache = None
_onnx_tried = False


class LinearModel:
    """Linear/ridge model: raw = x . w + b, then an optional monotone calibration
    maps the raw score back to the 0-1 efficiency scale.

    When the model is trained on a rank-transformed target (to optimise the
    ranking metric, Spearman, rather than squared error), the raw linear output
    is in rank units, not efficiency. ``calibration`` = (xs, ys) is a monotone
    non-decreasing lookup (isotonic fit on the training set) that maps raw -> the
    0-1 efficiency scale via ``np.interp``. Being monotone, it preserves the
    ranking exactly while restoring the efficiency scale the UI and the conformal
    intervals live on. Absent -> the raw score is already the efficiency estimate
    (the classic MSE-trained model), so behaviour is unchanged.
    """

    def __init__(self, weights: np.ndarray, intercept: float, meta: dict | None = None,
                 calibration: tuple[np.ndarray, np.ndarray] | None = None,
                 normconformal: dict | None = None):
        self.weights = np.asarray(weights, dtype=np.float64)
        self.intercept = float(intercept)
        self.meta = meta or {}
        self.calibration = None
        if calibration is not None:
            xs, ys = calibration
            self.calibration = (np.asarray(xs, dtype=np.float64), np.asarray(ys, dtype=np.float64))
        # Normalized (locally-adaptive) conformal: a lightweight sigma-model gives
        # each guide its own interval WIDTH (heteroscedastic), so uncertainty is
        # guide-specific -- and can drive ranking -- while marginal coverage is
        # still guaranteed. dict: {weights, intercept, floor, q:{q80,q90}}.
        self.normconformal = None
        if normconformal is not None:
            self.normconformal = {
                "weights": np.asarray(normconformal["weights"], dtype=np.float64),
                "intercept": float(normconformal["intercept"]),
                "floor": float(normconformal["floor"]),
                "q": {k: float(v) for k, v in normconformal["q"].items()},
            }

    def sigma(self, x: np.ndarray) -> float:
        """Per-guide uncertainty scale (floored), or None if not calibrated."""
        nc = self.normconformal
        if nc is None:
            return None
        return max(float(x @ nc["weights"] + nc["intercept"]), nc["floor"])

    def predict(self, guide: str, ngg_context: str | None = None,
                up: str = "", down: str = "") -> float:
        x = featurize(guide, ngg_context, up, down)
        if x.shape[0] != self.weights.shape[0]:
            raise ValueError("feature/weight dimension mismatch")
        y = float(x @ self.weights + self.intercept)
        if self.calibration is not None:
            y = float(np.interp(y, self.calibration[0], self.calibration[1]))
        return round(max(0.0, min(1.0, y)), 3)

    def to_json(self) -> dict:
        d = {"weights": self.weights.tolist(), "intercept": self.intercept, "meta": self.meta}
        if self.calibration is not None:
            d["calibration"] = {"x": self.calibration[0].tolist(), "y": self.calibration[1].tolist()}
        if self.normconformal is not None:
            nc = self.normconformal
            d["normconformal"] = {"weights": nc["weights"].tolist(), "intercept": nc["intercept"],
                                  "floor": nc["floor"], "q": nc["q"]}
        return d

    @classmethod
    def from_json(cls, data: dict) -> "LinearModel":
        c = data.get("calibration")
        calib = (c["x"], c["y"]) if c else None
        return cls(np.array(data["weights"], dtype=np.float64), data["intercept"],
                   data.get("meta", {}), calib, data.get("normconformal"))

    def save(self, path: str) -> None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(json.dumps(self.to_json()))


def _load_linear() -> "LinearModel | None":
    global _linear_cache
    if _linear_cache is not None:
        return _linear_cache
    for path in (_LINEAR_PATH, _DEFAULT_PATH):  # user-trained model wins over shipped default
        p = Path(path)
        if p.is_file():
            try:
                _linear_cache = LinearModel.from_json(json.loads(p.read_text()))
                break
            except Exception:
                _linear_cache = None
    return _linear_cache


def _load_oof() -> "LinearModel | None":
    """Load the out-of-frame (knockout-mode) model, if shipped."""
    global _oof_cache
    if _oof_cache is not False:
        return _oof_cache or None
    _oof_cache = None
    p = Path(_OOF_PATH)
    if p.is_file():
        try:
            _oof_cache = LinearModel.from_json(json.loads(p.read_text()))
        except Exception:
            _oof_cache = None
    return _oof_cache


def _goal_model(goal: str) -> "LinearModel | None":
    """Model for a design goal: 'knockout' -> out-of-frame model (if present)."""
    if goal == "knockout":
        m = _load_oof()
        if m is not None:
            return m
    return _load_linear()


def _load_onnx():
    global _onnx_cache, _onnx_tried
    if _onnx_tried:
        return _onnx_cache
    _onnx_tried = True
    if not Path(_ONNX_PATH).is_file():
        return None
    try:
        import onnxruntime as ort  # type: ignore
        _onnx_cache = ort.InferenceSession(_ONNX_PATH, providers=["CPUExecutionProvider"])
    except Exception:
        _onnx_cache = None
    return _onnx_cache


def available_backends() -> list[str]:
    # CRISPRscan (peer-reviewed, published weights) is always available for
    # 20-nt NGG guides with flanking context; applied in the design pipeline.
    backends = ["heuristic", "crisprscan"]
    if _load_linear() is not None:
        backends.append("linear")
    if _load_onnx() is not None:
        backends.append("onnx")
    return backends


def active_backend(backend: str = "auto") -> str:
    if backend != "auto":
        return backend
    if _load_onnx() is not None:
        return "onnx"
    if _load_linear() is not None:
        return "linear"
    return "heuristic"


def predict_on_target(guide: str, ngg_context: str | None = None, backend: str = "auto",
                      goal: str = "general", up: str = "", down: str = "") -> float:
    """On-target score. ``goal='knockout'`` uses the out-of-frame model (rank by
    predicted frameshift) when available; otherwise the general cutting model.
    ``up``/``down`` are the flanking sequence context (6 nt 5', PAM+6 nt 3')."""
    if goal == "knockout":
        m = _load_oof()
        if m is not None:
            try:
                return m.predict(guide, ngg_context, up, down)
            except Exception:
                pass
    chosen = active_backend(backend)
    if chosen == "onnx":
        sess = _load_onnx()
        if sess is not None:
            try:
                x = featurize(guide, ngg_context, up, down).astype(np.float32)[None, :]
                out = sess.run(None, {sess.get_inputs()[0].name: x})[0]
                return round(float(np.clip(np.ravel(out)[0], 0.0, 1.0)), 3)
            except Exception:
                pass  # fall through to next backend
    if chosen in ("onnx", "linear"):
        lm = _load_linear()
        if lm is not None:
            try:
                return lm.predict(guide, ngg_context, up, down)
            except Exception:
                pass
    return on_target_score(guide, ngg_context)


def predict_interval(guide: str, ngg_context: str | None = None, level: str = "q90",
                     goal: str = "general", up: str = "", down: str = "") -> dict | None:
    """Conformal prediction interval for the (goal-appropriate) on-target model.

    Returns {"point", "low", "high", "level", "coverage"} on the 0-1 scale, or
    None when the active model has no conformal calibration (e.g. heuristic).
    """
    lm = _goal_model(goal)
    if lm is None or active_backend() == "onnx":
        return None
    point = lm.predict(guide, ngg_context, up, down)
    # Prefer the normalized (per-guide width) interval when the model ships it;
    # fall back to the classic constant-width conformal quantile otherwise.
    nc = lm.normconformal
    if nc is not None and level in nc["q"]:
        sig = lm.sigma(featurize(guide, ngg_context, up, down))
        half = nc["q"][level] * sig
    else:
        conf = (lm.meta or {}).get("conformal") or {}
        if level not in conf:
            return None
        half = float(conf[level])
    lo, hi = _conf_interval(point, half)
    return {"point": point, "low": lo, "high": hi, "level": level,
            "coverage": int(level[1:]) / 100.0, "halfwidth": round(float(half), 3)}


def reset_caches() -> None:
    """Test/CLI helper to forget loaded models."""
    global _linear_cache, _oof_cache, _onnx_cache, _onnx_tried
    _linear_cache = None
    _oof_cache = False
    _onnx_cache = None
    _onnx_tried = False
