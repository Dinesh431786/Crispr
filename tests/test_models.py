import numpy as np
import pytest

from crispr_app import models
from crispr_app.benchmark import spearman
from crispr_app.features import featurize, featurize_many, n_features, feature_names
from crispr_app.models import LinearModel, active_backend, predict_on_target
from crispr_app.train import fit_ridge


def test_feature_vector_shape():
    v = featurize("GACGATCAGTCAGGATCACC", "A")
    assert v.shape == (n_features(),)
    assert len(feature_names()) == n_features()


def test_active_backend_defaults_to_heuristic(monkeypatch, tmp_path):
    # With no model files at all, auto resolution must fall back to the surrogate.
    models.reset_caches()
    monkeypatch.setattr(models, "_LINEAR_PATH", str(tmp_path / "none.json"))
    monkeypatch.setattr(models, "_DEFAULT_PATH", str(tmp_path / "none_default.json"))
    monkeypatch.setattr(models, "_ONNX_PATH", str(tmp_path / "none.onnx"))
    models.reset_caches()
    assert active_backend() == "heuristic"
    s = predict_on_target("GACGATCAGTCAGGATCACC")
    assert 0.0 <= s <= 1.0
    models.reset_caches()


def test_shipped_default_model_is_used():
    # The version-controlled default model should make 'linear' the active backend.
    models.reset_caches()
    assert active_backend() == "linear"
    assert 0.0 <= predict_on_target("GACGATCAGTCAGGATCACC") <= 1.0
    models.reset_caches()


def test_linear_model_save_load_roundtrip(tmp_path):
    w = np.arange(n_features(), dtype=float) * 0.0 + 0.01
    m = LinearModel(w, 0.2, {"k": 1})
    p = tmp_path / "m.json"
    m.save(str(p))
    import json
    m2 = LinearModel.from_json(json.loads(p.read_text()))
    assert m2.intercept == 0.2
    assert np.allclose(m2.weights, w)


def test_ridge_training_learns_gc_signal():
    # Synthetic data: efficiency is a monotonic function of GC content.
    rng = np.random.default_rng(0)
    bases = "ACGT"
    guides, y = [], []
    for _ in range(120):
        g = "".join(rng.choice(list(bases)) for _ in range(20))
        gc = (g.count("G") + g.count("C")) / 20.0
        guides.append(g)
        y.append(gc + rng.normal(0, 0.02))
    X = featurize_many(guides)
    y = np.array(y)
    model = fit_ridge(X, y, alpha=1.0)
    pred = np.clip(X @ model.weights + model.intercept, 0, 1)
    assert spearman(pred, y) > 0.8


def test_predict_with_linear_model(monkeypatch, tmp_path):
    # A trained linear model should be picked up and used by the registry.
    rng = np.random.default_rng(1)
    guides = ["".join(rng.choice(list("ACGT")) for _ in range(20)) for _ in range(60)]
    y = np.array([(g.count("G") + g.count("C")) / 20.0 for g in guides])
    model = fit_ridge(featurize_many(guides), y, alpha=1.0)
    path = tmp_path / "linear.json"
    model.save(str(path))

    models.reset_caches()
    monkeypatch.setattr(models, "_LINEAR_PATH", str(path))
    monkeypatch.setattr(models, "_ONNX_PATH", str(tmp_path / "none.onnx"))
    models.reset_caches()
    assert active_backend() == "linear"
    assert 0.0 <= predict_on_target(guides[0]) <= 1.0
    models.reset_caches()
