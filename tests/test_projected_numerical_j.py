from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


pytest.importorskip("scipy")

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_projected_importance_correction_matches_full_forms():
    full = _load(
        "projected_correction_full_reference",
        ROOT / "scripts/build_numerical_j_candidate_matrix.py",
    )
    projected = _load(
        "projected_correction_under_test",
        ROOT / "scripts/build_projected_numerical_j.py",
    )
    support_path = ROOT / "reproduction/212/support.json"
    support = full.builder.load_support(support_path)
    k, degree = 5, 2
    dimension = len(full.q.basis_indices(degree))
    rng = np.random.default_rng(212)
    projection, _ = np.linalg.qr(rng.normal(size=(dimension, 2)))
    zeros = np.zeros((dimension, dimension))
    _, full_i, full_j, _, _ = full.build_candidate_matrix(
        k=k,
        degree=degree,
        log2_n=8,
        seed=71201,
        batch_log2=4,
        support=support,
        base_i=zeros,
        base_j=zeros,
        include_i=True,
    )
    projected_support = projected.builder.load_support(support_path)
    projected_i, projected_j, _, _ = projected.projected_importance_correction(
        k=k,
        degree=degree,
        projection=projection,
        support=projected_support,
        log2_n=8,
        seed=71201,
        batch_log2=4,
        base_i=np.zeros((2, 2)),
        base_j=np.zeros((2, 2)),
    )
    np.testing.assert_allclose(projected_i, projection.T @ full_i @ projection)
    np.testing.assert_allclose(projected_j, projection.T @ full_j @ projection)
