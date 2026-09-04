from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "physical_parity_modulus_reach.py"


def _load_experiment():
    spec = importlib.util.spec_from_file_location("physical_parity_modulus_reach", EXPERIMENT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_reach_endpoint_reproduces_full_face_on_small_mesh():
    experiment = _load_experiment()
    upstream = experiment.parity.load_inputs()
    model = experiment.ReachModel(upstream, dimension=39, intervals=128)
    endpoint = 2 * upstream.RHO_STAR * upstream.OUTER_RADIUS
    reached = model.reach_face_matrix(endpoint)
    full = model.face_matrices()["Jfull"]
    assert np.max(np.abs(reached - full)) <= 2e-13 * np.max(np.abs(full))


def test_modulus_reach_record_replays():
    completed = subprocess.run(
        [sys.executable, ROOT / "scripts" / "check_physical_parity_modulus_reach.py"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    replay = json.loads(completed.stdout)
    assert replay["status"] == "checked-exploratory-record"
    assert replay["score_at_half"] < 1 < replay["score_at_0.52"]
    assert Fraction(1, 2) < replay["crossing_estimate"] < 0.52
