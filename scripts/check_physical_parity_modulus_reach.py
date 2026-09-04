#!/usr/bin/env python3
"""Cheap replay of the recorded k=39 physical modulus-reach scan."""

from __future__ import annotations

import argparse
import json
from math import isclose
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECORD = ROOT / "experiments" / "physical_parity_modulus_reach.json"
FULL_FACE_RECORD = ROOT / "experiments" / "physical_parity_viability.json"


def _row_by_theta(payload: dict, theta: str) -> dict:
    return next(row for row in payload["rows"] if row["theta"] == theta)


def _production(row: dict, method: str = "fixed_finest_vector_quadratic_last_three") -> float:
    return float(row[method]["production"])


def check(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.physical-parity-modulus-reach-analysis.v1":
        raise ValueError("unexpected modulus-reach result schema")
    if payload.get("dimension") != 39 or payload.get("meshes") != [1024, 2048, 4096, 8192]:
        raise ArithmeticError("unexpected dimension or mesh ladder")
    if payload.get("production_intervals") != 98304:
        raise ArithmeticError("unexpected production mesh")
    rows = payload["rows"]
    theta_values = [float(row["theta_float"]) for row in rows]
    if theta_values != sorted(set(theta_values)):
        raise ArithmeticError("theta scan is not strictly ordered")

    score_half = _production(_row_by_theta(payload, "1/2"))
    score_051 = _production(_row_by_theta(payload, "51/100"))
    score_052 = _production(_row_by_theta(payload, "13/25"))
    if not score_half < score_051 < 1 < score_052:
        raise ArithmeticError("recorded half/0.51/0.52 classification changed")

    crossing = payload["crossing"]["fixed_finest_vector"]
    lower = _row_by_theta(payload, crossing["lower_theta"])
    upper = _row_by_theta(payload, crossing["upper_theta"])
    lower_score, upper_score = _production(lower), _production(upper)
    if not lower_score < 1 < upper_score:
        raise ArithmeticError("recorded crossing bracket does not straddle one")
    estimate = lower["theta_float"] + (1 - lower_score) * (
        upper["theta_float"] - lower["theta_float"]
    ) / (upper_score - lower_score)
    if not isclose(estimate, crossing["linear_interpolation"], abs_tol=2e-15):
        raise ArithmeticError("crossing interpolation mismatch")

    endpoint = _production(_row_by_theta(payload, "2742997/5000000"))
    full_face = json.loads(FULL_FACE_RECORD.read_text())
    k39 = next(row for row in full_face["dimensions"] if row["dimension"] == 39)
    expected_endpoint = k39["quadratic_in_inverse_mesh"]["production_mesh_score"]
    if not isclose(endpoint, expected_endpoint, abs_tol=2e-10):
        raise ArithmeticError("full-face endpoint no longer matches the parity record")

    optimized = payload["crossing"]["separately_optimized_crosscheck"]
    if abs(optimized["linear_interpolation"] - estimate) > 5e-5:
        raise ArithmeticError("fixed-vector and optimized crossing estimates disagree")
    return {
        "schema": "primegaps.physical-parity-modulus-reach-replay.v1",
        "status": "checked-exploratory-record",
        "classification": "theta approximately 0.51 (strictly above 1/2)",
        "score_at_half": score_half,
        "score_at_0.51": score_051,
        "score_at_0.52": score_052,
        "crossing_bracket": [lower["theta_float"], upper["theta_float"]],
        "crossing_estimate": estimate,
        "full_face_endpoint_score": endpoint,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("record", type=Path, nargs="?", default=DEFAULT_RECORD)
    args = parser.parse_args()
    print(json.dumps(check(args.record), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
