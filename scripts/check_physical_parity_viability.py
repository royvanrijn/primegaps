#!/usr/bin/env python3
"""Cheap replay of the recorded quarter-rough parity viability calculation."""

from __future__ import annotations

import argparse
import json
from math import isclose
from pathlib import Path

from primegaps.parity import (
    parity_contributions,
    parity_error_budget,
    rough_factor_constants,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECORD = ROOT / "experiments/physical_parity_viability.json"


def check(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.physical-parity-production-extrapolation.v1":
        raise ValueError("unexpected parity viability result schema")
    constants = rough_factor_constants(payload["beta"])
    recorded_constants = payload["rough_factor_constants"]
    for field in (
        "prime",
        "semiprime",
        "triprime",
        "omega_choose_1",
        "omega_choose_2",
        "omega_choose_3",
        "signed_identity",
        "gross_signed_condition",
    ):
        if not isclose(getattr(constants, field), recorded_constants[field], abs_tol=2e-13):
            raise ArithmeticError(f"rough-factor constant mismatch: {field}")

    checked = {}
    for row in payload["dimensions"]:
        dimension = int(row["dimension"])
        score = float(row["quadratic_in_inverse_mesh"]["production_mesh_score"])
        terms = parity_contributions(score, constants)
        eta, relative = parity_error_budget(terms)
        recorded_terms = row["factorial_contributions_normalized_by_I"]
        comparisons = {
            "plus_omega": terms.plus_omega,
            "minus_2_choose_2": terms.minus_2_choose_2,
            "plus_3_choose_3": terms.plus_3_choose_3,
        }
        for field, value in comparisons.items():
            if not isclose(value, recorded_terms[field], abs_tol=3e-13):
                raise ArithmeticError(f"k={dimension} contribution mismatch: {field}")
        budget = row["parity_error_budget"]
        if not isclose(terms.signed_sum, row["signed_sum"], abs_tol=4e-13):
            raise ArithmeticError(f"k={dimension} signed sum mismatch")
        if not isclose(eta, budget["eta_max_for_abs_error_le_eta_I"], abs_tol=4e-13):
            raise ArithmeticError(f"k={dimension} absolute error budget mismatch")
        if not isclose(
            relative,
            budget["common_relative_error_budget_across_unsigned_terms"],
            abs_tol=4e-13,
        ):
            raise ArithmeticError(f"k={dimension} relative error budget mismatch")
        if score <= 1.0 or eta <= 0.0:
            raise ArithmeticError(f"k={dimension} does not cross one")
        checked[str(dimension)] = {
            "score": score,
            "eta_max": eta,
            "common_relative_budget": relative,
        }

    calibration = payload["calibration"]
    if abs(calibration["prediction_minus_published_lower"]) > 1e-5:
        raise ArithmeticError("the mesh calibration misses the published k=40 bound")
    return {
        "schema": "primegaps.physical-parity-replay.v1",
        "status": "checked-exploratory-record",
        "dimensions": checked,
        "calibration_error": calibration["prediction_minus_published_lower"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("record", type=Path, nargs="?", default=DEFAULT_RECORD)
    args = parser.parse_args()
    print(json.dumps(check(args.record), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
