#!/usr/bin/env python3
"""Analyze mesh convergence for the k=39 modulus-reach scan."""

from __future__ import annotations

import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "experiments" / "physical_parity_viability.py"
RHO_STAR = Fraction(2_624_989, 10_000_000)
PRODUCTION_INTERVALS = 98_304


def load_model_module():
    spec = importlib.util.spec_from_file_location("physical_parity_viability", SOURCE)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


parity = load_model_module()


def key(prefix: str, theta: str) -> str:
    return f"{prefix}_{theta.replace('/', '_')}"


def rayleigh(matrix, metric, vector):
    return float(vector @ matrix @ vector) / float(vector @ metric @ vector)


def quadratic_prediction(meshes, values, start=1):
    x = 1 / np.asarray(meshes[start:], dtype=float)
    coefficients = np.polyfit(x, np.asarray(values[start:], dtype=float), 2)
    return {
        "continuum": float(coefficients[-1]),
        "production": float(np.polyval(coefficients, 1 / PRODUCTION_INTERVALS)),
    }


def interpolate_crossing(rows, value_key):
    for left, right in zip(rows, rows[1:]):
        a, b = left[value_key], right[value_key]
        if a <= 1 < b:
            x0, x1 = left["theta_float"], right["theta_float"]
            estimate = x0 + (1 - a) * (x1 - x0) / (b - a)
            return {
                "lower_theta": left["theta"],
                "lower_score": a,
                "upper_theta": right["theta"],
                "upper_score": b,
                "linear_interpolation": estimate,
            }
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, nargs="+", required=True)
    parser.add_argument("--matrices", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if len(args.results) != len(args.matrices):
        raise ValueError("one matrix archive per result is required")
    results = [json.loads(path.read_text()) for path in args.results]
    archives = [np.load(path) for path in args.matrices]
    meshes = [int(row["intervals"]) for row in results]
    if meshes != sorted(meshes) or len(meshes) < 4:
        raise ValueError("four increasingly fine meshes are required")
    theta_rows = results[-1]["rows"]
    theta_labels = [row["theta"] for row in theta_rows]
    for result in results:
        if [row["theta"] for row in result["rows"]] != theta_labels:
            raise ValueError("theta grids differ")

    analyzed = []
    for theta_row in theta_rows:
        theta = theta_row["theta"]
        finest_vector = archives[-1][key("vector", theta)]
        fixed_scores = []
        optimized_scores = []
        for result, archive in zip(results, archives, strict=True):
            fixed_scores.append(
                float(RHO_STAR)
                * rayleigh(archive[key("theta", theta)], archive["I"], finest_vector)
            )
            optimized_scores.append(
                next(row["score"] for row in result["rows"] if row["theta"] == theta)
            )
        fixed_fit = quadratic_prediction(meshes, fixed_scores)
        fixed_fit_all = quadratic_prediction(meshes, fixed_scores, start=0)
        optimized_fit = quadratic_prediction(meshes, optimized_scores)
        finest_matrix, finest_metric = (
            archives[-1][key("theta", theta)],
            archives[-1]["I"],
        )
        sensitivity = []
        for cutoff in (1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13):
            eig = parity.spectral_generalized_maximum(
                finest_matrix, finest_metric, cutoff
            )
            sensitivity.append(
                {
                    "cutoff": cutoff,
                    "retained_dimension": eig["retained_dimension"],
                    "score": float(RHO_STAR) * float(eig["quotient"]),
                }
            )
        analyzed.append(
            {
                "theta": theta,
                "theta_float": theta_row["theta_float"],
                "fixed_finest_vector_scores": dict(zip(map(str, meshes), fixed_scores)),
                "optimized_mesh_scores": dict(zip(map(str, meshes), optimized_scores)),
                "fixed_finest_vector_quadratic_last_three": fixed_fit,
                "fixed_finest_vector_quadratic_all_four": fixed_fit_all,
                "optimized_quadratic_last_three": optimized_fit,
                "finest_mesh_cutoff_sensitivity": sensitivity,
            }
        )

    fixed_crossing = interpolate_crossing(
        [
            {
                **row,
                "production_score": row["fixed_finest_vector_quadratic_last_three"][
                    "production"
                ],
            }
            for row in analyzed
        ],
        "production_score",
    )
    optimized_crossing = interpolate_crossing(
        [
            {
                **row,
                "production_score": row["optimized_quadratic_last_three"]["production"],
            }
            for row in analyzed
        ],
        "production_score",
    )
    payload = {
        "schema": "primegaps.physical-parity-modulus-reach-analysis.v1",
        "status": "exploratory-float64-mesh-extrapolation-not-a-certificate",
        "definition": results[-1]["definition"],
        "dimension": 39,
        "meshes": meshes,
        "production_intervals": PRODUCTION_INTERVALS,
        "rho_star": str(RHO_STAR),
        "rows": analyzed,
        "crossing": {
            "fixed_finest_vector": fixed_crossing,
            "separately_optimized_crosscheck": optimized_crossing,
        },
        "inputs": [
            {
                "result_path": str(result_path),
                "result_sha256": parity.file_hash(result_path),
                "matrix_path": str(matrix_path),
                "matrix_sha256": parity.file_hash(matrix_path),
            }
            for result_path, matrix_path in zip(
                args.results, args.matrices, strict=True
            )
        ],
    }
    parity.atomic_json(args.output, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
