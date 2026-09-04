#!/usr/bin/env python3
"""Extrapolate the parity-rough physical trial to the production mesh."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import tempfile

import numpy as np


RHO_STAR = Fraction(2_624_989, 10_000_000)
PRODUCTION_INTERVALS = 98_304
PUBLISHED_I_UPPER = Fraction(23_685_317_890, 10**24)
PUBLISHED_J_LOWER = Fraction(90_248_755_123, 10**24)


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def quadratic_prediction(meshes, values, *, start: int) -> tuple[float, float]:
    x = 1 / np.asarray(meshes[start:], dtype=float)
    coefficients = np.polyfit(x, np.asarray(values[start:], dtype=float), 2)
    return float(coefficients[-1]), float(np.polyval(coefficients, 1 / PRODUCTION_INTERVALS))


def rayleigh(matrix, metric, vector) -> float:
    return float(vector @ matrix @ vector) / float(vector @ metric @ vector)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, nargs="+", required=True)
    parser.add_argument("--matrices", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if len(args.results) != len(args.matrices):
        raise ValueError("one matrix archive is required per result")
    loaded = [json.loads(path.read_text()) for path in args.results]
    meshes = [int(row["dimensions"][0]["intervals"]) for row in loaded]
    if meshes != sorted(meshes) or len(meshes) < 4:
        raise ValueError("at least four increasingly fine meshes are required")
    archives = [np.load(path) for path in args.matrices]
    constants = loaded[-1]["model"]["rough_factor_constants"]
    dimensions = []
    for dimension in (39, 38):
        finest_vector = archives[-1][f"k{dimension}_selected_vector"]
        scores = []
        optimized_scores = []
        for row, archive in zip(loaded, archives):
            metric = archive[f"k{dimension}_I"]
            matrix = archive[f"k{dimension}_Jfull"]
            scores.append(float(RHO_STAR) * rayleigh(matrix, metric, finest_vector))
            record = next(item for item in row["dimensions"] if item["dimension"] == dimension)
            optimized_scores.append(record["full_face_optimization"]["score"])
        continuum, production = quadratic_prediction(meshes, scores, start=1)
        continuum_all, production_all = quadratic_prediction(meshes, scores, start=0)
        optimized_continuum, optimized_production = quadratic_prediction(
            meshes, optimized_scores, start=1
        )
        contributions = {
            "plus_omega": production * constants["omega_choose_1"],
            "minus_2_choose_2": -2 * production * constants["omega_choose_2"],
            "plus_3_choose_3": 3 * production * constants["omega_choose_3"],
        }
        gross = sum(abs(value) for value in contributions.values())
        eta = production - 1
        dimensions.append(
            {
                "dimension": dimension,
                "finest_mesh_vector_scores": dict(zip(map(str, meshes), scores)),
                "separately_optimized_mesh_scores": dict(zip(map(str, meshes), optimized_scores)),
                "quadratic_in_inverse_mesh": {
                    "fit_meshes": meshes[1:],
                    "continuum_score": continuum,
                    "production_mesh_score": production,
                    "all_mesh_production_score_crosscheck": production_all,
                    "all_mesh_continuum_score_crosscheck": continuum_all,
                    "separately_optimized_production_score_crosscheck": optimized_production,
                    "separately_optimized_continuum_score_crosscheck": optimized_continuum,
                },
                "factorial_contributions_normalized_by_I": contributions,
                "signed_sum": sum(contributions.values()),
                "gross_absolute_contribution": gross,
                "parity_error_budget": {
                    "eta_max_for_abs_error_le_eta_I": eta,
                    "common_relative_error_budget_across_unsigned_terms": eta / gross,
                },
            }
        )

    hybrid_scores = []
    for row in loaded:
        k40 = next(item for item in row["dimensions"] if item["dimension"] == 40)
        hybrid_scores.append(k40["fixed_published_k40_coefficients"]["hybrid_score"])
    hybrid_continuum, hybrid_production = quadratic_prediction(meshes, hybrid_scores, start=1)
    published_lower = float(RHO_STAR * PUBLISHED_J_LOWER / PUBLISHED_I_UPPER)
    payload = {
        "schema": "primegaps.physical-parity-production-extrapolation.v1",
        "status": "exploratory-float64-mesh-extrapolation-not-a-certificate",
        "production_intervals": PRODUCTION_INTERVALS,
        "beta": loaded[-1]["model"]["beta"],
        "rough_factor_constants": constants,
        "calibration": {
            "mesh_scores": dict(zip(map(str, meshes), hybrid_scores)),
            "continuum_score": hybrid_continuum,
            "predicted_production_score": hybrid_production,
            "published_rigorous_lower_endpoint_score": published_lower,
            "prediction_minus_published_lower": hybrid_production - published_lower,
        },
        "dimensions": dimensions,
        "inputs": [
            {
                "result": str(result),
                "result_sha256": file_hash(result),
                "matrices": str(matrix),
                "matrices_sha256": file_hash(matrix),
            }
            for result, matrix in zip(args.results, args.matrices)
        ],
    }
    atomic_json(args.output, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
