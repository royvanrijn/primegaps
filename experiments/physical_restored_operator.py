#!/usr/bin/env python3
"""Build the signed 77-dimensional PrimeGaps186 cap-restoration operator.

The calculation retains the exact rational hybrid coefficients and optimizes

    rho_star * (J0 + (a+b) Jplus + b Jtail) / I.

Every source-cover loss in the complete certificate is a positive-semidefinite
quadratic form that is subtracted from this numerator.  The operator computed
here is therefore an optimistic spectral upper screen for the full restored
certificate at fixed geometry.  This is float64 mesh extrapolation, not a
rigorous eigenvalue enclosure.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

import numpy as np


F = Fraction
ROOT = Path(__file__).resolve().parents[1]
MODEL_SOURCE = ROOT / "experiments" / "physical_parity_viability.py"
DEFAULT_INPUT = ROOT / "reproduction" / "186" / "physical-parity-input.json"
DEFAULT_OUTPUT = Path(__file__).with_suffix(".json")
PRODUCTION_INTERVALS = 98_304
PUBLISHED_I_UPPER = F(23_685_317_890, 10**24)
PUBLISHED_J_LOWER = F(90_248_755_123, 10**24)


def load_model_module():
    spec = importlib.util.spec_from_file_location(
        "physical_parity_viability_for_restoration", MODEL_SOURCE
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {MODEL_SOURCE}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


model_module = load_model_module()


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


def quadratic_prediction(
    meshes: list[int], values: list[float], *, start: int
) -> dict[str, object]:
    x = 1 / np.asarray(meshes[start:], dtype=float)
    coefficients = np.polyfit(x, np.asarray(values[start:], dtype=float), 2)
    return {
        "fit_meshes": meshes[start:],
        "coefficients_descending": [float(value) for value in coefficients],
        "continuum": float(coefficients[-1]),
        "production": float(
            np.polyval(coefficients, 1 / PRODUCTION_INTERVALS)
        ),
    }


def restored_matrix(inputs, matrices: dict[str, np.ndarray]) -> np.ndarray:
    hybrid = inputs.DERIVED_INPUTS["hybrid"]
    a, b = F(hybrid["a"]), F(hybrid["b"])
    return (
        matrices["J0"]
        + float(a + b) * matrices["Jplus"]
        + float(b) * matrices["Jtail"]
    )


def build_result(
    *, dimensions: list[int], meshes: list[int], input_path: Path
) -> dict[str, object]:
    if meshes != sorted(set(meshes)) or len(meshes) < 4:
        raise ValueError("at least four distinct increasing meshes are required")
    if meshes[-1] > 8192:
        raise ValueError(
            "this float64 model is validated only through the direct-convolution 8192 mesh"
        )
    inputs = model_module.load_inputs(input_path)
    rho_star = float(inputs.RHO_STAR)
    stored: dict[int, list[dict[str, np.ndarray]]] = {
        dimension: [] for dimension in dimensions
    }
    optimized: dict[int, list[dict[str, object]]] = {
        dimension: [] for dimension in dimensions
    }
    published_vector = np.asarray(
        [
            F(value, 10**10)
            for row in inputs.CAP_COEFFICIENTS
            for value in row
        ],
        dtype=float,
    )

    for mesh in meshes:
        for dimension in dimensions:
            print(f"BUILD dimension={dimension} intervals={mesh}", flush=True)
            cap_model = model_module.FloatPhysicalCapModel(
                inputs, dimension=dimension, intervals=mesh
            )
            matrices = cap_model.matrices()
            matrices["restored"] = restored_matrix(inputs, matrices)
            run = model_module.spectral_generalized_maximum(
                matrices["restored"], matrices["I"], 1e-11
            )
            full = model_module.spectral_generalized_maximum(
                matrices["Jfull"], matrices["I"], 1e-11
            )
            optimized[dimension].append(
                {
                    "mesh": mesh,
                    "restored_score": rho_star * float(run["quotient"]),
                    "full_face_score": rho_star * float(full["quotient"]),
                    "retained_dimension": run["retained_dimension"],
                    "projected_relative_residual": run[
                        "projected_relative_residual"
                    ],
                    "published_vector_score": rho_star
                    * model_module.rayleigh(
                        matrices["restored"], matrices["I"], published_vector
                    ),
                }
            )
            matrices["restored_vector"] = np.asarray(run["vector"])
            stored[dimension].append(matrices)

    dimensions_result = []
    for dimension in dimensions:
        finest_vector = stored[dimension][-1]["restored_vector"]
        fixed_scores = []
        full_at_restored_vector = []
        for matrices in stored[dimension]:
            fixed_scores.append(
                rho_star
                * model_module.rayleigh(
                    matrices["restored"], matrices["I"], finest_vector
                )
            )
            full_at_restored_vector.append(
                rho_star
                * model_module.rayleigh(
                    matrices["Jfull"], matrices["I"], finest_vector
                )
            )
        separate_scores = [
            float(row["restored_score"]) for row in optimized[dimension]
        ]
        full_optimized = [
            float(row["full_face_score"]) for row in optimized[dimension]
        ]
        primary = quadratic_prediction(meshes, fixed_scores, start=1)
        separate = quadratic_prediction(meshes, separate_scores, start=1)
        all_mesh = quadratic_prediction(meshes, fixed_scores, start=0)
        full = quadratic_prediction(meshes, full_optimized, start=1)
        dimensions_result.append(
            {
                "dimension": dimension,
                "mesh_runs": optimized[dimension],
                "finest_restored_vector_scores": dict(
                    zip(map(str, meshes), fixed_scores, strict=True)
                ),
                "full_face_scores_at_finest_restored_vector": dict(
                    zip(map(str, meshes), full_at_restored_vector, strict=True)
                ),
                "restored_quadratic_inverse_mesh_fit": primary,
                "separately_optimized_fit_crosscheck": separate,
                "all_mesh_fit_crosscheck": all_mesh,
                "full_face_optimized_fit": full,
                "production_deficit_from_one": 1 - float(primary["production"]),
                "classification": (
                    "below_one_before_source_losses"
                    if primary["production"] < 1
                    else "above_one_before_source_losses"
                ),
            }
        )

    k40 = next((row for row in dimensions_result if row["dimension"] == 40), None)
    calibration = None
    if k40 is not None:
        published_lower = float(
            inputs.RHO_STAR * PUBLISHED_J_LOWER / PUBLISHED_I_UPPER
        )
        predicted = float(
            k40["separately_optimized_fit_crosscheck"]["production"]
        )
        calibration = {
            "published_rigorous_fixed_vector_lower_endpoint": published_lower,
            "predicted_optimized_production_score": predicted,
            "difference": predicted - published_lower,
            "note": (
                "The prediction is an optimized value and should be at least the "
                "published fixed-vector value up to mesh-extrapolation error."
            ),
        }

    hybrid = inputs.DERIVED_INPUTS["hybrid"]
    return {
        "schema": "primegaps.physical-restored-operator.v1",
        "status": "exploratory-float64-signed-operator-upper-screen",
        "operator": {
            "basis_dimension": 77,
            "formula": "J0 + (a+b) Jplus + b Jtail",
            "a": str(hybrid["a"]),
            "b": str(hybrid["b"]),
            "a_plus_b": str(F(hybrid["a"]) + F(hybrid["b"])),
            "rho_star": str(inputs.RHO_STAR),
            "source_loss_dominance": (
                "The complete fixed-geometry certificate subtracts 97 "
                "positive-semidefinite source-cover forms, so its largest "
                "generalized eigenvalue cannot exceed this signed cap value."
            ),
        },
        "dimensions": dimensions_result,
        "meshes": meshes,
        "production_intervals": PRODUCTION_INTERVALS,
        "calibration": calibration,
        "inputs": {
            "physical_input": str(input_path.resolve().relative_to(ROOT)),
            "physical_input_sha256": file_hash(input_path),
            "matrix_builder": str(MODEL_SOURCE.relative_to(ROOT)),
            "matrix_builder_sha256": file_hash(MODEL_SOURCE),
            "upstream": inputs.UPSTREAM,
        },
        "conclusion": (
            "At frozen geometry the exact signed hybrid coefficients erase the "
            "idealized full-face crossing for k=39 and k=38 before any source "
            "loss is charged. Trial-only optimization of the complete restored "
            "certificate is therefore a numerical no-go; geometry or hybrid/source "
            "parameters must move."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimensions", type=int, nargs="+", default=[40, 39, 38])
    parser.add_argument(
        "--meshes", type=int, nargs="+", default=[1024, 2048, 4096, 8192]
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build_result(
        dimensions=args.dimensions,
        meshes=args.meshes,
        input_path=args.input,
    )
    result["implementation_sha256"] = file_hash(Path(__file__))
    atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
