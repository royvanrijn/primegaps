#!/usr/bin/env python3
"""Optimize averaged candidate-space J with projected Davidson enrichment."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy import linalg


HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in (HERE, *HERE.parents) if (parent / "pyproject.toml").is_file())
SOURCE = ROOT / "scripts/build_projected_numerical_j.py"
spec = importlib.util.spec_from_file_location("projected_j_basis_builder", SOURCE)
projected_builder = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = projected_builder
spec.loader.exec_module(projected_builder)


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def metric_append(projection, direction, metric, tolerance=1e-12):
    candidate = np.asarray(direction, dtype=float).copy()
    gram = projection.T @ metric @ projection
    for _ in range(2):
        coefficients = linalg.solve(gram, projection.T @ metric @ candidate, assume_a="sym")
        candidate -= projection @ coefficients
    norm_squared = float(candidate @ metric @ candidate)
    reference = float(direction @ metric @ direction)
    if norm_squared <= tolerance * max(reference, 1.0):
        return projection, False, norm_squared
    candidate /= np.sqrt(norm_squared)
    return np.column_stack((projection, candidate)), True, norm_squared


def solve_projected(matrix_j, matrix_i, projection):
    projected_i = projection.T @ matrix_i @ projection
    projected_j = projection.T @ matrix_j @ projection
    values, vectors = linalg.eigh(
        (projected_j + projected_j.T) / 2,
        (projected_i + projected_i.T) / 2,
        subset_by_index=[projection.shape[1] - 1, projection.shape[1] - 1],
    )
    coefficient = projection @ vectors[:, 0]
    coefficient /= np.sqrt(float(coefficient @ matrix_i @ coefficient))
    value = float(coefficient @ matrix_j @ coefficient)
    residual = matrix_j @ coefficient - value * (matrix_i @ coefficient)
    relative_residual = float(
        np.linalg.norm(residual)
        / (np.linalg.norm(matrix_j @ coefficient) + abs(value) * np.linalg.norm(matrix_i @ coefficient))
    )
    return value, coefficient, residual, relative_residual


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-matrix", type=Path, required=True)
    parser.add_argument("--metric-matrix", type=Path, required=True)
    parser.add_argument("--j-matrices", type=Path, nargs="+", required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--k", type=int, default=48)
    parser.add_argument("--outer", default="2097/8000")
    parser.add_argument("--top", type=int, default=64)
    parser.add_argument("--spectral-cutoff", type=float, default=1e-13)
    parser.add_argument("--davidson-rounds", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    basis, _source_i, _source_j, projection, _candidate, diagnostics = (
        projected_builder.build_projection(
            args.source_matrix, args.degree, args.k, args.top,
            args.spectral_cutoff, args.candidate,
            projected_builder.Fraction(args.outer),
        )
    )
    with np.load(args.metric_matrix, allow_pickle=False) as payload:
        matrix_i = np.asarray(payload["m1"], dtype=float)
    matrices_j = []
    for path in args.j_matrices:
        with np.load(path, allow_pickle=False) as payload:
            matrices_j.append(np.asarray(payload["j"], dtype=float))
    if any(matrix.shape != matrix_i.shape for matrix in matrices_j):
        raise ValueError("I/J matrix dimensions differ")
    matrix_j = np.mean(matrices_j, axis=0)

    rounds = []
    coefficient = None
    for round_number in range(args.davidson_rounds + 1):
        value, coefficient, residual, relative_residual = solve_projected(
            matrix_j, matrix_i, projection
        )
        replicate_values = [
            float(coefficient @ replicate @ coefficient)
            for replicate in matrices_j
        ]
        record = {
            "round": round_number,
            "projected_dimension": projection.shape[1],
            "unrestricted_I_normalized_quotient": value,
            "relative_residual": relative_residual,
            "replicate_values": replicate_values,
            "replicate_standard_deviation": float(np.std(replicate_values, ddof=1))
            if len(replicate_values) > 1 else None,
        }
        rounds.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
        if round_number == args.davidson_rounds:
            break
        diagonal = np.diag(matrix_j) - value * np.diag(matrix_i)
        floor = max(np.max(np.abs(diagonal)) * 1e-8, np.finfo(float).tiny)
        denominator = np.where(
            np.abs(diagonal) >= floor,
            diagonal,
            np.where(diagonal < 0, -floor, floor),
        )
        direction = residual / denominator
        projection, added, norm_squared = metric_append(
            projection, direction, matrix_i
        )
        record["enrichment_added"] = added
        record["enrichment_metric_norm_squared"] = norm_squared
        if not added:
            break

    output_npz = args.output.with_suffix(".npz")
    output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_npz, projection=projection, coefficient=coefficient,
        mean_j=matrix_j,
    )
    payload = {
        "schema": "primegaps-projected-J-davidson-v1",
        "status": "numerical-discovery; I is unrestricted until exact candidate replay",
        "k": args.k, "degree": args.degree,
        "initial_top_unrestricted_directions": args.top,
        "initial_projection_diagnostics": diagnostics,
        "rounds": rounds,
        "output_npz": str(output_npz),
        "output_npz_sha256": file_hash(output_npz),
        "inputs": {
            "source_matrix": str(args.source_matrix),
            "source_matrix_sha256": file_hash(args.source_matrix),
            "metric_matrix": str(args.metric_matrix),
            "metric_matrix_sha256": file_hash(args.metric_matrix),
            "candidate": str(args.candidate),
            "candidate_sha256": file_hash(args.candidate),
            "j_matrices": [str(path) for path in args.j_matrices],
            "j_matrix_sha256": [file_hash(path) for path in args.j_matrices],
        },
    }
    atomic_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
