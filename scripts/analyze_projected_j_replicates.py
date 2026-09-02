#!/usr/bin/env python3
"""Cross-validate nested projected-J subspaces across independent replicates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np
from scipy import linalg


def solve(matrix_j, matrix_i):
    value, vector = linalg.eigh(
        (matrix_j + matrix_j.T) / 2,
        (matrix_i + matrix_i.T) / 2,
        subset_by_index=[len(matrix_i) - 1, len(matrix_i) - 1],
    )
    return float(value[0]), vector[:, 0]


def atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrices", type=Path, nargs="+", required=True)
    parser.add_argument("--prefixes", default="1,2,4,8,16,32,64")
    parser.add_argument("--certificate-column", type=int, default=-1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    matrices_i = []
    matrices_j = []
    projections = []
    for path in args.matrices:
        with np.load(path, allow_pickle=False) as payload:
            matrices_i.append(np.asarray(payload["matrix_i"], dtype=float))
            matrices_j.append(np.asarray(payload["matrix_j"], dtype=float))
            projections.append(np.asarray(payload["projection"], dtype=float))
    dimension = matrices_i[0].shape[0]
    certificate = args.certificate_column % dimension
    if any(matrix.shape != (dimension, dimension) for matrix in matrices_i + matrices_j):
        raise ValueError("projected matrix dimensions differ")
    if any(
        projection.shape != projections[0].shape
        or not np.allclose(projection, projections[0], rtol=0.0, atol=0.0)
        for projection in projections[1:]
    ):
        raise ValueError("replicates use different projection matrices")
    reference_i = np.mean(matrices_i, axis=0)
    results = []
    full_vectors = {}
    for prefix in (int(value) for value in args.prefixes.split(",")):
        columns = list(range(prefix))
        if certificate not in columns:
            columns.append(certificate)
        selection = np.ix_(columns, columns)
        current_i = reference_i[selection]
        current_j = [matrix[selection] for matrix in matrices_j]
        training_value, vector = solve(np.mean(current_j, axis=0), current_i)
        replicate_values = [
            float(vector @ matrix @ vector / (vector @ current_i @ vector))
            for matrix in current_j
        ]
        held_out = []
        for index, held_matrix in enumerate(current_j):
            training = np.mean(
                [matrix for other, matrix in enumerate(current_j) if other != index],
                axis=0,
            )
            _value, held_vector = solve(training, current_i)
            held_out.append(float(
                held_vector @ held_matrix @ held_vector
                / (held_vector @ current_i @ held_vector)
            ))
        results.append({
            "unrestricted_prefix": prefix,
            "dimension_with_certificate": len(columns),
            "training_value": training_value,
            "replicate_values": replicate_values,
            "replicate_standard_deviation": float(np.std(replicate_values, ddof=1)),
            "leave_one_out_values": held_out,
            "leave_one_out_mean": float(np.mean(held_out)),
            "leave_one_out_minimum": float(np.min(held_out)),
            "leave_one_out_standard_deviation": float(np.std(held_out, ddof=1)),
        })
        full_vectors[f"prefix_{prefix}"] = projections[0][:, columns] @ vector
    vector_output = args.output.with_suffix(".npz")
    np.savez_compressed(vector_output, **full_vectors)
    payload = {
        "schema": "primegaps-projected-J-replicate-cross-validation-v1",
        "matrix_paths": [str(path) for path in args.matrices],
        "replicate_count": len(args.matrices),
        "vector_output": str(vector_output),
        "results": results,
    }
    atomic_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
