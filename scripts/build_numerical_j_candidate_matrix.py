#!/usr/bin/env python3
"""Accumulate numerical J directly in candidate space, one SYRK per batch."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import tempfile
import time

import numpy as np

from primegaps.fast_exact.j_block import (
    MarginalMap,
    accumulate_candidate_gram,
    accumulate_gram_difference,
    candidate_feature_values,
)


HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in (HERE, *HERE.parents) if (parent / "pyproject.toml").is_file())
SOURCE = ROOT / "scripts/build_numerical_j_block_operator.py"
spec = importlib.util.spec_from_file_location("candidate_j_feature_builder", SOURCE)
builder = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = builder
spec.loader.exec_module(builder)
q = builder.qmc_verifier


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


def load_base(path: Path, dimension: int, *, include_i: bool):
    with np.load(path, allow_pickle=False) as payload:
        key = "m2_unrestricted" if "m2_unrestricted" in payload else "j"
        matrix_j = np.asarray(payload[key], dtype=float)
        matrix_i = np.asarray(payload["m1"], dtype=float) if include_i else None
    if matrix_j.shape != (dimension, dimension):
        raise ValueError("unrestricted base matrix has the wrong dimension")
    if include_i and matrix_i.shape != matrix_j.shape:
        raise ValueError("unrestricted I base matrix has the wrong dimension")
    return matrix_i, matrix_j


def support_acceptance(points, support):
    total = points.sum(axis=1)
    big = points > support["delta"]
    count = big.sum(axis=1)
    big_sum = np.where(big, points, 0.0).sum(axis=1)
    accepted = np.zeros(len(points), dtype=bool)
    for band_index, (lower, upper, caps) in enumerate(support["bands"]):
        in_band = (total >= lower) & (
            total <= upper if band_index + 1 == len(support["bands"]) else total < upper
        )
        current = count == 0
        active = (count > 0) & (count <= len(caps))
        indices = np.flatnonzero(active)
        cap_array = np.asarray(caps)
        current[indices] = big_sum[indices] <= cap_array[count[indices] - 1]
        accepted |= in_band & current
    return accepted


def build_candidate_matrix(
    *, k, degree, log2_n, seed, batch_log2, support,
    base_i=None, base_j=None, include_i=False,
):
    q.U, q.R, q.DELTA = support["U"], support["R"], support["delta"]
    basis = tuple(q.basis_indices(degree))
    parts = q.all_partitions(degree // 2)
    marginal_map = MarginalMap.from_basis(basis)
    coordinate_scale = k / (q.U * q.U)
    scale = math.exp(
        2 * math.log(k) + (k - 1) * math.log(q.R) - k * math.log(q.U)
    )
    correction_mode = base_j is not None
    if include_i and correction_mode != (base_i is not None):
        raise ValueError("I and J must use the same direct/correction mode")
    matrix_i = None if base_i is None else np.array(base_i, copy=True)
    matrix_j = None if base_j is None else np.array(base_j, copy=True)
    component_count = len(q.IMPORTANCE_COMPONENTS)
    component_bits = int(round(math.log2(component_count)))
    if 1 << component_bits != component_count or log2_n < component_bits:
        raise ValueError("importance component count must fit sample size")
    component_log2 = log2_n - component_bits
    batch_n = 1 << min(batch_log2, component_log2)
    batches = 1 << max(0, component_log2 - batch_log2)
    sample_n = 1 << log2_n
    started = time.perf_counter()
    completed = 0
    changed_rows = 0
    changed_i_rows = 0
    feature_seconds = 0.0
    assembly_seconds = 0.0
    syrk_seconds = 0.0
    for component_index, component in enumerate(q.IMPORTANCE_COMPONENTS):
        component_seed = seed + 10_007 * component_index + 1_000_003
        points = q.importance_simplex_points(
            k - 1, q.R, component_log2, component_seed, component
        )
        weights = q.importance_weights(points, q.R)
        if include_i:
            i_points = q.importance_simplex_points(
                k, q.U, component_log2, component_seed - 1_000_003, component
            )
            i_weights = q.importance_weights(i_points, q.U)
        for batch in range(batches):
            selection = slice(batch * batch_n, (batch + 1) * batch_n)
            if include_i:
                current_i = i_points[selection]
                current_i_weights = i_weights[selection] / sample_n
                all_i = q.features(
                    current_i, basis, parts,
                    degree, k, coordinate_scale,
                )
                accepted_i = support_acceptance(current_i, support)
                if correction_mode:
                    legal_i = np.array(all_i, copy=True)
                    legal_i[~accepted_i] = 0.0
                    matrix_i, changed_i = accumulate_gram_difference(
                        legal_i, all_i, current_i_weights, gram=matrix_i
                    )
                    changed_i_rows += changed_i
                elif np.any(accepted_i):
                    matrix_i = accumulate_candidate_gram(
                        all_i[accepted_i], current_i_weights[accepted_i],
                        gram=matrix_i,
                    )
                    changed_i_rows += int(accepted_i.sum())
            batch_points = points[selection]
            batch_weights = scale * weights[selection] / sample_n
            tick = time.perf_counter()
            legal_blocks = builder.evaluated_blocks(
                batch_points, marginal_map, degree, k, coordinate_scale,
                legal=True, support=support,
            )
            feature_seconds += time.perf_counter() - tick
            tick = time.perf_counter()
            legal = candidate_feature_values(marginal_map, legal_blocks)
            assembly_seconds += time.perf_counter() - tick
            tick = time.perf_counter()
            if not correction_mode:
                matrix_j = accumulate_candidate_gram(
                    legal, batch_weights, gram=matrix_j
                )
                changed_rows += len(legal)
            else:
                full_blocks = builder.evaluated_blocks(
                    batch_points, marginal_map, degree, k, coordinate_scale,
                    legal=False, support=support,
                )
                full = candidate_feature_values(marginal_map, full_blocks)
                matrix_j, changed = accumulate_gram_difference(
                    legal, full, batch_weights, gram=matrix_j
                )
                changed_rows += changed
            syrk_seconds += time.perf_counter() - tick
            completed += 1
            print(json.dumps({
                "completed_batches": completed,
                "total_batches": component_count * batches,
                "elapsed_seconds": time.perf_counter() - started,
                "changed_rows": changed_rows,
                "changed_i_rows": changed_i_rows,
            }), flush=True)
    timings = {
        "elapsed_seconds": time.perf_counter() - started,
        "feature_seconds": feature_seconds,
        "candidate_assembly_seconds": assembly_seconds,
        "gram_update_seconds": syrk_seconds,
    }
    counts = {
        "candidate_dimension": len(basis),
        "candidate_upper_triangle_entries": len(basis) * (len(basis) + 1) // 2,
        "marginal_feature_dimension": sum(map(len, marginal_map.feature_keys.values())),
        "signature_pair_upper_triangle_entries": sum(
            len(marginal_map.feature_keys[left]) * len(marginal_map.feature_keys[right])
            for left_index, left in enumerate(marginal_map.feature_keys)
            for right in tuple(marginal_map.feature_keys)[left_index:]
        ),
        "changed_rows": changed_rows,
        "changed_i_rows": changed_i_rows,
        "sample_count": sample_n,
    }
    return basis, matrix_i, matrix_j, timings, counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--k", type=int, default=48)
    parser.add_argument("--log2-n", type=int, default=15)
    parser.add_argument("--batch-log2", type=int, default=9)
    parser.add_argument("--seed", type=int, default=55101)
    parser.add_argument("--support-config", type=Path, required=True)
    parser.add_argument(
        "--unrestricted-base", type=Path,
        help="endpoint-compatible full unrestricted J; enables correction mode",
    )
    parser.add_argument(
        "--include-i", action="store_true",
        help="also accumulate the legal I matrix in candidate space",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    support = builder.load_support(args.support_config)
    dimension = len(q.basis_indices(args.degree))
    if args.unrestricted_base:
        base_i, base_j = load_base(
            args.unrestricted_base, dimension, include_i=args.include_i
        )
    else:
        base_i = base_j = None
    basis, matrix_i, matrix_j, timings, counts = build_candidate_matrix(
        k=args.k, degree=args.degree, log2_n=args.log2_n, seed=args.seed,
        batch_log2=args.batch_log2, support=support,
        base_i=base_i, base_j=base_j, include_i=args.include_i,
    )
    matrix_path = args.output.with_suffix(".npz")
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {"j": matrix_j}
    if matrix_i is not None:
        arrays["i"] = matrix_i
    np.savez_compressed(matrix_path, **arrays)
    metadata = {
        "schema": "primegaps-numerical-J-candidate-matrix-v1",
        "method": (
            "importance-QMC legal-minus-unrestricted candidate-space cross update"
            if base_j is not None
            else "importance-QMC direct legal candidate-space Gram update"
        ),
        "k": args.k, "degree": args.degree,
        "log2_n": args.log2_n, "seed": args.seed,
        "batch_log2": args.batch_log2,
        "includes_i": args.include_i,
        "counts": counts, "timings": timings,
        "matrix_output": str(matrix_path),
        "matrix_sha256": file_hash(matrix_path),
        "support_config": str(args.support_config),
        "support_config_sha256": file_hash(args.support_config),
        "unrestricted_base": str(args.unrestricted_base) if args.unrestricted_base else None,
        "unrestricted_base_sha256": file_hash(args.unrestricted_base) if args.unrestricted_base else None,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }
    atomic_json(args.output, metadata)
    print(json.dumps(metadata, sort_keys=True))


if __name__ == "__main__":
    main()
