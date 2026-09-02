#!/usr/bin/env python3
"""Control-variate verification of fixed D=21 candidates at the B cutoffs."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy import linalg, special


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("qmc_verifier", HERE / "qmc_verifier.py")
q = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(q)


def candidates(k: int, degree: int, tolerances: list[float]):
    m1 = q.full_simplex_reference_gram(k, degree)
    j = q.unrestricted_marginal_gram(k, degree)
    scale = math.exp(2 * math.log(k) + (k - 1) * math.log(q.R) - k * math.log(q.U))
    m2 = scale * j
    diagonal_scale = 1.0 / np.sqrt(np.diag(m1))
    correlation = m1 * diagonal_scale[:, None] * diagonal_scale[None, :]
    eigenvalues, eigenvectors = linalg.eigh((correlation + correlation.T) * 0.5)
    raw_vectors = []
    records = []
    for tolerance in tolerances:
        keep = eigenvalues > tolerance * eigenvalues[-1]
        transform = (
            diagonal_scale[:, None]
            * eigenvectors[:, keep]
            / np.sqrt(eigenvalues[keep])[None, :]
        )
        reduced_m2 = transform.T @ m2 @ transform
        value, vector = linalg.eigh(
            (reduced_m2 + reduced_m2.T) * 0.5,
            subset_by_index=[reduced_m2.shape[0] - 1, reduced_m2.shape[0] - 1],
        )
        raw = transform @ vector[:, 0]
        denominator = float(raw @ m1 @ raw)
        numerator = float(raw @ m2 @ raw)
        raw_vectors.append(raw)
        records.append({
            "tolerance": tolerance,
            "retained_dimension": int(keep.sum()),
            "unrestricted_quotient": numerator / denominator,
            "unrestricted_denominator": denominator,
            "unrestricted_numerator": numerator,
        })
    return m1, m2, np.column_stack(raw_vectors), records


def correction_replicate(
    k: int,
    degree: int,
    vectors: np.ndarray,
    unrestricted_denominators: np.ndarray,
    unrestricted_numerators: np.ndarray,
    log2_n: int,
    seed: int,
    batch_log2: int,
):
    basis = q.basis_indices(degree)
    parts = q.all_partitions(degree // 2)
    coordinate_scale = k / (q.U * q.U)
    batch_n = 1 << min(batch_log2, log2_n)
    batches = 1 << max(0, log2_n - batch_log2)
    delta_denominator = np.zeros(vectors.shape[1])
    delta_numerator = np.zeros(vectors.shape[1])
    weighted_acceptance = 0.0
    proposal_acceptance = 0
    sample_n = 1 << log2_n
    scale = math.exp(2 * math.log(k) + (k - 1) * math.log(q.R) - k * math.log(q.U))
    component_deltas = []

    def stratum_shape(dimension: int, m: int) -> float:
        if m == 0:
            return 1.0
        target_big_mass = 0.68
        return target_big_mass * (dimension + 1 - m) / (m * (1.0 - target_big_mass))

    def fixed_component_weights(points: np.ndarray, radius: float, m: int, shape: float):
        if m == 0:
            return np.ones(points.shape[0])
        dimension = points.shape[1]
        y = np.maximum(points[:, :m] / radius, np.finfo(float).tiny)
        log_ratio = (
            special.gammaln(dimension + 1 + m * (shape - 1.0))
            - m * special.gammaln(shape)
            - special.gammaln(dimension + 1)
            + (shape - 1.0) * np.log(y).sum(axis=1)
        )
        return np.exp(-log_ratio)

    # M1 correction: partition the rejected set by the exact number m of big
    # coordinates.  Symmetry lets us tilt the first m and multiply by C(k,m).
    max_big_x = int(q.U // q.DELTA)
    for m in range(1, max_big_x + 1):
        shape = stratum_shape(k, m)
        component_seed = seed + 10_007 * m
        x = q.importance_simplex_points(k, q.U, log2_n, component_seed, (m, shape))
        wx = fixed_component_weights(x, q.U, m, shape) * math.comb(k, m)
        component_den = np.zeros(vectors.shape[1])
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            xb = x[sl]
            big = xb > q.DELTA
            exact_cell = big[:, :m].all(axis=1) & (~big[:, m:]).all(axis=1)
            accept = q.support_acceptance(xb)
            rejected = exact_cell & ~accept
            proposal_acceptance += int(rejected.sum())
            if np.any(rejected):
                f = q.features(
                    xb[rejected], basis, parts, degree, k, coordinate_scale
                ) @ vectors
                component_den -= np.sum(wx[sl][rejected, None] * f * f, axis=0)
        delta_denominator += component_den
        component_deltas.append({
            "kind": "M1",
            "m": m,
            "shape": shape,
            "delta_mean": (component_den / sample_n).tolist(),
        })

    # M2 correction: the same exact-count stratification in the k-1 marginal
    # variables.  The last-coordinate integral itself remains analytic Gauss.
    dimension_u = k - 1
    max_big_u = int(q.R // q.DELTA)
    for m in range(0, max_big_u + 1):
        shape = stratum_shape(dimension_u, m)
        component_seed = seed + 500_009 + 10_007 * m
        u = q.importance_simplex_points(
            dimension_u, q.R, log2_n, component_seed, (m, shape)
        )
        wu = fixed_component_weights(u, q.R, m, shape) * math.comb(dimension_u, m)
        component_num = np.zeros(vectors.shape[1])
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            ub = u[sl]
            big = ub > q.DELTA
            exact_cell = big[:, :m].all(axis=1) & (~big[:, m:]).all(axis=1)
            if not np.any(exact_cell):
                continue
            cell_points = ub[exact_cell]
            h_support = q.marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale
            ) @ vectors
            h_full = q.unrestricted_marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale
            ) @ vectors
            component_num += scale * np.sum(
                wu[sl][exact_cell, None] * (h_support * h_support - h_full * h_full), axis=0
            )
        delta_numerator += component_num
        component_deltas.append({
            "kind": "M2",
            "m": m,
            "shape": shape,
            "delta_mean": (component_num / sample_n).tolist(),
        })

    delta_denominator /= sample_n
    delta_numerator /= sample_n
    denominators = unrestricted_denominators + delta_denominator
    numerators = unrestricted_numerators + delta_numerator
    quotients = numerators / denominators
    return {
        "seed": seed,
        "log2_n": log2_n,
        "rejected_proposal_points": proposal_acceptance,
        "denominator_corrections": delta_denominator.tolist(),
        "numerator_corrections": delta_numerator.tolist(),
        "quotients": quotients.tolist(),
        "deficits_from_1": (1.0 - quotients).tolist(),
        "stratum_diagnostics": component_deltas,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--degree", type=int, default=21)
    parser.add_argument("--tolerances", default="1e-12,3e-13,1e-13")
    parser.add_argument("--log2-n", type=int, default=18)
    parser.add_argument("--seeds", default="49101,49102,49103,49104")
    parser.add_argument("--batch-log2", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    tolerances = [float(x) for x in args.tolerances.split(",")]
    m1, m2, vectors, records = candidates(args.k, args.degree, tolerances)
    den = np.array([x["unrestricted_denominator"] for x in records])
    num = np.array([x["unrestricted_numerator"] for x in records])
    replicates = [
        correction_replicate(
            args.k, args.degree, vectors, den, num, args.log2_n, int(seed), args.batch_log2
        )
        for seed in args.seeds.split(",") if seed
    ]
    quotient_array = np.array([x["quotients"] for x in replicates])
    result = {
        "method": "analytic unrestricted matrices plus importance-QMC cutoff corrections",
        "k": args.k,
        "degree": args.degree,
        "basis_size": len(m1),
        "parameters": {
            "epsilon": 0.0075,
            "delta": q.DELTA,
            "A": [-0.0075, 0.253],
            "B": [0.15, 0.15, "0.17 for m>=3"],
            "c1": 0,
            "c2": 0,
        },
        "candidates": records,
        "replicates": replicates,
        "quotient_means": quotient_array.mean(axis=0).tolist(),
        "quotient_standard_errors": (
            quotient_array.std(axis=0, ddof=1) / math.sqrt(len(replicates))
        ).tolist() if len(replicates) > 1 else None,
        "quotient_mins": quotient_array.min(axis=0).tolist(),
        "quotient_maxs": quotient_array.max(axis=0).tolist(),
        "elapsed_seconds": time.time() - started,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    np.savez_compressed(
        args.output.with_suffix(".npz"), m1=m1, m2_unrestricted=m2, vectors=vectors
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

