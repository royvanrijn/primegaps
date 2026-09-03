#!/usr/bin/env python3
"""Run the small localized sparse-PSD/sum-of-squares sieve screen.

This is an opt-in numerical experiment.  It builds I and kJ by scrambled Sobol
integration, uses the production distribution oracle for the cell-pair mask,
and solves both the rank-one clique problem and the full sparse SDP.  Run it in
an environment providing SciPy and CVXOPT (the ``sdp`` optional dependency).
"""

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

import numpy as np
from scipy import linalg

from primegaps.distribution import Minorant, RegionCell, constraint_failures, is_certified
from primegaps.sos import best_rank_one, maximal_cliques, solve_sparse_psd


ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / "pyproject.toml").is_file())
QMC_PATH = ROOT / "reproduction/240/numerical/qmc_verifier.py"
QMC_SPEC = importlib.util.spec_from_file_location("sparse_sos_qmc", QMC_PATH)
QMC = importlib.util.module_from_spec(QMC_SPEC)
assert QMC_SPEC.loader is not None
sys.modules[QMC_SPEC.name] = QMC
QMC_SPEC.loader.exec_module(QMC)

EPSILON = 0.0075
DELTA = 0.028
SUPPORT_MAX = 0.253
MINORANT = Minorant("0.38", "0.4", "0.4")

# (A_j, exact large-coordinate count, B_jm, local lower total-mass cutoff).
# The production oracle gives the induced cycle 0--1--2--3--0, with 0--2 and
# 1--3 forbidden.  Each local support is a subset of the named Xi cell.
CELL_DATA = (
    (0.253, 3, 0.156, 0.2485),
    (0.250, 3, 0.196, 0.2455),
    (0.249, 4, 0.229, 0.2445),
    (0.251, 3, 0.176, 0.2465),
)


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=path.name, delete=False) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def cells() -> tuple[RegionCell, ...]:
    return tuple(
        RegionCell(
            str(a), count, str(cap), str(DELTA), support_max=str(SUPPORT_MAX),
            label=f"A={a},m={count},B={cap}",
        )
        for a, count, cap, _ in CELL_DATA
    )


def oracle_mask() -> tuple[np.ndarray, list[list[dict[str, object]]]]:
    regions = cells()
    mask = np.empty((len(regions), len(regions)), dtype=bool)
    failures: list[list[dict[str, object]]] = []
    for left, first in enumerate(regions):
        row = []
        for right, second in enumerate(regions):
            mask[left, right] = bool(is_certified(first, second, MINORANT))
            row.append(
                {
                    "certified": bool(mask[left, right]),
                    "failures": [item.constraint_id for item in constraint_failures(first, second, MINORANT)],
                }
            )
        failures.append(row)
    return mask, failures


def shifted_stratum_weight(dimension: int, count: int, radius: float) -> float:
    residual = radius - count * DELTA
    return math.comb(dimension, count) * (residual / radius) ** dimension


def marginal_features(
    points: np.ndarray,
    cell: tuple[float, int, float, float],
    basis: list[tuple[tuple[int, ...], int]],
    partitions: list[tuple[int, ...]],
    degree: int,
    k: int,
    coordinate_scale: float,
    outer_radius: float,
) -> np.ndarray:
    a_upper, target_count, cap, lower = cell
    upper = a_upper + EPSILON
    total = points.sum(axis=1)
    large = points > DELTA
    count = large.sum(axis=1)
    large_sum = np.where(large, points, 0.0).sum(axis=1)
    maximum_power = max(map(sum, partitions), default=0)
    integrated = np.zeros((maximum_power + 1, degree + 1, len(points)))

    small_lo = np.maximum(0.0, lower - total)
    small_hi = np.minimum.reduce((np.full(len(points), DELTA), upper - total))
    small_hi = np.maximum(small_hi, small_lo)
    small_active = (count == target_count) & (large_sum <= cap)
    small_hi = np.where(small_active, small_hi, small_lo)
    integrated += QMC.integrated_jacobi_moments(
        outer_radius - total, small_lo, small_hi, maximum_power, degree, k
    )

    large_lo = np.maximum(DELTA, lower - total)
    large_hi = np.minimum(upper - total, cap - large_sum)
    large_hi = np.maximum(large_hi, large_lo)
    large_hi = np.where(count + 1 == target_count, large_hi, large_lo)
    integrated += QMC.integrated_jacobi_moments(
        outer_radius - total, large_lo, large_hi, maximum_power, degree, k
    )

    symmetric = QMC.monomial_symmetric_values(points, partitions, coordinate_scale)
    columns = []
    for partition, radial_degree in basis:
        value = symmetric[partition] * integrated[0, radial_degree]
        for power in set(partition):
            value += (
                coordinate_scale**power
                * symmetric[QMC.remove_one(partition, power)]
                * integrated[power, radial_degree]
            )
        columns.append(value)
    return np.column_stack(columns)


def build_full_forms(
    *, k: int, degree: int, log2_samples: int, seed: int
) -> tuple[np.ndarray, np.ndarray, int]:
    basis = QMC.basis_indices(degree)
    partitions = QMC.all_partitions(degree // 2)
    local_dimension = len(basis)
    dimension = len(CELL_DATA) * local_dimension
    outer_radius = max(a + EPSILON for a, _, _, _ in CELL_DATA)
    marginal_radius = max(a - EPSILON for a, _, _, _ in CELL_DATA)
    QMC.U = outer_radius
    coordinate_scale = k / outer_radius**2
    denominator = np.zeros((dimension, dimension))
    numerator = np.zeros((dimension, dimension))

    for count in sorted({item[1] for item in CELL_DATA}):
        residual = outer_radius - count * DELTA
        points = QMC.simplex_points(k, residual, log2_samples, seed + 10_007 * count)
        points[:, :count] += DELTA
        exact_cell = (points[:, count:] <= DELTA).all(axis=1)
        total = points.sum(axis=1)
        large_sum = points[:, :count].sum(axis=1)
        raw = QMC.features(points, basis, partitions, degree, k, coordinate_scale)
        localized = np.zeros((len(points), dimension))
        for cell_index, (a_upper, target_count, cap, lower) in enumerate(CELL_DATA):
            if target_count != count:
                continue
            active = (
                exact_cell
                & (total >= lower)
                & (total < a_upper + EPSILON)
                & (large_sum <= cap)
            )
            block = slice(cell_index * local_dimension, (cell_index + 1) * local_dimension)
            localized[active, block] = raw[active]
        denominator += (
            shifted_stratum_weight(k, count, outer_radius)
            * localized.T @ localized
            / len(points)
        )

    cutoffs = np.array([item[0] - EPSILON for item in CELL_DATA])
    common_counts = sorted(
        {item[1] for item in CELL_DATA} | {item[1] - 1 for item in CELL_DATA}
    )
    for count in common_counts:
        residual = marginal_radius - count * DELTA
        points = QMC.simplex_points(
            k - 1, residual, log2_samples, seed + 500_009 + 10_007 * count
        )
        points[:, :count] += DELTA
        exact_cell = (points[:, count:] <= DELTA).all(axis=1)
        total = points.sum(axis=1)
        marginals = []
        for cell in CELL_DATA:
            values = marginal_features(
                points, cell, basis, partitions, degree, k, coordinate_scale, outer_radius
            )
            values[~exact_cell] = 0.0
            marginals.append(values)
        weight = shifted_stratum_weight(k - 1, count, marginal_radius) / len(points)
        for left in range(len(CELL_DATA)):
            left_block = slice(left * local_dimension, (left + 1) * local_dimension)
            for right in range(left, len(CELL_DATA)):
                right_block = slice(right * local_dimension, (right + 1) * local_dimension)
                active = total <= max(cutoffs[left], cutoffs[right])
                value = weight * marginals[left][active].T @ marginals[right][active]
                numerator[left_block, right_block] += value
                if left != right:
                    numerator[right_block, left_block] += value.T

    volume_ratio = k * marginal_radius ** (k - 1) / outer_radius**k
    objective = k * volume_ratio * numerator
    return denominator, objective, local_dimension


def local_mode_compression(
    denominator: np.ndarray,
    objective: np.ndarray,
    local_dimension: int,
    modes: int,
    relative_cutoff: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict[str, object]]]:
    transform = np.zeros((len(denominator), len(CELL_DATA) * modes))
    diagnostics = []
    for cell_index in range(len(CELL_DATA)):
        block = slice(cell_index * local_dimension, (cell_index + 1) * local_dimension)
        values, vectors = linalg.eigh(denominator[block, block])
        keep = values > values[-1] * relative_cutoff
        whitening = vectors[:, keep] / np.sqrt(values[keep])[None, :]
        local_objective = whitening.T @ objective[block, block] @ whitening
        local_values, local_vectors = linalg.eigh((local_objective + local_objective.T) / 2)
        if len(local_values) < modes:
            raise np.linalg.LinAlgError("too few stable local modes")
        transform[block, cell_index * modes : (cell_index + 1) * modes] = (
            whitening @ local_vectors[:, -modes:]
        )
        diagnostics.append(
            {
                "stable_rank": int(np.sum(keep)),
                "largest_local_value": float(local_values[-1]),
                "smallest_retained_mode_value": float(local_values[-modes]),
            }
        )
    compressed_i = transform.T @ denominator @ transform
    compressed_objective = transform.T @ objective @ transform
    return (
        (compressed_i + compressed_i.T) / 2,
        (compressed_objective + compressed_objective.T) / 2,
        transform,
        diagnostics,
    )


def run_one(args: argparse.Namespace, k: int, mask: np.ndarray) -> dict[str, object]:
    raw_i, raw_objective, local_dimension = build_full_forms(
        k=k, degree=args.degree, log2_samples=args.log2_samples, seed=args.seed
    )
    matrix_i, matrix_objective, transform, local_diagnostics = local_mode_compression(
        raw_i, raw_objective, local_dimension, args.modes, args.relative_cutoff
    )
    groups = np.repeat(np.arange(len(CELL_DATA)), args.modes)
    rank_one = best_rank_one(matrix_i, matrix_objective, mask, groups)
    psd = solve_sparse_psd(
        matrix_i, matrix_objective, mask, groups,
        tolerance=args.solver_tolerance, factor_cutoff=args.factor_cutoff,
    )
    output_base = args.output_directory / f"k{k}-d{args.degree}"
    matrix_path = output_base.with_suffix(".npz")
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        matrix_path,
        matrix_i=matrix_i,
        matrix_objective=matrix_objective,
        local_transform=transform,
        rank_one_vector=rank_one.vector,
        psd_matrix=psd.matrix,
        sos_factor=psd.factor,
        group_of_component=groups,
        allowed_mask=mask,
    )
    return {
        "k": k,
        "degree": args.degree,
        "localized_component_count": int(len(matrix_i)),
        "support_region_count": len(CELL_DATA),
        "modes_per_region": args.modes,
        "raw_local_basis_dimension": local_dimension,
        "log2_samples_per_stratum": args.log2_samples,
        "seed": args.seed,
        "local_diagnostics": local_diagnostics,
        "maximal_cliques": [list(item) for item in maximal_cliques(mask)],
        "rank_one": {"value": rank_one.value, "clique": list(rank_one.clique)},
        "psd": {
            "value": psd.value,
            "rank": psd.rank,
            "normalization": psd.normalization,
            "forbidden_max_abs": psd.forbidden_max_abs,
            "minimum_eigenvalue": psd.minimum_eigenvalue,
            "relative_gap": psd.relative_gap,
            "status": psd.status,
        },
        "relative_advantage": psd.value / rank_one.value - 1.0,
        "exceeds_2_2_percent": psd.value > 1.022 * rank_one.value,
        "artifact": str(matrix_path),
        "artifact_sha256": file_hash(matrix_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, nargs="+", default=[47, 46])
    parser.add_argument("--degree", type=int, default=7)
    parser.add_argument("--modes", type=int, default=15)
    parser.add_argument("--log2-samples", type=int, default=16)
    parser.add_argument("--seed", type=int, default=260903)
    parser.add_argument("--relative-cutoff", type=float, default=1e-10)
    parser.add_argument("--solver-tolerance", type=float, default=1e-8)
    parser.add_argument("--factor-cutoff", type=float, default=1e-7)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    if not 7 <= args.degree <= 12:
        parser.error("degree must be between 7 and 12")
    if not 1 <= args.modes * len(CELL_DATA) <= 100:
        parser.error("the compressed experiment must contain 1--100 components")

    mask, oracle_diagnostics = oracle_mask()
    if not np.all(np.diag(mask)):
        raise RuntimeError("the chosen bank contains an uncertified diagonal")
    results = [run_one(args, k, mask) for k in args.k]
    payload = {
        "schema": "primegaps-sparse-sos-screen-v1",
        "method": "stratified scrambled-Sobol I/kJ, local D7--D12 compression, CVXOPT SDP",
        "claim_scope": "numerical support-cell mask screen; not coefficient-level cancellation or an exact sieve certificate",
        "cells": [
            {
                "a_upper": a,
                "large_count": count,
                "large_sum_bound": cap,
                "local_total_lower": lower,
                "local_total_upper": a + EPSILON,
            }
            for a, count, cap, lower in CELL_DATA
        ],
        "allowed_mask": mask.astype(int).tolist(),
        "oracle_diagnostics": oracle_diagnostics,
        "results": results,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "cvxopt": __import__("cvxopt").__version__,
            "platform": platform.platform(),
        },
        "inputs": {"qmc_source": str(QMC_PATH), "qmc_source_sha256": file_hash(QMC_PATH)},
    }
    atomic_json(args.summary, payload)
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
