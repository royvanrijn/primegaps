#!/usr/bin/env python3
"""Build a stable numerical J block operator from evaluated Jacobi features."""

from __future__ import annotations

import argparse
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

from primegaps.fast_exact.j_block import (
    JBlockOperator,
    MarginalMap,
    accumulate_feature_gram_blocks,
    factorized_feature_values,
    save_block_operator,
)


HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in (HERE, *HERE.parents) if (parent / "pyproject.toml").is_file())
QMC_PATH = ROOT / "reproduction/240/numerical/qmc_verifier.py"
spec = importlib.util.spec_from_file_location("j_block_qmc_verifier", QMC_PATH)
qmc_verifier = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = qmc_verifier
spec.loader.exec_module(qmc_verifier)


def _as_float(value):
    """Parse JSON numbers or exact ``numerator/denominator`` strings."""
    if isinstance(value, str):
        return float(Fraction(value))
    return float(value)


def default_support():
    maximum_large = int(qmc_verifier.U // qmc_verifier.DELTA)
    caps = (0.15, 0.15) + (0.17,) * max(0, maximum_large - 2)
    return {
        "U": float(qmc_verifier.U),
        "R": float(qmc_verifier.R),
        "delta": float(qmc_verifier.DELTA),
        "bands": ((0.0, float(qmc_verifier.U), caps),),
    }


def load_support(path):
    """Load the total-sum bands and large-coordinate caps for legal J."""
    if path is None:
        return default_support()
    payload = json.loads(path.read_text())
    bands = tuple(
        (
            _as_float(row["lower"]),
            _as_float(row["upper"]),
            tuple(_as_float(value) for value in row["large_caps"]),
        )
        for row in payload["bands"]
    )
    support = {
        "U": _as_float(payload["U"]),
        "R": _as_float(payload["R"]),
        "delta": _as_float(payload["delta"]),
        "bands": bands,
    }
    if not bands or bands[0][0] != 0.0:
        raise ValueError("support bands must start at total sum zero")
    for index, (lower, upper, caps) in enumerate(bands):
        if lower < 0 or upper <= lower or not caps:
            raise ValueError("invalid support band")
        if index and not math.isclose(lower, bands[index - 1][1]):
            raise ValueError("support bands must be contiguous")
        if any(value <= 0 for value in caps):
            raise ValueError("large-coordinate caps must be positive")
    if not math.isclose(bands[-1][1], support["U"]):
        raise ValueError("last support band must end at U")
    if not (0 < support["R"] <= support["U"] and support["delta"] > 0):
        raise ValueError("invalid U, R, or delta")
    return support


def integrated_jacobi_moments(total_room, lo, hi, maximum_power, degree, k):
    """Gauss-integrate every required power/Jacobi product through D=30+."""
    point_count = total_room.size
    result = np.zeros((maximum_power + 1, degree + 1, point_count))
    half = np.maximum((hi - lo) * 0.5, 0.0)
    middle = (hi + lo) * 0.5
    order = max(24, degree + 1)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    exponents = 2 * np.arange(maximum_power + 1)
    for node, weight in zip(nodes, weights):
        coordinate = middle + half * node
        powers = coordinate[None, :] ** exponents[:, None]
        radial_coordinate = np.clip(
            (total_room - coordinate) / qmc_verifier.U,
            0.0,
            1.0,
        )
        jacobi = qmc_verifier.eval_jacobi_basis(degree, k, radial_coordinate)
        result += np.einsum(
            "rn,nb,n->rbn",
            powers,
            jacobi,
            half * weight,
            optimize=True,
        )
    return result


def radial_values(points, degree, k, *, legal, support):
    """Evaluate integrated last-coordinate Jacobi/power factors."""
    total = points.sum(axis=1)
    room = qmc_verifier.U - total
    maximum_power = degree // 2
    if not legal:
        return integrated_jacobi_moments(
            room,
            np.zeros_like(room),
            room,
            maximum_power,
            degree,
            k,
        )

    big_mask = points > support["delta"]
    count = big_mask.sum(axis=1)
    big_sum = np.where(big_mask, points, 0.0).sum(axis=1)
    result = np.zeros((maximum_power + 1, degree + 1, len(points)))
    for lower, upper, caps in support["bands"]:
        caps = np.asarray(caps)
        band_lo = lower - total
        band_hi = upper - total

        current_ok = count == 0
        active = (count > 0) & (count <= len(caps))
        indices = np.flatnonzero(active)
        current_ok[indices] = big_sum[indices] <= caps[count[indices] - 1]
        small_lo = np.maximum(0.0, band_lo)
        small_hi = np.where(
            current_ok,
            np.minimum.reduce((
                np.full(len(points), support["delta"]), room, band_hi,
            )),
            small_lo,
        )
        small_hi = np.maximum(small_hi, small_lo)

        new_count = count + 1
        big_lo = np.maximum(support["delta"], band_lo)
        big_hi = np.array(big_lo, copy=True)
        active = new_count <= len(caps)
        indices = np.flatnonzero(active)
        big_hi[indices] = np.minimum.reduce((
            room[indices],
            band_hi[indices],
            caps[new_count[indices] - 1] - big_sum[indices],
        ))
        big_hi = np.maximum(big_hi, big_lo)

        result += integrated_jacobi_moments(
            room, small_lo, small_hi, maximum_power, degree, k
        )
        result += integrated_jacobi_moments(
            room, big_lo, big_hi, maximum_power, degree, k
        )
    return result


def evaluated_blocks(
    points, marginal_map, degree, k, coordinate_scale, *, legal, support=None
):
    support = default_support() if support is None else support
    signature_values = qmc_verifier.monomial_symmetric_values(
        points,
        qmc_verifier.all_partitions(degree // 2),
        coordinate_scale,
    )
    return factorized_feature_values(
        marginal_map,
        signature_values,
        radial_values(points, degree, k, legal=legal, support=support),
        power_scale=coordinate_scale,
    )


def unrestricted_blocks(marginal_map, degree, k):
    """Stable exact-Gauss unrestricted J blocks in the evaluated basis."""
    dimension = k - 1
    coordinate_scale = k / (qmc_verifier.U * qmc_verifier.U)
    maximum_power = degree // 2
    flat_size = (maximum_power + 1) * (degree + 1)
    radial = {}
    for angular_degree in range(degree + 1):
        beta_dimension = dimension + 2 * angular_degree
        nodes, weights = qmc_verifier.special.roots_jacobi(
            degree + 2, 0, beta_dimension - 1
        )
        total = qmc_verifier.R * (1.0 + nodes) * 0.5
        room = qmc_verifier.U - total
        integrated = integrated_jacobi_moments(
            room,
            np.zeros_like(room),
            room,
            maximum_power,
            degree,
            k,
        )
        for power in range(maximum_power + 1):
            integrated[power] *= coordinate_scale**power
        flat = integrated.reshape(flat_size, -1)
        radial[angular_degree] = (
            flat * (weights / weights.sum())
        ) @ flat.T

    blocks = {}
    signatures = tuple(marginal_map.feature_keys)
    for left_index, left in enumerate(signatures):
        left_keys = marginal_map.feature_keys[left]
        for right in signatures[left_index:]:
            right_keys = marginal_map.feature_keys[right]
            angular = qmc_verifier.monomial_angular_coefficient(
                dimension,
                coordinate_scale,
                qmc_verifier.R,
                left,
                right,
            )
            radial_matrix = radial[sum(left) + sum(right)]
            block = np.empty((len(left_keys), len(right_keys)))
            for row, (left_power, left_degree) in enumerate(left_keys):
                left_flat = left_power * (degree + 1) + left_degree
                for column, (right_power, right_degree) in enumerate(right_keys):
                    right_flat = right_power * (degree + 1) + right_degree
                    block[row, column] = (
                        angular * radial_matrix[left_flat, right_flat]
                    )
            if left == right:
                block = (block + block.T) / 2
            blocks[(left, right)] = block
    return blocks


def build_operator(
    k, degree, log2_n, seed, batch_log2, *, control_variate, support=None
):
    support = default_support() if support is None else support
    qmc_verifier.U = support["U"]
    qmc_verifier.R = support["R"]
    qmc_verifier.DELTA = support["delta"]
    basis = qmc_verifier.basis_indices(degree)
    marginal_map = MarginalMap.from_basis(basis)
    coordinate_scale = k / (qmc_verifier.U * qmc_verifier.U)
    scale = math.exp(
        2 * math.log(k)
        + (k - 1) * math.log(qmc_verifier.R)
        - k * math.log(qmc_verifier.U)
    )
    blocks = (
        unrestricted_blocks(marginal_map, degree, k)
        if control_variate
        else None
    )
    if blocks is not None:
        for key in blocks:
            blocks[key] *= scale

    component_count = len(qmc_verifier.IMPORTANCE_COMPONENTS)
    component_bits = int(round(math.log2(component_count)))
    if 1 << component_bits != component_count or log2_n < component_bits:
        raise ValueError("importance component count must fit the requested sample")
    component_log2 = log2_n - component_bits
    batch_n = 1 << min(batch_log2, component_log2)
    batch_count = 1 << max(0, component_log2 - batch_log2)
    sample_n = 1 << log2_n
    started = time.perf_counter()
    completed_batches = 0
    for component_index, component in enumerate(
        qmc_verifier.IMPORTANCE_COMPONENTS
    ):
        component_seed = seed + 10_007 * component_index + 1_000_003
        points = qmc_verifier.importance_simplex_points(
            k - 1,
            qmc_verifier.R,
            component_log2,
            component_seed,
            component,
        )
        weights = qmc_verifier.importance_weights(points, qmc_verifier.R)
        for batch in range(batch_count):
            selection = slice(batch * batch_n, (batch + 1) * batch_n)
            batch_points = points[selection]
            batch_weights = scale * weights[selection] / sample_n
            legal = evaluated_blocks(
                batch_points,
                marginal_map,
                degree,
                k,
                coordinate_scale,
                legal=True,
                support=support,
            )
            blocks = accumulate_feature_gram_blocks(
                marginal_map,
                legal,
                batch_weights,
                blocks=blocks,
            )
            if control_variate:
                unrestricted = evaluated_blocks(
                    batch_points,
                    marginal_map,
                    degree,
                    k,
                    coordinate_scale,
                    legal=False,
                    support=support,
                )
                blocks = accumulate_feature_gram_blocks(
                    marginal_map,
                    unrestricted,
                    -batch_weights,
                    blocks=blocks,
                )
            completed_batches += 1
            print(json.dumps({
                "completed_batches": completed_batches,
                "total_batches": component_count * batch_count,
                "elapsed_seconds": time.perf_counter() - started,
            }), flush=True)
    metadata = {
        "builder": "direct-evaluated-Jacobi-features",
        "control_variate": control_variate,
        "degree": degree,
        "k": k,
        "log2_n": log2_n,
        "m2_scale": scale,
        "random_seed": seed,
        "sample_count": sample_n,
        "support": {
            "U": support["U"],
            "R": support["R"],
            "delta": support["delta"],
            "bands": [
                {
                    "lower": lower,
                    "upper": upper,
                    "large_caps": list(caps),
                }
                for lower, upper, caps in support["bands"]
            ],
        },
    }
    return JBlockOperator(marginal_map, blocks), metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=49)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--log2-n", type=int, default=15)
    parser.add_argument("--seed", type=int, default=49101)
    parser.add_argument("--batch-log2", type=int, default=9)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument(
        "--support-config",
        type=Path,
        help="JSON support bands; defaults to the historical one-band support",
    )
    parser.add_argument(
        "--direct",
        action="store_true",
        help="estimate all of J rather than only the cutoff correction",
    )
    args = parser.parse_args()
    if args.k < 2 or args.degree < 0 or args.log2_n < 3:
        raise ValueError("invalid k, degree, or sample size")
    support = load_support(args.support_config)
    operator, metadata = build_operator(
        args.k,
        args.degree,
        args.log2_n,
        args.seed,
        args.batch_log2,
        control_variate=not args.direct,
        support=support,
    )
    save_block_operator(args.output_directory, operator, metadata=metadata)
    print(json.dumps({
        "output_directory": str(args.output_directory),
        "blocks": len(operator.blocks),
        "basis_dimension": operator.shape[0],
        "metadata": metadata,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
