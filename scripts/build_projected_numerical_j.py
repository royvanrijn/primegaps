#!/usr/bin/env python3
"""Build legal I/J directly in a small candidate projection.

The integration loop never constructs signature-pair blocks.  For J it forms
``GQ`` from evaluated marginal feature blocks and the sparse marginal map, then
performs one projected rank-k update per batch.  The calibrated production path
uses the symmetric Dirichlet importance mixture; translated exact-large-count
strata remain available as an explicitly experimental sampler.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
import math
from math import comb, sqrt
from pathlib import Path
import sys
import tempfile
import time

import numpy as np
from scipy import linalg

from primegaps.fast_exact.j_block import (
    MarginalMap,
    accumulate_candidate_gram,
    accumulate_gram_difference,
    projected_feature_values,
)


HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in (HERE, *HERE.parents) if (parent / "pyproject.toml").is_file())
BUILDER_PATH = ROOT / "scripts/build_numerical_j_block_operator.py"
spec = importlib.util.spec_from_file_location("projected_j_feature_builder", BUILDER_PATH)
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


def jacobi_q_coefficients(degree: int, beta: int) -> tuple[int, ...]:
    answer = [0] * (degree + 1)
    for middle in range(degree + 1):
        base = (
            comb(degree, middle)
            * comb(degree + beta, degree - middle)
            * (-1) ** (degree - middle)
        )
        for offset in range(middle + 1):
            answer[degree - middle + offset] += (
                base * comb(middle, offset) * (-1) ** offset
            )
    return tuple(answer)


def recover_candidate(candidate, basis, *, k: int, outer: Fraction):
    """Recover the normalized Jacobi-basis vector by exact triangular solve."""
    expanded = {
        (tuple(value // 2 for value in term["signature"]), term["slack_power"]):
        Fraction(int(term["numerator"]), int(term["denominator"]))
        for term in candidate["terms"]
    }
    if any(value % 2 for term in candidate["terms"] for value in term["signature"]):
        raise ValueError("candidate has a non-even angular exponent")
    coordinate_scale = Fraction(k) / outer**2
    coefficients = {}
    for partition in sorted({partition for partition, _degree in basis}):
        maximum_degree = int(candidate["degree"]) - 2 * sum(partition)
        angular_scale = coordinate_scale ** sum(partition)
        right = [
            expanded.get((partition, slack), Fraction())
            * outer**slack / angular_scale
            for slack in range(maximum_degree + 1)
        ]
        normalized = [Fraction() for _ in range(maximum_degree + 1)]
        jacobi = [
            jacobi_q_coefficients(degree, k - 1)
            for degree in range(maximum_degree + 1)
        ]
        for slack in range(maximum_degree, -1, -1):
            remainder = right[slack] - sum(
                normalized[degree] * jacobi[degree][slack]
                for degree in range(slack + 1, maximum_degree + 1)
            )
            normalized[slack] = remainder / jacobi[slack][slack]
        for degree, value in enumerate(normalized):
            coefficients[(partition, degree)] = (
                float(value) / sqrt((2 * degree + k) / k)
            )
    return np.asarray([coefficients.get(key, 0.0) for key in basis])


def conditioned_top_projection(matrix_i, matrix_j, count, cutoff):
    """Return the top generalized eigenspace after spectral whitening."""
    if count < 0:
        raise ValueError("top projection count must be nonnegative")
    if count == 0:
        return np.zeros((len(matrix_i), 0)), np.asarray([]), 0
    diagonal = np.sqrt(np.diag(matrix_i))
    correlation_i = matrix_i / diagonal[:, None] / diagonal[None, :]
    values_i, vectors_i = linalg.eigh((correlation_i + correlation_i.T) / 2)
    keep = values_i > cutoff * values_i[-1]
    whitening = (
        vectors_i[:, keep] / np.sqrt(values_i[keep])[None, :]
    ) / diagonal[:, None]
    reduced = whitening.T @ matrix_j @ whitening
    count = min(int(count), reduced.shape[0])
    values, vectors = linalg.eigh(
        (reduced + reduced.T) / 2,
        subset_by_index=[reduced.shape[0] - count, reduced.shape[0] - 1],
    )
    return whitening @ vectors[:, ::-1], values[::-1], int(keep.sum())


def append_metric_direction(projection, direction, metric):
    """Append a direction literally and report its distance from the span.

    Keeping the original coefficient column is important for calibration:
    reconstructing a high-degree certificate from a normalized residual loses
    meaningful digits even when the metric error is tiny.
    """
    residual = np.asarray(direction, dtype=float).copy()
    residual -= projection @ (projection.T @ metric @ residual)
    residual -= projection @ (projection.T @ metric @ residual)
    norm_squared = float(residual @ metric @ residual)
    if norm_squared <= 1e-24:
        return projection, False, math.sqrt(max(norm_squared, 0.0))
    return np.column_stack((projection, direction)), True, math.sqrt(norm_squared)


def build_projection(matrix_path, degree, k, top, cutoff, candidate_path, outer):
    with np.load(matrix_path, allow_pickle=False) as payload:
        matrix_i = np.asarray(payload["m1"], dtype=float)
        matrix_j = np.asarray(
            payload["m2_unrestricted"] if "m2_unrestricted" in payload else payload["j"],
            dtype=float,
        )
    basis = tuple(q.basis_indices(degree))
    if matrix_i.shape != (len(basis), len(basis)) or matrix_j.shape != matrix_i.shape:
        raise ValueError("source matrices do not match requested basis")
    projection, values, effective_rank = conditioned_top_projection(
        matrix_i, matrix_j, top, cutoff
    )
    candidate_added = False
    candidate_residual_norm = None
    embedded = None
    if candidate_path is not None:
        candidate = json.loads(candidate_path.read_text())
        source_basis = tuple(q.basis_indices(int(candidate["degree"])))
        source = recover_candidate(
            candidate, source_basis, k=k, outer=Fraction(outer)
        )
        position = {key: index for index, key in enumerate(basis)}
        embedded = np.zeros(len(basis))
        for key, value in zip(source_basis, source):
            if key not in position:
                raise ValueError("candidate degree exceeds projection degree")
            embedded[position[key]] = value
        projection, candidate_added, candidate_residual_norm = append_metric_direction(
            projection, embedded, matrix_i
        )
    metric_gram = projection.T @ matrix_i @ projection
    top_dimension = min(top, projection.shape[1])
    diagnostics = {
        "effective_i_rank": effective_rank,
        "source_top_eigenvalues": values.tolist(),
        "candidate_added": candidate_added,
        "candidate_metric_residual_norm": candidate_residual_norm,
        "metric_gram_condition": float(np.linalg.cond(metric_gram)),
        "top_space_metric_orthogonality_max_abs": float(
            np.max(np.abs(metric_gram[:top_dimension, :top_dimension] - np.eye(top_dimension)))
        ) if top_dimension else 0.0,
    }
    return basis, matrix_i, matrix_j, projection, embedded, diagnostics


def load_forms(path, expected_dimension):
    with np.load(path, allow_pickle=False) as payload:
        matrix_i = np.asarray(payload["m1"], dtype=float)
        matrix_j = np.asarray(
            payload["m2_unrestricted"] if "m2_unrestricted" in payload else payload["j"],
            dtype=float,
        )
    expected = (expected_dimension, expected_dimension)
    if matrix_i.shape != expected or matrix_j.shape != expected:
        raise ValueError(f"base matrices have shape {matrix_i.shape}/{matrix_j.shape}, expected {expected}")
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


def shifted_volume_weight(dimension, large_count, radius, delta):
    residual = radius - large_count * delta
    if residual < 0:
        return 0.0
    return math.comb(dimension, large_count) * (residual / radius) ** dimension


def projected_stratified_forms(
    *, k, degree, projection, support, log2_n, seed, batch_log2,
    base_i=None, base_j=None,
):
    """Estimate direct legal forms or corrections to compatible base forms."""
    q.U, q.R, q.DELTA = support["U"], support["R"], support["delta"]
    basis = tuple(q.basis_indices(degree))
    marginal_map = MarginalMap.from_basis(basis)
    mapped_projection = marginal_map.forward_matrix(projection)
    coordinate_scale = k / (q.U * q.U)
    parts = q.all_partitions(degree // 2)
    sample_n = 1 << log2_n
    batch_n = 1 << min(batch_log2, log2_n)
    batches = 1 << max(0, log2_n - batch_log2)
    correction_mode = base_i is not None or base_j is not None
    if correction_mode != (base_i is not None and base_j is not None):
        raise ValueError("both unrestricted base forms are required for correction mode")
    matrix_i = np.array(base_i, copy=True) if correction_mode else None
    matrix_j = np.array(base_j, copy=True) if correction_mode else None
    strata = []
    started = time.perf_counter()

    first_m = 1 if correction_mode else 0
    for m in range(first_m, int(q.U // q.DELTA) + 1):
        residual = q.U - m * q.DELTA
        points = q.simplex_points(k, residual, log2_n, seed + 10_007 * m)
        if m:
            points[:, :m] += q.DELTA
        volume = shifted_volume_weight(k, m, q.U, q.DELTA)
        selected = 0
        for batch in range(batches):
            selection = slice(batch * batch_n, (batch + 1) * batch_n)
            current = points[selection]
            exact_cell = (current[:, m:] <= q.DELTA).all(axis=1)
            accepted = support_acceptance(current, support)
            active = exact_cell & ((~accepted) if correction_mode else accepted)
            selected += int(active.sum())
            if not np.any(active):
                continue
            values = q.features(
                current[active], basis, parts, degree, k, coordinate_scale
            ) @ projection
            weight = volume / sample_n * (-1.0 if correction_mode else 1.0)
            matrix_i = accumulate_candidate_gram(
                values, np.full(len(values), weight), gram=matrix_i
            )
        strata.append({
            "form": "I", "m": m, "selected": selected,
            "volume_weight": volume,
        })
        print(json.dumps({
            "form": "I", "m": m, "selected": selected,
            "elapsed_seconds": time.perf_counter() - started,
        }), flush=True)

    scale = math.exp(2 * math.log(k) + (k - 1) * math.log(q.R) - k * math.log(q.U))
    for m in range(0, int(q.R // q.DELTA) + 1):
        residual = q.R - m * q.DELTA
        points = q.simplex_points(k - 1, residual, log2_n, seed + 500_009 + 10_007 * m)
        if m:
            points[:, :m] += q.DELTA
        volume = shifted_volume_weight(k - 1, m, q.R, q.DELTA)
        exact_count = 0
        changed_rows = 0
        for batch in range(batches):
            selection = slice(batch * batch_n, (batch + 1) * batch_n)
            current = points[selection]
            exact_cell = (current[:, m:] <= q.DELTA).all(axis=1)
            exact_count += int(exact_cell.sum())
            if not np.any(exact_cell):
                continue
            cell = current[exact_cell]
            legal_blocks = builder.evaluated_blocks(
                cell, marginal_map, degree, k, coordinate_scale,
                legal=True, support=support,
            )
            legal = projected_feature_values(
                marginal_map, legal_blocks, projection,
                mapped_projection=mapped_projection,
            )
            weights = np.full(len(cell), scale * volume / sample_n)
            if correction_mode:
                full_blocks = builder.evaluated_blocks(
                    cell, marginal_map, degree, k, coordinate_scale,
                    legal=False, support=support,
                )
                full = projected_feature_values(
                    marginal_map, full_blocks, projection,
                    mapped_projection=mapped_projection,
                )
                matrix_j, changed = accumulate_gram_difference(
                    legal, full, weights, gram=matrix_j
                )
                changed_rows += changed
            else:
                matrix_j = accumulate_candidate_gram(legal, weights, gram=matrix_j)
                changed_rows += len(cell)
        strata.append({
            "form": "J", "m": m, "exact_count": exact_count,
            "changed_rows": changed_rows, "volume_weight": volume,
        })
        print(json.dumps({
            "form": "J", "m": m, "exact_count": exact_count,
            "changed_rows": changed_rows,
            "elapsed_seconds": time.perf_counter() - started,
        }), flush=True)
    return matrix_i, matrix_j, strata, time.perf_counter() - started


def projected_importance_correction(
    *, k, degree, projection, support, log2_n, seed, batch_log2,
    base_i, base_j, control_u=None, control_r=None,
):
    """Stream the calibrated importance correction directly through ``GQ``."""
    q.U, q.R, q.DELTA = support["U"], support["R"], support["delta"]
    basis = tuple(q.basis_indices(degree))
    parts = q.all_partitions(degree // 2)
    marginal_map = MarginalMap.from_basis(basis)
    mapped_projection = marginal_map.forward_matrix(projection)
    coordinate_scale = k / (q.U * q.U)
    safe_scale = math.exp(
        2 * math.log(k) + (k - 1) * math.log(q.R) - k * math.log(q.U)
    )
    coupled_endpoint = control_u is not None or control_r is not None
    if coupled_endpoint != (control_u is not None and control_r is not None):
        raise ValueError("control U and R must be supplied together")
    control_u = q.U if control_u is None else float(control_u)
    control_r = q.R if control_r is None else float(control_r)
    control_scale = math.exp(
        2 * math.log(k) + (k - 1) * math.log(control_r) - k * math.log(control_u)
    )
    component_count = len(q.IMPORTANCE_COMPONENTS)
    component_bits = int(round(math.log2(component_count)))
    if 1 << component_bits != component_count or log2_n < component_bits:
        raise ValueError("importance component count must fit sample size")
    component_log2 = log2_n - component_bits
    batch_n = 1 << min(batch_log2, component_log2)
    batches = 1 << max(0, component_log2 - batch_log2)
    sample_n = 1 << log2_n
    matrix_i = np.array(base_i, copy=True)
    matrix_j = np.array(base_j, copy=True)
    changed_rows = 0
    changed_i_rows = 0
    started = time.perf_counter()
    completed = 0
    for component_index, component in enumerate(q.IMPORTANCE_COMPONENTS):
        component_seed = seed + 10_007 * component_index + 1_000_003
        points = q.importance_simplex_points(
            k - 1, q.R, component_log2, component_seed, component
        )
        weights = q.importance_weights(points, q.R) / sample_n
        i_points = q.importance_simplex_points(
            k, q.U, component_log2, component_seed - 1_000_003, component
        )
        i_weights = q.importance_weights(i_points, q.U) / sample_n
        for batch in range(batches):
            selection = slice(batch * batch_n, (batch + 1) * batch_n)
            current_i = i_points[selection]
            full_i_points = (
                current_i * (control_u / support["U"])
                if coupled_endpoint else current_i
            )
            try:
                if coupled_endpoint:
                    q.U = control_u
                full_i = q.features(
                    full_i_points, basis, parts, degree, k,
                    k / (control_u * control_u),
                ) @ projection
            finally:
                q.U = support["U"]
            legal_i = q.features(
                current_i, basis, parts, degree, k, coordinate_scale
            ) @ projection
            legal_i[~support_acceptance(current_i, support)] = 0.0
            matrix_i, changed_i = accumulate_gram_difference(
                legal_i, full_i, i_weights[selection], gram=matrix_i
            )
            changed_i_rows += changed_i
            current = points[selection]
            legal_blocks = builder.evaluated_blocks(
                current, marginal_map, degree, k, coordinate_scale,
                legal=True, support=support,
            )
            legal = projected_feature_values(
                marginal_map, legal_blocks, projection,
                mapped_projection=mapped_projection,
            )
            if coupled_endpoint:
                control_points = current * (control_r / support["R"])
                try:
                    q.U, q.R = control_u, control_r
                    control_coordinate_scale = k / (control_u * control_u)
                    full_blocks = builder.evaluated_blocks(
                        control_points, marginal_map, degree, k,
                        control_coordinate_scale, legal=False, support=support,
                    )
                finally:
                    q.U, q.R = support["U"], support["R"]
            else:
                full_blocks = builder.evaluated_blocks(
                    current, marginal_map, degree, k, coordinate_scale,
                    legal=False, support=support,
                )
            full = projected_feature_values(
                marginal_map, full_blocks, projection,
                mapped_projection=mapped_projection,
            )
            matrix_j, changed = accumulate_gram_difference(
                math.sqrt(safe_scale) * legal,
                math.sqrt(control_scale) * full,
                weights[selection], gram=matrix_j
            )
            changed_rows += changed
            completed += 1
            print(json.dumps({
                "completed_batches": completed,
                "total_batches": component_count * batches,
                "changed_rows": changed_rows,
                "changed_i_rows": changed_i_rows,
                "elapsed_seconds": time.perf_counter() - started,
            }), flush=True)
    common = {
        "sampler": "symmetric Dirichlet importance mixture",
        "coupled_control_endpoint": coupled_endpoint,
        "control_U": control_u,
        "control_R": control_r,
    }
    return matrix_i, matrix_j, [
        {"form": "I", **common, "changed_rows": changed_i_rows},
        {"form": "J", **common, "changed_rows": changed_rows},
    ], time.perf_counter() - started


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-matrix", type=Path, required=True)
    parser.add_argument(
        "--base-matrix", type=Path,
        help="endpoint-compatible unrestricted I/J; defaults to source matrix",
    )
    parser.add_argument("--support-config", type=Path, required=True)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--k", type=int, default=48)
    parser.add_argument("--top", type=int, default=64)
    parser.add_argument("--spectral-cutoff", type=float, default=1e-13)
    parser.add_argument("--log2-n", type=int, default=10)
    parser.add_argument("--batch-log2", type=int, default=8)
    parser.add_argument("--seed", type=int, default=54101)
    parser.add_argument("--correction", action="store_true")
    parser.add_argument(
        "--sampler", choices=("translated", "importance"),
        default="importance",
    )
    parser.add_argument(
        "--control-u", type=builder._as_float,
        help="U of a different unrestricted source-matrix endpoint",
    )
    parser.add_argument(
        "--control-r", type=builder._as_float,
        help="R of a different unrestricted source-matrix endpoint",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    support = builder.load_support(args.support_config)
    q.U, q.R, q.DELTA = support["U"], support["R"], support["delta"]
    basis, source_i, source_j, projection, candidate_vector, projection_diagnostics = build_projection(
        args.source_matrix, args.degree, args.k, args.top,
        args.spectral_cutoff, args.candidate, Fraction(str(support["U"])),
    )
    base_path = args.base_matrix or args.source_matrix
    if args.correction:
        unrestricted_i, unrestricted_j = load_forms(base_path, len(basis))
        base_i = projection.T @ unrestricted_i @ projection
        base_j = projection.T @ unrestricted_j @ projection
    else:
        base_i = base_j = None
    if args.sampler == "importance":
        if not args.correction:
            raise ValueError("importance projected mode requires --correction")
        matrix_i, matrix_j, strata, elapsed = projected_importance_correction(
            k=args.k, degree=args.degree, projection=projection, support=support,
            log2_n=args.log2_n, seed=args.seed, batch_log2=args.batch_log2,
            base_i=base_i, base_j=base_j,
            control_u=args.control_u, control_r=args.control_r,
        )
    else:
        matrix_i, matrix_j, strata, elapsed = projected_stratified_forms(
            k=args.k, degree=args.degree, projection=projection, support=support,
            log2_n=args.log2_n, seed=args.seed, batch_log2=args.batch_log2,
            base_i=base_i, base_j=base_j,
        )
    values, vectors = linalg.eigh(
        (matrix_j + matrix_j.T) / 2,
        (matrix_i + matrix_i.T) / 2,
        subset_by_index=[projection.shape[1] - 1, projection.shape[1] - 1],
    )
    fixed_candidate = None
    if candidate_vector is not None:
        candidate_coordinates = np.zeros(projection.shape[1])
        candidate_coordinates[-1] = 1.0
        fixed_i = float(candidate_coordinates @ matrix_i @ candidate_coordinates)
        fixed_j = float(candidate_coordinates @ matrix_j @ candidate_coordinates)
        fixed_candidate = {
            "projected_I": fixed_i,
            "projected_kJ": fixed_j,
            "projected_quotient": fixed_j / fixed_i,
            "source_span_reconstruction_relative_error": 0.0,
        }
    output_npz = args.output.with_suffix(".npz")
    output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_npz, projection=projection, matrix_i=matrix_i, matrix_j=matrix_j,
        ritz_vector=vectors[:, 0], candidate_vector=projection @ vectors[:, 0],
    )
    payload = {
        "schema": "primegaps-projected-legal-forms-v1",
        "status": (
            "numerical discovery; cross-validate optimized directions and use exact I"
        ),
        "method": (
            "projected importance legal-minus-unrestricted correction"
            if args.sampler == "importance" else (
                "translated exact-m legal-minus-unrestricted correction"
                if args.correction else "translated exact-m direct legal integration"
            )
        ),
        "k": args.k, "degree": args.degree,
        "candidate_dimension": len(basis),
        "projected_dimension": projection.shape[1],
        "top_unrestricted_directions": args.top,
        "spectral_cutoff": args.spectral_cutoff,
        "log2_n_per_stratum": args.log2_n,
        "seed": args.seed,
        "control_endpoint": {
            "U": args.control_u, "R": args.control_r,
        } if args.control_u is not None else None,
        "largest_projected_quotient": float(values[0]),
        "fixed_candidate": fixed_candidate,
        "elapsed_seconds": elapsed,
        "projection": projection_diagnostics,
        "strata": strata,
        "inputs": {
            "source_matrix": str(args.source_matrix),
            "source_matrix_sha256": file_hash(args.source_matrix),
            "base_matrix": str(base_path) if args.correction else None,
            "base_matrix_sha256": file_hash(base_path) if args.correction else None,
            "support_config": str(args.support_config),
            "support_config_sha256": file_hash(args.support_config),
            "candidate": str(args.candidate) if args.candidate else None,
            "candidate_sha256": file_hash(args.candidate) if args.candidate else None,
        },
        "matrix_output": str(output_npz),
    }
    atomic_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
