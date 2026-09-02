#!/usr/bin/env python3
"""Common-random-number D=21 screening for one-band support geometries.

This is a numerical screening layer, not an exact certificate.  It keeps the
published Harman/distribution parameters fixed, rejects candidates not proved
by the repository's distribution oracle, and evaluates a bank of already
computed D=21 coefficient vectors using analytic unrestricted matrices plus
importance-QMC corrections at the delta/B boundaries.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from primegaps.distribution import Minorant, cells_from_support, is_certified
from primegaps.support import SupportParameters


HERE = Path(__file__).resolve().parent
QMC_PATH = HERE / "qmc_verifier.py"
spec = importlib.util.spec_from_file_location("qmc_verifier", QMC_PATH)
q = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(q)

PUBLISHED_MINORANT = Minorant("0.38", "0.4", "0.4")


def support_from_geometry(geometry: dict[str, float]) -> SupportParameters:
    delta = float(geometry["delta"])
    outer = float(geometry["outer"])
    marginal = float(geometry["marginal"])
    epsilon = (outer - marginal) / 2.0
    width = math.floor(1.0 / delta)
    if "bands" in geometry:
        bands = geometry["bands"]
        A = [-epsilon]
        rows = []
        for band in bands:
            A.append(float(band["a_upper"]))
            b12 = float(band["b12"])
            b3 = float(band["b3plus"])
            rows.append((b12, b12) + (b3,) * max(0, width - 2))
        expected_a1 = (outer + marginal) / 2.0
        if not math.isclose(A[-1], expected_a1, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("last band endpoint disagrees with outer/marginal")
        support = SupportParameters(delta, epsilon, tuple(A), tuple(rows))
    else:
        a1 = (outer + marginal) / 2.0
        b12 = float(geometry["b12"])
        b3 = float(geometry["b3plus"])
        row = (b12, b12) + (b3,) * max(0, width - 2)
        support = SupportParameters(delta, epsilon, (-epsilon, a1), (row,))
    support.validate()
    return support


def analytic_certificate(geometry: dict[str, float]) -> tuple[bool, str | None]:
    try:
        support = support_from_geometry(geometry)
    except ValueError as exc:
        return False, f"geometry: {exc}"
    cells = cells_from_support(support)
    for left in cells:
        for right in cells:
            certificate = is_certified(left, right, PUBLISHED_MINORANT)
            if not certificate:
                return False, f"{left.label} x {right.label}: {certificate.reason}"
    return True, None


def support_acceptance(points: np.ndarray, support: SupportParameters) -> np.ndarray:
    big = points > support.delta
    count = big.sum(axis=1)
    big_sum = np.where(big, points, 0.0).sum(axis=1)
    total = points.sum(axis=1)
    result = np.zeros(points.shape[0], dtype=bool)
    for band_index, row_values in enumerate(support.B):
        row = np.asarray(row_values)
        lower = support.A[band_index] + support.epsilon
        upper = support.A[band_index + 1] + support.epsilon
        in_band = (total >= lower) & (total <= upper if band_index + 1 == len(support.B) else total < upper)
        accepted = count == 0
        active = (count > 0) & (count <= len(row))
        indices = np.flatnonzero(active)
        accepted[indices] = big_sum[indices] <= row[count[indices] - 1]
        result |= in_band & accepted
    return result


def marginal_features(
    points: np.ndarray,
    basis,
    parts,
    degree: int,
    k: int,
    coordinate_scale: float,
    support: SupportParameters,
) -> np.ndarray:
    n = points.shape[0]
    total = points.sum(axis=1)
    room = q.U - total
    big_mask = points > support.delta
    count = big_mask.sum(axis=1)
    big_sum = np.where(big_mask, points, 0.0).sum(axis=1)
    max_r = max(map(sum, parts), default=0)
    integrated = np.zeros((max_r + 1, degree + 1, n))
    for band_index, row_values in enumerate(support.B):
        row = np.asarray(row_values)
        band_lo = support.A[band_index] + support.epsilon - total
        band_hi = support.A[band_index + 1] + support.epsilon - total

        current_ok = count == 0
        active = (count > 0) & (count <= len(row))
        indices = np.flatnonzero(active)
        current_ok[indices] = big_sum[indices] <= row[count[indices] - 1]
        small_lo = np.maximum(0.0, band_lo)
        small_hi = np.where(
            current_ok,
            np.minimum.reduce((np.full(n, support.delta), room, band_hi)),
            small_lo,
        )
        small_hi = np.maximum(small_hi, small_lo)

        new_count = count + 1
        big_lo = np.maximum(support.delta, band_lo)
        big_hi = np.array(big_lo, copy=True)
        active = new_count <= len(row)
        indices = np.flatnonzero(active)
        big_hi[indices] = np.minimum.reduce((
            room[indices],
            band_hi[indices],
            row[new_count[indices] - 1] - big_sum[indices],
        ))
        big_hi = np.maximum(big_hi, big_lo)

        integrated += q.integrated_jacobi_moments(
            room, small_lo, small_hi, max_r, degree, k
        )
        integrated += q.integrated_jacobi_moments(
            room, big_lo, big_hi, max_r, degree, k
        )
    mvals = q.monomial_symmetric_values(points, parts, coordinate_scale)
    columns = []
    for lam, b in basis:
        value = mvals[lam] * integrated[0, b]
        for r in set(lam):
            value = value + coordinate_scale**r * mvals[q.remove_one(lam, r)] * integrated[r, b]
        columns.append(value)
    return np.column_stack(columns)


def score_replicate(
    *,
    k: int,
    degree: int,
    geometry: dict[str, float],
    vectors: np.ndarray,
    base_denominators: np.ndarray,
    base_numerators: np.ndarray,
    log2_n: int,
    seed: int,
    batch_log2: int,
    direct: bool,
) -> dict[str, object]:
    q.U = float(geometry["outer"])
    q.R = float(geometry["marginal"])
    q.DELTA = float(geometry["delta"])
    support = support_from_geometry(geometry)
    basis = q.basis_indices(degree)
    parts = q.all_partitions(degree // 2)
    coordinate_scale = k / (q.U * q.U)
    component_bits = int(round(math.log2(len(q.IMPORTANCE_COMPONENTS))))
    component_log2 = log2_n - component_bits
    batch_n = 1 << min(batch_log2, component_log2)
    batches = 1 << max(0, component_log2 - batch_log2)
    sample_n = 1 << log2_n
    delta_denominator = np.zeros(vectors.shape[1])
    delta_numerator = np.zeros(vectors.shape[1])
    direct_denominator = np.zeros(vectors.shape[1])
    direct_numerator = np.zeros(vectors.shape[1])
    weighted_acceptance = 0.0
    proposal_acceptance = 0
    scale = math.exp(2 * math.log(k) + (k - 1) * math.log(q.R) - k * math.log(q.U))

    for component_index, component in enumerate(q.IMPORTANCE_COMPONENTS):
        component_seed = seed + 10_007 * component_index
        x = q.importance_simplex_points(k, q.U, component_log2, component_seed, component)
        u = q.importance_simplex_points(k - 1, q.R, component_log2, component_seed + 1_000_003, component)
        wx = q.importance_weights(x, q.U)
        wu = q.importance_weights(u, q.R)
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            xb = x[sl]
            accept = support_acceptance(xb, support)
            proposal_acceptance += int(accept.sum())
            weighted_acceptance += float(wx[sl][accept].sum())
            if direct and np.any(accept):
                f_accept = q.features(
                    xb[accept], basis, parts, degree, k, coordinate_scale
                ) @ vectors
                direct_denominator += np.sum(
                    wx[sl][accept, None] * f_accept * f_accept, axis=0
                )
            rejected = ~accept
            if np.any(rejected):
                f = q.features(xb[rejected], basis, parts, degree, k, coordinate_scale) @ vectors
                delta_denominator -= np.sum(wx[sl][rejected, None] * f * f, axis=0)

            ub = u[sl]
            h_support = marginal_features(
                ub, basis, parts, degree, k, coordinate_scale, support
            ) @ vectors
            h_full = q.unrestricted_marginal_features(
                ub, basis, parts, degree, k, coordinate_scale
            ) @ vectors
            if direct:
                direct_numerator += scale * np.sum(
                    wu[sl, None] * h_support * h_support, axis=0
                )
            delta_numerator += scale * np.sum(
                wu[sl, None] * (h_support * h_support - h_full * h_full), axis=0
            )

    delta_denominator /= sample_n
    delta_numerator /= sample_n
    direct_denominator /= sample_n
    direct_numerator /= sample_n
    quotients = (base_numerators + delta_numerator) / (base_denominators + delta_denominator)
    result = {
        "seed": seed,
        "quotients": quotients.tolist(),
        "denominator_corrections": delta_denominator.tolist(),
        "numerator_corrections": delta_numerator.tolist(),
        "weighted_acceptance": weighted_acceptance / sample_n,
        "proposal_acceptance": proposal_acceptance / sample_n,
    }
    if direct:
        direct_quotients = direct_numerator / direct_denominator
        result.update({
            "direct_quotients": direct_quotients.tolist(),
            "direct_denominators": direct_denominator.tolist(),
            "direct_numerators": direct_numerator.tolist(),
        })
    return result


def score_replicate_stratified(
    *,
    k: int,
    degree: int,
    geometry: dict[str, float],
    vectors: np.ndarray,
    base_denominators: np.ndarray,
    base_numerators: np.ndarray,
    log2_n: int,
    seed: int,
    batch_log2: int,
) -> dict[str, object]:
    """Boundary correction partitioned by the exact large-coordinate count."""
    q.U = float(geometry["outer"])
    q.R = float(geometry["marginal"])
    q.DELTA = float(geometry["delta"])
    support = support_from_geometry(geometry)
    basis = q.basis_indices(degree)
    parts = q.all_partitions(degree // 2)
    coordinate_scale = k / (q.U * q.U)
    batch_n = 1 << min(batch_log2, log2_n)
    batches = 1 << max(0, log2_n - batch_log2)
    sample_n = 1 << log2_n
    scale = math.exp(2 * math.log(k) + (k - 1) * math.log(q.R) - k * math.log(q.U))
    delta_denominator = np.zeros(vectors.shape[1])
    delta_numerator = np.zeros(vectors.shape[1])
    strata = []

    def stratum_shape(dimension: int, m: int, target_big_mass: float) -> float:
        if m == 0:
            return 1.0
        target_big_mass = min(0.9, max(0.05, target_big_mass))
        return target_big_mass * (dimension + 1 - m) / (m * (1.0 - target_big_mass))

    def fixed_component_weights(points: np.ndarray, radius: float, m: int, shape: float):
        if m == 0:
            return np.ones(points.shape[0])
        dimension = points.shape[1]
        y = np.maximum(points[:, :m] / radius, np.finfo(float).tiny)
        log_ratio = (
            q.special.gammaln(dimension + 1 + m * (shape - 1.0))
            - m * q.special.gammaln(shape)
            - q.special.gammaln(dimension + 1)
            + (shape - 1.0) * np.log(y).sum(axis=1)
        )
        return np.exp(-log_ratio)

    for m in range(1, int(q.U // q.DELTA) + 1):
        shape = stratum_shape(k, m, 0.68)
        component_seed = seed + 10_007 * m
        x = q.importance_simplex_points(k, q.U, log2_n, component_seed, (m, shape))
        wx = fixed_component_weights(x, q.U, m, shape) * math.comb(k, m)
        component = np.zeros(vectors.shape[1])
        rejected_count = 0
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            xb = x[sl]
            big = xb > q.DELTA
            exact_cell = big[:, :m].all(axis=1) & (~big[:, m:]).all(axis=1)
            rejected = exact_cell & ~support_acceptance(xb, support)
            rejected_count += int(rejected.sum())
            if np.any(rejected):
                f = q.features(xb[rejected], basis, parts, degree, k, coordinate_scale) @ vectors
                component -= np.sum(wx[sl][rejected, None] * f * f, axis=0)
        delta_denominator += component
        strata.append({"kind": "M1", "m": m, "shape": shape, "selected": rejected_count,
                       "delta_mean": (component / sample_n).tolist()})

    dimension_u = k - 1
    for m in range(0, int(q.R // q.DELTA) + 1):
        shape = stratum_shape(dimension_u, m, 0.68)
        component_seed = seed + 500_009 + 10_007 * m
        u = q.importance_simplex_points(dimension_u, q.R, log2_n, component_seed, (m, shape))
        wu = fixed_component_weights(u, q.R, m, shape) * math.comb(dimension_u, m)
        component = np.zeros(vectors.shape[1])
        selected_count = 0
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            ub = u[sl]
            big = ub > q.DELTA
            exact_cell = big[:, :m].all(axis=1) & (~big[:, m:]).all(axis=1)
            selected_count += int(exact_cell.sum())
            if not np.any(exact_cell):
                continue
            cell_points = ub[exact_cell]
            h_support = marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale, support
            ) @ vectors
            h_full = q.unrestricted_marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale
            ) @ vectors
            component += scale * np.sum(
                wu[sl][exact_cell, None] * (h_support * h_support - h_full * h_full), axis=0
            )
        delta_numerator += component
        strata.append({"kind": "M2", "m": m, "shape": shape, "selected": selected_count,
                       "delta_mean": (component / sample_n).tolist()})

    delta_denominator /= sample_n
    delta_numerator /= sample_n
    quotients = (base_numerators + delta_numerator) / (base_denominators + delta_denominator)
    return {
        "seed": seed,
        "quotients": quotients.tolist(),
        "denominator_corrections": delta_denominator.tolist(),
        "numerator_corrections": delta_numerator.tolist(),
        "strata": strata,
    }


def evaluate(
    *,
    k: int,
    degree: int,
    geometry: dict[str, float],
    matrix_path: Path,
    log2_n: int,
    seeds: list[int],
    batch_log2: int,
    direct: bool,
) -> dict[str, object]:
    certified, reason = analytic_certificate(geometry)
    if not certified:
        return {"geometry": geometry, "analytic_feasible": False, "reason": reason}
    with np.load(matrix_path, allow_pickle=False) as payload:
        m1 = payload["m1"]
        m2 = payload["m2_unrestricted"]
        vectors = payload["vectors"]
    base_denominators = np.einsum("ik,ij,jk->k", vectors, m1, vectors, optimize=True)
    base_numerators = np.einsum("ik,ij,jk->k", vectors, m2, vectors, optimize=True)
    replicates = [
        score_replicate_stratified(
            k=k,
            degree=degree,
            geometry=geometry,
            vectors=vectors,
            base_denominators=base_denominators,
            base_numerators=base_numerators,
            log2_n=log2_n,
            seed=seed,
            batch_log2=batch_log2,
        )
        for seed in seeds
    ]
    values = np.asarray([replicate["quotients"] for replicate in replicates])
    means = values.mean(axis=0)
    errors = values.std(axis=0, ddof=1) / math.sqrt(len(values)) if len(values) > 1 else np.full(values.shape[1], np.nan)
    best = int(np.argmax(means))
    result = {
        "geometry": geometry,
        "analytic_feasible": True,
        "candidate_count": vectors.shape[1],
        "quotient_means": means.tolist(),
        "quotient_standard_errors": errors.tolist(),
        "best_candidate": best,
        "score": float(means[best]),
        "score_standard_error": float(errors[best]),
        "replicates": replicates,
    }
    if direct and "direct_quotients" in replicates[0]:
        direct_values = np.asarray([replicate["direct_quotients"] for replicate in replicates])
        direct_means = direct_values.mean(axis=0)
        direct_errors = direct_values.std(axis=0, ddof=1) / math.sqrt(len(direct_values)) if len(direct_values) > 1 else np.full(direct_values.shape[1], np.nan)
        direct_best = int(np.argmax(direct_means))
        result.update({
            "direct_quotient_means": direct_means.tolist(),
            "direct_quotient_standard_errors": direct_errors.tolist(),
            "direct_best_candidate": direct_best,
            "direct_score": float(direct_means[direct_best]),
            "direct_score_standard_error": float(direct_errors[direct_best]),
        })
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=49)
    parser.add_argument("--degree", type=int, default=21)
    parser.add_argument("--log2-n", type=int, default=14)
    parser.add_argument("--seeds", default="49101,49102")
    parser.add_argument("--batch-log2", type=int, default=9)
    parser.add_argument("--direct", action="store_true")
    args = parser.parse_args()
    candidates = json.loads(args.candidates.read_text())
    results = [
        evaluate(
            k=args.k,
            degree=args.degree,
            geometry=candidate,
            matrix_path=args.matrix,
            log2_n=args.log2_n,
            seeds=[int(seed) for seed in args.seeds.split(",") if seed],
            batch_log2=args.batch_log2,
            direct=args.direct,
        )
        for candidate in candidates
    ]
    payload = {
        "schema": "primegaps-support-geometry-screen-v1",
        "method": "fixed D=21 vector bank; analytic unrestricted forms plus exact-large-count stratified QMC boundary corrections",
        "k": args.k,
        "degree": args.degree,
        "log2_n": args.log2_n,
        "seeds": [int(seed) for seed in args.seeds.split(",") if seed],
        "results": results,
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for result in results:
        summary = {key: result.get(key) for key in ("geometry", "analytic_feasible", "score", "score_standard_error", "direct_score", "direct_score_standard_error", "reason")}
        print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
