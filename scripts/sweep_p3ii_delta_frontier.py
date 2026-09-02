#!/usr/bin/env python3
"""Numerical A_max -> lambda_48 frontier for the P3.II.delta relaxation.

This is an explicit screening computation, not a proof certificate.  At every
endpoint it rebuilds the degree-21 unrestricted vector bank, then evaluates the
fixed two-band support with exact-large-count stratified randomized QMC.  The
same seeds are used at every endpoint so later analysis can retain the paired
covariance of curve differences.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import platform
import time
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

from primegaps.distribution import Minorant, support_constraint_failures


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / ".research/work/p3ii-delta-frontier/delta-frontier-root"
ARCHIVE = ROOT / "reproduction/240/numerical"
QMC_PATH = ARCHIVE / "qmc_verifier.py"
CONTROL_PATH = ARCHIVE / "control_variate_verifier.py"
GEOMETRY_PATH = ARCHIVE / "geometry_screen.py"

K = 48
DEGREE = 21
DELTA = Fraction(7, 250)
EPSILON = Fraction(17, 2000)
THEOREM_EPSILON = Fraction(1, 10**10)
MINORANT = Minorant("0.38", "0.4", "0.4")
TOLERANCES = [3e-12, 1e-12, 3e-13, 1e-13]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_geometry_screen():
    """Load the archived scorer while redirecting its moved qmc dependency."""
    original = importlib.util.spec_from_file_location

    def redirected(name, location, *args, **kwargs):
        path = Path(location)
        if path.name == "qmc_verifier.py" and not path.exists():
            location = QMC_PATH
        return original(name, location, *args, **kwargs)

    importlib.util.spec_from_file_location = redirected
    try:
        return load_module("p3ii_delta_geometry", GEOMETRY_PATH)
    finally:
        importlib.util.spec_from_file_location = original


def exact_constraint_bounds() -> dict[str, Fraction]:
    e = THEOREM_EPSILON
    d = DELTA
    x1 = Fraction(19, 50)
    x2 = x3 = Fraction(2, 5)
    return {
        "P3.II.delta.branch2": (x2 / 4 + Fraction(11, 16) - 2 * e - d) / 3,
        "P3.II.delta.branch1": (
            x2 / 10 + Fraction(4, 5) - 2 * e - d
        ) * Fraction(5, 16),
        "P3.II.range": (Fraction(19, 2) - 13 * d + 100 * e) / 36,
        "P3.I.branch1": (x1 + Fraction(2, 3) - 2 * e - d) / 4,
        "P3.I.branch2": (Fraction(9, 7) - 2 * e - d) * Fraction(7, 34),
        "P3.III": (
            Fraction(11, 8) - Fraction(9, 8) * x3 - 2 * e - d
        ) * Fraction(2, 7),
    }


def geometry_for(endpoint: Decimal) -> dict[str, object]:
    a = float(endpoint)
    epsilon = float(EPSILON)
    return {
        "delta": float(DELTA),
        "outer": a + epsilon,
        "marginal": a - epsilon,
        "bands": [
            {"a_upper": 0.23, "b12": 0.18, "b3plus": 0.20},
            {"a_upper": a, "b12": 0.15, "b3plus": 0.17},
        ],
    }


def curve_grid(interior_count: int) -> list[Decimal]:
    """Dense frontier grid plus the previously reported counterfactual points."""
    if interior_count < 2:
        raise ValueError("interior_count must be at least 2")
    getcontext().prec = 30
    bounds = exact_constraint_bounds()
    start = Decimal(bounds["P3.II.delta.branch2"].numerator) / Decimal(
        bounds["P3.II.delta.branch2"].denominator
    )
    stop = Decimal(bounds["P3.II.range"].numerator) / Decimal(
        bounds["P3.II.range"].denominator
    )
    points = {Decimal("0.253"), start, stop}
    for index in range(interior_count + 1):
        points.add(start + (stop - start) * Decimal(index) / Decimal(interior_count))
    points.update(
        Decimal(value)
        for value in (
            "0.2531",
            "0.2532",
            "0.2533",
            "0.2534",
            "0.2535",
            "0.2536",
            "0.2537",
            "0.25375",
        )
    )
    return sorted(point for point in points if point <= stop)


def write_checkpoint(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def shifted_stratum_volume_weight(
    dimension: int, large_count: int, radius: float, delta: float
) -> float:
    """Normalized volume of all permutations of one translated simplex."""
    residual = radius - large_count * delta
    if dimension < 1 or not 0 <= large_count <= dimension or residual < 0:
        return 0.0
    return math.comb(dimension, large_count) * (residual / radius) ** dimension


def score_replicate_stabilized(
    geometry_module,
    *,
    k: int,
    degree: int,
    geometry: dict[str, object],
    vectors: np.ndarray,
    base_denominators: np.ndarray,
    base_numerators: np.ndarray,
    log2_n: int,
    seed: int,
    batch_log2: int,
) -> dict[str, object]:
    """Stratified correction with bounded-balance weights on the I side.

    The legacy single Dirichlet tilt has unbounded importance weights and its
    m=2 denominator correction produces a long right tail in lambda.  Here the
    proposal is an equal mixture containing the uniform simplex law.  The
    balance weight is therefore bounded by the number of mixture components.
    The J-side estimator is retained unchanged; its observed variance is much
    smaller and its marginal-feature evaluation is the expensive part.
    """
    q = geometry_module.q
    q.U = float(geometry["outer"])
    q.R = float(geometry["marginal"])
    q.DELTA = float(geometry["delta"])
    support = geometry_module.support_from_geometry(geometry)
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

    def shape_for_mass(dimension: int, m: int, mass: float) -> float:
        if m == 0:
            return 1.0
        return mass * (dimension + 1 - m) / (m * (1.0 - mass))

    def density_ratio(points: np.ndarray, radius: float, m: int, shape: float):
        """Tilted Dirichlet density divided by the uniform-simplex density."""
        if shape == 1.0 or m == 0:
            return np.ones(points.shape[0])
        dimension = points.shape[1]
        y = np.maximum(points[:, :m] / radius, np.finfo(float).tiny)
        log_ratio = (
            q.special.gammaln(dimension + 1 + m * (shape - 1.0))
            - m * q.special.gammaln(shape)
            - q.special.gammaln(dimension + 1)
            + (shape - 1.0) * np.log(y).sum(axis=1)
        )
        return np.exp(log_ratio)

    for m in range(1, int(q.U // q.DELTA) + 1):
        shapes = tuple(
            dict.fromkeys(
                [
                    1.0,
                    shape_for_mass(k, m, 0.45),
                    shape_for_mass(k, m, 0.50),
                    shape_for_mass(k, m, 0.54),
                    shape_for_mass(k, m, 0.58),
                    shape_for_mass(k, m, 0.62),
                    shape_for_mass(k, m, 0.66),
                    shape_for_mass(k, m, 0.70),
                ]
            )
        )
        component = np.zeros(vectors.shape[1])
        rejected_counts = []
        maximum_weight = 0.0
        for component_index, shape in enumerate(shapes):
            component_seed = seed + 10_007 * m + 1_000_003 * component_index
            x = q.importance_simplex_points(k, q.U, log2_n, component_seed, (m, shape))
            mixture_ratio = np.mean(
                [density_ratio(x, q.U, m, other) for other in shapes], axis=0
            )
            weights = math.comb(k, m) / mixture_ratio
            maximum_weight = max(maximum_weight, float(weights.max()))
            rejected_count = 0
            for batch in range(batches):
                sl = slice(batch * batch_n, (batch + 1) * batch_n)
                xb = x[sl]
                big = xb > q.DELTA
                exact_cell = big[:, :m].all(axis=1) & (~big[:, m:]).all(axis=1)
                rejected = exact_cell & ~geometry_module.support_acceptance(xb, support)
                rejected_count += int(rejected.sum())
                if np.any(rejected):
                    f = q.features(
                        xb[rejected], basis, parts, degree, k, coordinate_scale
                    ) @ vectors
                    component -= np.sum(
                        weights[sl][rejected, None] * f * f, axis=0
                    )
            rejected_counts.append(rejected_count)
        component /= len(shapes) * sample_n
        delta_denominator += component
        strata.append(
            {
                "kind": "M1-balanced-mixture",
                "m": m,
                "shapes": list(shapes),
                "rejected_by_component": rejected_counts,
                "maximum_weight": maximum_weight,
                "delta_mean": component.tolist(),
            }
        )

    # Retain the legacy exact-count J correction.  This avoids multiplying the
    # costly analytic marginal evaluation by the four I-side mixture components.
    dimension_u = k - 1
    for m in range(0, int(q.R // q.DELTA) + 1):
        shape = shape_for_mass(dimension_u, m, 0.68)
        component_seed = seed + 500_009 + 10_007 * m
        u = q.importance_simplex_points(
            dimension_u, q.R, log2_n, component_seed, (m, shape)
        )
        inverse_ratio = 1.0 / density_ratio(u, q.R, m, shape)
        weights = inverse_ratio * math.comb(dimension_u, m)
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
            h_support = geometry_module.marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale, support
            ) @ vectors
            h_full = q.unrestricted_marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale
            ) @ vectors
            component += scale * np.sum(
                weights[sl][exact_cell, None]
                * (h_support * h_support - h_full * h_full),
                axis=0,
            )
        component /= sample_n
        delta_numerator += component
        strata.append(
            {
                "kind": "M2",
                "m": m,
                "shape": shape,
                "selected": selected_count,
                "delta_mean": component.tolist(),
            }
        )

    quotients = (base_numerators + delta_numerator) / (
        base_denominators + delta_denominator
    )
    return {
        "seed": seed,
        "quotients": quotients.tolist(),
        "denominator_corrections": delta_denominator.tolist(),
        "numerator_corrections": delta_numerator.tolist(),
        "strata": strata,
    }


def score_replicate_shifted(
    geometry_module,
    *,
    k: int,
    degree: int,
    geometry: dict[str, object],
    vectors: np.ndarray,
    base_denominators: np.ndarray,
    base_numerators: np.ndarray,
    log2_n: int,
    seed: int,
    batch_log2: int,
) -> dict[str, object]:
    """Exact-count QMC using translated simplices and constant weights.

    In the stratum with exactly m large coordinates, translate the nominated m
    coordinates by delta and sample uniformly from the residual simplex.  The
    normalized-volume multiplier is then C(d,m)((R-m*delta)/R)^d.  This avoids
    the unbounded Dirichlet importance weights of the legacy estimator.
    """
    q = geometry_module.q
    q.U = float(geometry["outer"])
    q.R = float(geometry["marginal"])
    q.DELTA = float(geometry["delta"])
    support = geometry_module.support_from_geometry(geometry)
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

    for m in range(1, int(q.U // q.DELTA) + 1):
        residual = q.U - m * q.DELTA
        x = q.simplex_points(k, residual, log2_n, seed + 10_007 * m)
        x[:, :m] += q.DELTA
        volume_weight = shifted_stratum_volume_weight(k, m, q.U, q.DELTA)
        component = np.zeros(vectors.shape[1])
        exact_count = 0
        rejected_count = 0
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            xb = x[sl]
            exact_cell = (xb[:, m:] <= q.DELTA).all(axis=1)
            rejected = exact_cell & ~geometry_module.support_acceptance(xb, support)
            exact_count += int(exact_cell.sum())
            rejected_count += int(rejected.sum())
            if np.any(rejected):
                f = q.features(
                    xb[rejected], basis, parts, degree, k, coordinate_scale
                ) @ vectors
                component -= volume_weight * np.sum(f * f, axis=0)
        component /= sample_n
        delta_denominator += component
        strata.append(
            {
                "kind": "M1-shifted-simplex",
                "m": m,
                "residual_radius": residual,
                "volume_weight": volume_weight,
                "exact_count": exact_count,
                "rejected_count": rejected_count,
                "delta_mean": component.tolist(),
            }
        )

    dimension_u = k - 1
    for m in range(0, int(q.R // q.DELTA) + 1):
        residual = q.R - m * q.DELTA
        u = q.simplex_points(
            dimension_u, residual, log2_n, seed + 500_009 + 10_007 * m
        )
        if m:
            u[:, :m] += q.DELTA
        volume_weight = shifted_stratum_volume_weight(
            dimension_u, m, q.R, q.DELTA
        )
        component = np.zeros(vectors.shape[1])
        exact_count = 0
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            ub = u[sl]
            exact_cell = (ub[:, m:] <= q.DELTA).all(axis=1)
            exact_count += int(exact_cell.sum())
            if not np.any(exact_cell):
                continue
            cell_points = ub[exact_cell]
            h_support = geometry_module.marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale, support
            ) @ vectors
            h_full = q.unrestricted_marginal_features(
                cell_points, basis, parts, degree, k, coordinate_scale
            ) @ vectors
            component += scale * volume_weight * np.sum(
                h_support * h_support - h_full * h_full, axis=0
            )
        component /= sample_n
        delta_numerator += component
        strata.append(
            {
                "kind": "M2-shifted-simplex",
                "m": m,
                "residual_radius": residual,
                "volume_weight": volume_weight,
                "exact_count": exact_count,
                "delta_mean": component.tolist(),
            }
        )

    quotients = (base_numerators + delta_numerator) / (
        base_denominators + delta_denominator
    )
    return {
        "seed": seed,
        "quotients": quotients.tolist(),
        "denominator_corrections": delta_denominator.tolist(),
        "numerator_corrections": delta_numerator.tolist(),
        "strata": strata,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interior-count", type=int, default=16)
    parser.add_argument("--log2-n", type=int, default=11)
    parser.add_argument(
        "--seeds",
        default="48301,48302,48303,48304,48305,48306,48307,48308",
    )
    parser.add_argument("--batch-log2", type=int, default=9)
    parser.add_argument(
        "--balanced-m1",
        action="store_true",
        help="use bounded balance-heuristic mixture weights for I corrections",
    )
    parser.add_argument(
        "--shifted-strata",
        action="store_true",
        help="sample each exact large-count cell via a translated simplex",
    )
    parser.add_argument("--output", type=Path, default=WORK / "curve-raw.json")
    parser.add_argument(
        "--only",
        help="comma-separated decimal endpoints; overrides the standard grid",
    )
    parser.add_argument("--stop-after", type=int)
    args = parser.parse_args()
    if args.balanced_m1 and args.shifted_strata:
        parser.error("choose at most one alternative stratum estimator")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    seeds = [int(value) for value in args.seeds.split(",") if value]
    endpoints = (
        [Decimal(value) for value in args.only.split(",") if value]
        if args.only
        else curve_grid(args.interior_count)
    )
    geometry = load_geometry_screen()
    control = load_module("p3ii_delta_control", CONTROL_PATH)

    if args.output.exists():
        payload = json.loads(args.output.read_text())
        completed = {item["A_max_decimal"] for item in payload["results"]}
    else:
        payload = {
            "schema": "primegaps.p3ii-delta-frontier-raw.v1",
            "status": "numerical-screening-not-certificate",
            "target": "P3.II.delta",
            "response": "lambda_48",
            "k": K,
            "degree": DEGREE,
            "delta": str(DELTA),
            "support_epsilon": str(EPSILON),
            "minorant": {"xi1": "0.38", "xi2": "0.4", "xi3": "0.4"},
            "fixed_support_family": {
                "inner_band": {"A_upper": "0.23", "B12": "0.18", "B3plus": "0.20"},
                "outer_band": {"B12": "0.15", "B3plus": "0.17"},
            },
            "constraint_bounds": {
                name: {
                    "fraction": str(value),
                    "decimal": format(float(value), ".15g"),
                }
                for name, value in exact_constraint_bounds().items()
            },
            "method": (
                "endpoint-specific degree-21 unrestricted vector banks; "
                + (
                    "translated-simplex exact-count randomized-QMC corrections"
                    if args.shifted_strata
                    else (
                        "bounded-balance I mixture and exact-large-count randomized-QMC corrections"
                        if args.balanced_m1
                        else "exact-large-count stratified randomized QMC boundary corrections"
                    )
                )
            ),
            "uncertainty": "standard error across common-seed randomized-QMC replicates",
            "tolerances": TOLERANCES,
            "log2_n_per_stratum": args.log2_n,
            "seeds": seeds,
            "environment": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "numpy": np.__version__,
            },
            "results": [],
        }
        completed = set()

    computed_this_run = 0
    for endpoint in endpoints:
        endpoint_text = format(endpoint, "f")
        if endpoint_text in completed:
            continue
        started = time.time()
        candidate = geometry_for(endpoint)
        support = geometry.support_from_geometry(candidate)
        failures = support_constraint_failures(support, MINORANT)
        failure_ids = sorted({item.constraint_id for item in failures})

        control.q.U = float(candidate["outer"])
        control.q.R = float(candidate["marginal"])
        control.q.DELTA = float(candidate["delta"])
        m1, m2, vectors, bank_records = control.candidates(
            K, DEGREE, TOLERANCES
        )
        base_denominators = np.einsum(
            "ik,ij,jk->k", vectors, m1, vectors, optimize=True
        )
        base_numerators = np.einsum(
            "ik,ij,jk->k", vectors, m2, vectors, optimize=True
        )

        replicate_function = (
            score_replicate_shifted
            if args.shifted_strata
            else (
                score_replicate_stabilized
                if args.balanced_m1
                else geometry.score_replicate_stratified
            )
        )
        replicates = [
            replicate_function(
                **(
                    {"geometry_module": geometry}
                    if args.balanced_m1 or args.shifted_strata
                    else {}
                ),
                k=K,
                degree=DEGREE,
                geometry=candidate,
                vectors=vectors,
                base_denominators=base_denominators,
                base_numerators=base_numerators,
                log2_n=args.log2_n,
                seed=seed,
                batch_log2=args.batch_log2,
            )
            for seed in seeds
        ]
        quotient_array = np.asarray([item["quotients"] for item in replicates])
        means = quotient_array.mean(axis=0)
        errors = (
            quotient_array.std(axis=0, ddof=1) / math.sqrt(len(replicates))
            if len(replicates) > 1
            else np.full(quotient_array.shape[1], np.nan)
        )
        best = int(np.argmax(means))
        result = {
            "A_max_decimal": endpoint_text,
            "A_max": float(endpoint),
            "geometry": candidate,
            "failed_constraints": failure_ids,
            "bank_records": bank_records,
            "quotient_means": means.tolist(),
            "quotient_standard_errors": errors.tolist(),
            "best_candidate": best,
            "best_tolerance": TOLERANCES[best],
            "lambda_48": float(means[best]),
            "lambda_48_standard_error": float(errors[best]),
            "lambda_48_replicates": quotient_array[:, best].tolist(),
            "replicates": replicates,
            "elapsed_seconds": time.time() - started,
        }
        payload["results"].append(result)
        payload["results"].sort(key=lambda item: item["A_max"])
        write_checkpoint(args.output, payload)
        completed.add(endpoint_text)
        computed_this_run += 1
        print(
            json.dumps(
                {
                    "A_max": endpoint_text,
                    "failed_constraints": failure_ids,
                    "lambda_48": result["lambda_48"],
                    "standard_error": result["lambda_48_standard_error"],
                    "best_tolerance": result["best_tolerance"],
                    "elapsed_seconds": result["elapsed_seconds"],
                },
                sort_keys=True,
            ),
            flush=True,
        )
        if args.stop_after is not None and computed_this_run >= args.stop_after:
            break


if __name__ == "__main__":
    main()
