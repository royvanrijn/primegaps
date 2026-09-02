#!/usr/bin/env python3
"""Independent QMC reproduction of Stadlmann's k=49, D=21 quotient.

This intentionally does not implement the paper's unpublished exact recurrence.
It integrates the published I and J forms directly with scrambled Sobol points,
using the same finite-dimensional polynomial space in a better-conditioned basis.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import linalg, special
from scipy.stats import qmc


U = 0.253 + 0.0075       # outer total bound A_1 + epsilon
R = 0.253 - 0.0075       # J marginal bound A_1 - epsilon
DELTA = 0.028

# Equal-weight symmetric Dirichlet mixture.  Component (r,a) tilts r randomly
# chosen coordinates from shape 1 to shape a.  Sampling the first r coordinates
# is equivalent for the symmetric integrands; weights use the full mixture.
IMPORTANCE_COMPONENTS = ((0, 1.0), (1, 20.0), (1, 70.0), (2, 15.0),
                         (2, 40.0), (3, 10.0), (3, 25.0), (4, 20.0))


def partitions(n: int, max_part: int | None = None):
    if n == 0:
        yield ()
        return
    max_part = n if max_part is None else min(max_part, n)
    for first in range(max_part, 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def all_partitions(max_degree: int) -> list[tuple[int, ...]]:
    out = [()]
    for n in range(1, max_degree + 1):
        out.extend(partitions(n))
    return sorted(out, key=lambda p: (len(p), sum(p), p))


def basis_indices(degree: int) -> list[tuple[tuple[int, ...], int]]:
    out = []
    for a in range(degree // 2 + 1):
        for lam in partitions(a):
            for b in range(degree - 2 * a + 1):
                out.append((lam, b))
    return out


def remove_one(p: tuple[int, ...], value: int) -> tuple[int, ...]:
    q = list(p)
    q.remove(value)
    return tuple(q)


def monomial_symmetric_values(
    coords: np.ndarray,
    parts: list[tuple[int, ...]],
    coordinate_scale: float,
) -> dict[tuple[int, ...], np.ndarray]:
    """Evaluate m_lambda(coordinate_scale * coords_i**2) for all parts.

    The recurrence p_r m_lambda is used only as an evaluation identity; it is
    independent of the integration recurrence sketched in the target paper.
    """
    n = coords.shape[0]
    max_degree = max(map(sum, parts), default=0)
    x = coordinate_scale * coords * coords
    power = {r: np.sum(x**r, axis=1) for r in range(1, max_degree + 1)}
    values: dict[tuple[int, ...], np.ndarray] = {(): np.ones(n)}
    for mu in parts:
        if not mu:
            continue
        r = mu[-1]
        lam = remove_one(mu, r)
        value = power[r] * values[lam]
        for q in set(lam):
            merged_list = list(remove_one(lam, q)) + [q + r]
            merged = tuple(sorted(merged_list, reverse=True))
            coefficient = merged.count(q + r)
            value = value - coefficient * values[merged]
        values[mu] = value / mu.count(r)
    return values


def eval_jacobi_basis(max_degree: int, k: int, q: np.ndarray) -> np.ndarray:
    """Stable orthonormal Jacobi basis for k(1-q)^(k-1)dq on [0,1]."""
    out = np.empty((q.size, max_degree + 1))
    x = 1.0 - 2.0 * q
    for b in range(max_degree + 1):
        # E[P_b^(0,k-1)(1-2q)^2] = k/(2b+k).
        out[:, b] = special.eval_jacobi(b, 0, k - 1, x) * math.sqrt((2 * b + k) / k)
    return out


def simplex_points(k: int, radius: float, log2_n: int, seed: int) -> np.ndarray:
    """Uniform points on {x_i>=0, sum x_i<=radius}, from scrambled Sobol."""
    engine = qmc.Sobol(d=k + 1, scramble=True, seed=seed)
    uniform = engine.random_base2(log2_n)
    exponential = -np.log(np.maximum(uniform, np.finfo(float).tiny))
    return radius * exponential[:, :k] / exponential.sum(axis=1, keepdims=True)


def importance_simplex_points(
    k: int, radius: float, log2_n: int, seed: int, component: tuple[int, float]
) -> np.ndarray:
    """Sample one fixed representative of a symmetric Dirichlet component."""
    r, shape = component
    engine = qmc.Sobol(d=k + 1, scramble=True, seed=seed)
    uniform = engine.random_base2(log2_n)
    gamma = -np.log1p(-np.minimum(uniform, 1.0 - np.finfo(float).eps))
    if r:
        gamma[:, :r] = special.gammaincinv(shape, uniform[:, :r])
    return radius * gamma[:, :k] / gamma.sum(axis=1, keepdims=True)


def log_elementary_symmetric(log_z: np.ndarray, degree: int) -> np.ndarray:
    """log e_degree(exp(log_z_i)), row by row."""
    n = log_z.shape[0]
    state = np.full((degree + 1, n), -np.inf)
    state[0] = 0.0
    for i in range(log_z.shape[1]):
        upper = min(degree, i + 1)
        for j in range(upper, 0, -1):
            state[j] = np.logaddexp(state[j], state[j - 1] + log_z[:, i])
    return state[degree]


def importance_weights(points: np.ndarray, radius: float) -> np.ndarray:
    """Uniform Dirichlet density divided by the symmetric proposal mixture."""
    k = points.shape[1]
    y = np.maximum(points / radius, np.finfo(float).tiny)
    log_ratios = []
    for r, shape in IMPORTANCE_COMPONENTS:
        if r == 0:
            log_ratios.append(np.zeros(points.shape[0]))
            continue
        log_constant = (
            special.gammaln(k + 1 + r * (shape - 1.0))
            - r * special.gammaln(shape)
            - special.gammaln(k + 1)
            - math.log(special.comb(k, r, exact=True))
        )
        log_e = log_elementary_symmetric((shape - 1.0) * np.log(y), r)
        log_ratios.append(log_constant + log_e)
    stacked = np.vstack(log_ratios)
    log_mixture_ratio = special.logsumexp(stacked, axis=0) - math.log(len(IMPORTANCE_COMPONENTS))
    return np.exp(-log_mixture_ratio)


def monomial_assignment_counts(
    k: int, lam: tuple[int, ...], mu: tuple[int, ...]
) -> dict[tuple[int, ...], int]:
    """Counts after fixing one m_lam monomial and assigning an m_mu monomial."""
    remaining = Counter(mu)
    totals: Counter[tuple[int, ...]] = Counter()

    def recurse(slot: int, current: list[int]) -> None:
        if slot == len(lam):
            rest = []
            divisor = 1
            count_rest = 0
            for value, count in remaining.items():
                rest.extend([value] * count)
                divisor *= math.factorial(count)
                count_rest += count
            placements = math.factorial(k - len(lam)) // math.factorial(k - len(lam) - count_rest)
            placements //= divisor
            nu = tuple(sorted(current + rest, reverse=True))
            totals[nu] += placements
            return
        recurse(slot + 1, current + [lam[slot]])
        for value in list(remaining):
            if remaining[value] == 0:
                continue
            remaining[value] -= 1
            recurse(slot + 1, current + [lam[slot] + value])
            remaining[value] += 1

    recurse(0, [])
    return dict(totals)


def monomial_term_count(k: int, lam: tuple[int, ...]) -> int:
    count = math.factorial(k) // math.factorial(k - len(lam))
    for multiplicity in Counter(lam).values():
        count //= math.factorial(multiplicity)
    return count


def full_simplex_reference_gram(k: int, degree: int) -> np.ndarray:
    """Analytic Gram under the full U-simplex, normalized by its volume."""
    basis = basis_indices(degree)
    size = len(basis)
    coordinate_scale = k / (U * U)
    jac_at_degree: dict[int, np.ndarray] = {}
    for a in range(degree + 1):
        kk = k + 2 * a
        nodes, weights = special.roots_jacobi(degree + 1, 0, kk - 1)
        q = (1.0 - nodes) * 0.5
        values = eval_jacobi_basis(degree, k, q)
        # Normalize weights to the Beta(1,kk) probability measure, then /kk
        # to recover the unnormalized radial integral.
        jac_at_degree[a] = (values.T * (weights / weights.sum())) @ values / kk

    overlap_cache: dict[tuple[tuple[int, ...], tuple[int, ...]], dict[tuple[int, ...], int]] = {}
    out = np.empty((size, size))
    for i, (lam, b) in enumerate(basis):
        nlam = monomial_term_count(k, lam)
        for j in range(i + 1):
            mu, d = basis[j]
            key = (lam, mu)
            counts = overlap_cache.get(key)
            if counts is None:
                counts = monomial_assignment_counts(k, lam, mu)
                overlap_cache[key] = counts
            a = sum(lam) + sum(mu)
            radial = jac_at_degree[a][b, d]
            total = 0.0
            common = (
                (coordinate_scale * U * U) ** a
                * math.exp(special.gammaln(k + 1) - special.gammaln(k + 2 * a))
                * radial
                * nlam
            )
            for nu, assignment_count in counts.items():
                factorial_product = math.prod(math.factorial(2 * value) for value in nu)
                total += common * assignment_count * factorial_product
            out[i, j] = out[j, i] = total
    return out


def whiten_transform(reference: np.ndarray) -> np.ndarray:
    """C with C.T @ reference @ C = identity."""
    chol = linalg.cholesky((reference + reference.T) * 0.5, lower=True)
    return linalg.solve_triangular(chol.T, np.eye(chol.shape[0]), lower=False)


def monomial_angular_coefficient(
    dimension: int,
    coordinate_scale: float,
    radius: float,
    eta: tuple[int, ...],
    theta: tuple[int, ...],
) -> float:
    """Angular/scale factor for E[m_eta m_theta f(sum)] on a simplex."""
    a = sum(eta) + sum(theta)
    counts = monomial_assignment_counts(dimension, eta, theta)
    combinatorial = 0
    for nu, assignment_count in counts.items():
        combinatorial += assignment_count * math.prod(math.factorial(2 * value) for value in nu)
    combinatorial *= monomial_term_count(dimension, eta)
    return (
        (coordinate_scale * radius * radius) ** a
        * math.exp(special.gammaln(dimension + 1) - special.gammaln(dimension + 2 * a + 1))
        * combinatorial
    )


def unrestricted_marginal_gram(k: int, degree: int) -> np.ndarray:
    """Exact-Gauss Gram of full last-coordinate marginals, before B cutoffs.

    The remaining k-1 coordinates are restricted only by sum(u)<=R.  All
    integrands are polynomials and the chosen Gauss orders are exact in exact
    arithmetic.
    """
    dimension = k - 1
    basis = basis_indices(degree)
    coordinate_scale = k / (U * U)
    max_r = degree // 2
    radial: dict[int, np.ndarray] = {}
    flat_size = (max_r + 1) * (degree + 1)
    for a in range(degree + 1):
        kk = dimension + 2 * a
        nodes, weights = special.roots_jacobi(degree + 2, 0, kk - 1)
        s = R * (1.0 + nodes) * 0.5
        room = U - s
        integrated = integrated_jacobi_moments(
            room,
            np.zeros_like(room),
            room,
            max_r,
            degree,
            k,
        )
        for r in range(max_r + 1):
            integrated[r] *= coordinate_scale**r
        flat = integrated.reshape(flat_size, -1)
        radial[a] = (flat * (weights / weights.sum())) @ flat.T

    angular_cache: dict[tuple[tuple[int, ...], tuple[int, ...]], float] = {}
    size = len(basis)
    out = np.empty((size, size))
    for i, (lam, b) in enumerate(basis):
        left_terms = [(lam, 0)] + [(remove_one(lam, r), r) for r in set(lam)]
        for j in range(i + 1):
            mu, d = basis[j]
            right_terms = [(mu, 0)] + [(remove_one(mu, r), r) for r in set(mu)]
            value = 0.0
            for eta, r in left_terms:
                for theta, s in right_terms:
                    key = (eta, theta)
                    angular = angular_cache.get(key)
                    if angular is None:
                        angular = monomial_angular_coefficient(
                            dimension, coordinate_scale, R, eta, theta
                        )
                        angular_cache[key] = angular
                    a = sum(eta) + sum(theta)
                    value += angular * radial[a][r * (degree + 1) + b, s * (degree + 1) + d]
            out[i, j] = out[j, i] = value
    return out


def support_acceptance(points: np.ndarray) -> np.ndarray:
    big = points > DELTA
    count = big.sum(axis=1)
    big_sum = np.where(big, points, 0.0).sum(axis=1)
    limit = np.where(count <= 2, 0.15, 0.17)
    return (points.sum(axis=1) <= U) & ((count == 0) | (big_sum <= limit))


def features(
    points: np.ndarray,
    basis: list[tuple[tuple[int, ...], int]],
    parts: list[tuple[int, ...]],
    degree: int,
    k: int,
    coordinate_scale: float,
) -> np.ndarray:
    mvals = monomial_symmetric_values(points, parts, coordinate_scale)
    q = (U - points.sum(axis=1)) / U
    jvals = eval_jacobi_basis(degree, k, q)
    return np.column_stack([mvals[lam] * jvals[:, b] for lam, b in basis])


def integrated_jacobi_moments(
    total_room: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    max_r: int,
    max_b: int,
    k: int,
) -> np.ndarray:
    """Gauss-integrate t^(2r) times the Jacobi basis over per-row intervals.

    With 24 nodes this is exact in exact arithmetic for every polynomial here:
    the largest possible degree is 2*10+21=41.
    """
    n = total_room.size
    out = np.zeros((max_r + 1, max_b + 1, n))
    half = np.maximum((hi - lo) * 0.5, 0.0)
    middle = (hi + lo) * 0.5
    nodes, weights = np.polynomial.legendre.leggauss(24)
    exponents = 2 * np.arange(max_r + 1)
    for node, weight in zip(nodes, weights):
        t = middle + half * node
        powers = t[None, :] ** exponents[:, None]
        q = np.clip((total_room - t) / U, 0.0, 1.0)
        jvals = eval_jacobi_basis(max_b, k, q)
        out += np.einsum("rn,nb,n->rbn", powers, jvals, half * weight, optimize=True)
    return out


def marginal_features(
    points: np.ndarray,
    basis: list[tuple[tuple[int, ...], int]],
    parts: list[tuple[int, ...]],
    degree: int,
    k: int,
    coordinate_scale: float,
) -> np.ndarray:
    """Compute h_i(u)=integral F_i(u,t)dt exactly in the last coordinate."""
    n = points.shape[0]
    total = points.sum(axis=1)
    room = U - total
    big_mask = points > DELTA
    count = big_mask.sum(axis=1)
    big_sum = np.where(big_mask, points, 0.0).sum(axis=1)

    old_limit = np.where(count <= 2, 0.15, 0.17)
    small_ok = (count == 0) | (big_sum <= old_limit)
    small_lo = np.zeros(n)
    small_hi = np.where(small_ok, np.minimum(DELTA, room), 0.0)

    new_count = count + 1
    new_limit = np.where(new_count <= 2, 0.15, 0.17)
    big_lo = np.full(n, DELTA)
    big_hi = np.minimum(room, new_limit - big_sum)
    big_hi = np.maximum(big_hi, 0.0)

    max_r = max(map(sum, parts), default=0)
    max_b = degree
    integrated = integrated_jacobi_moments(room, small_lo, small_hi, max_r, max_b, k)
    integrated += integrated_jacobi_moments(room, big_lo, big_hi, max_r, max_b, k)
    mvals = monomial_symmetric_values(points, parts, coordinate_scale)
    columns = []
    for lam, b in basis:
        value = mvals[lam] * integrated[0, b]
        for r in set(lam):
            value = value + (
                coordinate_scale**r
                * mvals[remove_one(lam, r)]
                * integrated[r, b]
            )
        columns.append(value)
    return np.column_stack(columns)


def unrestricted_marginal_features(
    points: np.ndarray,
    basis: list[tuple[tuple[int, ...], int]],
    parts: list[tuple[int, ...]],
    degree: int,
    k: int,
    coordinate_scale: float,
) -> np.ndarray:
    """Last-coordinate marginals with only the outer total bound imposed."""
    room = U - points.sum(axis=1)
    max_r = max(map(sum, parts), default=0)
    integrated = integrated_jacobi_moments(
        room, np.zeros_like(room), room, max_r, degree, k
    )
    mvals = monomial_symmetric_values(points, parts, coordinate_scale)
    columns = []
    for lam, b in basis:
        value = mvals[lam] * integrated[0, b]
        for r in set(lam):
            value = value + (
                coordinate_scale**r
                * mvals[remove_one(lam, r)]
                * integrated[r, b]
            )
        columns.append(value)
    return np.column_stack(columns)


def accumulate_grams(
    k: int,
    degree: int,
    log2_n: int,
    seed: int,
    batch_log2: int,
) -> tuple[np.ndarray, np.ndarray, dict, np.ndarray]:
    basis = basis_indices(degree)
    parts = all_partitions(degree // 2)
    coordinate_scale = k / (U * U)
    reference = full_simplex_reference_gram(k, degree)
    transform = whiten_transform(reference)
    size = len(basis)
    g1 = np.zeros((size, size))
    g2 = np.zeros((size, size))
    accepted = 0
    component_bits = int(round(math.log2(len(IMPORTANCE_COMPONENTS))))
    if (1 << component_bits) != len(IMPORTANCE_COMPONENTS) or log2_n < component_bits:
        raise ValueError("importance component count must be a power of two and fit sample size")
    component_log2 = log2_n - component_bits
    batch_n = 1 << min(batch_log2, component_log2)
    batches = 1 << max(0, component_log2 - batch_log2)
    weighted_acceptance = 0.0
    for component_index, component in enumerate(IMPORTANCE_COMPONENTS):
        component_seed = seed + 10_007 * component_index
        x = importance_simplex_points(k, U, component_log2, component_seed, component)
        u = importance_simplex_points(
            k - 1, R, component_log2, component_seed + 1_000_003, component
        )
        wx = importance_weights(x, U)
        wu = importance_weights(u, R)
        for batch in range(batches):
            sl = slice(batch * batch_n, (batch + 1) * batch_n)
            xb = x[sl]
            accept = support_acceptance(xb)
            accepted += int(accept.sum())
            weighted_acceptance += float(wx[sl][accept].sum())
            a = features(xb[accept], basis, parts, degree, k, coordinate_scale) @ transform
            weighted_a = a * wx[sl][accept, None]
            g1 += a.T @ weighted_a
            h = marginal_features(u[sl], basis, parts, degree, k, coordinate_scale) @ transform
            weighted_h = h * wu[sl, None]
            g2 += h.T @ weighted_h
    sample_n = 1 << log2_n
    g1 /= sample_n
    g2 /= sample_n
    scale = math.exp(2 * math.log(k) + (k - 1) * math.log(R) - k * math.log(U))
    metadata = {
        "k": k,
        "degree": degree,
        "basis_size": size,
        "log2_n": log2_n,
        "seed": seed,
        "proposal_acceptance": accepted / sample_n,
        "acceptance": weighted_acceptance / sample_n,
        "m2_scale": scale,
        "importance_components": IMPORTANCE_COMPONENTS,
    }
    return g1, scale * g2, metadata, transform


def solve(g1: np.ndarray, g2: np.ndarray) -> tuple[float, np.ndarray, dict]:
    g1 = (g1 + g1.T) * 0.5
    g2 = (g2 + g2.T) * 0.5
    m1_eigs = linalg.eigvalsh(g1, check_finite=False)
    condition = float(m1_eigs[-1] / m1_eigs[0])
    values, vectors = linalg.eigh(
        g2,
        g1,
        subset_by_index=[g1.shape[0] - 1, g1.shape[0] - 1],
        check_finite=False,
        driver="gvx",
    )
    vector = vectors[:, 0]
    return float(values[0]), vector, {
        "m1_min_eigenvalue": float(m1_eigs[0]),
        "m1_max_eigenvalue": float(m1_eigs[-1]),
        "m1_condition": condition,
    }


def validate_vector(
    k: int,
    degree: int,
    raw_vector: np.ndarray,
    log2_n: int,
    seeds: list[int],
    batch_log2: int,
) -> list[dict]:
    basis = basis_indices(degree)
    parts = all_partitions(degree // 2)
    coordinate_scale = k / (U * U)
    scale = math.exp(2 * math.log(k) + (k - 1) * math.log(R) - k * math.log(U))
    results = []
    for seed in seeds:
        sample_n = 1 << log2_n
        component_bits = int(round(math.log2(len(IMPORTANCE_COMPONENTS))))
        component_log2 = log2_n - component_bits
        batch_n = 1 << min(batch_log2, component_log2)
        batches = 1 << max(0, component_log2 - batch_log2)
        denominator = 0.0
        numerator = 0.0
        accepted = 0
        weighted_acceptance = 0.0
        for component_index, component in enumerate(IMPORTANCE_COMPONENTS):
            component_seed = seed + 10_007 * component_index
            x = importance_simplex_points(k, U, component_log2, component_seed, component)
            u = importance_simplex_points(
                k - 1, R, component_log2, component_seed + 1_000_003, component
            )
            wx = importance_weights(x, U)
            wu = importance_weights(u, R)
            for batch in range(batches):
                sl = slice(batch * batch_n, (batch + 1) * batch_n)
                xb = x[sl]
                accept = support_acceptance(xb)
                accepted += int(accept.sum())
                weighted_acceptance += float(wx[sl][accept].sum())
                f = features(xb[accept], basis, parts, degree, k, coordinate_scale) @ raw_vector
                denominator += float(np.dot(wx[sl][accept], f * f))
                h = marginal_features(u[sl], basis, parts, degree, k, coordinate_scale) @ raw_vector
                numerator += float(np.dot(wu[sl], h * h))
        quotient = scale * numerator / denominator
        results.append({
            "seed": seed,
            "quotient": quotient,
            "deficit_from_1": 1.0 - quotient,
            "proposal_acceptance": accepted / sample_n,
            "acceptance": weighted_acceptance / sample_n,
        })
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--degree", type=int, default=21)
    parser.add_argument("--train-log2", type=int, default=16)
    parser.add_argument("--train-seed", type=int, default=240049)
    parser.add_argument("--validate-log2", type=int, default=18)
    parser.add_argument("--validate-seeds", default="49001,49002,49003,49004")
    parser.add_argument("--batch-log2", type=int, default=12)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    started = time.time()
    g1, g2, metadata, transform = accumulate_grams(
        args.k, args.degree, args.train_log2, args.train_seed, args.batch_log2
    )
    train_value, vector, diagnostics = solve(g1, g2)
    raw_vector = transform @ vector
    validation = validate_vector(
        args.k,
        args.degree,
        raw_vector,
        args.validate_log2,
        [int(x) for x in args.validate_seeds.split(",") if x],
        args.batch_log2,
    )
    quotients = np.array([x["quotient"] for x in validation])
    result = {
        **metadata,
        **diagnostics,
        "train_generalized_eigenvalue": train_value,
        "validate_log2": args.validate_log2,
        "validation": validation,
        "validation_mean": float(quotients.mean()),
        "validation_standard_error": float(quotients.std(ddof=1) / math.sqrt(len(quotients)))
        if len(quotients) > 1 else None,
        "validation_min": float(quotients.min()),
        "validation_max": float(quotients.max()),
        "elapsed_seconds": time.time() - started,
        "parameters": {
            "epsilon": 0.0075,
            "delta": DELTA,
            "A": [-0.0075, 0.253],
            "B_1_1": 0.15,
            "B_1_2": 0.15,
            "B_1_m_ge_3": 0.17,
            "c1": 0,
            "c2": 0,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output.with_suffix(".npz"),
        g1=g1,
        g2=g2,
        whitened_vector=vector,
        raw_vector=raw_vector,
        transform=transform,
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

