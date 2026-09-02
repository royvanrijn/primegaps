#!/usr/bin/env python3
"""Exact symmetry-compressed large/small status densities.

This is a standalone research prototype.  It is deliberately outside tracked
source until its normalization and scaling have been reviewed.

For one canonical monomial ``prod_i t_i**g_i``, split each coordinate at
``delta``.  A large coordinate and an inclusion--exclusion-shifted small
coordinate are written ``t = delta + y``.  After collapsing the coordinates
in each status to their group sums ``x`` (large) and ``z`` (small), the
Dirichlet density only needs

    (r, h, P, Q)

where r is the number of large coordinates, h is the number of shifted-small
coordinates, and P/Q are the residual power totals in the two groups.

The returned coefficient already contains:

* every labeled large/small status choice for the canonical monomial;
* the inclusion--exclusion sign;
* every binomial term from ``(delta+y)**g``;
* the products of residual factorials required by the group-sum density.

It does *not* contain the final density denominators
``(P+r-1)! (Q+k-r-1)!`` or the monomial-symmetric orbit size.  Keeping those
outside makes the result reusable by both I and J assemblers.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import product
from math import comb, factorial, floor
import os
from time import perf_counter
from typing import Iterable


State = tuple[int, int, int, int]
Density = dict[State, Fraction]


def _fraction(value: int | str | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


def shifted_terms(exponent: int, delta: Fraction) -> tuple[tuple[int, Fraction], ...]:
    """Return ``(residual_power, binomial*delta_power*residual_factorial)``."""
    if exponent < 0:
        raise ValueError("exponents must be non-negative")
    g_factorial = factorial(exponent)
    return tuple(
        (power, Fraction(g_factorial, factorial(exponent - power)) * delta ** (exponent - power))
        for power in range(exponent + 1)
    )


def orbit_status_densities(
    signature: Iterable[int],
    *,
    k: int,
    delta: int | str | Fraction,
    max_large: int,
    max_offset_count: int,
) -> Density:
    """Build the exact truncated ``(r,h,P,Q)`` density generating function.

    ``signature`` contains only the positive exponents; zeroes are padded to
    length ``k``.  Terms with ``r > max_large`` or
    ``r+h > max_offset_count`` are discarded as geometrically empty for the
    caller's total cap.
    """
    exponents = tuple(int(value) for value in signature)
    if any(value <= 0 for value in exponents):
        raise ValueError("signature must contain positive exponents only")
    if len(exponents) > k:
        raise ValueError("signature is longer than k")
    if max_large < 0 or max_offset_count < 0:
        raise ValueError("truncation bounds must be non-negative")
    max_large = min(max_large, k, max_offset_count)
    delta_value = _fraction(delta)
    if not 0 < delta_value < 1:
        raise ValueError("delta must be in (0,1)")

    # Positive exponents first keeps the later zero-coordinate transitions
    # particularly cheap.  Coordinate labels are irrelevant to the product.
    padded = exponents + (0,) * (k - len(exponents))
    shifted_cache = {
        exponent: shifted_terms(exponent, delta_value) for exponent in set(padded)
    }
    density: Density = {(0, 0, 0, 0): Fraction(1)}
    for exponent in padded:
        shifted = shifted_cache[exponent]
        unshifted_small = factorial(exponent)
        next_density: defaultdict[State, Fraction] = defaultdict(Fraction)
        for (large_count, shifted_small_count, large_power, small_power), coefficient in density.items():
            # Small and not selected by inclusion--exclusion.
            next_density[
                (large_count, shifted_small_count, large_power, small_power + exponent)
            ] += coefficient * unshifted_small

            # Large: every such coordinate is shifted by delta.
            if large_count < max_large and large_count + shifted_small_count < max_offset_count:
                for power, multiplier in shifted:
                    next_density[
                        (large_count + 1, shifted_small_count, large_power + power, small_power)
                    ] += coefficient * multiplier

            # Small face removed by inclusion--exclusion: same shift, negative sign.
            if large_count + shifted_small_count < max_offset_count:
                for power, multiplier in shifted:
                    next_density[
                        (large_count, shifted_small_count + 1, large_power, small_power + power)
                    ] -= coefficient * multiplier
        density = {state: value for state, value in next_density.items() if value}
    return density


def monomial_symmetric_orbit_size(signature: Iterable[int], k: int) -> int:
    """Number of distinct permutations of ``signature`` padded by zeroes."""
    signature = tuple(signature)
    if len(signature) > k:
        raise ValueError("signature is longer than k")
    counts: dict[int, int] = defaultdict(int)
    counts[0] = k - len(signature)
    for exponent in signature:
        counts[int(exponent)] += 1
    result = factorial(k)
    for multiplicity in counts.values():
        result //= factorial(multiplicity)
    return result


def normalized_density_terms(density: Density, k: int):
    """Yield group-sum powers and exact density coefficients.

    The large group is absent when r=0.  The small group is absent when r=k.
    Their corresponding power is reported as ``None``.
    """
    for (r, h, p, q), coefficient in density.items():
        if r == 0:
            if p:
                raise AssertionError("an absent large group cannot carry power")
            x_power = None
            large_denominator = 1
        else:
            x_power = p + r - 1
            large_denominator = factorial(x_power)
        if r == k:
            if q:
                raise AssertionError("an absent small group cannot carry power")
            z_power = None
            small_denominator = 1
        else:
            z_power = q + (k - r) - 1
            small_denominator = factorial(z_power)
        yield (r, h, x_power, z_power), coefficient / (
            large_denominator * small_denominator
        )


def _explicit_status_densities(
    signature: tuple[int, ...],
    *,
    k: int,
    delta: Fraction,
    max_large: int,
    max_offset_count: int,
) -> Density:
    """Slow labeled-coordinate enumeration used only as an independent check."""
    exponents = signature + (0,) * (k - len(signature))
    result: defaultdict[State, Fraction] = defaultdict(Fraction)

    # status 0 = unshifted small, 1 = large, 2 = shifted small face.
    for statuses in product(range(3), repeat=k):
        r = statuses.count(1)
        h = statuses.count(2)
        if r > max_large or r + h > max_offset_count:
            continue
        choices = []
        for exponent, status in zip(exponents, statuses):
            if status == 0:
                choices.append(((exponent, Fraction(factorial(exponent))),))
            else:
                choices.append(shifted_terms(exponent, delta))
        for selected in product(*choices):
            p = sum(power for (power, _), status in zip(selected, statuses) if status == 1)
            q = sum(power for (power, _), status in zip(selected, statuses) if status != 1)
            coefficient = Fraction(-1 if h % 2 else 1)
            for _, multiplier in selected:
                coefficient *= multiplier
            result[(r, h, p, q)] += coefficient
    return {state: value for state, value in result.items() if value}


def _two_group_monomial_integral(
    x_power: int | None,
    z_power: int | None,
    height: Fraction,
    large_cap: Fraction | None,
) -> Fraction:
    """Integrate ``x^X z^Z`` over x+z<=height and x<=large_cap."""
    if height <= 0:
        return Fraction(0)
    if x_power is None and z_power is None:
        return Fraction(1)
    if x_power is None:
        assert z_power is not None
        return height ** (z_power + 1) / (z_power + 1)
    maximum = height if large_cap is None else min(height, large_cap)
    if maximum <= 0:
        return Fraction(0)
    if z_power is None:
        return maximum ** (x_power + 1) / (x_power + 1)
    answer = Fraction(0)
    n = z_power + 1
    for j in range(n + 1):
        answer += (
            Fraction((-1) ** j * comb(n, j), n * (x_power + j + 1))
            * height ** (n - j)
            * maximum ** (x_power + j + 1)
        )
    return answer


def compressed_one_band_monomial_integral(
    signature: tuple[int, ...],
    *,
    k: int,
    delta: Fraction,
    total_cap: Fraction,
    large_caps: tuple[Fraction, ...],
) -> Fraction:
    """Exact canonical-monomial I integral for one support band.

    This small wrapper exists for validation.  The production assembler should
    reuse the density polynomial directly and include radial slack powers.
    """
    maximum_offset = floor(total_cap / delta)
    maximum_large = 0
    for r, cap in enumerate(large_caps, start=1):
        if cap > r * delta:
            maximum_large = r
    maximum_large = min(maximum_large, k)
    density = orbit_status_densities(
        signature,
        k=k,
        delta=delta,
        max_large=maximum_large,
        max_offset_count=maximum_offset,
    )
    answer = Fraction(0)
    for (r, h, x_power, z_power), coefficient in normalized_density_terms(density, k):
        height = total_cap - (r + h) * delta
        if r == 0:
            residual_large_cap = None
        else:
            residual_large_cap = large_caps[r - 1] - r * delta
        answer += coefficient * _two_group_monomial_integral(
            x_power, z_power, height, residual_large_cap
        )
    return answer


def self_check() -> dict[str, object]:
    """Run explicit-DP and existing-reference exact checks."""
    delta = Fraction(1, 4)
    explicit_cases = (((), 3), ((1,), 3), ((2, 1), 3), ((2, 1, 1), 4))
    for signature, k in explicit_cases:
        dynamic = orbit_status_densities(
            signature,
            k=k,
            delta=delta,
            max_large=min(k, 2),
            max_offset_count=2,
        )
        explicit = _explicit_status_densities(
            signature,
            k=k,
            delta=delta,
            max_large=min(k, 2),
            max_offset_count=2,
        )
        if dynamic != explicit:
            missing = set(dynamic) ^ set(explicit)
            differing = {
                key for key in set(dynamic) & set(explicit) if dynamic[key] != explicit[key]
            }
            raise AssertionError(
                f"DP mismatch for {signature}, k={k}: missing={missing}, differing={differing}"
            )

    # Import the repo only for the second, end-to-end check.
    from primegaps.integrals import ExactSupportParameters, _integrate_support_monomial

    support = ExactSupportParameters.from_values(
        delta=delta,
        epsilon=Fraction(1, 20),
        A=(Fraction(-1, 20), Fraction(9, 20)),
        B=((Fraction(3, 10), Fraction(9, 20), Fraction(9, 20), Fraction(9, 20)),),
    )
    reference_cases = (
        ((), 2),
        ((1,), 2),
        ((2,), 3),
        ((2, 1), 3),
        ((3, 1, 1), 4),
    )
    comparisons = []
    for signature, k in reference_cases:
        exponents = signature + (0,) * (k - len(signature))
        compressed = compressed_one_band_monomial_integral(
            signature,
            k=k,
            delta=support.delta,
            total_cap=support.A[1] + support.epsilon,
            large_caps=support.B[0],
        )
        reference = _integrate_support_monomial(exponents, support)
        if compressed != reference:
            raise AssertionError(
                f"integral mismatch for {signature}, k={k}: {compressed} != {reference}"
            )
        comparisons.append(
            {"signature": signature, "k": k, "value": str(compressed)}
        )
    return {
        "explicit_density_cases": len(explicit_cases),
        "reference_integral_cases": comparisons,
    }


def _partitions(n: int, maximum: int | None = None):
    if n == 0:
        yield ()
        return
    maximum = n if maximum is None else min(maximum, n)
    for first in range(maximum, 0, -1):
        for rest in _partitions(n - first, first):
            yield (first,) + rest


def _benchmark_i_signature(signature: tuple[int, ...]) -> dict[str, object]:
    started = perf_counter()
    density = orbit_status_densities(
        signature,
        k=49,
        delta=Fraction(7, 250),
        max_large=6,
        max_offset_count=9,
    )
    return {
        "signature": signature,
        "states": len(density),
        "seconds": perf_counter() - started,
    }


def benchmark_d21(*, scan_all: bool = True, workers: int = 1) -> dict[str, object]:
    """Benchmark I/J truncations on degree-40 D21 product signatures."""
    signatures = [tuple(2 * part for part in lam) for lam in _partitions(20)]
    if not scan_all:
        signatures = [
            (40,),
            (20, 20),
            (10, 8, 6, 4, 4, 2, 2, 2, 2),
            (2,) * 20,
        ]
    started = perf_counter()
    if workers == 1:
        rows = [_benchmark_i_signature(signature) for signature in signatures]
    else:
        # chunksize=1 balances the very uneven signature costs.
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_benchmark_i_signature, signatures, chunksize=1))
    elapsed = perf_counter() - started
    slowest = max(rows, key=lambda row: row["seconds"])
    largest = max(rows, key=lambda row: row["states"])

    # The J common-coordinate truncation for the largest I state case.
    j_started = perf_counter()
    j_density = orbit_status_densities(
        largest["signature"],
        k=48,
        delta=Fraction(7, 250),
        max_large=6,
        max_offset_count=8,
    )
    j_row = {
        "signature": largest["signature"],
        "states": len(j_density),
        "seconds": perf_counter() - j_started,
    }
    return {
        "signatures_scanned": len(rows),
        "workers": workers,
        "total_seconds": elapsed,
        "slowest_i": slowest,
        "largest_i": largest,
        "corresponding_j": j_row,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--selected-only",
        action="store_true",
        help="benchmark four representative signatures instead of all 627 partitions of 20",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=min(16, os.cpu_count() or 1),
        help="processes for the exhaustive benchmark (default: min(16,cpu_count))",
    )
    args = parser.parse_args()
    print({"self_check": self_check()})
    print(
        {
            "benchmark": benchmark_d21(
                scan_all=not args.selected_only, workers=max(1, args.workers)
            )
        }
    )


if __name__ == "__main__":
    main()
