"""Exact I-side accelerators for the Stadlmann one-band verifier.

This prototype deliberately has no dependency on the frozen verifier. The
caller supplies the rational type and the final radial moment callback. It
therefore supports both stdlib Fraction reference checks and the production
gmpy2.mpq ring.
"""

from __future__ import annotations

from collections import defaultdict
from math import comb, factorial
from typing import Callable, Iterable


State = tuple[int, int, int, int]


def _coerce(value, rational):
    if hasattr(value, "numerator") and hasattr(value, "denominator"):
        return rational(int(value.numerator), int(value.denominator))
    return rational(value)


def shifted_terms(exponent: int, delta, rational) -> tuple[tuple[int, object], ...]:
    """Terms of (delta+y)^exponent including residual factorials."""
    if exponent < 0:
        raise ValueError("exponents must be non-negative")
    delta = _coerce(delta, rational)
    numerator = factorial(exponent)
    return tuple(
        (
            power,
            rational(numerator, factorial(exponent - power))
            * delta ** (exponent - power),
        )
        for power in range(exponent + 1)
    )


def positive_status_density(
    signature: Iterable[int],
    *,
    delta,
    max_large: int,
    max_offset_count: int,
    rational,
) -> dict[State, object]:
    """Process only positive exponents in the status-density recurrence.

    The result is independent of k. It can be reused for several nearby
    dimensions and then completed with apply_zero_block.
    """
    exponents = tuple(int(value) for value in signature)
    if any(value <= 0 for value in exponents):
        raise ValueError("signature must contain positive exponents only")
    if max_large < 0 or max_offset_count < 0:
        raise ValueError("truncation bounds must be non-negative")
    max_large = min(max_large, max_offset_count)
    delta = _coerce(delta, rational)
    if not rational(0) < delta < rational(1):
        raise ValueError("delta must be in (0,1)")

    shifted = {
        exponent: shifted_terms(exponent, delta, rational)
        for exponent in set(exponents)
    }
    density: dict[State, object] = {(0, 0, 0, 0): rational(1)}
    for exponent in exponents:
        next_density = defaultdict(rational)
        unshifted_small = factorial(exponent)
        for (large, shifted_small, p, q), coefficient in density.items():
            next_density[(large, shifted_small, p, q + exponent)] += (
                coefficient * unshifted_small
            )
            if large < max_large and large + shifted_small < max_offset_count:
                for power, multiplier in shifted[exponent]:
                    next_density[(large + 1, shifted_small, p + power, q)] += (
                        coefficient * multiplier
                    )
            if large + shifted_small < max_offset_count:
                for power, multiplier in shifted[exponent]:
                    next_density[(large, shifted_small + 1, p, q + power)] -= (
                        coefficient * multiplier
                    )
        density = {state: value for state, value in next_density.items() if value}
    return density


def apply_zero_block(
    positive_density: dict[State, object],
    zero_count: int,
    *,
    max_large: int,
    max_offset_count: int,
    rational,
) -> dict[State, object]:
    """Apply all zero-exponent coordinate transitions in one convolution.

    Choosing a large zeroes and b inclusion-exclusion faces contributes
    binom(z,a) binom(z-a,b) (-1)^b. This is exactly the multinomial
    expansion of the repeated zero-coordinate transition.
    """
    if zero_count < 0:
        raise ValueError("zero_count must be non-negative")
    answer = defaultdict(rational)
    for (large, shifted_small, p, q), coefficient in positive_density.items():
        maximum_large_zeroes = min(zero_count, max_large - large)
        for added_large in range(maximum_large_zeroes + 1):
            remaining = zero_count - added_large
            maximum_shifted = min(
                remaining,
                max_offset_count - large - shifted_small - added_large,
            )
            if maximum_shifted < 0:
                continue
            choose_large = comb(zero_count, added_large)
            for added_shifted in range(maximum_shifted + 1):
                multiplier = (
                    choose_large
                    * comb(remaining, added_shifted)
                    * (-1) ** added_shifted
                )
                answer[
                    (
                        large + added_large,
                        shifted_small + added_shifted,
                        p,
                        q,
                    )
                ] += coefficient * multiplier
    return {state: value for state, value in answer.items() if value}


def orbit_status_densities(
    signature: Iterable[int],
    *,
    k: int,
    delta,
    max_large: int,
    max_offset_count: int,
    rational,
    positive_cache: dict | None = None,
) -> dict[State, object]:
    """Closed-zero exact replacement for the coordinatewise density DP."""
    signature = tuple(int(value) for value in signature)
    if len(signature) > k:
        raise ValueError("signature is longer than k")
    max_large = min(max_large, k, max_offset_count)
    cache_key = (
        signature,
        str(delta),
        max_large,
        max_offset_count,
        rational.__module__,
        rational.__name__,
    )
    if positive_cache is not None and cache_key in positive_cache:
        positive = positive_cache[cache_key]
    else:
        positive = positive_status_density(
            signature,
            delta=delta,
            max_large=max_large,
            max_offset_count=max_offset_count,
            rational=rational,
        )
        if positive_cache is not None:
            positive_cache[cache_key] = positive
    return apply_zero_block(
        positive,
        k - len(signature),
        max_large=max_large,
        max_offset_count=max_offset_count,
        rational=rational,
    )


def normalized_density_terms(density: dict[State, object], k: int):
    """Yield normalized two-group density powers and coefficients."""
    for (large, shifted_small, p, q), coefficient in density.items():
        if large == 0:
            if p:
                raise AssertionError("absent large group carries power")
            large_power = None
            large_denominator = 1
        else:
            large_power = p + large - 1
            large_denominator = factorial(large_power)
        if large == k:
            if q:
                raise AssertionError("absent small group carries power")
            small_power = None
            small_denominator = 1
        else:
            small_power = q + (k - large) - 1
            small_denominator = factorial(small_power)
        yield (
            large,
            shifted_small,
            large_power,
            small_power,
        ), coefficient / (large_denominator * small_denominator)


def signature_moments(
    signature: tuple[int, ...],
    slacks: Iterable[int],
    *,
    k: int,
    delta,
    total_cap,
    large_caps: tuple,
    rational,
    orbit_size: Callable[[tuple[int, ...], int], int],
    radial_moment: Callable,
    positive_cache: dict | None = None,
):
    """Compute candidate-independent I moments for requested radial powers."""
    slacks = tuple(sorted(set(int(slack) for slack in slacks)))
    if any(slack < 0 for slack in slacks):
        raise ValueError("slack powers must be non-negative")
    density = orbit_status_densities(
        signature,
        k=k,
        delta=delta,
        max_large=len(large_caps),
        max_offset_count=int(total_cap // delta),
        rational=rational,
        positive_cache=positive_cache,
    )
    orbit = orbit_size(signature, k)
    values = {slack: rational(0) for slack in slacks}
    for (
        large,
        shifted_small,
        large_power,
        small_power,
    ), density_coefficient in normalized_density_terms(density, k):
        height = total_cap - (large + shifted_small) * delta
        large_cap = (
            None
            if large == 0
            else large_caps[large - 1] - large * delta
        )
        weight = orbit * density_coefficient
        for slack in slacks:
            values[slack] += weight * radial_moment(
                large_power,
                small_power,
                slack,
                height,
                large_cap,
            )
    return values


def signature_value(
    signature: tuple[int, ...],
    coefficient_by_slack: dict[int, object],
    *,
    k: int,
    delta,
    total_cap,
    large_caps: tuple,
    rational,
    orbit_size: Callable[[tuple[int, ...], int], int],
    radial_moment: Callable,
    positive_cache: dict | None = None,
):
    """Contract one I signature using the closed-zero density."""
    density = orbit_status_densities(
        signature,
        k=k,
        delta=delta,
        max_large=len(large_caps),
        max_offset_count=int(total_cap // delta),
        rational=rational,
        positive_cache=positive_cache,
    )
    orbit = orbit_size(signature, k)
    value = rational(0)
    for (
        large,
        shifted_small,
        large_power,
        small_power,
    ), density_coefficient in normalized_density_terms(density, k):
        height = total_cap - (large + shifted_small) * delta
        large_cap = (
            None
            if large == 0
            else large_caps[large - 1] - large * delta
        )
        for slack, candidate_coefficient in coefficient_by_slack.items():
            value += (
                orbit
                * density_coefficient
                * candidate_coefficient
                * radial_moment(
                    large_power,
                    small_power,
                    slack,
                    height,
                    large_cap,
                )
            )
    return value
