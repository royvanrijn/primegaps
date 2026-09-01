from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
from itertools import combinations, product


def beta_integer(a: int, b: int) -> Fraction:
    """B(a+1,b+1) for non-negative integer exponents."""
    if a < 0 or b < 0:
        raise ValueError("exponents must be non-negative")
    return Fraction(factorial(a) * factorial(b), factorial(a + b + 1))


def small_interval_coefficients(a: int, b: int) -> tuple[Fraction, ...]:
    """Exact coefficients of int_0^delta t^a (1-t)^b dt.

    This is the k=1 base case C_{m,i}(1,(a),b) in Section 5.2.1.  It is
    independent of m as long as 0 < delta < 1.
    """
    if a < 0 or b < 0:
        raise ValueError("exponents must be non-negative")
    coeff = [Fraction(0) for _ in range(a + b + 2)]
    for j in range(b + 1):
        power = a + j + 1
        coeff[power] += Fraction((-1) ** j * comb(b, j), power)
    return tuple(coeff)


def big_interval_coefficients(a: int, b: int) -> tuple[Fraction, ...]:
    """Exact coefficients of int_delta^1 t^a (1-t)^b dt.

    This is the k=1 base case D_{m,i}(1,(a),b).
    """
    coeff = list(-x for x in small_interval_coefficients(a, b))
    coeff[0] += beta_integer(a, b)
    return tuple(coeff)


def evaluate_polynomial(coeff: tuple[Fraction, ...], x: Fraction) -> Fraction:
    """Evaluate an exact coefficient vector using Horner's rule."""
    result = Fraction(0)
    for c in reversed(coeff):
        result = result * x + c
    return result


def _add_term(coeff: list[Fraction], power: int, value: Fraction) -> None:
    if power >= len(coeff):
        coeff.extend(Fraction(0) for _ in range(power + 1 - len(coeff)))
    coeff[power] += value


def _shifted_power_terms(a: int) -> tuple[tuple[int, int, int], ...]:
    """Terms ``(y_power, delta_power, multiplier)`` in ``(y+delta)^a``."""
    return tuple((q, a - q, comb(a, q)) for q in range(a + 1))


def small_cube_coefficients(
    m: int, exponents: tuple[int, ...], b: int
) -> tuple[Fraction, ...]:
    """Return Stadlmann's exact ``C_{m,i}(k,a,b)`` coefficient vector.

    The vector represents the integral over

    ``{0 <= t_j <= delta, sum(t) <= 1}``

    on the chamber ``1/(m+1) < delta <= 1/m``.  It is the fully expanded
    inclusion--exclusion recurrence omitted after Section 5.2.1 of
    arXiv:2608.31126v1.  Values on chamber boundaries agree by continuity.
    """
    if m < 1:
        raise ValueError("m must be positive")
    if b < 0 or any(a < 0 for a in exponents):
        raise ValueError("exponents must be non-negative")
    k = len(exponents)
    degree = k + sum(exponents) + b
    if k == 0:
        return (Fraction(1),) + (Fraction(0),) * degree

    coeff = [Fraction(0) for _ in range(degree + 1)]
    indices = range(k)
    # Inclusion--exclusion removes t_j >= delta.  Only |S| <= m can
    # contribute in the open chamber; the omitted faces have zero volume.
    for size in range(min(k, m) + 1):
        for shifted in combinations(indices, size):
            shifted_set = set(shifted)
            choices = [
                _shifted_power_terms(exponents[j]) if j in shifted_set
                else ((exponents[j], 0, 1),)
                for j in indices
            ]
            sign = -1 if size % 2 else 1
            for selected in product(*choices):
                powers = tuple(term[0] for term in selected)
                delta_power = sum(term[1] for term in selected)
                multiplier = sign
                for term in selected:
                    multiplier *= term[2]
                simplex_power = b + k + sum(powers)
                scale = Fraction(
                    multiplier * factorial(b) * _factorial_product(powers),
                    factorial(simplex_power),
                )
                # (1-size*delta)^simplex_power
                for j in range(simplex_power + 1):
                    _add_term(
                        coeff,
                        delta_power + j,
                        scale * comb(simplex_power, j) * (-size) ** j,
                    )
    return tuple(coeff)


def large_simplex_coefficients(
    m: int, exponents: tuple[int, ...], b: int
) -> tuple[Fraction, ...]:
    """Return Stadlmann's exact ``D_{m,i}(k,a,b)`` coefficient vector.

    The vector represents the integral over

    ``{delta <= t_j <= 1, sum(t) <= 1}``

    on ``1/(m+1) < delta <= 1/m``.  The upper bounds ``t_j <= 1`` are
    redundant.  Shifting every coordinate by ``delta`` leaves a simplex;
    expanding that Dirichlet integral gives the recurrence in closed form.
    """
    if m < 1:
        raise ValueError("m must be positive")
    if b < 0 or any(a < 0 for a in exponents):
        raise ValueError("exponents must be non-negative")
    k = len(exponents)
    degree = k + sum(exponents) + b
    if k == 0:
        return (Fraction(1),) + (Fraction(0),) * degree
    if k > m:
        return (Fraction(0),) * (degree + 1)

    coeff = [Fraction(0) for _ in range(degree + 1)]
    choices = [_shifted_power_terms(a) for a in exponents]
    for selected in product(*choices):
        powers = tuple(term[0] for term in selected)
        delta_power = sum(term[1] for term in selected)
        multiplier = 1
        for term in selected:
            multiplier *= term[2]
        simplex_power = b + k + sum(powers)
        scale = Fraction(
            multiplier * factorial(b) * _factorial_product(powers),
            factorial(simplex_power),
        )
        # (1-k*delta)^simplex_power
        for j in range(simplex_power + 1):
            _add_term(
                coeff,
                delta_power + j,
                scale * comb(simplex_power, j) * (-k) ** j,
            )
    return tuple(coeff)


def _factorial_product(values: tuple[int, ...]) -> int:
    result = 1
    for value in values:
        result *= factorial(value)
    return result
