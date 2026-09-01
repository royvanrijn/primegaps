from __future__ import annotations

from fractions import Fraction
from math import comb, factorial


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
