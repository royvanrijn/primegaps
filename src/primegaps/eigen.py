from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import Iterable

import numpy as np


def largest_generalized_eigenpair(m1: np.ndarray, m2: np.ndarray) -> tuple[float, np.ndarray]:
    """Largest eigenpair of M2 c = lambda M1 c for symmetric positive-definite M1.

    This is the floating-point search stage used only to locate a promising c.
    Exact verification should be done with `exact_rayleigh_quotient`.
    """
    m1 = np.asarray(m1, dtype=float)
    m2 = np.asarray(m2, dtype=float)
    if m1.shape != m2.shape or m1.ndim != 2 or m1.shape[0] != m1.shape[1]:
        raise ValueError("M1 and M2 must be square matrices of equal shape")
    l = np.linalg.cholesky(m1)
    # A = L^{-1} M2 L^{-T}, preserving symmetry up to roundoff.
    x = np.linalg.solve(l, m2)
    a = np.linalg.solve(l, x.T).T
    a = (a + a.T) * 0.5
    values, vectors = np.linalg.eigh(a)
    idx = int(np.argmax(values))
    y = vectors[:, idx]
    c = np.linalg.solve(l.T, y)
    return float(values[idx]), c


def rationalize_vector(vector: Iterable[float], max_denominator: int = 1_000_000) -> tuple[Fraction, ...]:
    """Convert a floating eigenvector into a primitive rational direction."""
    fractions = tuple(Fraction(float(v)).limit_denominator(max_denominator) for v in vector)
    return fractions


def exact_quadratic_form(matrix: list[list[Fraction | int]], vector: tuple[Fraction, ...]) -> Fraction:
    n = len(vector)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("matrix/vector dimensions differ")
    return sum(
        vector[i] * Fraction(matrix[i][j]) * vector[j]
        for i in range(n)
        for j in range(n)
    )


def exact_rayleigh_quotient(
    m1: list[list[Fraction | int]],
    m2: list[list[Fraction | int]],
    vector: tuple[Fraction, ...],
) -> Fraction:
    denominator = exact_quadratic_form(m1, vector)
    if denominator == 0:
        raise ZeroDivisionError("c M1 c^T is zero")
    return exact_quadratic_form(m2, vector) / denominator
