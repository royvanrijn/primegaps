from fractions import Fraction

import numpy as np

from primegaps.coefficients import (
    beta_integer,
    big_interval_coefficients,
    evaluate_polynomial,
    small_interval_coefficients,
)
from primegaps.eigen import exact_rayleigh_quotient, largest_generalized_eigenpair


def test_small_big_interval_partition_beta():
    delta = Fraction(7, 250)  # published 0.028
    for a in range(5):
        for b in range(5):
            small = evaluate_polynomial(small_interval_coefficients(a, b), delta)
            big = evaluate_polynomial(big_interval_coefficients(a, b), delta)
            assert small + big == beta_integer(a, b)


def test_known_small_interval_integral():
    # int_0^d t(1-t)^2 dt = d^2/2 - 2d^3/3 + d^4/4
    assert small_interval_coefficients(1, 2) == (
        Fraction(0), Fraction(0), Fraction(1, 2), Fraction(-2, 3), Fraction(1, 4)
    )


def test_generalized_eigen_and_exact_certificate():
    m1 = np.array([[2.0, 0.0], [0.0, 1.0]])
    m2 = np.array([[4.0, 0.0], [0.0, 1.0]])
    value, vector = largest_generalized_eigenpair(m1, m2)
    assert abs(value - 2.0) < 1e-12
    q = exact_rayleigh_quotient(
        [[2, 0], [0, 1]], [[4, 0], [0, 1]], (Fraction(1), Fraction(0))
    )
    assert q == 2
