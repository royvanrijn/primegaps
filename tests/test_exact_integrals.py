from fractions import Fraction

import numpy as np

from primegaps.coefficients import (
    evaluate_polynomial,
    large_simplex_coefficients,
    small_cube_coefficients,
)
from primegaps.integrals import (
    ExactSupportParameters,
    exact_ijk_matrices,
    monomial,
)


def _gauss_integral(function, lower: float, upper: float, order: int = 24) -> float:
    if upper <= lower:
        return 0.0
    nodes, weights = np.polynomial.legendre.leggauss(order)
    points = (upper - lower) * 0.5 * nodes + (upper + lower) * 0.5
    return float((upper - lower) * 0.5 * np.dot(weights, [function(x) for x in points]))


def _value(polynomial, point) -> float:
    return sum(
        float(coefficient) * np.prod([x**power for x, power in zip(point, exponents)])
        for exponents, coefficient in polynomial.items()
    )


def _full_simplex_support() -> ExactSupportParameters:
    # T_2 is exactly {t1,t2 >= 0, t1+t2 <= 1/2}; B never cuts it.
    return ExactSupportParameters.from_values(
        delta=Fraction(1, 4),
        epsilon=Fraction(1, 10),
        A=(Fraction(-1, 10), Fraction(2, 5)),
        B=((Fraction(1, 2),) * 4,),
    )


def test_exact_ijk_matrix_against_tiny_numerical_quadrature():
    support = _full_simplex_support()
    basis = [
        monomial((0, 0)),
        {(1, 0): Fraction(1), (0, 1): Fraction(1)},
    ]
    matrices = exact_ijk_matrices(basis, support)
    total = 0.5
    common_cut = 0.3

    for row, left in enumerate(basis):
        for column, right in enumerate(basis):
            numeric_i = _gauss_integral(
                lambda x: _gauss_integral(
                    lambda y: _value(left, (x, y)) * _value(right, (x, y)),
                    0.0,
                    total - x,
                ),
                0.0,
                total,
            )
            numeric_j = _gauss_integral(
                lambda common: _gauss_integral(
                    lambda x: _gauss_integral(
                        lambda y: _value(left, (common, x))
                        * _value(right, (common, y)),
                        0.0,
                        total - common,
                    ),
                    0.0,
                    total - common,
                ),
                0.0,
                common_cut,
            )
            numeric_k = _gauss_integral(
                lambda common: _gauss_integral(
                    lambda x: _value(left, (common, x))
                    * _value(right, (common, x)),
                    0.0,
                    total - common,
                ),
                common_cut,
                total,
            )
            assert abs(float(matrices.I[row][column]) - numeric_i) < 2e-12
            assert abs(float(matrices.J[row][column]) - numeric_j) < 2e-12
            assert abs(float(matrices.K[row][column]) - numeric_k) < 2e-12


def test_cd_coefficients_against_tiny_numerical_quadrature():
    exponents = (1, 2)
    b = 1
    integrand = lambda x, y: x * y**2 * (1.0 - x - y)

    small_delta = Fraction(3, 5)  # chamber m=1; the simplex cuts the square
    small_numeric = _gauss_integral(
        lambda x: _gauss_integral(lambda y: integrand(x, y), 0.0, 0.6),
        0.0,
        0.4,
    ) + _gauss_integral(
        lambda x: _gauss_integral(lambda y: integrand(x, y), 0.0, 1.0 - x),
        0.4,
        0.6,
    )
    small_exact = evaluate_polynomial(
        small_cube_coefficients(1, exponents, b), small_delta
    )

    large_delta = Fraction(21, 100)  # chamber m=4
    large_numeric = _gauss_integral(
        lambda x: _gauss_integral(
            lambda y: integrand(x, y), float(large_delta), 1.0 - x
        ),
        float(large_delta),
        1.0 - float(large_delta),
    )
    large_exact = evaluate_polynomial(
        large_simplex_coefficients(4, exponents, b), large_delta
    )

    assert abs(float(small_exact) - small_numeric) < 2e-13
    assert abs(float(large_exact) - large_numeric) < 2e-13


def test_exact_i_reproduces_dirichlet_identity():
    support = _full_simplex_support()
    basis = [monomial((2, 1)), monomial((0, 0))]
    matrices = exact_ijk_matrices(basis, support)
    # int_{t1+t2<=1/2} t1^2 t2 = 2! 1! / 5! * (1/2)^5
    assert matrices.I[0][1] == Fraction(1, 1920)


def test_constant_full_simplex_intermediate_identities():
    support = _full_simplex_support()
    matrices = exact_ijk_matrices([monomial((0, 0))], support)
    assert matrices.I == ((Fraction(1, 8),),)
    assert matrices.J == ((Fraction(39, 1000),),)
    assert matrices.K == ((Fraction(1, 50),),)


def test_nontrivial_large_coordinate_cut_identity():
    support = ExactSupportParameters.from_values(
        delta=Fraction(1, 4),
        epsilon=Fraction(1, 20),
        A=(Fraction(-1, 20), Fraction(9, 20)),
        B=((Fraction(3, 10), Fraction(9, 20), Fraction(9, 20), Fraction(9, 20)),),
    )
    matrices = exact_ijk_matrices([monomial((0, 0))], support)
    # I: [0,1/4]^2 plus two strips with area int_{1/4}^{3/10}(1/2-t)dt.
    assert matrices.I == ((Fraction(17, 200),),)
    # J: slice length is 3/10 up to common=1/5, then 1/2-common up to 3/10.
    assert matrices.J == ((Fraction(73, 3000),),)
    assert matrices.K == ((Fraction(0),),)
