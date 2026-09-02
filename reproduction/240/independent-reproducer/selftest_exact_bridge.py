#!/usr/bin/env python3
"""Independent low-dimensional checks for the compressed exact bridge."""

from __future__ import annotations

from fractions import Fraction
from itertools import product

from primegaps.integrals import ExactSupportParameters, exact_i_entry, exact_j_entry
from primegaps.symmetric import monomial_symmetric_polynomial

import exact_symmetric_verifier as verifier


def support() -> ExactSupportParameters:
    return ExactSupportParameters.from_values(
        delta=verifier.DELTA,
        epsilon=Fraction(3, 400),
        A=(-Fraction(3, 400), Fraction(253, 1000)),
        B=((Fraction(3, 20), Fraction(3, 20)) + (Fraction(17, 100),) * 33,),
    )


def multiply(left, right):
    answer = {}
    for a, ca in left.items():
        for b, cb in right.items():
            exponent = tuple(x + y for x, y in zip(a, b))
            answer[exponent] = answer.get(exponent, Fraction(0)) + ca * cb
    return {key: value for key, value in answer.items() if value}


def radial(dimension: int, slack: int):
    answer = {}
    # (U-sum t)^slack, expanded by repeated multiplication.
    one = {(0,) * dimension: verifier.U}
    for coordinate in range(dimension):
        exponent = [0] * dimension
        exponent[coordinate] = 1
        one[tuple(exponent)] = Fraction(-1)
    answer[(0,) * dimension] = Fraction(1)
    for _ in range(slack):
        answer = multiply(answer, one)
    return answer


def expand(terms, dimension: int):
    answer = {}
    for term in terms:
        shape = monomial_symmetric_polynomial(term.signature, dimension)
        piece = multiply(shape, radial(dimension, term.slack))
        for exponent, coefficient in piece.items():
            answer[exponent] = answer.get(exponent, Fraction(0)) + term.coefficient * coefficient
    return {key: value for key, value in answer.items() if value}


def main() -> None:
    cases = (
        (2, (verifier.Term((), 0, Fraction(2)), verifier.Term((2,), 1, Fraction(-3, 2)))),
        (3, (verifier.Term((), 2, Fraction(5, 3)), verifier.Term((2,), 0, Fraction(-1)),
             verifier.Term((2, 2), 1, Fraction(2, 5)))),
    )
    parameters = support()
    for dimension, terms in cases:
        polynomial = expand(terms, dimension)
        compressed_i = verifier.exact_i_grouped(terms, dimension)
        compressed_j = verifier.exact_j(terms, dimension)
        assert compressed_i == exact_i_entry(polynomial, polynomial, parameters)
        assert compressed_j == exact_j_entry(polynomial, polynomial, parameters, dimension)

    # Product structure constants are independently checked by labeled orbits.
    for dimension, left, right in ((3, (2,), (2,)), (4, (4, 2), (2, 2))):
        expected = multiply(
            monomial_symmetric_polynomial(left, dimension),
            monomial_symmetric_polynomial(right, dimension),
        )
        actual = {}
        for signature, coefficient in verifier.product_signatures(dimension, left, right):
            for exponent in monomial_symmetric_polynomial(signature, dimension):
                actual[exponent] = actual.get(exponent, Fraction(0)) + coefficient
        assert actual == expected
    print("exact bridge self-test passed")


if __name__ == "__main__":
    main()
