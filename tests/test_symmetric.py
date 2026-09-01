from fractions import Fraction
from itertools import product
from math import factorial

from primegaps.integrals import ExactSupportParameters, exact_ijk_matrices
from primegaps.symmetric import (
    FactorialMomentTable,
    SparseSymmetricTerm,
    SymmetricBasisTerm,
    assemble_symmetric_simplex_matrices,
    canonical_signature,
    enlarged_simplex_pairing,
    evaluate_sparse_symmetric_certificate,
    factorial_moment,
    load_sparse_symmetric_terms,
    monomial_symmetric_exponents,
    monomial_symmetric_polynomial,
    simplex_marginal_pairing,
)


def _brute_factorial_moment(dimension, left, right):
    return sum(
        sum(
            _product(factorial(x + y) for x, y in zip(a, b))
            for b in monomial_symmetric_exponents(right, dimension)
        )
        for a in monomial_symmetric_exponents(left, dimension)
    )


def _product(values):
    result = 1
    for value in values:
        result *= value
    return result


def test_canonical_signature_and_duplicate_orbit_expansion():
    assert canonical_signature((0, 2, 1, 2, 0)) == (2, 2, 1)
    orbit = tuple(monomial_symmetric_exponents((2, 2), 3))
    assert orbit == ((2, 2, 0), (2, 0, 2), (0, 2, 2))


def test_overlap_and_nilpotent_factorial_moments_match_brute_force():
    signatures = [(), (1,), (2,), (2, 1), (2, 2)]
    table = FactorialMomentTable(signatures, (2, 3, 4))
    for dimension in (2, 3, 4):
        for left, right in product(signatures, repeat=2):
            if max(len(left), len(right)) > dimension:
                continue
            expected = _brute_factorial_moment(dimension, left, right)
            assert factorial_moment(dimension, left, right) == expected
            assert table.get(dimension, left, right) == expected


def _full_simplex_support():
    return ExactSupportParameters.from_values(
        delta=Fraction(1, 4),
        epsilon=Fraction(1, 10),
        A=(Fraction(-1, 10), Fraction(2, 5)),
        B=((Fraction(1, 2),) * 4,),
    )


def test_symmetric_matrix_assembly_matches_independent_reference_engine():
    basis = [
        SymmetricBasisTerm(0, (2,)),
        SymmetricBasisTerm(0, (1, 1)),
        SymmetricBasisTerm(0, ()),
    ]
    assembled = assemble_symmetric_simplex_matrices(
        basis, dimension=2, inner_radius=Fraction(1, 2), outer_radius=Fraction(3, 10)
    )
    expanded = [monomial_symmetric_polynomial(term.signature, 2) for term in basis]
    reference = exact_ijk_matrices(expanded, _full_simplex_support())
    assert assembled.mass == reference.I
    assert assembled.marginal == reference.J


def test_sparse_cleared_contraction_matches_direct_rational_pairings():
    terms = (
        SparseSymmetricTerm(0, (), 2),
        SparseSymmetricTerm(1, (), -1),
        SparseSymmetricTerm(0, (2,), 1),
    )
    dimension = 4
    epsilon_denominator = 3
    degree_bound = 2
    evaluation = evaluate_sparse_symmetric_certificate(
        terms,
        dimension,
        epsilon_denominator,
        degree_bound=degree_bound,
    )
    epsilon = Fraction(1, epsilon_denominator)
    inner = 1 + epsilon
    outer = 1 - epsilon
    mass = sum(
        left.coefficient
        * right.coefficient
        * enlarged_simplex_pairing(
            dimension,
            inner,
            left.slack_power,
            left.signature,
            right.slack_power,
            right.signature,
        )
        for left, right in product(terms, repeat=2)
    )
    marginal = sum(
        left.coefficient
        * right.coefficient
        * simplex_marginal_pairing(
            dimension,
            inner,
            outer,
            left.slack_power,
            left.signature,
            right.slack_power,
            right.signature,
        )
        for left, right in product(terms, repeat=2)
    )
    common_denominator = (
        epsilon_denominator ** (2 * dimension + 1)
        * factorial(2 * dimension + 1)
        * factorial(degree_bound + 1) ** 2
    )
    assert evaluation.mass == common_denominator * mass
    assert evaluation.marginal == common_denominator * marginal


def test_sparse_term_loader_canonicalizes_and_combines_duplicates(tmp_path):
    path = tmp_path / "terms.json"
    path.write_text("[[1,[0,2],3],[1,[2],-1],[0,[],5]]\n")
    assert load_sparse_symmetric_terms(path) == (
        SparseSymmetricTerm(0, (), 5),
        SparseSymmetricTerm(1, (2,), 2),
    )
