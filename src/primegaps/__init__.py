"""Exploratory reconstruction of Stadlmann's bounded-gap optimization."""

from .coefficients import large_simplex_coefficients, small_cube_coefficients
from .integrals import (
    ExactSupportParameters,
    IntegralMatrices,
    exact_i_entry,
    exact_ijk_matrices,
    exact_j_entry,
    exact_k_entry,
    monomial,
)
from .distribution import (
    DistributionCertificate,
    Minorant,
    RegionCell,
    cells_from_support,
    is_certified,
)
from .symmetric import (
    FactorialMomentTable,
    SparseCertificateEvaluation,
    SparseSymmetricTerm,
    SymmetricBasisTerm,
    SymmetricSimplexMatrices,
    assemble_symmetric_simplex_matrices,
    enlarged_simplex_pairing,
    evaluate_sparse_symmetric_certificate,
    factorial_moment,
    load_sparse_symmetric_terms,
    monomial_symmetric_polynomial,
    simplex_marginal_pairing,
)

__all__ = [
    "DistributionCertificate",
    "Minorant",
    "RegionCell",
    "cells_from_support",
    "is_certified",
    "ExactSupportParameters",
    "FactorialMomentTable",
    "IntegralMatrices",
    "SparseCertificateEvaluation",
    "SparseSymmetricTerm",
    "SymmetricBasisTerm",
    "SymmetricSimplexMatrices",
    "assemble_symmetric_simplex_matrices",
    "enlarged_simplex_pairing",
    "exact_i_entry",
    "exact_ijk_matrices",
    "exact_j_entry",
    "exact_k_entry",
    "evaluate_sparse_symmetric_certificate",
    "factorial_moment",
    "large_simplex_coefficients",
    "load_sparse_symmetric_terms",
    "monomial",
    "monomial_symmetric_polynomial",
    "simplex_marginal_pairing",
    "small_cube_coefficients",
]
