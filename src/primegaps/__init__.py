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

__all__ = [
    "DistributionCertificate",
    "Minorant",
    "RegionCell",
    "cells_from_support",
    "is_certified",
    "ExactSupportParameters",
    "IntegralMatrices",
    "exact_i_entry",
    "exact_ijk_matrices",
    "exact_j_entry",
    "exact_k_entry",
    "large_simplex_coefficients",
    "monomial",
    "small_cube_coefficients",
]
