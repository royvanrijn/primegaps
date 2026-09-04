"""Degree-three rough-almost-prime identities and viability bookkeeping.

This module deliberately separates the elementary arithmetic calculation from
the physical-fragment numerical experiment.  If ``m`` is of size ``x`` and
``P^-(m) > x**beta`` with ``beta > 1/4``, then (asymptotically, for fixed
``beta``) ``Omega(m) <= 3``.  On that range the signed factorial-moment
polynomial below is an exact prime indicator.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb, log

import numpy as np


@dataclass(frozen=True)
class RoughFactorConstants:
    """Leading ``x/log(x)`` constants for beta-rough 1-, 2-, and 3-almost-primes."""

    beta: float
    prime: float
    semiprime: float
    triprime: float
    omega_choose_1: float
    omega_choose_2: float
    omega_choose_3: float
    signed_identity: float
    gross_signed_condition: float


@dataclass(frozen=True)
class ParityContributions:
    """Signed degree-three contributions, normalized by the physical mass I."""

    plus_omega: float
    minus_2_choose_2: float
    plus_3_choose_3: float

    @property
    def signed_sum(self) -> float:
        return self.plus_omega + self.minus_2_choose_2 + self.plus_3_choose_3

    @property
    def gross_absolute(self) -> float:
        return (
            abs(self.plus_omega)
            + abs(self.minus_2_choose_2)
            + abs(self.plus_3_choose_3)
        )


def _validate_omega(omega: int) -> None:
    if not isinstance(omega, int) or isinstance(omega, bool) or not 0 <= omega <= 3:
        raise ValueError("omega must be an integer in {0, 1, 2, 3}")


def degree_three_prime_indicator(omega: int) -> int:
    """Return ``1_{omega=1}`` from factorial moments when ``0 <= omega <= 3``."""

    _validate_omega(omega)
    return omega - 2 * comb(omega, 2) + 3 * comb(omega, 3)


def liouville_third_moment_prime_indicator(omega: int) -> int:
    """Return the equivalent parity form ``(1-lambda)/2 - C(omega, 3)``."""

    _validate_omega(omega)
    liouville = -1 if omega % 2 else 1
    return (1 - liouville) // 2 - comb(omega, 3)


def rough_factor_constants(
    beta: float, *, quadrature_order: int = 96
) -> RoughFactorConstants:
    """Compute the beta-rough factorial-count constants through degree three.

    The formula applies for ``1/4 < beta < 1/3``.  It is the scale-invariant
    prime-factor simplex calculation for integers in a fixed dyadic interval.
    The three-almost-prime integral is evaluated by tensor Gauss--Legendre
    quadrature; it is only used as numerical viability bookkeeping.
    """

    beta = float(beta)
    if not 0.25 < beta < 1.0 / 3.0:
        raise ValueError("beta must satisfy 1/4 < beta < 1/3")
    if not isinstance(quadrature_order, int) or quadrature_order < 8:
        raise ValueError("quadrature_order must be an integer at least 8")

    semiprime = log((1.0 - beta) / beta)
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)

    u_left = beta
    u_right = 1.0 - 2.0 * beta
    u_mid = 0.5 * (u_left + u_right)
    u_half_width = 0.5 * (u_right - u_left)
    triprime_integral = 0.0
    for node_u, weight_u in zip(nodes, weights, strict=True):
        u = u_mid + u_half_width * node_u
        v_left = beta
        v_right = 1.0 - beta - u
        v_mid = 0.5 * (v_left + v_right)
        v_half_width = 0.5 * (v_right - v_left)
        v = v_mid + v_half_width * nodes
        integrand = 1.0 / (u * v * (1.0 - u - v))
        triprime_integral += (
            weight_u * u_half_width * v_half_width * float(np.dot(weights, integrand))
        )
    triprime = triprime_integral / 6.0

    omega_choose_1 = 1.0 + 2.0 * semiprime + 3.0 * triprime
    omega_choose_2 = semiprime + 3.0 * triprime
    omega_choose_3 = triprime
    signed_identity = (
        omega_choose_1 - 2.0 * omega_choose_2 + 3.0 * omega_choose_3
    )
    gross_signed_condition = (
        omega_choose_1 + 2.0 * omega_choose_2 + 3.0 * omega_choose_3
    )
    return RoughFactorConstants(
        beta=beta,
        prime=1.0,
        semiprime=semiprime,
        triprime=triprime,
        omega_choose_1=omega_choose_1,
        omega_choose_2=omega_choose_2,
        omega_choose_3=omega_choose_3,
        signed_identity=signed_identity,
        gross_signed_condition=gross_signed_condition,
    )


def parity_contributions(
    physical_score: float, constants: RoughFactorConstants
) -> ParityContributions:
    """Scale the three factorial-count terms by a full-face physical score."""

    score = float(physical_score)
    return ParityContributions(
        plus_omega=score * constants.omega_choose_1,
        minus_2_choose_2=-2.0 * score * constants.omega_choose_2,
        plus_3_choose_3=3.0 * score * constants.omega_choose_3,
    )


def parity_error_budget(contributions: ParityContributions) -> tuple[float, float]:
    """Return the absolute-I and common unsigned-relative budgets for crossing 1."""

    eta = contributions.signed_sum - 1.0
    return eta, eta / contributions.gross_absolute
