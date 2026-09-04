"""Rough-almost-prime identities and viability bookkeeping.

This module deliberately separates the elementary arithmetic calculation from
the physical-fragment numerical experiment.  If ``m`` is of size ``x`` and
``P^-(m) > x**beta`` with ``beta > 1/4``, then (asymptotically, for fixed
``beta``) ``Omega(m) <= 3``.  The bound improves to two when ``beta > 1/3``.
On those ranges the signed factorial-moment polynomials below are exact prime
indicators.
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

    @property
    def detector_degree(self) -> int:
        """Degree of the exact detector in the applicable open beta regime."""

        return 3 if self.beta <= 1.0 / 3.0 else 2

    @property
    def rough_carrier(self) -> float:
        """Leading constant for all beta-rough integers in the dyadic range."""

        return self.prime + self.semiprime + self.triprime


@dataclass(frozen=True)
class ParityContributions:
    """Signed degree-two/three contributions, normalized by physical mass I."""

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


def _validate_omega(omega: int, maximum: int) -> None:
    if (
        not isinstance(omega, int)
        or isinstance(omega, bool)
        or not 0 <= omega <= maximum
    ):
        allowed = ", ".join(map(str, range(maximum + 1)))
        raise ValueError(f"omega must be an integer in {{{allowed}}}")


def degree_two_prime_indicator(omega: int) -> int:
    """Return ``1_{omega=1}`` from factorial moments when ``0 <= omega <= 2``."""

    _validate_omega(omega, 2)
    return omega - 2 * comb(omega, 2)


def degree_three_prime_indicator(omega: int) -> int:
    """Return ``1_{omega=1}`` from factorial moments when ``0 <= omega <= 3``."""

    _validate_omega(omega, 3)
    return omega - 2 * comb(omega, 2) + 3 * comb(omega, 3)


def liouville_second_moment_prime_indicator(omega: int) -> int:
    """Return ``(1-lambda)/2``, exact for ``0 <= omega <= 2``."""

    _validate_omega(omega, 2)
    liouville = -1 if omega % 2 else 1
    return (1 - liouville) // 2


def liouville_third_moment_prime_indicator(omega: int) -> int:
    """Return the equivalent parity form ``(1-lambda)/2 - C(omega, 3)``."""

    _validate_omega(omega, 3)
    liouville = -1 if omega % 2 else 1
    return (1 - liouville) // 2 - comb(omega, 3)


def rough_factor_constants(
    beta: float, *, quadrature_order: int = 96
) -> RoughFactorConstants:
    """Compute beta-rough factorial-count constants in the degree-2/3 regimes.

    The formula applies for ``1/4 < beta < 1/2``.  For ``beta <= 1/3`` it uses
    the scale-invariant three-factor simplex calculation; at the boundary the
    triprime leading constant is zero, although the pointwise degree-two bound
    requires strict ``beta > 1/3``.  Above that boundary the triprime state is
    absent and the exact detector has degree two.  The three-almost-prime
    integral is evaluated by tensor Gauss--Legendre quadrature and is used only
    as numerical viability bookkeeping.
    """

    beta = float(beta)
    if not 0.25 < beta < 0.5:
        raise ValueError("beta must satisfy 1/4 < beta < 1/2")
    if not isinstance(quadrature_order, int) or quadrature_order < 8:
        raise ValueError("quadrature_order must be an integer at least 8")

    semiprime = log((1.0 - beta) / beta)
    triprime_integral = 0.0
    if beta < 1.0 / 3.0:
        nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
        u_left = beta
        u_right = 1.0 - 2.0 * beta
        u_mid = 0.5 * (u_left + u_right)
        u_half_width = 0.5 * (u_right - u_left)
        for node_u, weight_u in zip(nodes, weights, strict=True):
            u = u_mid + u_half_width * node_u
            v_left = beta
            v_right = 1.0 - beta - u
            v_mid = 0.5 * (v_left + v_right)
            v_half_width = 0.5 * (v_right - v_left)
            v = v_mid + v_half_width * nodes
            integrand = 1.0 / (u * v * (1.0 - u - v))
            triprime_integral += (
                weight_u
                * u_half_width
                * v_half_width
                * float(np.dot(weights, integrand))
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
