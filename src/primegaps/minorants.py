"""Harman/Buchstab prime-minorant parameter and mass calculations.

This module owns no support search and no GPY test-function optimization.  It
describes the convolution regimes a candidate minorant asks an external
equidistribution checker to cover.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log

import numpy as np
from numpy.polynomial.legendre import leggauss


@dataclass(frozen=True)
class HarmanRegimes:
    """The three convolution classes in Stadlmann's Definition 5."""

    xi1: float
    xi2: float
    xi3: float

    @property
    def type_i_smooth_gamma(self) -> tuple[float, float]:
        return (self.xi1, 1.0)

    @property
    def type_ii_gamma(self) -> tuple[float, float]:
        return (self.xi2, 1.0 - self.xi2)

    @property
    def type_iii_each_gamma(self) -> tuple[float, float]:
        return (1.0 - 2.0 * self.xi3, self.xi3)

    @property
    def type_iii_pair_sum_min(self) -> float:
        return 1.0 - self.xi3


def stadlmann_admissible(xi1: float, xi2: float, xi3: float) -> bool:
    """Check the strict inequalities of arXiv:2608.31126v1, Proposition 2."""
    return (
        0.0 < xi1 < 1.0
        and 0.0 < xi2 < 1.0
        and 0.0 < xi3 < 1.0
        and 2.0 * xi1 + 3.0 * xi2 < 2.0
        and xi2 <= xi3
        and xi1 + 9.0 * xi2 < 4.0
        and 2.0 * xi1 + xi2 > 1.0
        and 17.0 * xi2 < 7.0
    )


def stadlmann_xi1_interval(xi2: float) -> tuple[float, float]:
    """Open admissible xi1 interval when xi3 >= xi2 is available."""
    return (
        0.5 * (1.0 - xi2),
        min(1.0 - 1.5 * xi2, 4.0 - 9.0 * xi2),
    )


def regime_frontier(xi2: float, xi1_slack: float) -> HarmanRegimes:
    """Choose xi3=xi2 and xi1 just below its least-demanding upper edge.

    Increasing xi3 only enlarges the Type III class and changes no retained
    mass.  Increasing xi1 shrinks the required Type I class and changes no
    retained mass.  Thus xi3=xi2 and xi1 approaching the open upper endpoint
    describe the regime-demand frontier; ``xi1_slack`` keeps a strict witness.
    """
    lower, upper = stadlmann_xi1_interval(xi2)
    xi1 = upper - xi1_slack
    candidate = HarmanRegimes(xi1=xi1, xi2=xi2, xi3=xi2)
    if not lower < xi1 < upper or not stadlmann_admissible(xi1, xi2, xi2):
        raise ValueError("xi2/slack does not give an interior admissible witness")
    return candidate


@dataclass(frozen=True)
class DiscardVariant:
    """One choice of which nonnegative exceptional pieces to discard."""

    name: str
    retained_mass: float
    c2: int
    extra_required_regimes: tuple[str, ...]


def discard_variants(loss_a: float, loss_b: float) -> tuple[DiscardVariant, ...]:
    """Return the four conditional minorants from two exceptional pieces.

    ``A`` is the ordered five-prime piece, bounded pointwise by 4. ``B`` is
    the reversal five-prime piece, bounded pointwise by 20.  A retained piece
    must itself be equidistributed over the candidate's modulus family.
    """
    if loss_a < 0.0 or loss_b < 0.0:
        raise ValueError("discarded masses must be nonnegative")
    return (
        DiscardVariant("discard-A-and-B", 1.0 - loss_a - loss_b, 24, ()),
        DiscardVariant("discard-A-retain-B", 1.0 - loss_a, 4, ("exception-B",)),
        DiscardVariant("retain-A-discard-B", 1.0 - loss_b, 20, ("exception-A",)),
        DiscardVariant(
            "retain-A-and-B", 1.0, 0, ("exception-A", "exception-B")
        ),
    )


class _Quadrature:
    def __init__(self, order: int):
        if order < 1:
            raise ValueError("quadrature order must be positive")
        self.nodes, self.weights = leggauss(order)

    def interval(self, lower: float, upper: float) -> tuple[np.ndarray, np.ndarray]:
        if upper <= lower:
            return np.empty(0), np.empty(0)
        half = (upper - lower) / 2.0
        return lower + half * (self.nodes + 1.0), half * self.weights


def _reciprocal_pair_integral(lower: float, upper: float, total: float) -> float:
    if upper <= lower:
        return 0.0
    return (
        log(upper)
        - log(total - upper)
        - log(lower)
        + log(total - lower)
    ) / total


def stadlmann_loss_components(xi2: float, order: int = 96) -> tuple[float, float]:
    """Numerically integrate the two losses in Proposition 2.

    The innermost integral is evaluated analytically and the remaining three
    dimensions use deterministic Gauss--Legendre quadrature.  The result is a
    numerical estimate, not a rigorous enclosure.
    """
    if xi2 <= 0.4:
        return (0.0, 0.0)
    if not xi2 < 7.0 / 17.0:
        raise ValueError("the proposition requires 17*xi2 < 7")
    q = _Quadrature(order)
    lower = 1.0 - 2.0 * xi2

    loss_a = 0.0
    alpha1s, w1s = q.interval(lower, 3.0 * xi2 - 1.0)
    for alpha1, w1 in zip(alpha1s, w1s, strict=True):
        alpha2s, w2s = q.interval(lower, min(alpha1, xi2 - alpha1))
        for alpha2, w2 in zip(alpha2s, w2s, strict=True):
            alpha3_lo = max(lower, (1.0 - xi2 - alpha2) / 2.0)
            alpha3s, w3s = q.interval(alpha3_lo, alpha2)
            for alpha3, w3 in zip(alpha3s, w3s, strict=True):
                lo4 = max(lower, 1.0 - xi2 - alpha2 - alpha3)
                hi4 = min(alpha3, (1.0 - alpha1 - alpha2 - alpha3) / 2.0)
                total4 = 1.0 - alpha1 - alpha2 - alpha3
                inner = _reciprocal_pair_integral(lo4, hi4, total4)
                loss_a += w1 * w2 * w3 * inner / (alpha1 * alpha2 * alpha3)

    upper = 8.0 * xi2 - 3.0
    loss_b = 0.0
    alpha2s, w2s = q.interval(lower, upper)
    for alpha2, w2 in zip(alpha2s, w2s, strict=True):
        alpha3s, w3s = q.interval(lower, min(upper, alpha2, xi2 - alpha2))
        for alpha3, w3 in zip(alpha3s, w3s, strict=True):
            alpha4s, w4s = q.interval(alpha3, min(upper, xi2 - alpha2))
            for alpha4, w4 in zip(alpha4s, w4s, strict=True):
                lo5 = max(lower, 1.0 - xi2 - alpha3 - alpha4)
                hi5 = min(upper, (1.0 - alpha2 - alpha3 - alpha4) / 2.0)
                total5 = 1.0 - alpha2 - alpha3 - alpha4
                inner = _reciprocal_pair_integral(lo5, hi5, total5)
                loss_b += w2 * w3 * w4 * inner / (alpha2 * alpha3 * alpha4)
    return (loss_a, loss_b)


def baker_irving_parameters(eta: float) -> tuple[HarmanRegimes, float, float]:
    """Return equivalent regimes, sift exponent, and modulus exponent theta."""
    if not 0.0 < eta < 22.0 / 3295.0:
        raise ValueError("Baker--Irving requires 0 < eta < 22/3295")
    regimes = HarmanRegimes(
        xi1=199.0 / 600.0 + 119.0 * eta / 240.0,
        xi2=0.4 + eta,
        xi3=0.4 + eta,
    )
    beta = 0.2 - 2.0 * eta
    theta = 0.5 + 7.0 / 300.0 + 17.0 * eta / 120.0
    return regimes, beta, theta


def baker_irving_base_loss(eta: float, order: int = 96) -> float:
    """Compute I(E(eta)); the published total loss is 6*I."""
    regimes, _, _ = baker_irving_parameters(eta)
    loss_a, _ = stadlmann_loss_components(regimes.xi2, order)
    return loss_a
