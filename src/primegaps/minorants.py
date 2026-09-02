"""Harman/Buchstab prime-minorant parameter and mass calculations.

This module owns no support search and no GPY test-function optimization.  It
describes the convolution regimes a candidate minorant asks an external
equidistribution checker to cover.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
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


@dataclass(frozen=True)
class RationalInterval:
    """Closed interval with exact rational endpoints."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if self.upper < self.lower:
            raise ValueError("interval upper endpoint is below its lower endpoint")


def log_fraction_enclosure(
    value: Fraction, terms: int = 12
) -> RationalInterval:
    """Rigorously enclose ``log(value)`` using the atanh series.

    All arithmetic, including the geometric remainder bound, is rational.
    ``terms`` is the largest included zero-based series index.
    """
    if value <= 0:
        raise ValueError("logarithm requires a positive value")
    if terms < 0:
        raise ValueError("terms must be nonnegative")
    if value < 1:
        reflected = log_fraction_enclosure(1 / value, terms)
        return RationalInterval(-reflected.upper, -reflected.lower)
    if value == 1:
        return RationalInterval(Fraction(0), Fraction(0))

    z = (value - 1) / (value + 1)
    z_squared = z * z
    power = z
    partial = Fraction(0)
    for index in range(terms + 1):
        partial += power / (2 * index + 1)
        power *= z_squared
    partial *= 2
    remainder = 2 * power / ((2 * terms + 3) * (1 - z_squared))
    return RationalInterval(partial, partial + remainder)


def type_iic_gamma_cutoff(
    support_max: Fraction, delta: Fraction, epsilon: Fraction
) -> Fraction:
    """Upper gamma endpoint of Proposition 3's Case IIc."""
    omega = support_max - Fraction(1, 4)
    return Fraction(1, 3) + 8 * omega + Fraction(7, 3) * delta + 3 * epsilon


def type_iic_middle_high_loss_enclosure(
    gamma_upper: Fraction,
    *,
    xi2: Fraction = Fraction(2, 5),
    parts: int = 128,
    log_terms: int = 12,
) -> RationalInterval:
    r"""Enclose a mandatory positive Type-IIc Buchstab-branch mass.

    This is the high-sum slice of the positive middle term in the Buchstab
    identity, with ``gamma=1-alpha1-alpha2`` and
    ``gamma/2 <= alpha2 <= (1-gamma)/2``.  Here the Buchstab argument lies in
    ``[1,2]``, so its value is exact and the normalized mass is

    ``integral log((2-3g)/g)/(g(1-g)) dg`` from ``xi2`` to ``gamma_upper``.

    The result is a rigorous enclosure obtained from monotone rational
    rectangles and :func:`log_fraction_enclosure`.  It is only a subset of the
    mass that a complete Type-IIc deletion would remove.
    """
    if parts < 1:
        raise ValueError("parts must be positive")
    if log_terms < 0:
        raise ValueError("log_terms must be nonnegative")
    if gamma_upper <= xi2:
        return RationalInterval(Fraction(0), Fraction(0))
    if xi2 < Fraction(2, 5) or gamma_upper >= Fraction(1, 2):
        raise ValueError("enclosure requires 2/5 <= xi2 < gamma_upper < 1/2")

    width = (gamma_upper - xi2) / parts
    lower = Fraction(0)
    upper = Fraction(0)
    for index in range(parts):
        left = xi2 + index * width
        right = left + width
        # (2-3g)/g decreases and g(1-g) increases on [2/5,1/2].
        log_lower = log_fraction_enclosure(
            (2 - 3 * right) / right, log_terms
        ).lower
        log_upper = log_fraction_enclosure(
            (2 - 3 * left) / left, log_terms
        ).upper
        lower += width * log_lower / (right * (1 - right))
        upper += width * log_upper / (left * (1 - left))
    return RationalInterval(lower, upper)


@dataclass(frozen=True)
class OptimisticMinorantScreen:
    """Cheap no-K rejection using a rigorous mandatory-loss enclosure."""

    loss: RationalInterval
    retained_mass_upper: Fraction
    raw_score_cap: Fraction
    optimistic_score_upper: Fraction
    required_raw_score_lower: Fraction
    survives: bool


def optimistic_no_k_screen(
    loss: RationalInterval, raw_score_cap: Fraction
) -> OptimisticMinorantScreen:
    """Bound ``rho * max_F(kJ/I)`` before any candidate-specific I/J/K work.

    For a proof-bearing rejection, ``raw_score_cap`` itself must be a valid
    upper bound for the chosen function space.  A measured score may still be
    supplied for a discovery-stage gate, but does not become rigorous merely
    by passing through this function.
    """
    if raw_score_cap <= 0:
        raise ValueError("raw score cap must be positive")
    if loss.lower < 0 or loss.lower >= 1:
        raise ValueError("loss lower bound must lie in [0,1)")
    retained_mass_upper = 1 - loss.lower
    required = 1 / retained_mass_upper
    optimistic = retained_mass_upper * raw_score_cap
    return OptimisticMinorantScreen(
        loss=loss,
        retained_mass_upper=retained_mass_upper,
        raw_score_cap=raw_score_cap,
        optimistic_score_upper=optimistic,
        required_raw_score_lower=required,
        survives=optimistic >= 1,
    )
