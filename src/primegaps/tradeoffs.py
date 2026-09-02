"""Coupled prime-minorant, support-unlock, and sieve-score tradeoffs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np

from .eigen import GeneralizedEigenResult, solve_generalized_eigenproblem


@dataclass(frozen=True)
class MinorantCandidate:
    """The objective-relevant data for one decomposition choice."""

    candidate_id: str
    retained_mass: float
    pointwise_negative_bound_c2: float
    theorem_status: str
    extra_required_regimes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("candidate_id must be non-empty")
        if not 0.0 <= self.retained_mass <= 1.0:
            raise ValueError("retained_mass must lie in [0,1]")
        if self.pointwise_negative_bound_c2 < 0.0:
            raise ValueError("pointwise_negative_bound_c2 must be nonnegative")
        if self.theorem_status not in {"theorem-backed", "conditional"}:
            raise ValueError("theorem_status must be theorem-backed or conditional")

    @property
    def retained_mass_loss(self) -> float:
        return 1.0 - self.retained_mass


def minorant_objective_matrix(
    J: np.ndarray,
    K: np.ndarray,
    *,
    k: int,
    candidate: MinorantCandidate,
) -> np.ndarray:
    """Build ``k * (retained_mass * J - c2 * K)``.

    ``J`` and ``K`` must use the unscaled integral convention from Section 2.
    """
    if k <= 0:
        raise ValueError("k must be positive")
    left = np.asarray(J, dtype=float)
    right = np.asarray(K, dtype=float)
    if left.shape != right.shape or left.ndim != 2 or left.shape[0] != left.shape[1]:
        raise ValueError("J and K must be square matrices with equal shape")
    return k * (
        candidate.retained_mass * left
        - candidate.pointwise_negative_bound_c2 * right
    )


def optimize_minorant_score(
    I: np.ndarray,
    J: np.ndarray,
    K: np.ndarray,
    *,
    k: int,
    candidate: MinorantCandidate,
    method: str = "dense",
) -> GeneralizedEigenResult:
    """Optimize the complete minorant score over a supplied matrix basis."""
    objective = minorant_objective_matrix(J, K, k=k, candidate=candidate)
    return solve_generalized_eigenproblem(
        np.asarray(I, dtype=float),
        objective,
        method=method,
        exploit_blocks=False,
    )


@dataclass(frozen=True)
class TradeoffMeasurement:
    """One decomposition's optimized support/score tradeoff."""

    candidate_id: str
    retained_mass_loss: float
    newly_unlocked_interactions: int
    support_measure_gain: float
    baseline_score: float
    optimized_score: float
    delta_score: float
    score_standard_error: float
    theorem_status: str
    best_support_id: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def pareto_frontier(
    measurements: Iterable[TradeoffMeasurement],
) -> tuple[TradeoffMeasurement, ...]:
    """Return points not dominated on loss, unlock, support measure, and score."""
    items = tuple(measurements)

    def dominates(left: TradeoffMeasurement, right: TradeoffMeasurement) -> bool:
        weak = (
            left.retained_mass_loss <= right.retained_mass_loss
            and left.newly_unlocked_interactions >= right.newly_unlocked_interactions
            and left.support_measure_gain >= right.support_measure_gain
            and left.delta_score >= right.delta_score
        )
        strict = (
            left.retained_mass_loss < right.retained_mass_loss
            or left.newly_unlocked_interactions > right.newly_unlocked_interactions
            or left.support_measure_gain > right.support_measure_gain
            or left.delta_score > right.delta_score
        )
        return weak and strict

    result = [
        item for item in items
        if not any(other is not item and dominates(other, item) for other in items)
    ]
    return tuple(sorted(result, key=lambda item: (item.retained_mass_loss, -item.delta_score, item.candidate_id)))
