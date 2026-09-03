"""Minimum simultaneous analytic-relaxation experiments.

The numerical scores are supplied by an external variational engine. This
module performs only cheap replay: exact rational slack accounting and a
weighted selection over the supplied, already-scored support candidates.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite
from typing import Iterable, Mapping

from .distribution import (
    ANALYTIC_CONSTRAINT_IDS,
    AnalyticSlack,
    Minorant,
    support_constraint_slacks,
)
from .shadow_prices import ScoredSupport


Q = Fraction


@dataclass(frozen=True)
class BreakthroughDiagnostic:
    candidate_id: str
    score: float
    score_standard_error: float
    score_gate: float
    reaches_target: bool
    weighted_cost: Q
    slacks: tuple[AnalyticSlack, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "candidate_id": self.candidate_id,
            "score": self.score,
            "score_standard_error": self.score_standard_error,
            "score_gate": self.score_gate,
            "reaches_target": self.reaches_target,
            "weighted_cost": str(self.weighted_cost),
            "weighted_cost_float": float(self.weighted_cost),
            "slacks": [item.as_dict() for item in self.slacks],
        }


@dataclass(frozen=True)
class MinimumBreakthroughExperiment:
    target_score: float
    score_standard_error_multiplier: float
    weights: tuple[tuple[str, Q], ...]
    optimum_candidate_id: str | None
    optimum_weighted_cost: Q | None
    candidates: tuple[BreakthroughDiagnostic, ...]
    caveats: tuple[str, ...] = (
        "The optimum is over the supplied scored candidates, not a continuous global optimum.",
        "Local slacks relax the implemented sufficient witness family, not every possible Proposition 3 partition.",
        "Strict analytic inequalities are costed against their closure; zero slack at equality is not a certificate.",
        "Numerical score gates are not rigorous variational certificates.",
    )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema": "primegaps.minimum-breakthrough.v1",
            "target_score": self.target_score,
            "score_standard_error_multiplier": self.score_standard_error_multiplier,
            "weights": {key: str(value) for key, value in self.weights},
            "optimum_candidate_id": self.optimum_candidate_id,
            "optimum_weighted_cost": (
                str(self.optimum_weighted_cost)
                if self.optimum_weighted_cost is not None else None
            ),
            "optimum_weighted_cost_float": (
                float(self.optimum_weighted_cost)
                if self.optimum_weighted_cost is not None else None
            ),
            "candidates": [item.as_dict() for item in self.candidates],
            "caveats": list(self.caveats),
        }


def minimum_breakthrough(
    candidates: Iterable[ScoredSupport],
    minorant: Minorant,
    weights: Mapping[str, int | float | str | Q] | None = None,
    *,
    target_score: float = 1.0,
    score_standard_error_multiplier: float = 0.0,
) -> MinimumBreakthroughExperiment:
    """Minimize weighted simultaneous theorem slack over scored supports."""
    items = tuple(candidates)
    if not items:
        raise ValueError("at least one scored support is required")
    if len({item.candidate_id for item in items}) != len(items):
        raise ValueError("candidate_id values must be unique")
    if not isfinite(target_score):
        raise ValueError("target_score must be finite")
    if score_standard_error_multiplier < 0 or not isfinite(score_standard_error_multiplier):
        raise ValueError(
            "score_standard_error_multiplier must be finite and nonnegative"
        )

    supplied = {} if weights is None else dict(weights)
    unknown = sorted(set(supplied) - set(ANALYTIC_CONSTRAINT_IDS))
    if unknown:
        raise ValueError(f"unknown analytic constraints: {', '.join(unknown)}")
    parsed_weights = {
        identifier: Q(str(supplied.get(identifier, 1)))
        for identifier in ANALYTIC_CONSTRAINT_IDS
    }
    if any(value < 0 for value in parsed_weights.values()):
        raise ValueError("constraint weights must be nonnegative")

    diagnostics = []
    for candidate in items:
        slacks = support_constraint_slacks(candidate.support, minorant)
        cost = sum(parsed_weights[item.constraint_id] * item.slack for item in slacks)
        gate = candidate.score - (
            score_standard_error_multiplier * candidate.score_standard_error
        )
        diagnostics.append(
            BreakthroughDiagnostic(
                candidate.candidate_id,
                candidate.score,
                candidate.score_standard_error,
                gate,
                gate >= target_score,
                cost,
                slacks,
            )
        )
    eligible = [item for item in diagnostics if item.reaches_target]
    optimum = min(
        eligible,
        key=lambda item: (item.weighted_cost, -item.score_gate, item.candidate_id),
        default=None,
    )
    return MinimumBreakthroughExperiment(
        target_score,
        score_standard_error_multiplier,
        tuple(
            (identifier, parsed_weights[identifier])
            for identifier in ANALYTIC_CONSTRAINT_IDS
        ),
        optimum.candidate_id if optimum is not None else None,
        optimum.weighted_cost if optimum is not None else None,
        tuple(sorted(diagnostics, key=lambda item: item.candidate_id)),
    )
