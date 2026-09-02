"""Finite constraint-relaxation experiments for analytic bottleneck ranking.

The quantities here are counterfactual score gains, not proof certificates and
not infinitesimal KKT multipliers. The theorem oracle itself is never weakened:
we diagnose which stable constraint IDs reject a support, then allow supports
whose only failures have the one ID under study.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import sqrt
from typing import Iterable

from .distribution import (
    ANALYTIC_CONSTRAINT_IDS,
    Minorant,
    support_constraint_failures,
)
from .support import SupportParameters


@dataclass(frozen=True)
class ScoredSupport:
    """One objective evaluation, whether theorem-feasible or not."""

    candidate_id: str
    support: SupportParameters
    score: float
    score_standard_error: float = 0.0

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("candidate_id must be non-empty")
        if self.score_standard_error < 0.0:
            raise ValueError("score_standard_error must be nonnegative")
        self.support.validate()


@dataclass(frozen=True)
class CandidateDiagnostic:
    candidate_id: str
    score: float
    score_standard_error: float
    failed_constraints: tuple[str, ...]


@dataclass(frozen=True)
class ConstraintRelaxation:
    constraint_id: str
    status: str
    rank: int | None
    baseline_candidate_id: str
    relaxed_candidate_id: str
    baseline_score: float
    relaxed_score: float
    delta_score: float
    delta_standard_error_independent: float
    newly_admitted_candidates: int
    changed_optimizer: bool


@dataclass(frozen=True)
class ShadowPriceExperiment:
    baseline_candidate_id: str
    baseline_score: float
    baseline_score_standard_error: float
    constraints: tuple[ConstraintRelaxation, ...]
    candidates: tuple[CandidateDiagnostic, ...]
    caveats: tuple[str, ...] = (
        "Relaxed optima are counterfactual numerical screens, not theorem-backed certificates.",
        "Delta scores are finite one-constraint gains over the supplied search set, not KKT derivatives.",
        "Local P3 failures mean the implemented sufficient witness search failed; they are not impossibility proofs.",
        "The reported combined standard error assumes independent score estimates; use paired errors when available.",
    )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema": "primegaps.constraint-shadow-prices.v1",
            "baseline_candidate_id": self.baseline_candidate_id,
            "baseline_score": self.baseline_score,
            "baseline_score_standard_error": self.baseline_score_standard_error,
            "constraints": [asdict(item) for item in self.constraints],
            "candidates": [asdict(item) for item in self.candidates],
            "caveats": list(self.caveats),
        }


def _best(items: Iterable[ScoredSupport]) -> ScoredSupport:
    try:
        return max(items, key=lambda item: (item.score, item.candidate_id))
    except ValueError as exc:
        raise ValueError("no theorem-feasible baseline candidate was supplied") from exc


def rank_constraint_relaxations(
    candidates: Iterable[ScoredSupport],
    minorant: Minorant,
    constraint_ids: Iterable[str] = ANALYTIC_CONSTRAINT_IDS,
) -> ShadowPriceExperiment:
    """Re-optimize a scored support grid after waiving each constraint once.

    The objective evaluations may come from an expensive external engine. This
    replay is deliberately cheap: for each constraint C_i, it recomputes the
    argmax over candidates with no failures other than C_i.
    """
    candidate_items = tuple(candidates)
    if not candidate_items:
        raise ValueError("at least one scored support is required")
    if len({item.candidate_id for item in candidate_items}) != len(candidate_items):
        raise ValueError("candidate_id values must be unique")

    requested = tuple(constraint_ids)
    unknown = sorted(set(requested) - set(ANALYTIC_CONSTRAINT_IDS))
    if unknown:
        raise ValueError(f"unknown analytic constraints: {', '.join(unknown)}")
    if len(set(requested)) != len(requested):
        raise ValueError("constraint_ids must not contain duplicates")

    failures_by_id: dict[str, frozenset[str]] = {}
    diagnostics: list[CandidateDiagnostic] = []
    for candidate in candidate_items:
        failures = support_constraint_failures(candidate.support, minorant)
        identifiers = frozenset(failure.constraint_id for failure in failures)
        failures_by_id[candidate.candidate_id] = identifiers
        diagnostics.append(
            CandidateDiagnostic(
                candidate.candidate_id,
                candidate.score,
                candidate.score_standard_error,
                tuple(sorted(identifiers)),
            )
        )

    baseline = _best(
        candidate for candidate in candidate_items if not failures_by_id[candidate.candidate_id]
    )
    provisional: list[ConstraintRelaxation] = []
    for identifier in requested:
        allowed = tuple(
            candidate
            for candidate in candidate_items
            if failures_by_id[candidate.candidate_id] <= {identifier}
        )
        relaxed = _best(allowed)
        newly_admitted = sum(
            bool(failures_by_id[candidate.candidate_id])
            for candidate in allowed
        )
        changed = relaxed.candidate_id != baseline.candidate_id
        combined_error = (
            sqrt(baseline.score_standard_error**2 + relaxed.score_standard_error**2)
            if changed
            else 0.0
        )
        provisional.append(
            ConstraintRelaxation(
                constraint_id=identifier,
                status="pending",
                rank=0,
                baseline_candidate_id=baseline.candidate_id,
                relaxed_candidate_id=relaxed.candidate_id,
                baseline_score=baseline.score,
                relaxed_score=relaxed.score,
                delta_score=relaxed.score - baseline.score,
                delta_standard_error_independent=combined_error,
                newly_admitted_candidates=newly_admitted,
                changed_optimizer=changed,
            )
        )

    measured = sorted(
        (item for item in provisional if item.newly_admitted_candidates),
        key=lambda item: (-item.delta_score, item.constraint_id),
    )
    ranks = {item.constraint_id: rank for rank, item in enumerate(measured, 1)}
    ordered = measured + sorted(
        (item for item in provisional if not item.newly_admitted_candidates),
        key=lambda item: item.constraint_id,
    )
    ranked = tuple(
        ConstraintRelaxation(
            **{
                **asdict(item),
                "status": "measured" if item.newly_admitted_candidates else "unprobed",
                "rank": ranks.get(item.constraint_id),
            }
        )
        for item in ordered
    )
    return ShadowPriceExperiment(
        baseline.candidate_id,
        baseline.score,
        baseline.score_standard_error,
        ranked,
        tuple(sorted(diagnostics, key=lambda item: item.candidate_id)),
    )
