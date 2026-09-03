"""Executable distribution-region certificates from Propositions 2 and 3.

The terminology and equation labels in this module refer to arXiv:2608.31126v1.
All arithmetic used to accept a certificate is rational.  In particular, a
floating-point rounding accident can never turn an uncertified cell pair into a
certified one.

``is_certified`` is deliberately a *sound separation oracle*: ``True`` means
that the returned finite witness proves every tuple in the continuous Xi cell
is covered.  ``False`` means that the implemented theorem witnesses did not
prove coverage; it is not a claim that no more elaborate partition can exist.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from typing import Iterable


Q = Fraction
PAPER_EPSILON = Q(1, 10**10)


@dataclass(frozen=True)
class AnalyticConstraint:
    """Stable identifier and provenance for one relaxable theorem condition.

    These identifiers are also used by the numerical shadow-price experiment.
    Relaxing one of them is a counterfactual calculation and never produces a
    distribution certificate.
    """

    identifier: str
    description: str
    source: str
    diagnostic_kind: str = "necessary-inequality"


ANALYTIC_CONSTRAINTS = (
    AnalyticConstraint("P2.1", "2*xi1 + 3*xi2 < 2", "Proposition 2"),
    AnalyticConstraint("P2.2", "xi2 <= xi3", "Proposition 2"),
    AnalyticConstraint("P2.3", "xi1 + 9*xi2 < 4", "Proposition 2"),
    AnalyticConstraint("P2.4", "2*xi1 + xi2 > 1", "Proposition 2"),
    AnalyticConstraint("P2.5", "17*xi2 < 7", "Proposition 2"),
    AnalyticConstraint("P2.domain.xi1.lower", "xi1 > 0", "Proposition 2"),
    AnalyticConstraint("P2.domain.xi1.upper", "xi1 < 1", "Proposition 2"),
    AnalyticConstraint("P2.domain.xi2.lower", "xi2 > 0", "Proposition 2"),
    AnalyticConstraint("P2.domain.xi2.upper", "xi2 < 1", "Proposition 2"),
    AnalyticConstraint("P2.domain.xi3.lower", "xi3 > 0", "Proposition 2"),
    AnalyticConstraint("P2.domain.xi3.upper", "xi3 < 1", "Proposition 2"),
    AnalyticConstraint("P3.I", "global Type I inequality", "Proposition 3 (I)"),
    AnalyticConstraint("P3.II.range", "global Type II range inequality", "Proposition 3 (II)"),
    AnalyticConstraint("P3.II.delta", "global Type II delta inequality", "Proposition 3 (II)"),
    AnalyticConstraint("P3.III", "global Type III inequality", "Proposition 3 (III)"),
    AnalyticConstraint(
        "P3.local.A",
        "universal Type I partition witness",
        "Proposition 3 (A)",
        "sufficient-witness-search",
    ),
    AnalyticConstraint(
        "P3.local.B",
        "universal Type IIa partition witness",
        "Proposition 3 (B)",
        "sufficient-witness-search",
    ),
    AnalyticConstraint(
        "P3.local.C",
        "universal Type IIb partition witness",
        "Proposition 3 (C)",
        "sufficient-witness-search",
    ),
    AnalyticConstraint(
        "P3.local.D",
        "universal Type IIc partition witness",
        "Proposition 3 (D)",
        "sufficient-witness-search",
    ),
    AnalyticConstraint(
        "P3.local.E",
        "universal Type III partition witness",
        "Proposition 3 (E)",
        "sufficient-witness-search",
    ),
)
ANALYTIC_CONSTRAINT_IDS = tuple(item.identifier for item in ANALYTIC_CONSTRAINTS)


def _q(value: int | float | str | Q) -> Q:
    """Convert public decimal inputs without importing their binary error."""
    if isinstance(value, Q):
        return value
    if isinstance(value, float):
        return Q(str(value))
    return Q(value)


def _show(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class RegionCell:
    """One ``(j,m)`` support cell from Definition 1 / Definition 4.

    ``a_upper`` is ``A_j`` (not the epsilon-enlarged endpoint), ``large_count``
    is ``m``, and ``large_sum_bound`` is ``B[j,m]``.  ``support_max`` must be
    the common ``A_n`` of the full support when this cell is part of a larger
    support; it defaults to ``a_upper`` for a one-band support.
    """

    a_upper: Q | int | float | str
    large_count: int
    large_sum_bound: Q | int | float | str
    delta: Q | int | float | str
    support_max: Q | int | float | str | None = None
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "a_upper", _q(self.a_upper))
        object.__setattr__(self, "large_sum_bound", _q(self.large_sum_bound))
        object.__setattr__(self, "delta", _q(self.delta))
        maximum = self.a_upper if self.support_max is None else _q(self.support_max)
        object.__setattr__(self, "support_max", maximum)
        if self.large_count < 0:
            raise ValueError("large_count must be nonnegative")
        if self.delta <= 0:
            raise ValueError("delta must be positive")
        if self.large_sum_bound < 0:
            raise ValueError("large_sum_bound must be nonnegative")
        if self.large_count == 0 and self.large_sum_bound != 0:
            raise ValueError("the m=0 cell must use large_sum_bound=0")
        if maximum < self.a_upper:
            raise ValueError("support_max must be at least a_upper")

    @property
    def is_empty(self) -> bool:
        """Whether its Xi block is empty (endpoints only strengthen this test)."""
        return self.large_count * self.delta > self.large_sum_bound


@dataclass(frozen=True)
class Minorant:
    """Proposition 2 parameters defining the Harman decomposition/minorant.

    When ``xi2 <= 2/5``, Proposition 2 says the constructed minorant is exactly
    the prime indicator.  Larger legal ``xi2`` values select the genuine Harman
    minorant and have ``c2=24``.
    """

    xi1: Q | int | float | str
    xi2: Q | int | float | str
    xi3: Q | int | float | str

    def __post_init__(self) -> None:
        object.__setattr__(self, "xi1", _q(self.xi1))
        object.__setattr__(self, "xi2", _q(self.xi2))
        object.__setattr__(self, "xi3", _q(self.xi3))

    @property
    def kind(self) -> str:
        return "prime-indicator" if self.xi2 <= Q(2, 5) else "harman-minorant"

    @property
    def c2(self) -> int:
        return 0 if self.xi2 <= Q(2, 5) else 24


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    statement: str
    source: str


@dataclass(frozen=True)
class ConstraintFailure:
    """One failed analytic check, for diagnostics and counterfactual screens."""

    constraint_id: str
    left_cell: str
    right_cell: str
    detail: str
    diagnostic_kind: str


@dataclass(frozen=True)
class AnalyticSlack:
    """Least additive exponent slack for one registered analytic condition.

    Global and Proposition 2 inequalities use their closed-boundary shortfall.
    For a local Proposition 3 condition, the slack is added to every displayed
    bin capacity and minimized over the implemented universal partition
    witnesses. It therefore measures this executable sufficient witness
    family, not every possible use of Proposition 3.
    """

    constraint_id: str
    slack: Q
    diagnostic_kind: str
    worst_left_cell: str | None = None
    worst_right_cell: str | None = None

    def as_dict(self) -> dict[str, object]:
        return {
            "constraint_id": self.constraint_id,
            "slack": _show(self.slack),
            "slack_float": float(self.slack),
            "diagnostic_kind": self.diagnostic_kind,
            "worst_left_cell": self.worst_left_cell,
            "worst_right_cell": self.worst_right_cell,
        }


@dataclass(frozen=True)
class PartitionWitness:
    condition: str
    capacities: tuple[Q, ...]
    primary_bin: int | None
    side_bin_order: tuple[int, ...]
    side_counts_a: tuple[int, ...]
    side_counts_b: tuple[int, ...]
    worst_bin_sums: tuple[Q, ...]
    explanation: str

    def as_dict(self) -> dict[str, object]:
        return {
            "condition": self.condition,
            "capacities": [_show(x) for x in self.capacities],
            "primary_bin": self.primary_bin,
            "side_bin_order": list(self.side_bin_order),
            "side_counts_a": list(self.side_counts_a),
            "side_counts_b": list(self.side_counts_b),
            "worst_bin_sums": [_show(x) for x in self.worst_bin_sums],
            "explanation": self.explanation,
        }


@dataclass(frozen=True)
class DistributionCertificate:
    certified: bool
    theorem: str | None
    reason: str
    modulus_exponent_bound: Q
    omega: Q
    minorant_kind: str
    checks: tuple[Check, ...] = ()
    partitions: tuple[PartitionWitness, ...] = ()
    caveats: tuple[str, ...] = ()

    def __bool__(self) -> bool:
        return self.certified

    def as_dict(self) -> dict[str, object]:
        return {
            "certified": self.certified,
            "theorem": self.theorem,
            "reason": self.reason,
            "modulus_exponent_bound": _show(self.modulus_exponent_bound),
            "omega": _show(self.omega),
            "minorant_kind": self.minorant_kind,
            "checks": [
                {
                    "name": check.name,
                    "passed": check.passed,
                    "statement": check.statement,
                    "source": check.source,
                }
                for check in self.checks
            ],
            "partitions": [witness.as_dict() for witness in self.partitions],
            "caveats": list(self.caveats),
        }


def _comparison(name: str, left: Q, relation: str, right: Q, source: str) -> Check:
    operations = {
        "<": left < right,
        "<=": left <= right,
        ">": left > right,
        ">=": left >= right,
    }
    passed = operations[relation]
    return Check(name, passed, f"{_show(left)} {relation} {_show(right)}", source)


def _proposition_2_checks(minorant: Minorant) -> tuple[Check, ...]:
    x1, x2, x3 = minorant.xi1, minorant.xi2, minorant.xi3
    source = "Proposition 2 (prop:harman)"
    checks = [
        _comparison("P2.1", 2 * x1 + 3 * x2, "<", Q(2), source),
        _comparison("P2.2", x2, "<=", x3, source),
        _comparison("P2.3", x1 + 9 * x2, "<", Q(4), source),
        _comparison("P2.4", 2 * x1 + x2, ">", Q(1), source),
        _comparison("P2.5", 17 * x2, "<", Q(7), source),
    ]
    for index, value in enumerate((x1, x2, x3), 1):
        checks.append(_comparison(f"P2.domain.xi{index}.lower", value, ">", Q(0), source))
        checks.append(_comparison(f"P2.domain.xi{index}.upper", value, "<", Q(1), source))
    return tuple(checks)


def _proposition_3_global_checks(
    delta: Q, support_max: Q, minorant: Minorant
) -> tuple[Check, ...]:
    e = PAPER_EPSILON
    x1, x2, x3 = minorant.xi1, minorant.xi2, minorant.xi3
    source = "Proposition 3 (prop:tupleconditions), global conditions (I)-(III)"
    return (
        _comparison(
            "P3.I",
            min(x1 - 4 * support_max + Q(2, 3), Q(9, 7) - Q(34, 7) * support_max) - 2 * e,
            ">",
            delta,
            source,
        ),
        _comparison("P3.II.range", Q(19, 2) - 36 * support_max - 13 * delta + 100 * e, ">=", Q(0), source),
        _comparison(
            "P3.II.delta",
            min(x2 / 10 - 16 * support_max / 5 + Q(4, 5), x2 / 4 + Q(11, 16) - 3 * support_max) - 2 * e,
            ">=",
            delta,
            source,
        ),
        _comparison(
            "P3.III",
            Q(11, 8) - Q(7, 2) * support_max - Q(9, 8) * x3 - 2 * e,
            ">",
            delta,
            source,
        ),
    )


def _count_vectors(count: int, bins: int) -> Iterable[tuple[int, ...]]:
    for values in product(range(count + 1), repeat=bins):
        if sum(values) <= count:
            yield values


def _group_bin_bounds(
    count: int,
    total_bound: Q,
    delta: Q,
    side_counts: tuple[int, ...],
) -> tuple[tuple[Q, ...], Q]:
    """Bounds from successively taking blocks of the smallest remaining items."""
    if count == 0:
        return (tuple(Q(0) for _ in side_counts), Q(0))
    assigned = 0
    side_bounds: list[Q] = []
    for take in side_counts:
        remaining = count - assigned
        if take == 0:
            side_bounds.append(Q(0))
        else:
            # The ``take`` smallest of ``remaining`` items carry at most this
            # fraction of their remaining total mass.
            remaining_mass = total_bound - assigned * delta
            side_bounds.append(take * remaining_mass / remaining)
        assigned += take
    primary = Q(0) if assigned == count else total_bound - assigned * delta
    return tuple(side_bounds), primary


def _partition_witness(
    condition: str,
    capacities: tuple[Q, ...],
    a: RegionCell,
    b: RegionCell,
) -> PartitionWitness | None:
    """Find a universal order-statistic partition witness for a Xi product.

    Within each of Xi's two blocks, sort the items.  Some consecutive blocks of
    the smallest remaining items go to side bins and all remaining items go to
    one primary bin.  The average of the smallest ``r`` items is at most the
    average of the remaining block, yielding exact upper bounds valid for every
    point of Xi.
    """
    if any(capacity < 0 for capacity in capacities):
        return None
    bin_count = len(capacities)
    vectors_a = tuple(_count_vectors(a.large_count, bin_count - 1))
    vectors_b = tuple(_count_vectors(b.large_count, bin_count - 1))
    for primary in range(bin_count):
        side_bins = tuple(index for index in range(bin_count) if index != primary)
        for order in permutations(side_bins):
            for counts_a in vectors_a:
                bounds_a, primary_a = _group_bin_bounds(
                    a.large_count, a.large_sum_bound, a.delta, counts_a
                )
                for counts_b in vectors_b:
                    bounds_b, primary_b = _group_bin_bounds(
                        b.large_count, b.large_sum_bound, b.delta, counts_b
                    )
                    worst = [Q(0)] * bin_count
                    worst[primary] = primary_a + primary_b
                    for position, target in enumerate(order):
                        worst[target] = bounds_a[position] + bounds_b[position]
                    if all(mass <= capacity for mass, capacity in zip(worst, capacities)):
                        explanation = (
                            "For each Xi block, sort the large-factor exponents; send "
                            "the stated consecutive blocks of smallest remaining items "
                            "to the side bins and the remainder to the primary bin. "
                            "The displayed worst-bin sums are exact rational upper bounds."
                        )
                        return PartitionWitness(
                            condition,
                            capacities,
                            primary,
                            order,
                            counts_a,
                            counts_b,
                            tuple(worst),
                            explanation,
                        )
    return None


def _partition_uniform_slack(
    capacities: tuple[Q, ...],
    a: RegionCell,
    b: RegionCell,
) -> Q:
    """Least common additive capacity slack in the implemented witness family."""
    return _partition_uniform_slack_values(
        capacities,
        a.large_count,
        a.large_sum_bound,
        a.delta,
        b.large_count,
        b.large_sum_bound,
    )


@lru_cache(maxsize=None)
def _partition_uniform_slack_values(
    capacities: tuple[Q, ...],
    count_a: int,
    total_a: Q,
    delta: Q,
    count_b: int,
    total_b: Q,
) -> Q:
    """Cached value-only implementation, independent of cell labels/endpoints."""
    bin_count = len(capacities)
    options_a = _group_bound_options(
        count_a, total_a, delta, bin_count - 1
    )
    options_b = _group_bound_options(
        count_b, total_b, delta, bin_count - 1
    )
    best: Q | None = None
    for primary in range(bin_count):
        side_bins = tuple(index for index in range(bin_count) if index != primary)
        for order in permutations(side_bins):
            for bounds_a, primary_a in options_a:
                for bounds_b, primary_b in options_b:
                    worst = [Q(0)] * bin_count
                    worst[primary] = primary_a + primary_b
                    for position, target in enumerate(order):
                        worst[target] = bounds_a[position] + bounds_b[position]
                    required = max(
                        (mass - capacity for mass, capacity in zip(worst, capacities)),
                        default=Q(0),
                    )
                    required = max(Q(0), required)
                    if best is None or required < best:
                        best = required
    if best is None:  # pragma: no cover - every condition has at least one bin
        raise ValueError("local condition has no partition assignments")
    return best


@lru_cache(maxsize=None)
def _group_bound_options(
    count: int,
    total_bound: Q,
    delta: Q,
    side_bin_count: int,
) -> tuple[tuple[tuple[Q, ...], Q], ...]:
    return tuple(
        _group_bin_bounds(count, total_bound, delta, counts)
        for counts in _count_vectors(count, side_bin_count)
    )


def _local_conditions(
    a: RegionCell, b: RegionCell, minorant: Minorant, omega: Q
) -> tuple[tuple[str, tuple[Q, ...]], ...]:
    e = PAPER_EPSILON
    d = a.delta
    x1, x2, x3 = minorant.xi1, minorant.xi2, minorant.xi3
    conditions: list[tuple[str, tuple[Q, ...]]] = [
        ("A / Type I", (x1 - 2 * e, Q(1, 6) - 4 * omega - 2 * e)),
        (
            "B / Type IIa",
            (Q(2, 5) + Q(24, 5) * omega + Q(7, 5) * d - 2 * e,
             Q(1, 14) - Q(24, 7) * omega - 2 * e),
        ),
        (
            "C / Type IIb",
            (Q(1, 3) + 8 * omega + Q(7, 3) * d - 4 * e,
             Q(1, 10) - Q(34, 5) * omega - Q(7, 5) * d - 4 * e,
             Q(1, 35) + Q(22, 35) * omega + Q(3, 5) * d - 4 * e),
        ),
    ]

    gamma_low = x2 - e
    gamma_high = Q(1, 3) + 8 * omega + Q(7, 3) * d + 3 * e
    if gamma_low <= gamma_high:
        # The negative-omega modulus subrange is BV.  On 0 <= omega_0 <= omega,
        # bins 3 and 4 may be empty, while these componentwise minima for bins 1
        # and 2 hold uniformly over every gamma and omega_0 in Proposition 3(D).
        conditions.append(
            (
                "D / Type IIc (BV for omega_0 <= 0; uniform positive range)",
                (gamma_low - 2 * d - 8 * omega - e,
                 Q(1, 2) - gamma_high - 2 * omega - e,
                 d - e,
                 Q(0)),
            )
        )
    conditions.append(
        (
            "E / Type III",
            (Q(1) - 6 * omega - Q(3, 2) * x3 - 2 * e,
             Q(5, 2) * omega + Q(3, 8) * x3 - 2 * e),
        )
    )
    return tuple(conditions)


def is_certified(
    region_a: RegionCell,
    region_b: RegionCell,
    minorant: Minorant,
) -> DistributionCertificate:
    """Certify a pair of support cells by BV or Propositions 2 and 3.

    The result is truthy exactly when ``result.certified`` is true.  Call
    ``result.as_dict()`` for a stable, serializable theorem certificate.
    """
    modulus_bound = region_a.a_upper + region_b.a_upper
    omega = modulus_bound / 2 - Q(1, 4)
    caveat = (
        "A negative result is conservative: no implemented universal "
        "order-statistic witness was found, but a different partition may exist.",
    )
    if region_a.delta != region_b.delta:
        return DistributionCertificate(False, None, "cell deltas differ", modulus_bound, omega, minorant.kind, caveats=caveat)
    if region_a.support_max != region_b.support_max:
        return DistributionCertificate(False, None, "cells do not name the same full-support A_n", modulus_bound, omega, minorant.kind, caveats=caveat)

    p2_checks = _proposition_2_checks(minorant)
    if not all(check.passed for check in p2_checks):
        return DistributionCertificate(
            False, None, "Proposition 2 does not construct the requested minorant", modulus_bound, omega,
            minorant.kind, checks=p2_checks, caveats=caveat,
        )

    if region_a.is_empty or region_b.is_empty:
        return DistributionCertificate(
            True, "vacuous (empty Xi cell)", "At least one cell has m*delta > B, so it produces no modulus factors.",
            modulus_bound, omega, minorant.kind, checks=p2_checks,
        )

    if modulus_bound <= Q(1, 2):
        return DistributionCertificate(
            True,
            "Bombieri–Vinogradov + Proposition 2",
            "The cell-pair modulus exponent is at most 1/2; BV supplies the distribution input for every Proposition 2 term.",
            modulus_bound,
            omega,
            minorant.kind,
            checks=p2_checks,
        )

    global_checks = _proposition_3_global_checks(region_a.delta, region_a.support_max, minorant)
    checks = p2_checks + global_checks
    if not all(check.passed for check in global_checks):
        return DistributionCertificate(
            False, None, "A global Type I/II/III hypothesis of Proposition 3 fails", modulus_bound, omega,
            minorant.kind, checks=checks, caveats=caveat,
        )

    # Proposition 3 imposes the Xi partition hypotheses only when m+m' > 0.
    if region_a.large_count + region_b.large_count == 0:
        return DistributionCertificate(
            True,
            "Propositions 2 and 3 (fully smooth cell pair)",
            "The global Type I/II/III inequalities hold and Proposition 3 has no local Xi condition when m+m'=0.",
            modulus_bound,
            omega,
            minorant.kind,
            checks=checks,
        )

    witnesses: list[PartitionWitness] = []
    for condition, capacities in _local_conditions(region_a, region_b, minorant, omega):
        witness = _partition_witness(condition, capacities, region_a, region_b)
        if witness is None:
            capacity_text = ", ".join(_show(value) for value in capacities)
            return DistributionCertificate(
                False,
                None,
                f"No universal partition witness found for {condition}; capacities are ({capacity_text})",
                modulus_bound,
                omega,
                minorant.kind,
                checks=checks,
                partitions=tuple(witnesses),
                caveats=caveat,
            )
        witnesses.append(witness)

    return DistributionCertificate(
        True,
        "Propositions 2 and 3 (Type I/IIa/IIb/IIc/III)",
        "Every global inequality holds and each continuous Xi cell has an explicit universal partition witness.",
        modulus_bound,
        omega,
        minorant.kind,
        checks=checks,
        partitions=tuple(witnesses),
        caveats=(
            "For Type IIc, moduli with omega_0 <= 0 are assigned to Bombieri–Vinogradov; Proposition 3(D)'s factorization is checked on 0 <= omega_0 <= omega.",
        ),
    )


def constraint_failures(
    region_a: RegionCell,
    region_b: RegionCell,
    minorant: Minorant,
) -> tuple[ConstraintFailure, ...]:
    """Report every failed relaxable check without weakening ``is_certified``.

    This is a diagnostic API for counterfactual optimization. In particular,
    absence of a local partition witness is only a failure of the implemented
    sufficient witness family, not a proof that Proposition 3 cannot apply.
    Structural mismatches such as unequal deltas are rejected with
    ``ValueError`` because they are not analytic constraints one may relax.
    """
    if region_a.delta != region_b.delta:
        raise ValueError("cell deltas differ")
    if region_a.support_max != region_b.support_max:
        raise ValueError("cells do not name the same full-support A_n")

    labels = (region_a.label or "unlabelled", region_b.label or "unlabelled")
    kinds = {item.identifier: item.diagnostic_kind for item in ANALYTIC_CONSTRAINTS}
    failures: list[ConstraintFailure] = []

    def add(identifier: str, detail: str) -> None:
        failures.append(ConstraintFailure(identifier, *labels, detail, kinds[identifier]))

    for check in _proposition_2_checks(minorant):
        if not check.passed:
            add(check.name, check.statement)

    # Empty Xi cells impose no modulus-distribution demand. BV cells require
    # only a valid Proposition 2 minorant.
    if region_a.is_empty or region_b.is_empty:
        return tuple(failures)
    modulus_bound = region_a.a_upper + region_b.a_upper
    if modulus_bound <= Q(1, 2):
        return tuple(failures)

    for check in _proposition_3_global_checks(region_a.delta, region_a.support_max, minorant):
        if not check.passed:
            add(check.name, check.statement)

    if region_a.large_count + region_b.large_count == 0:
        return tuple(failures)
    omega = modulus_bound / 2 - Q(1, 4)
    for condition, capacities in _local_conditions(region_a, region_b, minorant, omega):
        if _partition_witness(condition, capacities, region_a, region_b) is None:
            identifier = f"P3.local.{condition[0]}"
            detail = "capacities=" + ",".join(_show(value) for value in capacities)
            add(identifier, detail)
    return tuple(failures)


def support_constraint_failures(
    parameters: object,
    minorant: Minorant,
) -> tuple[ConstraintFailure, ...]:
    """Aggregate analytic failures over every Cartesian pair of support cells."""
    cells = cells_from_support(parameters)
    failures: list[ConstraintFailure] = []
    for left in cells:
        for right in cells:
            failures.extend(constraint_failures(left, right, minorant))
    return tuple(failures)


def support_constraint_slacks(
    parameters: object,
    minorant: Minorant,
) -> tuple[AnalyticSlack, ...]:
    """Return simultaneous additive slacks for every analytic constraint.

    The maximum shortfall over all nonvacuous support-cell pairs is used for a
    theorem condition, because one relaxed statement must cover the complete
    Cartesian support. Strict inequalities are measured against their closure:
    equality has zero numerical slack but is still not a certificate.
    """
    cells = cells_from_support(parameters)
    registry = {item.identifier: item for item in ANALYTIC_CONSTRAINTS}
    values = {identifier: Q(0) for identifier in ANALYTIC_CONSTRAINT_IDS}
    locations: dict[str, tuple[str | None, str | None]] = {
        identifier: (None, None) for identifier in ANALYTIC_CONSTRAINT_IDS
    }

    def update(identifier: str, slack: Q, left=None, right=None) -> None:
        slack = max(Q(0), slack)
        if slack > values[identifier]:
            values[identifier] = slack
            locations[identifier] = (left, right)

    x1, x2, x3 = minorant.xi1, minorant.xi2, minorant.xi3
    update("P2.1", 2 * x1 + 3 * x2 - 2)
    update("P2.2", x2 - x3)
    update("P2.3", x1 + 9 * x2 - 4)
    update("P2.4", 1 - 2 * x1 - x2)
    update("P2.5", 17 * x2 - 7)
    for index, value in enumerate((x1, x2, x3), 1):
        update(f"P2.domain.xi{index}.lower", -value)
        update(f"P2.domain.xi{index}.upper", value - 1)

    for left_index, left in enumerate(cells):
        for right in cells[left_index:]:
            if left.is_empty or right.is_empty:
                continue
            modulus_bound = left.a_upper + right.a_upper
            if modulus_bound <= Q(1, 2):
                continue
            for check in _proposition_3_global_checks(
                left.delta, left.support_max, minorant
            ):
                left_text, relation, right_text = check.statement.split()
                lhs, rhs = Q(left_text), Q(right_text)
                shortfall = rhs - lhs if relation in (">", ">=") else lhs - rhs
                update(check.name, shortfall, left.label, right.label)

            if left.large_count + right.large_count == 0:
                continue
            omega = modulus_bound / 2 - Q(1, 4)
            for condition, capacities in _local_conditions(left, right, minorant, omega):
                identifier = f"P3.local.{condition[0]}"
                update(
                    identifier,
                    _partition_uniform_slack(capacities, left, right),
                    left.label,
                    right.label,
                )

    return tuple(
        AnalyticSlack(
            identifier,
            values[identifier],
            registry[identifier].diagnostic_kind,
            *locations[identifier],
        )
        for identifier in ANALYTIC_CONSTRAINT_IDS
    )


def cells_from_support(parameters: object) -> tuple[RegionCell, ...]:
    """Expand a validated ``SupportParameters`` object into all ``(j,m)`` cells."""
    parameters.validate()
    cells: list[RegionCell] = []
    support_max = parameters.A[-1]
    for row_index, row in enumerate(parameters.B):
        upper = parameters.A[row_index + 1]
        cells.append(RegionCell(upper, 0, 0, parameters.delta, support_max, f"j={row_index + 1},m=0"))
        for count, bound in enumerate(row, 1):
            cells.append(
                RegionCell(upper, count, bound, parameters.delta, support_max, f"j={row_index + 1},m={count}")
            )
    return tuple(cells)
