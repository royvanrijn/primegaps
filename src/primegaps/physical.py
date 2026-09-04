"""Exact finite predicates for the PrimeGaps186 physical source geometry.

The numerical integrals live elsewhere.  This module deliberately handles only
the rational source ladders and point-configuration predicates, so comparisons
at strict endpoints do not depend on binary floating point.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Iterable


F = Fraction


@dataclass(frozen=True)
class PhysicalSourceRow:
    ladder: str
    index: int
    order: int
    activation: F
    outer_core: F
    inner_core: F
    outer_threshold: F
    inner_threshold: F
    owner_plateau: F | None


@dataclass(frozen=True)
class FiveHalvesGroup:
    name: str
    role: str
    rows: tuple[PhysicalSourceRow, ...]
    activation: F
    threshold: F
    radial_lower: F
    radial_upper: F
    hard_cap: F
    split: F


def _source_ladder(
    name: str,
    inner_radius: F,
    epsilon: F,
    limit: F,
    affine_rows: tuple[tuple[F, F], ...],
    outer_radius: F,
    rho: F,
    gap: F,
) -> tuple[PhysicalSourceRow, ...]:
    previous = F(0)
    error = rho * (outer_radius + inner_radius) - F(1, 2)
    rows: list[PhysicalSourceRow] = []
    for index in range(100):
        order = min(index // 12 + 1, 3)
        intercept, slope = affine_rows[order - 1]
        omega = min(
            limit,
            (intercept - epsilon - error + 2 * previous - gap) / slope,
        )
        delta = intercept - slope * omega - epsilon
        band = (F(1, 2) + 2 * previous) / rho
        activation = delta / rho
        outer_core = band - inner_radius
        inner_core = band - outer_radius
        eta = (
            activation
            if order < 3
            else (activation + outer_radius + inner_radius - band) / 2
        )
        outer_threshold = outer_core + eta
        inner_threshold = inner_core + eta
        if not previous < omega <= limit or delta <= 0:
            raise ArithmeticError("invalid physical source ladder")
        rows.append(
            PhysicalSourceRow(
                ladder=name,
                index=index,
                order=order,
                activation=activation,
                outer_core=outer_core,
                inner_core=inner_core,
                outer_threshold=outer_threshold,
                inner_threshold=inner_threshold,
                owner_plateau=(
                    None if order < 3 else 23 * inner_threshold / 40
                ),
            )
        )
        if omega == limit:
            return tuple(rows)
        previous = omega
    raise ArithmeticError("physical source ladder did not terminate")


def primegaps186_source_data() -> tuple[
    tuple[PhysicalSourceRow, ...],
    tuple[PhysicalSourceRow, ...],
    tuple[FiveHalvesGroup, ...],
    F,
]:
    """Rebuild the pinned old/new ladders and the three order-5/2 groups."""

    gap, tau = F(1, 10**7), F(1, 10**10)
    rho = F(1, 4) + F(12499, 10**6)
    rho_star = rho - gap
    outer_radius = F(2742997, 10**7) / rho_star
    new_inner_radius = F(251, 1000) / rho_star
    old_inner_radius = 2 - F(3, 1000) - outer_radius
    sigma_old = F(100001, 10**6)
    sigma_new = F(1, 2) - F(40481, 100000) + tau
    old_affine = (
        ((1 - 5 * sigma_old) / 15, F(18, 5)),
        ((1 - 4 * sigma_old) / 16, F(7, 2)),
        (F(3, 80), F(3)),
    )
    new_affine = (
        ((1 - 5 * sigma_new) / 15, F(18, 5)),
        ((1 - 4 * sigma_new) / 16, F(7, 2)),
        ((1 - 2 * sigma_new) / 20, F(16, 5)),
    )
    old = _source_ladder(
        "old",
        old_inner_radius,
        F(1, 10**6),
        F(12499, 10**6),
        old_affine,
        outer_radius,
        rho,
        gap,
    )
    new = _source_ladder(
        "new",
        new_inner_radius,
        F(1, 10**7),
        F(253, 20000),
        new_affine,
        outer_radius,
        rho,
        gap,
    )
    if (len(old), len(new)) != (29, 43):
        raise ArithmeticError("PrimeGaps186 source row inventory changed")

    mesh = outer_radius / 98304
    specifications = (
        (
            "outer_h25",
            "outer",
            old[24:28] + new[24:39],
            "outer_threshold",
            "outer_core",
            98303,
            46580,
            19660,
        ),
        (
            "old_inner_h25",
            "inner",
            old[24:28],
            "inner_threshold",
            "inner_core",
            89563,
            35265,
            17912,
        ),
        (
            "new_inner_h25",
            "inner",
            new[24:39],
            "inner_threshold",
            "inner_core",
            89953,
            35419,
            17990,
        ),
    )
    groups = []
    for (
        name,
        role,
        rows,
        threshold_field,
        core_field,
        radial_upper,
        hard_cap,
        split,
    ) in specifications:
        groups.append(
            FiveHalvesGroup(
                name=name,
                role=role,
                rows=tuple(rows),
                activation=min(row.activation for row in rows),
                threshold=min(getattr(row, threshold_field) for row in rows),
                radial_lower=min(getattr(row, core_field) for row in rows),
                radial_upper=radial_upper * mesh,
                hard_cap=hard_cap * mesh,
                split=split * mesh,
            )
        )
    return old, new, tuple(groups), mesh


def inclusive_obstruction(
    fragments: Iterable[F],
    activation: F,
    allocation: Callable[[F], F],
) -> F:
    """Evaluate the inclusive-prefix obstruction with exact rationals."""

    active = sorted((F(value) for value in fragments if F(value) > activation), reverse=True)
    result = F(0)
    for fragment in set(active):
        prefix = sum((value for value in active if value >= fragment), F(0))
        result = max(result, prefix + allocation(fragment))
    return result


def nonlargest_five_halves_obstruction(
    fragments: Iterable[F], activation: F
) -> F:
    """Evaluate H_{5/2} after excluding the largest-fragment witness."""

    active = sorted((F(value) for value in fragments if F(value) > activation), reverse=True)
    result = F(0)
    for index, fragment in enumerate(active):
        if index == 0:
            continue
        prefix = sum(active[: index + 1], F(0))
        result = max(result, prefix + F(3, 2) * fragment)
    return result


def order_three_row_failure(
    row: PhysicalSourceRow,
    role: str,
    fragments: Iterable[F],
    total_mass: F,
) -> bool:
    """Evaluate the order-three ownership obstruction after opposite-root guards.

    The separate largest-fragment and opposite-root cap comparisons are assumed,
    exactly as in PrimeGaps186 Lemma 1.4.  This function consequently tests the
    remaining nonlinear obstruction, not the already discharged allocation guard.
    """

    if row.order != 3 or row.owner_plateau is None:
        raise ValueError("the row is not an order-three source row")
    if role not in {"outer", "inner"}:
        raise ValueError("role must be 'outer' or 'inner'")
    plateau = row.owner_plateau
    phi_d = lambda value: min(F(3, 2) * value, plateau)
    phi_e = lambda value: 3 * value - phi_d(value)
    if role == "outer":
        return total_mass > row.outer_core and inclusive_obstruction(
            fragments, row.activation, phi_d
        ) > row.outer_threshold
    return total_mass > row.inner_core and inclusive_obstruction(
        fragments, row.activation, phi_e
    ) > row.inner_threshold


def exact_group_failure(
    group: FiveHalvesGroup, fragments: Iterable[F], total_mass: F
) -> bool:
    fragments = tuple(map(F, fragments))
    total_mass = F(total_mass)
    if not group.radial_lower < total_mass <= group.radial_upper:
        return False
    if max(fragments, default=F(0)) > group.hard_cap:
        return False
    return any(
        order_three_row_failure(row, group.role, fragments, total_mass)
        for row in group.rows
    )


def grouped_five_halves_failure(
    group: FiveHalvesGroup, fragments: Iterable[F], total_mass: F
) -> bool:
    fragments = tuple(map(F, fragments))
    total_mass = F(total_mass)
    if not group.radial_lower < total_mass <= group.radial_upper:
        return False
    if max(fragments, default=F(0)) > group.hard_cap:
        return False
    return (
        nonlargest_five_halves_obstruction(fragments, group.activation)
        > group.threshold
    )
