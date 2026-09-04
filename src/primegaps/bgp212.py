"""Exact arithmetic model for the preliminary BGP212 support.

This module transcribes the public 2026-09-03 draft *A New Bound for Small
Gaps Between Primes*.  It deliberately separates three levels of evidence:

* the exact Table 3 datum;
* the five symbolic modulus classes and the continuous packing obligations;
* Appendix B's reported slack ledger.

The paper does not publish the 455 polyhedral certificate trees.  Consequently
``packing_problem()`` encodes their quantified statement and expected root
count, but does not claim that those roots have been independently replayed.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement
from typing import Mapping


Q = Fraction


def _q(value: int | str | Q) -> Q:
    return value if isinstance(value, Q) else Q(value)


@dataclass(frozen=True)
class AffineForm:
    """An exact affine expression in named logarithmic exponents."""

    constant: Q = Q(0)
    terms: tuple[tuple[str, Q], ...] = ()

    def __post_init__(self) -> None:
        combined: dict[str, Q] = {}
        for name, coefficient in self.terms:
            combined[name] = combined.get(name, Q(0)) + _q(coefficient)
        object.__setattr__(self, "constant", _q(self.constant))
        object.__setattr__(
            self,
            "terms",
            tuple(sorted((name, value) for name, value in combined.items() if value)),
        )

    def evaluate(self, values: Mapping[str, int | str | Q]) -> Q:
        missing = set(self.variables) - set(values)
        if missing:
            raise KeyError(f"missing affine variables: {sorted(missing)}")
        return self.constant + sum(
            (coefficient * _q(values[name]) for name, coefficient in self.terms),
            Q(0),
        )

    @property
    def variables(self) -> tuple[str, ...]:
        return tuple(name for name, _ in self.terms)

    def as_dict(self) -> dict[str, object]:
        return {
            "constant": str(self.constant),
            "terms": {name: str(coefficient) for name, coefficient in self.terms},
        }


def affine(constant: int | str | Q = 0, **terms: int | str | Q) -> AffineForm:
    return AffineForm(_q(constant), tuple((name, _q(value)) for name, value in terms.items()))


@dataclass(frozen=True)
class LinearConstraint:
    """A named affine comparison against zero."""

    name: str
    form: AffineForm
    relation: str

    def __post_init__(self) -> None:
        if self.relation not in {"<", "<=", ">", ">=", "="}:
            raise ValueError(f"unsupported relation {self.relation!r}")

    def holds(self, values: Mapping[str, int | str | Q]) -> bool:
        value = self.form.evaluate(values)
        return {
            "<": value < 0,
            "<=": value <= 0,
            ">": value > 0,
            ">=": value >= 0,
            "=": value == 0,
        }[self.relation]

    def as_dict(self) -> dict[str, object]:
        return {"name": self.name, "form": self.form.as_dict(), "relation": self.relation}


@dataclass(frozen=True)
class FactorInterval:
    """A factor and its open logarithmic exponent interval."""

    factor: str
    divides: str
    lower: AffineForm
    upper: AffineForm

    def width(self, values: Mapping[str, int | str | Q]) -> Q:
        return self.upper.evaluate(values) - self.lower.evaluate(values)

    def as_dict(self) -> dict[str, object]:
        return {
            "factor": self.factor,
            "divides": self.divides,
            "bounds": "open",
            "lower": self.lower.as_dict(),
            "upper": self.upper.as_dict(),
        }


@dataclass(frozen=True)
class ModulusCase:
    name: str
    domain: tuple[LinearConstraint, ...]
    factor_intervals: tuple[FactorInterval, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "domain": [item.as_dict() for item in self.domain],
            "factor_intervals": [item.as_dict() for item in self.factor_intervals],
        }


@dataclass(frozen=True)
class ModulusClass:
    """One of the five exact Section 5 equidistribution classes."""

    identifier: str
    source: str
    convolution: str
    modulus_upper_exponent: AffineForm
    cases: tuple[ModulusCase, ...]
    analytic_walls: tuple[LinearConstraint, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "identifier": self.identifier,
            "source": self.source,
            "convolution": self.convolution,
            "modulus_upper_exponent": self.modulus_upper_exponent.as_dict(),
            "cases": [case.as_dict() for case in self.cases],
            "analytic_walls": [wall.as_dict() for wall in self.analytic_walls],
        }


def modulus_classes() -> tuple[ModulusClass, ...]:
    """Return exact symbolic transcriptions of Lemmas 5.3--5.7."""

    q_upper = affine("1/2", omega=2)
    gamma_le_half = LinearConstraint("gamma <= 1/2", affine("-1/2", gamma=1), "<=")
    gamma_gt_half = LinearConstraint("gamma > 1/2", affine("-1/2", gamma=1), ">")
    gamma_le_top = LinearConstraint(
        "gamma <= 1/2 + 2 omega + epsilon",
        affine("-1/2", gamma=1, omega=-2, epsilon=-1),
        "<=",
    )
    gamma_gt_top = LinearConstraint(
        "gamma > 1/2 + 2 omega + epsilon",
        affine("-1/2", gamma=1, omega=-2, epsilon=-1),
        ">",
    )

    standard_r = FactorInterval(
        "r",
        "d",
        affine(0, gamma=1, delta=-1, epsilon=-3),
        affine(0, gamma=1, epsilon=-3),
    )
    complementary_r = FactorInterval(
        "r",
        "d",
        affine(1, gamma=-1, delta=-1, epsilon=-3),
        affine(1, gamma=-1, epsilon=-3),
    )

    return (
        ModulusClass(
            "D_I",
            "Lemma 5.3 (Type I)",
            "alpha * beta, beta smooth; alpha Siegel-Walfisz above gamma=1/2",
            q_upper,
            (
                ModulusCase("gamma-at-or-below-half", (gamma_le_half,), (standard_r,)),
                ModulusCase(
                    "gamma-just-above-half",
                    (gamma_gt_half, gamma_le_top),
                    (complementary_r,),
                ),
                ModulusCase("gamma-above-modulus-range", (gamma_gt_top,), ()),
            ),
            (
                LinearConstraint(
                    "3 gamma - 12 omega - 3 delta > 1",
                    affine(-1, gamma=3, omega=-12, delta=-3),
                    ">",
                ),
                LinearConstraint(
                    "68 omega + 14 delta < 1",
                    affine(1, omega=-68, delta=-14),
                    ">",
                ),
            ),
        ),
        ModulusClass(
            "D_IIa",
            "Lemma 5.4 (Type IIa)",
            "alpha * beta, beta Siegel-Walfisz, gamma <= 1/2",
            q_upper,
            (ModulusCase("type-iia", (gamma_le_half,), (standard_r,)),),
            (
                LinearConstraint(
                    "24 omega + 7 delta - 5 gamma < -2",
                    affine(-2, gamma=5, omega=-24, delta=-7),
                    ">",
                ),
                LinearConstraint(
                    "8 omega + 3 delta - gamma < 0",
                    affine(0, gamma=1, omega=-8, delta=-3),
                    ">",
                ),
            ),
        ),
        ModulusClass(
            "D_IIb",
            "Lemma 5.5 (Type IIb)",
            "alpha * beta, beta Siegel-Walfisz, gamma <= 1/2",
            q_upper,
            (
                ModulusCase(
                    "type-iib",
                    (gamma_le_half,),
                    (
                        standard_r,
                        FactorInterval(
                            "u",
                            "d/r",
                            affine("1/2", gamma=-1, omega=-2, epsilon=-6, delta=-1),
                            affine("1/2", gamma=-1, omega=-2, epsilon=-6),
                        ),
                    ),
                ),
            ),
            (
                LinearConstraint(
                    "24 omega + 7 delta - 3 gamma < -1",
                    affine(-1, gamma=3, omega=-24, delta=-7),
                    ">",
                ),
                LinearConstraint(
                    "8 omega + 3 delta - gamma < 0",
                    affine(0, gamma=1, omega=-8, delta=-3),
                    ">",
                ),
            ),
        ),
        ModulusClass(
            "D_IIc",
            "Lemma 5.6 (Type IIc)",
            "alpha * beta, beta Siegel-Walfisz, gamma <= 1/2",
            q_upper,
            (
                ModulusCase(
                    "type-iic",
                    (gamma_le_half,),
                    (
                        standard_r,
                        FactorInterval(
                            "u",
                            "d/r",
                            affine(1, gamma=-1, epsilon=-6, delta=-1, theta=-1),
                            affine(1, gamma=-1, epsilon=-6, theta=-1),
                        ),
                        FactorInterval(
                            "d1",
                            "r",
                            affine(2, rho=2, gamma=-1, epsilon=-52, delta=-1, theta=-4),
                            affine(2, rho=2, gamma=-1, epsilon=-52, theta=-4),
                        ),
                    ),
                ),
            ),
            (
                LinearConstraint(
                    "8 omega + 4 delta + 2 gamma < 1",
                    affine(1, omega=-8, delta=-4, gamma=-2),
                    ">",
                ),
                LinearConstraint(
                    "32 omega + 10 delta - gamma < 0",
                    affine(0, gamma=1, omega=-32, delta=-10),
                    ">",
                ),
                LinearConstraint(
                    "48 omega + 16 delta - 4 gamma < -1",
                    affine(-1, gamma=4, omega=-48, delta=-16),
                    ">",
                ),
            ),
        ),
        ModulusClass(
            "D_III",
            "Lemma 5.7 (Type III)",
            "alpha * psi1 * psi2 * psi3",
            q_upper,
            (
                ModulusCase(
                    "type-iii",
                    (),
                    (
                        FactorInterval(
                            "r",
                            "d",
                            affine("1/3", delta="1/3", omega="-4/3"),
                            affine("1/3", delta="4/3", omega="-4/3"),
                        ),
                    ),
                ),
            ),
            (
                LinearConstraint(
                    "28 omega + 9 gamma + 8 delta < 4",
                    affine(4, omega=-28, gamma=-9, delta=-8),
                    ">",
                ),
            ),
        ),
    )


@dataclass(frozen=True)
class BGP212Parameters:
    k: int
    omega: Q
    a0: Q
    a1: Q
    support_epsilon: Q
    delta: Q
    rough_caps: tuple[Q, ...]
    xi1: Q
    xi2: Q
    xi3: Q
    analytic_epsilon: Q

    @property
    def total_cap(self) -> Q:
        return self.a1 + self.support_epsilon

    @property
    def marginal_cap(self) -> Q:
        return self.a1 - self.support_epsilon

    @property
    def maximum_nonempty_rough_count(self) -> int:
        return max(
            count
            for count, cap in enumerate(self.rough_caps, 1)
            if count * self.delta <= cap
        )

    def rough_cap(self, count: int) -> Q:
        if count == 0:
            return Q(0)
        if count < 0 or count > len(self.rough_caps):
            raise ValueError("rough count outside the Table 3 support")
        return self.rough_caps[count - 1]


def parameters() -> BGP212Parameters:
    explicit = tuple(Q(value, 5000) for value in (777, 794, 875, 917, 953, 983, 1016, 1042, 1063, 1081))
    delta = Q(41, 2500)
    return BGP212Parameters(
        k=45,
        omega=Q(7, 1000),
        a0=Q(-1, 125),
        a1=Q(257, 1000),
        support_epsilon=Q(1, 125),
        delta=delta,
        rough_caps=explicit + (Q(1081, 5000),) * (int(Q(1, 1) / delta) - len(explicit)),
        xi1=Q(19, 50),
        xi2=Q(2, 5),
        xi3=Q(2, 5),
        analytic_epsilon=Q(1, 10**10),
    )


@dataclass(frozen=True)
class ParameterInterval:
    variable: str
    lower: Q
    upper: Q


@dataclass(frozen=True)
class PackingPartitionRequirement:
    name: str
    capacities: tuple[AffineForm, ...]


@dataclass(frozen=True)
class PackingCondition:
    identifier: str
    source: str
    parameter_domains: tuple[ParameterInterval, ...]
    partition_requirements: tuple[PackingPartitionRequirement, ...]


@dataclass(frozen=True)
class RoughProfile:
    left_count: int
    right_count: int
    delta: Q
    left_cap: Q
    right_cap: Q


@dataclass(frozen=True)
class PackingProblem:
    """The quantified Proposition 7.8(A)--(E) continuum problem."""

    delta: Q
    rough_caps: tuple[Q, ...]
    maximum_count: int
    conditions: tuple[PackingCondition, ...]
    reported_successful_roots: int

    @property
    def ordered_positive_count_pairs(self) -> tuple[tuple[int, int], ...]:
        return tuple(combinations_with_replacement(range(1, self.maximum_count + 1), 2))

    @property
    def expected_root_count(self) -> int:
        return len(self.ordered_positive_count_pairs) * len(self.conditions)

    def rough_profile(self, left_count: int, right_count: int) -> RoughProfile:
        if not 0 <= left_count <= self.maximum_count:
            raise ValueError("left rough count is outside the nonempty domain")
        if not 0 <= right_count <= self.maximum_count:
            raise ValueError("right rough count is outside the nonempty domain")
        left = Q(0) if left_count == 0 else self.rough_caps[left_count - 1]
        right = Q(0) if right_count == 0 else self.rough_caps[right_count - 1]
        return RoughProfile(left_count, right_count, self.delta, left, right)


def packing_problem(p: BGP212Parameters | None = None) -> PackingProblem:
    """Encode the full continuous packing statement of Lemma 7.10.

    A condition means: for every rough profile in the corresponding Xi
    polytope, and for every parameter in ``parameter_domains``, there exists a
    set partition satisfying each capacity list.  Condition A has two separate
    existential partitions.  Condition D's partition may depend on both gamma
    and omega0.
    """

    p = parameters() if p is None else p
    e, d, w = p.analytic_epsilon, p.delta, p.omega

    def fixed(value: Q) -> AffineForm:
        return affine(value)

    conditions = (
        PackingCondition(
            "A",
            "Proposition 7.8(A), Type I",
            (),
            (
                PackingPartitionRequirement(
                    "gamma-at-or-below-half",
                    (fixed(p.xi1 - 2 * e), fixed(Q(1, 6) - 4 * w - 2 * e)),
                ),
                PackingPartitionRequirement(
                    "gamma-just-above-half",
                    (fixed(Q(1, 2) - 2 * w - 2 * e), fixed(Q(1, 14) - Q(34, 7) * w - 2 * e)),
                ),
            ),
        ),
        PackingCondition(
            "B",
            "Proposition 7.8(B), Type IIa",
            (),
            (
                PackingPartitionRequirement(
                    "type-iia",
                    (
                        fixed(Q(2, 5) + Q(24, 5) * w + Q(7, 5) * d - 2 * e),
                        fixed(Q(1, 14) - Q(24, 7) * w - 2 * e),
                    ),
                ),
            ),
        ),
        PackingCondition(
            "C",
            "Proposition 7.8(C), Type IIb",
            (),
            (
                PackingPartitionRequirement(
                    "type-iib",
                    (
                        fixed(Q(1, 3) + 8 * w + Q(7, 3) * d - 4 * e),
                        fixed(Q(1, 10) - Q(34, 5) * w - Q(7, 5) * d - 4 * e),
                        fixed(2 * w + d - 4 * e),
                    ),
                ),
            ),
        ),
        PackingCondition(
            "D",
            "Proposition 7.8(D), Type IIc",
            (
                ParameterInterval("gamma", p.xi2 - e, Q(1, 3) + 8 * w + Q(7, 3) * d + 3 * e),
                ParameterInterval("omega0", Q(0), w),
            ),
            (
                PackingPartitionRequirement(
                    "type-iic",
                    (
                        affine(-2 * d - e, gamma=1, omega0=-8),
                        affine(Q(1, 2) - e, gamma=-1, omega0=-2),
                        affine(d - e, omega0=4),
                        affine(0, omega0=8),
                    ),
                ),
            ),
        ),
        PackingCondition(
            "E",
            "Proposition 7.8(E), Type III",
            (),
            (
                PackingPartitionRequirement(
                    "type-iii",
                    (
                        fixed(1 - 6 * w - Q(3, 2) * p.xi3 - 2 * e),
                        fixed(Q(5, 2) * w + Q(3, 8) * p.xi3 - 2 * e),
                    ),
                ),
            ),
        ),
    )
    return PackingProblem(
        p.delta,
        p.rough_caps,
        p.maximum_nonempty_rough_count,
        conditions,
        reported_successful_roots=455,
    )


@dataclass(frozen=True)
class Table6Row:
    identifier: str
    condition: str
    left: Q
    right: Q
    strict: bool
    use: str
    verification: str = "recomputed"

    @property
    def slack(self) -> Q:
        return self.right - self.left


def reported_table6_rows() -> tuple[Table6Row, ...]:
    """Appendix B values exactly as printed, including its apparent typo."""

    rows = (
        ("support.delta", "Support: delta > 0", 0, "41/2500", True, "support"),
        ("support.B1", "Support: B1 > delta", "41/2500", "777/5000", True, "support"),
        ("support.monotone", "Support: B_m <= B_(m+1)", 0, "17/5000", False, "monotonicity"),
        ("support.heredity", "Support: B_(m+1) <= B_m + delta", "81/5000", "82/5000", False, "heredity"),
        ("packing.first-rung", "First-rung cap wall", "777/5000", "778/5000", False, "arithmetic certificate"),
        ("prime.xi2", "Prime decomposition: xi2 <= 2/5", "2/5", "2/5", False, "direct prime"),
        ("prime.xi3", "Prime decomposition: xi3 >= xi2", "2/5", "2/5", False, "Type III split"),
        ("prime.harman1", "Prime decomposition: 2 xi1 + 3 xi2 < 2", "49/25", 2, True, "Harman decomposition"),
        ("prime.harman2", "Prime decomposition: xi1 + 9 xi2 < 4", "199/50", 4, True, "Harman decomposition"),
        ("prime.harman3", "Prime decomposition: 2 xi1 + xi2 > 1", 1, "29/25", True, "Harman decomposition"),
        ("prime.harman4", "Prime decomposition: 17 xi2 < 7", "34/5", 7, True, "Harman decomposition"),
        ("global.I.1", "Type I admissibility, first wall", "41/2500", "279999997/15000000000", True, "Prop. 7.8(I)"),
        ("global.I.2", "Type I admissibility, second wall", "41/2500", "1309999993/35000000000", True, "Prop. 7.8(I)"),
        ("global.II.range", "Type II admissibility, first wall", 0, "69599997/2000000000", False, "Prop. 7.8(II)"),
        ("global.II.cap1", "Type II admissibility, first cap wall", "41/2500", "87999999/5000000000", False, "Prop. 7.8(II)"),
        ("global.II.cap2", "Type II admissibility, second cap wall", "41/2500", "82499999/5000000000", False, "Prop. 7.8(II)"),
        ("global.III", "Type III admissibility", "41/2500", "127499999/5000000000", True, "Prop. 7.8(III)"),
        ("global.II.dominant", "Dominant Type II wall", "3937000001/5000000000", "63/80", False, "arithmetic certificate"),
        ("packing.roots", "Exact packing roots A-E", 455, 455, False, "all continuous packing cells"),
        ("transition.bin3", "Half-level transition, third bin", "11/400000000000", "109/1000000000000", False, "transition"),
        ("transition.bin4", "Half-level transition, fourth bin", "11/200000000000", "11/200000000000", False, "transition"),
    )
    return tuple(Table6Row(key, condition, _q(left), _q(right), strict, use) for key, condition, left, right, strict, use in rows)


def recomputed_table6_rows(p: BGP212Parameters | None = None) -> tuple[Table6Row, ...]:
    """Recompute Appendix B from Table 3 and Proposition 7.8."""

    p = parameters() if p is None else p
    e, d, a = p.analytic_epsilon, p.delta, p.a1
    cap_steps = tuple(right - left for left, right in zip(p.rough_caps, p.rough_caps[1:]))
    positive_steps = tuple(step for step in cap_steps if step > 0)
    maximum_step = max(cap_steps)
    transition_kappa = Q(11, 80) * e
    problem = packing_problem(p)

    values: dict[str, tuple[Q, Q, str]] = {
        "support.delta": (Q(0), d, "recomputed"),
        "support.B1": (d, p.rough_caps[0], "recomputed"),
        "support.monotone": (Q(0), min(positive_steps), "representative tight positive step"),
        "support.heredity": (maximum_step, d, "recomputed worst step"),
        "packing.first-rung": (p.rough_caps[0], Q(778, 5000), "reported certificate threshold; tree unavailable"),
        "prime.xi2": (p.xi2, Q(2, 5), "recomputed"),
        "prime.xi3": (p.xi2, p.xi3, "recomputed"),
        "prime.harman1": (2 * p.xi1 + 3 * p.xi2, Q(2), "recomputed"),
        "prime.harman2": (p.xi1 + 9 * p.xi2, Q(4), "recomputed"),
        "prime.harman3": (Q(1), 2 * p.xi1 + p.xi2, "recomputed with comparison reversed"),
        "prime.harman4": (17 * p.xi2, Q(7), "recomputed"),
        "global.I.1": (d, p.xi1 - 4 * a + Q(2, 3) - 2 * e, "recomputed"),
        "global.I.2": (d, Q(9, 7) - Q(34, 7) * a - 2 * e, "recomputed"),
        "global.II.range": (Q(0), Q(19, 2) - 36 * a - 13 * d - 9 * e, "recomputed"),
        "global.II.cap1": (d, p.xi2 / 10 - Q(32, 10) * a + Q(8, 10) - 2 * e, "recomputed"),
        "global.II.cap2": (d, p.xi2 / 4 + Q(11, 16) - 3 * a - 2 * e, "recomputed"),
        "global.III": (d, Q(11, 8) - Q(7, 2) * a - Q(9, 8) * p.xi3 - 2 * e, "recomputed"),
        "global.II.dominant": (3 * a + d + 2 * e, Q(63, 80), "recomputed"),
        "packing.roots": (Q(problem.expected_root_count), Q(problem.reported_successful_roots), "count only; certificate trees unavailable"),
        "transition.bin3": (2 * transition_kappa, Q(109, 100) * e, "recomputed"),
        "transition.bin4": (4 * transition_kappa, Q(55, 100) * e, "recomputed"),
    }
    reported = {row.identifier: row for row in reported_table6_rows()}
    return tuple(
        Table6Row(
            row.identifier,
            row.condition,
            values[row.identifier][0],
            values[row.identifier][1],
            row.strict,
            row.use,
            values[row.identifier][2],
        )
        for row in reported.values()
    )


def table6_source_discrepancies(p: BGP212Parameters | None = None) -> tuple[dict[str, str], ...]:
    """Return exact mismatches between printed Appendix B and recomputation."""

    reported = {row.identifier: row for row in reported_table6_rows()}
    recomputed = {row.identifier: row for row in recomputed_table6_rows(p)}
    discrepancies = []
    for identifier, source_row in reported.items():
        replay_row = recomputed[identifier]
        if source_row.left != replay_row.left or source_row.right != replay_row.right:
            discrepancies.append(
                {
                    "identifier": identifier,
                    "reported_left": str(source_row.left),
                    "reported_right": str(source_row.right),
                    "recomputed_left": str(replay_row.left),
                    "recomputed_right": str(replay_row.right),
                    "right_difference": str(replay_row.right - source_row.right),
                }
            )
    return tuple(discrepancies)


def section_9_stale_datum_discrepancy(p: BGP212Parameters | None = None) -> dict[str, str]:
    """Record the stale component values printed in Section 9.2.

    The draft writes ``3*(513/2000) + 179/10000`` where Table 3 has
    ``A1=257/1000`` and ``delta=41/2500``.  Both pairs have the same combined
    value 3937/5000, so the subsequent dominant-wall slack is unaffected.
    """

    p = parameters() if p is None else p
    printed_a1 = Q(513, 2000)
    printed_delta = Q(179, 10000)
    return {
        "table3_A1": str(p.a1),
        "section9_printed_A1": str(printed_a1),
        "table3_delta": str(p.delta),
        "section9_printed_delta": str(printed_delta),
        "table3_combination": str(3 * p.a1 + p.delta),
        "section9_printed_combination": str(3 * printed_a1 + printed_delta),
    }


def section_9_support_prose_discrepancy(
    p: BGP212Parameters | None = None,
) -> dict[str, object]:
    """Record the old cap values retained in the proof of Lemma 9.1.

    The lemma statement uses the symbolic cap function ``C`` and is compatible
    with Table 3.  Its concluding prose instead substitutes the earlier
    ``(3/20, 3/20, 17/100, ...)`` physical cap sequence.  Those displayed
    values are not the Table 3 sequence used by the new support.
    """

    p = parameters() if p is None else p
    actual = tuple(4 * cap for cap in p.rough_caps[:3])
    return {
        "location": "proof of Lemma 9.1",
        "printed_rescaled_caps": ["31/50", "31/50", "17/25"],
        "table3_rescaled_caps_first_three": [str(value) for value in actual],
        "symbolic_lemma_statement_uses_table3_C": True,
    }
