#!/usr/bin/env python3
"""Stage-1 viability check for the PrimeGaps186 depth-3 source cover.

The experiment requested a comparison with a barrier that includes the largest
fragment.  The published outer order-5/2 cover applies only after the largest
witness has been excluded.  This script extracts the exact rational parameters,
implements both predicates literally, and stops before estimation when the
one-fragment boundary check exposes that semantic mismatch.

The renewal discretization is retained here so that the requested reference
calculation has one auditable implementation if the event is corrected later.
It is deliberately not called by ``main`` while the pointwise majorization gate
fails.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Callable, Iterable, Sequence


F = Fraction
UPSTREAM_REVISION = "61340d0b74163003b32756bb16e91d9209a5e330"
DEFAULT_OUTPUT = Path(__file__).with_suffix(".json")


@dataclass(frozen=True)
class SourceParameters:
    n: int
    xi: F
    z: F
    threshold: F
    mu: F
    split: F
    radial_lower: F
    radial_upper: F
    radial_cell_first: int
    radial_cell_last: int
    grid_step: F


@dataclass(frozen=True)
class CoverEvaluation:
    total: float
    components: dict[str, float]


def _ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def _new_ladder() -> list[dict[str, F | int]]:
    """Reproduce the exact new-source ladder used by upstream ``build_inputs``."""

    gap = F(1, 10**7)
    tau = F(1, 10**10)
    rho = F(1, 4) + F(12499, 10**6)
    rho_star = rho - gap
    outer_radius = F(2742997, 10**7) / rho_star
    inner_radius = F(251, 1000) / rho_star
    epsilon = F(1, 10**7)
    limit = F(253, 20000)
    sigma = F(1, 2) - F(40481, 100000) + tau
    affine = (
        ((1 - 5 * sigma) / 15, F(18, 5)),
        ((1 - 4 * sigma) / 16, F(7, 2)),
        ((1 - 2 * sigma) / 20, F(16, 5)),
    )
    error = rho * (outer_radius + inner_radius) - F(1, 2)
    previous = F(0)
    rows: list[dict[str, F | int]] = []
    for index in range(100):
        order = min(index // 12 + 1, 3)
        intercept, slope = affine[order - 1]
        omega = min(limit, (intercept - epsilon - error + 2 * previous - gap) / slope)
        delta = intercept - slope * omega - epsilon
        band = (F(1, 2) + 2 * previous) / rho
        xi = delta / rho
        outer_core = band - inner_radius
        eta = xi if order < 3 else (xi + outer_radius + inner_radius - band) / 2
        rows.append(
            {
                "index": index,
                "order": order,
                "activation": xi,
                "outer_core": outer_core,
                "outer_threshold": outer_core + eta,
            }
        )
        if omega == limit:
            return rows
        previous = omega
    raise ArithmeticError("new source ladder did not reach its terminal value")


def extracted_parameters() -> SourceParameters:
    """Extract group G2 = outer order-5/2 and its sole aligned cap slice."""

    gap = F(1, 10**7)
    rho = F(1, 4) + F(12499, 10**6)
    rho_star = rho - gap
    outer_radius = F(2742997, 10**7) / rho_star
    advance = gap / rho
    step = outer_radius / 98304
    new = _new_ladder()
    answer = SourceParameters(
        n=40,
        xi=F(new[38]["activation"]),
        z=46580 * step,
        threshold=outer_radius + advance / 2,
        mu=F(5, 2),
        split=19660 * step,
        radial_lower=F(new[24]["outer_core"]),
        radial_upper=98303 * step,
        radial_cell_first=95599,
        radial_cell_last=98263,
        grid_step=step,
    )

    # Exact values emitted by the pinned Python source and duplicated in Lean.
    assert answer.xi == F(
        1038826867984921151804142858423732482601802307,
        30904730932085735018956267409864291787032494080000,
    )
    assert answer.z == F(31942200065, 64511729664)
    assert answer.threshold == F(14400682015049, 13781139750220)
    assert answer.split == F(13481830255, 64511729664)
    assert answer.radial_lower == F(
        16430591763736936545249922448197799591,
        16161921199408696007503616565983946000,
    )
    assert answer.radial_upper == F(269644834091, 258046918656)
    return answer


def low_boundaries(parameters: SourceParameters) -> tuple[F, ...]:
    """The 23 endpoints defining the 22 G2 low-fragment components."""

    a = tuple(F(1, 20) * F(6, 5) ** j for j in range(10))
    result = (
        tuple(parameters.xi * 2**j for j in range(9))
        + (F(1, 100), F(3, 200), F(9, 400), F(27, 800))
        + a[:8]
        + ((a[7] + parameters.split) / 2, parameters.split)
    )
    assert len(result) == 23 and all(x < y for x, y in zip(result, result[1:]))
    return result


def rank_boundaries(parameters: SourceParameters) -> tuple[F, ...]:
    """The 13 endpoints defining the 12 G2 rank-two components."""

    fractions = tuple(F(j, 16) for j in range(9)) + (
        F(5, 8),
        F(3, 4),
        F(7, 8),
        F(1),
    )
    q0 = parameters.threshold / (parameters.mu + 1)
    result = tuple(q0 + t * (parameters.z - q0) for t in fractions)
    assert len(result) == 13 and all(x < y for x, y in zip(result, result[1:]))
    return result


def activated_points(points: Iterable[F | float], parameters: SourceParameters) -> list[F]:
    """Apply the source convention ``xi < p <= z`` and sort decreasingly."""

    active = [F(p) for p in points if parameters.xi < F(p) <= parameters.z]
    active.sort(reverse=True)
    return active


def first_violation(
    points: Iterable[F | float], parameters: SourceParameters
) -> dict[str, object] | None:
    """Evaluate the requested inclusive barrier, including its largest point.

    Equal fragment sizes share the full inclusive prefix, as required by the
    numerical note's equation (1.12).  For a simple Poisson realization this is
    the same as scanning decreasing order with an ordinary prefix sum.
    """

    active = activated_points(points, parameters)
    prefix = F(0)
    index = 0
    while index < len(active):
        value = active[index]
        end = index
        while end < len(active) and active[end] == value:
            end += 1
        prefix += (end - index) * value
        barrier = prefix + (parameters.mu - 1) * value
        if barrier > parameters.threshold:
            return {
                "index": index + 1,
                "fragment": value,
                "fragments_at_or_above": end,
                "inclusive_prefix_mass": prefix,
                "barrier_value": barrier,
                "excess": barrier - parameters.threshold,
            }
        index = end
    return None


def bad_event(points: Iterable[F | float], parameters: SourceParameters) -> bool:
    return first_violation(points, parameters) is not None


def current_cover(
    points: Iterable[F | float], parameters: SourceParameters
) -> CoverEvaluation:
    """Evaluate the exact G2 low/rank-two/third-factorial cover.

    The fixed radial mask is assumed to equal one on the selected slice.  Its
    aligned cap is already ``parameters.z``.  This is the point-configuration
    specialization of ``physicalSourceCover`` in the pinned Lean source.
    """

    active = activated_points(points, parameters)
    components: dict[str, float] = {}

    low = low_boundaries(parameters)
    for index, (left, right) in enumerate(zip(low, low[1:]), start=1):
        count = sum(left < p <= right for p in active)
        if not count:
            value = 0.0
        else:
            theta = _ceil_fraction(F(7, 1) / right)
            mass_above_left = sum((p for p in active if p > left), F(0))
            exponent = theta * (
                mass_above_left
                + (parameters.mu - 1) * right
                - parameters.threshold
            )
            try:
                value = count * math.exp(float(exponent))
            except OverflowError:
                value = math.inf
        components[f"low.L{index}"] = value

    rank = rank_boundaries(parameters)
    maximum = active[0] if active else None
    for index, (left, right) in enumerate(zip(rank, rank[1:]), start=1):
        value = 0
        if maximum is not None:
            for q in active:
                if not (left < q <= right) or any(r > q for r in active):
                    continue
                secondary_cutoff = (parameters.threshold - q) / parameters.mu
                if sum(secondary_cutoff < r <= q for r in active) >= 2:
                    value += 1
        components[f"rank_two.P{index}"] = float(value)

    high_count = sum(p > parameters.split for p in active)
    components["third_factorial.H"] = float(math.comb(high_count, 3))
    return CoverEvaluation(total=sum(components.values()), components=components)


def renewal_dp(
    parameters: SourceParameters,
    *,
    u_bins: int,
    s_bins: int,
    terminal_weight: Callable[["object"], "object"] | None = None,
) -> float:
    """Coarse product-grid discretization of the recurrence in the request.

    The u grid and s grid are linear, with trapezoidal quadrature in u and
    linear interpolation in s.  This computes Q_w(0,z), not 1-Q(0,z).
    NumPy is imported lazily because the Stage-1 stop requires no numerical
    dependency.
    """

    if u_bins < 2 or s_bins < 2:
        raise ValueError("both grids need at least two bins")
    import numpy as np

    weight = terminal_weight or (lambda value: np.ones_like(value))
    xi, z, threshold, mu = map(
        float, (parameters.xi, parameters.z, parameters.threshold, parameters.mu)
    )
    u = np.linspace(xi, z, u_bins + 1, dtype=np.float64)
    s = np.linspace(0.0, threshold, s_bins + 1, dtype=np.float64)
    ds = s[1] - s[0]
    table = np.empty((s_bins + 1, u_bins + 1), dtype=np.float64)
    u_power = u**parameters.n
    xi_power = xi**parameters.n

    def cumulative_at(values: "object", cumulative: "object", endpoint: float) -> float:
        if endpoint <= xi:
            return 0.0
        if endpoint >= z:
            return float(cumulative[-1])
        right = int(np.searchsorted(u, endpoint, side="right"))
        left = right - 1
        width = endpoint - u[left]
        fraction = width / (u[right] - u[left])
        endpoint_value = values[left] + fraction * (values[right] - values[left])
        return float(cumulative[left] + width * (values[left] + endpoint_value) / 2)

    for si in range(s_bins, -1, -1):
        state = s[si]
        if si == s_bins:
            table[si, :] = xi_power * float(np.asarray(weight(np.asarray(state)))) / u_power
            continue
        targets = state + u
        lower = np.floor(targets / ds).astype(np.int64)
        lower = np.clip(lower, si + 1, s_bins)
        upper = np.minimum(lower + 1, s_bins)
        fraction = np.clip((targets - s[lower]) / ds, 0.0, 1.0)
        columns = np.arange(u_bins + 1)
        continuation = (
            (1.0 - fraction) * table[lower, columns]
            + fraction * table[upper, columns]
        )
        integrand = parameters.n * u ** (parameters.n - 1) * continuation
        cumulative = np.zeros_like(u)
        cumulative[1:] = np.cumsum(
            (u[1:] - u[:-1]) * (integrand[1:] + integrand[:-1]) / 2
        )
        allowed = (threshold - state) / mu
        if allowed <= xi:
            integral = np.zeros_like(u)
        else:
            capped = cumulative_at(integrand, cumulative, min(allowed, z))
            integral = np.where(u <= allowed, cumulative, capped)
        base = xi_power * float(np.asarray(weight(np.asarray(state))))
        table[si, :] = np.clip((base + integral) / u_power, 0.0, np.inf)
    return float(table[0, -1])


def _fraction_json(value: F) -> dict[str, object]:
    return {"fraction": str(value), "decimal": float(value)}


def _jsonify(value: object) -> object:
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonify(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonify(item) for item in value]
    return value


def boundary_checks(parameters: SourceParameters) -> list[dict[str, object]]:
    """Exercise the source endpoints and the decisive one-largest state."""

    tiny = min(parameters.xi, parameters.z - parameters.split) / 10**6
    barely_safe_largest = parameters.threshold / parameters.mu - 2 * tiny
    barely_violating_second = (
        (parameters.threshold - barely_safe_largest) / parameters.mu + tiny
    )
    cases: list[tuple[str, Sequence[F]]] = [
        ("empty", ()),
        ("at_activation_excluded", (parameters.xi,)),
        ("just_above_activation", (parameters.xi + tiny,)),
        ("at_split_low_inclusive_high_exclusive", (parameters.split,)),
        ("just_above_split", (parameters.split + tiny,)),
        ("one_large_mid", ((parameters.threshold / parameters.mu + parameters.z) / 2,)),
        ("one_large_at_cap", (parameters.z,)),
        (
            "two_large_barely_violating_at_second",
            (barely_safe_largest, barely_violating_second),
        ),
        (
            "three_above_split",
            (parameters.z, parameters.split + 2 * tiny, parameters.split + tiny),
        ),
    ]
    answer = []
    for name, points in cases:
        violation = first_violation(points, parameters)
        cover = current_cover(points, parameters)
        nonzero = {key: value for key, value in cover.components.items() if value != 0.0}
        answer.append(
            {
                "name": name,
                "points": [str(point) for point in points],
                "bad": violation is not None,
                "first_violation": _jsonify(violation),
                "cover": cover.total,
                "nonzero_cover_components": nonzero,
                "majorizes": cover.total >= int(violation is not None),
            }
        )
    return answer


def build_result() -> dict[str, object]:
    parameters = extracted_parameters()
    checks = boundary_checks(parameters)
    failures = [row for row in checks if not row["majorizes"]]
    assert {row["name"] for row in failures} == {"one_large_mid", "one_large_at_cap"}

    one_cap = next(row for row in checks if row["name"] == "one_large_at_cap")
    assert one_cap["bad"] and one_cap["cover"] == 0.0
    assert parameters.mu * parameters.z > parameters.threshold

    return {
        "schema": "primegaps.direct-source-barrier-viability.v1",
        "status": "INCONCLUSIVE",
        "stage_reached": "stage_1_semantic_validation",
        "upstream": {
            "repository": "https://github.com/openai/PrimeGaps186",
            "revision": UPSTREAM_REVISION,
            "group": "G2 / outer_h25 / outer order-5/2",
            "radial_cap_slice": "sole group slice, radial cell sums 95599..98263, aligned cap 46580h",
            "source_locations": {
                "python_group_derivation": "prime_gap_186_certificate.py:200-261,342-425",
                "python_component_schedule": "prime_gap_186_certificate.py:264-339",
                "lean_group_and_endpoints": "PrimeGaps186.lean:736-893",
                "lean_radial_mask_and_cover": "PrimeGaps186.lean:927-997",
                "numerical_note": "short_gaps_numerics.pdf, equations (1.12), (1.23)-(1.26), (1.34)-(1.47); Lemmas 1.4-1.7",
            },
        },
        "parameters": {
            key: value if isinstance(value, int) else _fraction_json(value)
            for key, value in asdict(parameters).items()
        },
        "cover_inventory": {
            "low_components": 22,
            "rank_two_components": 12,
            "third_factorial_components": 1,
            "total_components": 35,
            "low_endpoint_convention": "(left,right]",
            "rank_endpoint_convention": "(left,right]",
            "activation_convention": "xi < p",
            "cap_convention": "p <= z",
            "high_convention": "p > split",
        },
        "stage_1": {
            "passed": False,
            "boundary_checks": checks,
            "counterexample": {
                "points": [str(parameters.z)],
                "bad_event_reason": "mu*z > T",
                "mu_z_minus_T": _fraction_json(parameters.mu * parameters.z - parameters.threshold),
                "cover_value": 0,
                "reason_cover_is_zero": (
                    "the point is above every low bin; rank-two requires a second point; "
                    "choose(N_(split,infinity),3)=choose(1,3)=0"
                ),
            },
            "source_semantics": (
                "The upstream cover majorizes nonlargest order-5/2 failures only, after separate "
                "largest-fragment cap arguments. The requested B includes the largest fragment."
            ),
        },
        "stages_not_run": {
            "monte_carlo": "not run because the required pointwise majorization failed",
            "renewal_dp_256": "not run because the required pointwise majorization failed",
            "renewal_dp_512": "not run because the required pointwise majorization failed",
            "stress_test": "not run because C_eta is not a majorant of the stated B_eta",
            "weighted_root_side": "not run",
        },
        "conclusion": {
            "decision": "INCONCLUSIVE",
            "statement": (
                "No E[C]/P(B) viability conclusion is valid for the stated event. "
                "Replace B by the source-faithful nonlargest event (or add the omitted "
                "largest-witness cover) before estimation."
            ),
            "best_next_step": (
                "Specify whether the intended event is max over nonlargest fragments j>=2 and "
                "whether the PPP law is conditioned on the selected radial slice; then derive "
                "and validate the modified renewal initial condition before MC."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build_result()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
