#!/usr/bin/env python3
"""Exact symmetry-compressed verifier for Stadlmann's one-band support.

This is an independent bridge between the low-dimensional exact slice kernel
and the published even-signature symmetric basis.  It evaluates a *fixed*
rational polynomial; it does not search or optimize support parameters.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import cache, lru_cache
from itertools import product
import json
from math import comb, factorial, isfinite
from pathlib import Path
from typing import Iterable

import primegaps.integrals as kernel


Signature = tuple[int, ...]
Poly2 = dict[tuple[int, int], Fraction]
U = Fraction(521, 2000)       # A_1 + epsilon = .2605
R = Fraction(491, 2000)       # A_1 - epsilon = .2455
DELTA = Fraction(7, 250)      # .028
B = (Fraction(3, 20), Fraction(3, 20)) + (Fraction(17, 100),) * 4


def canonical(values: Iterable[int]) -> Signature:
    return tuple(sorted((int(v) for v in values if v), reverse=True))


def orbit_count(dimension: int, signature: Signature) -> int:
    if len(signature) > dimension:
        return 0
    result = factorial(dimension) // factorial(dimension - len(signature))
    for multiplicity in Counter(signature).values():
        result //= factorial(multiplicity)
    return result


@cache
def product_signatures(
    dimension: int, left: Signature, right: Signature
) -> tuple[tuple[Signature, int], ...]:
    """Structure constants for m_left*m_right in the monomial basis."""
    left = canonical(left)
    right = canonical(right)
    remaining = Counter(right)
    representative_counts: Counter[Signature] = Counter()

    def visit(slot: int, accumulated: list[int]) -> None:
        if slot == len(left):
            tail: list[int] = []
            divisor = 1
            tail_count = 0
            for value, multiplicity in remaining.items():
                tail.extend([value] * multiplicity)
                divisor *= factorial(multiplicity)
                tail_count += multiplicity
            if tail_count > dimension - len(left):
                return
            placements = factorial(dimension - len(left)) // factorial(
                dimension - len(left) - tail_count
            )
            placements //= divisor
            representative_counts[canonical((*accumulated, *tail))] += placements
            return
        visit(slot + 1, accumulated + [left[slot]])
        for value in tuple(remaining):
            if not remaining[value]:
                continue
            remaining[value] -= 1
            visit(slot + 1, accumulated + [left[slot] + value])
            remaining[value] += 1

    visit(0, [])
    left_terms = orbit_count(dimension, left)
    answer = []
    for signature, representative_count in representative_counts.items():
        numerator = left_terms * representative_count
        denominator = orbit_count(dimension, signature)
        coefficient, remainder = divmod(numerator, denominator)
        if remainder:
            raise ArithmeticError("non-integral monomial-symmetric structure constant")
        answer.append((signature, coefficient))
    return tuple(sorted(answer))


def _allocations(counts: tuple[int, ...], limit: int):
    for values in product(*(range(count + 1) for count in counts)):
        if sum(values) <= limit:
            yield values


@cache
def _shifted_factorial_polynomial(exponent: int, count: int) -> dict[int, Fraction]:
    """Product of count copies of sum_p C(e,p) delta^(e-p) p! X^p."""
    result = {0: Fraction(1)}
    one = {
        power: Fraction(comb(exponent, power))
        * DELTA ** (exponent - power)
        * factorial(power)
        for power in range(exponent + 1)
    }
    for _ in range(count):
        next_result: dict[int, Fraction] = defaultdict(Fraction)
        for a, ca in result.items():
            for b, cb in one.items():
                next_result[a + b] += ca * cb
        result = dict(next_result)
    return result


def _multiply_univariate(
    left: dict[int, Fraction], right: dict[int, Fraction]
) -> dict[int, Fraction]:
    result: dict[int, Fraction] = defaultdict(Fraction)
    for a, ca in left.items():
        for b, cb in right.items():
            result[a + b] += ca * cb
    return dict(result)


@cache
def orbit_status_densities(
    dimension: int, signature: Signature, large_count: int
) -> tuple[tuple[tuple[int, int | None, int | None], Fraction], ...]:
    """Aggregate orbit, status, IE, shift, and group-density coefficients.

    Keys are (number of inclusion-exclusion shifts among small coordinates,
    large-sum power, small-sum power).  A missing group has power None.
    """
    signature = canonical(signature)
    if len(signature) > dimension or not 0 <= large_count <= dimension:
        return ()
    multiplicities = Counter(signature)
    exponents = tuple(sorted(multiplicities))
    counts = tuple(multiplicities[e] for e in exponents)
    nonzero_count = len(signature)
    answer: dict[tuple[int, int | None, int | None], Fraction] = defaultdict(Fraction)

    for large_alloc in _allocations(counts, large_count):
        large_nonzero = sum(large_alloc)
        small_alloc = tuple(c - a for c, a in zip(counts, large_alloc))
        small_nonzero = sum(small_alloc)
        if small_nonzero > dimension - large_count:
            continue
        large_placements = factorial(large_count) // factorial(
            large_count - large_nonzero
        )
        small_placements = factorial(dimension - large_count) // factorial(
            dimension - large_count - small_nonzero
        )
        for value in large_alloc:
            large_placements //= factorial(value)
        for value in small_alloc:
            small_placements //= factorial(value)
        placements = comb(dimension, large_count) * large_placements * small_placements

        large_poly = {0: Fraction(1)}
        for exponent, count in zip(exponents, large_alloc):
            large_poly = _multiply_univariate(
                large_poly, _shifted_factorial_polynomial(exponent, count)
            )

        zero_small = dimension - large_count - small_nonzero
        for shifted_alloc in _allocations(small_alloc, dimension - large_count):
            shifted_nonzero = sum(shifted_alloc)
            unshifted_alloc = tuple(
                count - shifted
                for count, shifted in zip(small_alloc, shifted_alloc)
            )
            shifted_choice = 1
            small_poly = {0: Fraction(1)}
            fixed_degree = 0
            fixed_factorials = 1
            for exponent, available, shifted, unshifted in zip(
                exponents, small_alloc, shifted_alloc, unshifted_alloc
            ):
                shifted_choice *= comb(available, shifted)
                small_poly = _multiply_univariate(
                    small_poly, _shifted_factorial_polynomial(exponent, shifted)
                )
                fixed_degree += exponent * unshifted
                fixed_factorials *= factorial(exponent) ** unshifted
            small_poly = {
                degree + fixed_degree: coefficient * fixed_factorials
                for degree, coefficient in small_poly.items()
            }
            for shifted_zeros in range(zero_small + 1):
                shifted_count = shifted_nonzero + shifted_zeros
                coefficient = Fraction(
                    placements
                    * shifted_choice
                    * comb(zero_small, shifted_zeros)
                    * (-1) ** shifted_count
                )
                for large_degree, large_coefficient in large_poly.items():
                    large_power = (
                        None
                        if large_count == 0
                        else large_degree + large_count - 1
                    )
                    large_density = (
                        Fraction(1)
                        if large_power is None
                        else large_coefficient / factorial(large_power)
                    )
                    for small_degree, small_coefficient in small_poly.items():
                        small_power = (
                            None
                            if dimension == large_count
                            else small_degree + dimension - large_count - 1
                        )
                        small_density = (
                            Fraction(1)
                            if small_power is None
                            else small_coefficient / factorial(small_power)
                        )
                        answer[(shifted_count, large_power, small_power)] += (
                            coefficient * large_density * small_density
                        )
    return tuple((key, value) for key, value in sorted(answer.items()) if value)


@cache
def _beta_segment(power: int, other: int, height: Fraction, maximum: Fraction) -> Fraction:
    if maximum <= 0 or height <= 0:
        return Fraction(0)
    maximum = min(maximum, height)
    if other == 0:
        return maximum ** (power + 1) / (power + 1)
    # Integral recurrence M[p,q] = H M[p,q-1] - M[p+1,q-1].
    # Caching fills one small rectangular moment table for each status cap,
    # avoiding a fresh alternating binomial sum for every orbit-density key.
    return height * _beta_segment(
        power, other - 1, height, maximum
    ) - _beta_segment(power + 1, other - 1, height, maximum)


def _radial_group_integral(
    large_power: int | None,
    small_power: int | None,
    slack: int,
    height: Fraction,
    large_cap: Fraction | None,
) -> Fraction:
    if height <= 0:
        return Fraction(0)
    if large_power is None and small_power is None:
        return height**slack
    if large_power is None:
        assert small_power is not None
        return (
            Fraction(factorial(small_power) * factorial(slack),
                     factorial(small_power + slack + 1))
            * height ** (small_power + slack + 1)
        )
    maximum = height if large_cap is None else min(height, large_cap)
    if small_power is None:
        return _beta_segment(large_power, slack, height, maximum)
    inner_factor = Fraction(
        factorial(small_power) * factorial(slack),
        factorial(small_power + slack + 1),
    )
    return inner_factor * _beta_segment(
        large_power, small_power + slack + 1, height, maximum
    )


@cache
def i_orbit(dimension: int, signature: Signature, slack: int) -> Fraction:
    answer = Fraction(0)
    for large_count in range(min(dimension, len(B)) + 1):
        large_cap = None if large_count == 0 else B[large_count - 1] - large_count * DELTA
        if large_cap is not None and large_cap <= 0:
            continue
        for (shifted_count, large_power, small_power), coefficient in orbit_status_densities(
            dimension, signature, large_count
        ):
            height = U - (large_count + shifted_count) * DELTA
            answer += coefficient * _radial_group_integral(
                large_power, small_power, slack, height, large_cap
            )
    return answer


@dataclass(frozen=True)
class RadialSlice:
    power: int
    slack: int
    large: bool
    lower: Fraction = Fraction(0)
    upper: Fraction = U
    support_limit: Fraction | None = None


@cache
def _affine_power(form, exponent: int) -> Poly2:
    return kernel._affine_power(form, exponent)


@cache
def _slice_polynomial(
    spec: RadialSlice,
    total_offset: Fraction,
    large_offset: Fraction,
    sample: tuple[Fraction, Fraction],
) -> Poly2:
    common = (Fraction(1), Fraction(1), total_offset)
    large_sum = (Fraction(1), Fraction(0), large_offset)
    lower_band = kernel._affine_sub((Fraction(0), Fraction(0), spec.lower), common)
    upper_band = kernel._affine_sub((Fraction(0), Fraction(0), spec.upper), common)
    if spec.large:
        lowers = ((Fraction(0), Fraction(0), DELTA), lower_band)
        uppers = [(Fraction(0), Fraction(0), Fraction(1)), upper_band]
        if spec.support_limit is not None:
            uppers.append(
                kernel._affine_sub(
                    (Fraction(0), Fraction(0), spec.support_limit), large_sum
                )
            )
    else:
        if spec.support_limit is not None and kernel._affine_value(large_sum, sample) > spec.support_limit:
            return {}
        lowers = ((Fraction(0), Fraction(0), Fraction(0)), lower_band)
        uppers = [(Fraction(0), Fraction(0), DELTA), upper_band]
    lower = max(lowers, key=lambda form: kernel._affine_value(form, sample))
    upper = min(uppers, key=lambda form: kernel._affine_value(form, sample))
    if kernel._affine_value(upper, sample) <= kernel._affine_value(lower, sample):
        return {}
    room = (-common[0], -common[1], U - common[2])
    answer: Poly2 = {}
    for j in range(spec.slack + 1):
        exponent = spec.power + j + 1
        endpoint = kernel._poly_add(
            _affine_power(upper, exponent),
            kernel._poly_scale(_affine_power(lower, exponent), -1),
        )
        term = kernel._poly_mul(
            _affine_power(room, spec.slack - j), endpoint
        )
        term = kernel._poly_scale(
            term, Fraction((-1) ** j * comb(spec.slack, j), exponent)
        )
        answer = kernel._poly_add(answer, term)
    return answer


@cache
def _slice_geometry(
    large_present: bool,
    small_present: bool,
    total_offset: Fraction,
    large_offset: Fraction,
    specs: tuple[RadialSlice, ...],
):
    """Pre-split the rational line arrangement shared by all orbit powers."""
    residual_upper = R - total_offset
    if residual_upper <= 0:
        return ("empty", ())
    lines = kernel._unique_lines(
        line
        for spec in specs
        for line in kernel._slice_lines(spec, total_offset, large_offset, DELTA)
    )
    if not large_present and not small_present:
        return ("point", ((Fraction(0), Fraction(0)),))
    if not large_present or not small_present:
        variable_is_x = large_present
        cuts = {Fraction(0), residual_upper}
        for a, b, c in lines:
            slope = a if variable_is_x else b
            if slope:
                root = -c / slope
                if 0 < root < residual_upper:
                    cuts.add(root)
        ordered = sorted(cuts)
        intervals = tuple(
            (start, end, ((start + end) / 2, Fraction(0)) if variable_is_x else
             (Fraction(0), (start + end) / 2))
            for start, end in zip(ordered, ordered[1:])
        )
        return ("xintervals" if variable_is_x else "zintervals", intervals)
    polygons = [[
        (Fraction(0), Fraction(0)),
        (residual_upper, Fraction(0)),
        (Fraction(0), residual_upper),
    ]]
    for line in lines:
        next_polygons = []
        for polygon in polygons:
            positive, negative = kernel._split_polygon(polygon, line)
            if positive:
                next_polygons.append(positive)
            if negative:
                next_polygons.append(negative)
        polygons = next_polygons
    return (
        "polygons",
        tuple((tuple(polygon), kernel._polygon_centroid(polygon)) for polygon in polygons),
    )


@cache
def _edge_moment(
    start: tuple[Fraction, Fraction],
    end: tuple[Fraction, Fraction],
    x_power: int,
    z_power: int,
) -> Fraction:
    """Integral_0^1 x(t)^p z(t)^q dt on one affine edge."""
    x0, z0 = start
    x1, z1 = end
    dx, dz = x1 - x0, z1 - z0
    if dx:
        if z_power == 0:
            return (x1 ** (x_power + 1) - x0 ** (x_power + 1)) / (
                (x_power + 1) * dx
            )
        return (
            (x1 ** (x_power + 1) * z1**z_power -
             x0 ** (x_power + 1) * z0**z_power)
            - z_power * dz * _edge_moment(start, end, x_power + 1, z_power - 1)
        ) / ((x_power + 1) * dx)
    if dz:
        if x_power == 0:
            return (z1 ** (z_power + 1) - z0 ** (z_power + 1)) / (
                (z_power + 1) * dz
            )
        return (
            (x1**x_power * z1 ** (z_power + 1) -
             x0**x_power * z0 ** (z_power + 1))
            - x_power * dx * _edge_moment(start, end, x_power - 1, z_power + 1)
        ) / ((z_power + 1) * dz)
    return x0**x_power * z0**z_power


@cache
def _polygon_monomial_moment(
    x_power: int, z_power: int, polygon: tuple[tuple[Fraction, Fraction], ...]
) -> Fraction:
    # Cone the oriented polygon edges to the origin.  Homogeneity supplies the
    # radial factor 1/(p+q+2); `_edge_moment` supplies the boundary moment.
    answer = Fraction(0)
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        determinant = start[0] * end[1] - end[0] * start[1]
        answer += determinant * _edge_moment(start, end, x_power, z_power)
    return answer / (x_power + z_power + 2)


def _integrate_polynomial_polygon(
    polynomial: Poly2, polygon: tuple[tuple[Fraction, Fraction], ...]
) -> Fraction:
    return sum(
        coefficient * _polygon_monomial_moment(x_power, z_power, polygon)
        for (x_power, z_power), coefficient in polynomial.items()
    )


def _integrate_intervals(
    variable_is_x: bool,
    power: int,
    lower: Fraction,
    upper: Fraction,
    lines,
    total_offset: Fraction,
    large_offset: Fraction,
    specs: tuple[RadialSlice, ...],
) -> Fraction:
    cuts = {lower, upper}
    for a, b, c in lines:
        slope = a if variable_is_x else b
        if slope:
            root = -c / slope
            if lower < root < upper:
                cuts.add(root)
    answer = Fraction(0)
    ordered = sorted(cuts)
    for start, end in zip(ordered, ordered[1:]):
        middle = (start + end) / 2
        sample = (middle, Fraction(0)) if variable_is_x else (Fraction(0), middle)
        polynomial: Poly2 = {
            (power, 0) if variable_is_x else (0, power): Fraction(1)
        }
        for spec in specs:
            polynomial = kernel._poly_mul(
                polynomial, _slice_polynomial(spec, total_offset, large_offset, sample)
            )
        for (x_degree, z_degree), coefficient in polynomial.items():
            if variable_is_x and z_degree == 0:
                degree = x_degree
            elif not variable_is_x and x_degree == 0:
                degree = z_degree
            else:
                continue
            answer += coefficient * (end ** (degree + 1) - start ** (degree + 1)) / (
                degree + 1
            )
    return answer


def integrate_slice_cells(
    large_power: int | None,
    small_power: int | None,
    total_offset: Fraction,
    large_offset: Fraction,
    specs: tuple[RadialSlice, ...],
) -> Fraction:
    kind, cells = _slice_geometry(
        large_power is not None,
        small_power is not None,
        total_offset,
        large_offset,
        specs,
    )
    if kind == "empty":
        return Fraction(0)
    if kind == "point":
        sample = cells[0]
        value = Fraction(1)
        for spec in specs:
            value *= _slice_polynomial(spec, total_offset, large_offset, sample).get(
                (0, 0), Fraction(0)
            )
        return value
    if kind in ("xintervals", "zintervals"):
        variable_is_x = kind == "xintervals"
        power = large_power if variable_is_x else small_power
        assert power is not None
        answer = Fraction(0)
        for start, end, sample in cells:
            polynomial: Poly2 = {
                (power, 0) if variable_is_x else (0, power): Fraction(1)
            }
            for spec in specs:
                polynomial = kernel._poly_mul(
                    polynomial,
                    _slice_polynomial(spec, total_offset, large_offset, sample),
                )
            for (x_degree, z_degree), coefficient in polynomial.items():
                if variable_is_x and z_degree == 0:
                    degree = x_degree
                elif not variable_is_x and x_degree == 0:
                    degree = z_degree
                else:
                    continue
                answer += coefficient * (
                    end ** (degree + 1) - start ** (degree + 1)
                ) / (degree + 1)
        return answer
    answer = Fraction(0)
    assert kind == "polygons"
    for polygon, sample in cells:
        polynomial: Poly2 = {(large_power, small_power): Fraction(1)}
        for spec in specs:
            polynomial = kernel._poly_mul(
                polynomial, _slice_polynomial(spec, total_offset, large_offset, sample)
            )
            if not polynomial:
                break
        answer += _integrate_polynomial_polygon(polynomial, polygon)
    return answer


def _support_limit(common_large: int, special_large: bool) -> Fraction | None:
    count = common_large + int(special_large)
    if count == 0:
        return None
    if count > len(B):
        raise ValueError("empty large-coordinate status")
    return B[count - 1]


@cache
def j_orbit(
    dimension: int,
    signature: Signature,
    left_power: int,
    left_slack: int,
    right_power: int,
    right_slack: int,
) -> Fraction:
    answer = Fraction(0)
    for large_count in range(min(dimension, len(B)) + 1):
        densities = orbit_status_densities(dimension, signature, large_count)
        for left_large in (False, True):
            if large_count + int(left_large) > len(B):
                continue
            for right_large in (False, True):
                if large_count + int(right_large) > len(B):
                    continue
                specs = (
                    RadialSlice(
                        left_power, left_slack, left_large,
                        support_limit=_support_limit(large_count, left_large)
                    ),
                    RadialSlice(
                        right_power, right_slack, right_large,
                        support_limit=_support_limit(large_count, right_large)
                    ),
                )
                for (shifted_count, large_power, small_power), coefficient in densities:
                    answer += coefficient * integrate_slice_cells(
                        large_power,
                        small_power,
                        (large_count + shifted_count) * DELTA,
                        large_count * DELTA,
                        specs,
                    )
    return answer


@dataclass(frozen=True)
class Term:
    signature: Signature
    slack: int
    coefficient: Fraction


def rational_terms_from_candidate(path: str | Path, dimension: int) -> tuple[Term, ...]:
    """Load and validate a recorded exact-decimal rational candidate."""
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != "primegaps-stadlmann-rational-candidate-v1":
        raise ValueError("unsupported rational candidate schema")
    if payload.get("k") != dimension or payload.get("degree") != 21:
        raise ValueError("candidate k/degree does not match this evaluation")
    if payload.get("basis_dimension") != 846:
        raise ValueError("candidate basis dimension is not the D=21 value 846")
    if payload.get("radial_coordinate") != "(U-sum(t)) with U=521/2000":
        raise ValueError("candidate radial-coordinate convention does not match")
    combined: dict[tuple[Signature, int], Fraction] = defaultdict(Fraction)
    for item in payload["terms"]:
        raw_signature = item["signature"]
        if not isinstance(raw_signature, list) or any(
            isinstance(x, bool) or not isinstance(x, int) or x <= 0 for x in raw_signature
        ):
            raise ValueError("candidate signatures must contain positive JSON integers")
        signature = canonical(raw_signature)
        if tuple(raw_signature) != signature:
            raise ValueError("candidate signature is not canonical")
        raw_slack = item["slack_power"]
        if isinstance(raw_slack, bool) or not isinstance(raw_slack, int) or raw_slack < 0:
            raise ValueError("candidate slack power must be a non-negative JSON integer")
        slack = raw_slack
        if sum(signature) + slack > 21 or any(x % 2 for x in signature):
            raise ValueError("candidate term violates the even-signature D=21 basis")
        if (signature, slack) in combined:
            raise ValueError("candidate contains a duplicate term")
        if (
            isinstance(item["numerator"], bool)
            or not isinstance(item["numerator"], int)
            or isinstance(item["denominator"], bool)
            or not isinstance(item["denominator"], int)
        ):
            raise ValueError("candidate numerator/denominator must be JSON integers")
        if item["denominator"] <= 0:
            raise ValueError("candidate denominator must be positive")
        combined[(signature, slack)] += Fraction(
            int(item["numerator"]), int(item["denominator"])
        )
    terms = tuple(
        Term(signature, slack, coefficient)
        for (signature, slack), coefficient in sorted(combined.items())
        if coefficient
    )
    if len(terms) != 846:
        raise ValueError(f"candidate has {len(terms)} nonzero terms, expected 846")
    return terms


def jacobi_q_coefficients(degree: int, beta: int) -> tuple[int, ...]:
    """Coefficients of P_degree^(0,beta)(1-2q), in ascending q powers."""
    answer = [0] * (degree + 1)
    for m in range(degree + 1):
        base = comb(degree, m) * comb(degree + beta, degree - m) * (-1) ** (
            degree - m
        )
        for h in range(m + 1):
            answer[degree - m + h] += base * comb(m, h) * (-1) ** h
    return tuple(answer)


def rational_terms_from_orthogonal_vector(
    vector: Iterable[float],
    basis: Iterable[tuple[tuple[int, ...], int]],
    dimension: int,
    significant_digits: int = 12,
) -> tuple[Term, ...]:
    """Convert the conditioned Jacobi basis to rational (U-sum)^b terms.

    Each resulting floating coefficient is rounded independently to the stated
    significant digits and then interpreted as an exact decimal rational.
    """
    vector = tuple(vector)
    basis = tuple(basis)
    if len(vector) != len(basis):
        raise ValueError("orthogonal vector and basis have different lengths")
    if any(not isfinite(float(value)) for value in vector):
        raise ValueError("orthogonal vector contains a non-finite coefficient")
    combined: dict[tuple[Signature, int], float] = defaultdict(float)
    for raw_coefficient, (partition, radial_degree) in zip(vector, basis):
        signature = tuple(2 * value for value in partition)
        beta = dimension - 1 + 2 * sum(signature)
        for slack, jacobi_coefficient in enumerate(
            jacobi_q_coefficients(radial_degree, beta)
        ):
            combined[(signature, slack)] += (
                float(raw_coefficient) * jacobi_coefficient / float(U) ** slack
            )
    terms = []
    for (signature, slack), coefficient in sorted(combined.items()):
        if coefficient:
            rounded = Fraction(format(coefficient, f".{significant_digits}e"))
            if rounded:
                terms.append(Term(signature, slack, rounded))
    return tuple(terms)


def aggregate_i_atoms(
    terms: tuple[Term, ...], dimension: int
) -> dict[tuple[Signature, int], Fraction]:
    atoms: dict[tuple[Signature, int], Fraction] = defaultdict(Fraction)
    for i, left in enumerate(terms):
        for j in range(i, len(terms)):
            right = terms[j]
            multiplier = 1 if i == j else 2
            coefficient = multiplier * left.coefficient * right.coefficient
            for signature, structure in product_signatures(
                dimension, left.signature, right.signature
            ):
                atoms[(signature, left.slack + right.slack)] += coefficient * structure
    return {key: value for key, value in atoms.items() if value}


def exact_i_grouped(terms: tuple[Term, ...], dimension: int) -> Fraction:
    return sum(
        coefficient * i_orbit(dimension, signature, slack)
        for (signature, slack), coefficient in aggregate_i_atoms(terms, dimension).items()
    )


def exact_i(terms: tuple[Term, ...], dimension: int) -> Fraction:
    answer = Fraction(0)
    for i, left in enumerate(terms):
        for j in range(i, len(terms)):
            right = terms[j]
            multiplier = 1 if i == j else 2
            for signature, structure in product_signatures(
                dimension, left.signature, right.signature
            ):
                answer += (
                    multiplier * left.coefficient * right.coefficient * structure
                    * i_orbit(dimension, signature, left.slack + right.slack)
                )
    return answer


def marginal_features(terms: tuple[Term, ...]) -> tuple[tuple[Signature, int, int, Fraction], ...]:
    combined: dict[tuple[Signature, int, int], Fraction] = defaultdict(Fraction)
    for term in terms:
        combined[(term.signature, 0, term.slack)] += term.coefficient
        for exponent in set(term.signature):
            erased = list(term.signature)
            erased.remove(exponent)
            combined[(canonical(erased), exponent, term.slack)] += term.coefficient
    return tuple((*key, value) for key, value in sorted(combined.items()) if value)


def exact_j(terms: tuple[Term, ...], dimension: int) -> Fraction:
    features = marginal_features(terms)
    answer = Fraction(0)
    for i, (left_signature, left_power, left_slack, left_coefficient) in enumerate(features):
        for j in range(i, len(features)):
            right_signature, right_power, right_slack, right_coefficient = features[j]
            multiplier = 1 if i == j else 2
            for signature, structure in product_signatures(
                dimension - 1, left_signature, right_signature
            ):
                answer += (
                    multiplier * left_coefficient * right_coefficient * structure
                    * j_orbit(
                        dimension - 1, signature,
                        left_power, left_slack, right_power, right_slack
                    )
                )
    return answer


def aggregate_j_atoms(
    terms: tuple[Term, ...], dimension: int
) -> dict[Signature, dict[tuple[int, int, int, int], Fraction]]:
    """Contract a fixed candidate to exact J atom weights, grouped by orbit."""
    features = marginal_features(terms)
    groups: dict[Signature, dict[tuple[int, int, int, int], Fraction]] = {}
    for i, (left_signature, left_power, left_slack, left_coefficient) in enumerate(features):
        for j in range(i, len(features)):
            right_signature, right_power, right_slack, right_coefficient = features[j]
            multiplier = 1 if i == j else 2
            coefficient = multiplier * left_coefficient * right_coefficient
            left_key = (left_power, left_slack)
            right_key = (right_power, right_slack)
            if right_key < left_key:
                left_key, right_key = right_key, left_key
            atom = (*left_key, *right_key)
            for signature, structure in product_signatures(
                dimension - 1, left_signature, right_signature
            ):
                bucket = groups.setdefault(signature, defaultdict(Fraction))
                bucket[atom] += coefficient * structure
    return {
        signature: {key: value for key, value in atoms.items() if value}
        for signature, atoms in groups.items()
    }


def grouped_marginal_coefficients(
    terms: tuple[Term, ...]
) -> dict[Signature, dict[tuple[int, int], Fraction]]:
    groups: dict[Signature, dict[tuple[int, int], Fraction]] = {}
    for signature, power, slack, coefficient in marginal_features(terms):
        groups.setdefault(signature, {})[(power, slack)] = coefficient
    return groups


def grouped_signature_pairs(
    signatures: Iterable[Signature], dimension: int
) -> dict[Signature, tuple[tuple[Signature, Signature, int], ...]]:
    ordered = tuple(sorted(signatures))
    answer: dict[Signature, list[tuple[Signature, Signature, int]]] = defaultdict(list)
    for i, left in enumerate(ordered):
        for j in range(i, len(ordered)):
            right = ordered[j]
            symmetry = 1 if i == j else 2
            for target, structure in product_signatures(dimension - 1, left, right):
                answer[target].append((left, right, symmetry * structure))
    return {target: tuple(pairs) for target, pairs in answer.items()}


def _linear_slice_polynomial(
    coefficients: dict[tuple[int, int], Fraction],
    large: bool,
    support_limit: Fraction | None,
    total_offset: Fraction,
    large_offset: Fraction,
    sample: tuple[Fraction, Fraction],
) -> Poly2:
    return _linear_slice_polynomial_cached(
        tuple(sorted(coefficients.items())),
        large,
        support_limit,
        total_offset,
        large_offset,
        sample,
    )


@lru_cache(maxsize=120_000)
def _linear_slice_polynomial_cached(
    coefficient_items: tuple[tuple[tuple[int, int], Fraction], ...],
    large: bool,
    support_limit: Fraction | None,
    total_offset: Fraction,
    large_offset: Fraction,
    sample: tuple[Fraction, Fraction],
) -> Poly2:
    coefficients = dict(coefficient_items)
    common = (Fraction(1), Fraction(1), total_offset)
    large_sum = (Fraction(1), Fraction(0), large_offset)
    upper_band = kernel._affine_sub((Fraction(0), Fraction(0), U), common)
    if large:
        lowers = ((Fraction(0), Fraction(0), DELTA),)
        uppers = [(Fraction(0), Fraction(0), Fraction(1)), upper_band]
        if support_limit is not None:
            uppers.append(
                kernel._affine_sub(
                    (Fraction(0), Fraction(0), support_limit), large_sum
                )
            )
    else:
        if support_limit is not None and kernel._affine_value(large_sum, sample) > support_limit:
            return {}
        lowers = ((Fraction(0), Fraction(0), Fraction(0)),)
        uppers = [(Fraction(0), Fraction(0), DELTA), upper_band]
    lower = max(lowers, key=lambda form: kernel._affine_value(form, sample))
    upper = min(uppers, key=lambda form: kernel._affine_value(form, sample))
    if kernel._affine_value(upper, sample) <= kernel._affine_value(lower, sample):
        return {}
    room = (-common[0], -common[1], U - common[2])
    transformed: dict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for (power, slack), coefficient in coefficients.items():
        for j in range(slack + 1):
            endpoint_power = power + j + 1
            transformed[(slack - j, endpoint_power)] += (
                coefficient * (-1) ** j * comb(slack, j) / endpoint_power
            )
    answer: Poly2 = {}
    for (room_power, endpoint_power), coefficient in transformed.items():
        endpoint = kernel._poly_add(
            _affine_product_power(room, room_power, upper, endpoint_power),
            kernel._poly_scale(
                _affine_product_power(room, room_power, lower, endpoint_power), -1
            ),
        )
        answer = kernel._poly_add(answer, kernel._poly_scale(endpoint, coefficient))
    return answer


@cache
def _affine_product_power(first, first_power: int, second, second_power: int) -> Poly2:
    return kernel._poly_mul(
        _affine_power(first, first_power), _affine_power(second, second_power)
    )


def _integrate_on_geometry(kind: str, cells, polynomial: Poly2) -> Fraction:
    if kind == "empty" or not polynomial:
        return Fraction(0)
    if kind == "point":
        return polynomial.get((0, 0), Fraction(0))
    if kind in ("xintervals", "zintervals"):
        variable_is_x = kind == "xintervals"
        answer = Fraction(0)
        for start, end, _sample in cells:
            for (x_power, z_power), coefficient in polynomial.items():
                if variable_is_x and z_power == 0:
                    degree = x_power
                elif not variable_is_x and x_power == 0:
                    degree = z_power
                else:
                    continue
                answer += coefficient * (
                    end ** (degree + 1) - start ** (degree + 1)
                ) / (degree + 1)
        return answer
    return sum(
        _integrate_polynomial_polygon(polynomial, polygon)
        for polygon, _sample in cells
    )


def exact_j_signature_group(
    signature: Signature,
    pair_contributions: tuple[tuple[Signature, Signature, int], ...],
    feature_groups: dict[Signature, dict[tuple[int, int], Fraction]],
    dimension: int,
) -> Fraction:
    """Evaluate all candidate contractions for one common product signature."""
    common_dimension = dimension - 1
    answer = Fraction(0)
    for large_count in range(min(common_dimension, len(B)) + 1):
        by_shift: dict[int, Poly2] = defaultdict(dict)
        for (shifted_count, large_power, small_power), coefficient in orbit_status_densities(
            common_dimension, signature, large_count
        ):
            if large_count + shifted_count > int(R // DELTA):
                continue
            exponent = (
                0 if large_power is None else large_power,
                0 if small_power is None else small_power,
            )
            bucket = by_shift.setdefault(shifted_count, {})
            bucket[exponent] = bucket.get(exponent, Fraction(0)) + coefficient
        for shifted_count, density_polynomial in by_shift.items():
            total_offset = (large_count + shifted_count) * DELTA
            large_offset = large_count * DELTA
            for left_large in (False, True):
                if large_count + int(left_large) > len(B):
                    continue
                left_limit = _support_limit(large_count, left_large)
                for right_large in (False, True):
                    if large_count + int(right_large) > len(B):
                        continue
                    right_limit = _support_limit(large_count, right_large)
                    # Any representative specs produce the universal line
                    # arrangement: powers/slacks do not affect its cuts.
                    geometry_specs = (
                        RadialSlice(0, 0, left_large, support_limit=left_limit),
                        RadialSlice(0, 0, right_large, support_limit=right_limit),
                    )
                    kind, cells = _slice_geometry(
                        large_count > 0,
                        common_dimension > large_count,
                        total_offset,
                        large_offset,
                        geometry_specs,
                    )
                    for cell in cells if kind != "point" else (cells[0],):
                        if kind == "polygons":
                            _polygon, sample = cell
                        elif kind in ("xintervals", "zintervals"):
                            _start, _end, sample = cell
                        else:
                            sample = cell
                        left_cache: dict[Signature, Poly2] = {}
                        right_cache: dict[Signature, Poly2] = {}
                        candidate_polynomial: Poly2 = {}
                        for left_signature, right_signature, structure in pair_contributions:
                            left_poly = left_cache.get(left_signature)
                            if left_poly is None:
                                left_poly = _linear_slice_polynomial(
                                    feature_groups[left_signature], left_large, left_limit,
                                    total_offset, large_offset, sample
                                )
                                left_cache[left_signature] = left_poly
                            right_poly = right_cache.get(right_signature)
                            if right_poly is None:
                                right_poly = _linear_slice_polynomial(
                                    feature_groups[right_signature], right_large, right_limit,
                                    total_offset, large_offset, sample
                                )
                                right_cache[right_signature] = right_poly
                            candidate_polynomial = kernel._poly_add(
                                candidate_polynomial,
                                kernel._poly_scale(
                                    kernel._poly_mul(left_poly, right_poly), structure
                                ),
                            )
                        integrand = kernel._poly_mul(
                            density_polynomial, candidate_polynomial
                        )
                        if kind == "polygons":
                            answer += _integrate_polynomial_polygon(integrand, cell[0])
                        elif kind in ("xintervals", "zintervals"):
                            # Integrate only this interval, not every cell.
                            start, end, _ = cell
                            variable_is_x = kind == "xintervals"
                            for (xp, zp), coefficient in integrand.items():
                                if variable_is_x and zp == 0:
                                    degree = xp
                                elif not variable_is_x and xp == 0:
                                    degree = zp
                                else:
                                    continue
                                answer += coefficient * (
                                    end ** (degree + 1) - start ** (degree + 1)
                                ) / (degree + 1)
                        else:
                            answer += integrand.get((0, 0), Fraction(0))
    return answer
