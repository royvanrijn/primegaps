"""Exact Section 5 integral machinery for Stadlmann's sieve support.

The public API consumes sparse polynomials and rational support parameters and
returns matrices whose entries are :class:`fractions.Fraction` objects.  It
does not perform an eigenvalue calculation or choose support parameters.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import comb, factorial, floor, gcd, lcm

from .support import SupportParameters

Rational = int | Fraction
Exponent = tuple[int, ...]
SparsePolynomial = Mapping[Exponent, Rational]
_Poly2 = dict[tuple[int, int], Fraction]
_Affine = tuple[Fraction, Fraction, Fraction]
_Point = tuple[Fraction, Fraction]


def _fraction(value: int | float | str | Fraction) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, float):
        return Fraction(str(value))
    return Fraction(value)


@dataclass(frozen=True)
class ExactSupportParameters:
    """An exact-rational version of Definition 1's support parameters."""

    delta: Fraction
    epsilon: Fraction
    A: tuple[Fraction, ...]
    B: tuple[tuple[Fraction, ...], ...]

    @classmethod
    def from_values(
        cls,
        *,
        delta: int | float | str | Fraction,
        epsilon: int | float | str | Fraction,
        A: Sequence[int | float | str | Fraction],
        B: Sequence[Sequence[int | float | str | Fraction]],
    ) -> ExactSupportParameters:
        result = cls(
            delta=_fraction(delta),
            epsilon=_fraction(epsilon),
            A=tuple(_fraction(x) for x in A),
            B=tuple(tuple(_fraction(x) for x in row) for row in B),
        )
        result.validate()
        return result

    @classmethod
    def from_support_parameters(cls, support: SupportParameters) -> ExactSupportParameters:
        """Convert decimal-looking floats through their strings, not binary values."""
        return cls.from_values(
            delta=support.delta,
            epsilon=support.epsilon,
            A=support.A,
            B=support.B,
        )

    def validate(self) -> None:
        if not 0 < self.delta < 1:
            raise ValueError("delta must lie strictly between zero and one")
        if self.epsilon < 0:
            raise ValueError("epsilon must be non-negative")
        n = len(self.A) - 1
        if n <= 0 or len(self.B) != n:
            raise ValueError("B must have len(A)-1 rows")
        if self.A[0] != -self.epsilon:
            raise ValueError("Definition 1 requires A[0] = -epsilon")
        if any(x >= y for x, y in zip(self.A, self.A[1:])):
            raise ValueError("A must be strictly increasing")
        if self.A[-1] + self.epsilon >= 1:
            raise ValueError("the exact engine requires the published total bound below 1")
        width = floor(Fraction(1, 1) / self.delta)
        for row in self.B:
            if len(row) != width:
                raise ValueError(f"B rows must have {width} entries")
            for index, value in enumerate(row):
                if value <= self.delta:
                    raise ValueError("B[j,m] must exceed delta")
                if index and not row[index - 1] <= value <= row[index - 1] + self.delta:
                    raise ValueError("B rows must be monotone and grow by at most delta")


@dataclass(frozen=True)
class IntegralMatrices:
    """Exact bilinear matrices for the printed Section 2 ``I``, ``J``, ``K``."""

    I: tuple[tuple[Fraction, ...], ...]
    J: tuple[tuple[Fraction, ...], ...]
    K: tuple[tuple[Fraction, ...], ...]


def monomial(exponents: Sequence[int], coefficient: Rational = 1) -> dict[Exponent, Fraction]:
    """Construct a one-term sparse polynomial accepted by the matrix API."""
    exponent = tuple(exponents)
    if any(value < 0 for value in exponent):
        raise ValueError("exponents must be non-negative")
    value = _fraction(coefficient)
    return {} if value == 0 else {exponent: value}


def exact_ijk_matrices(
    basis: Sequence[SparsePolynomial], support: ExactSupportParameters
) -> IntegralMatrices:
    """Return exact ``I``, ``J`` and ``K`` matrices for a polynomial basis.

    Polynomials are mappings from exponent tuples to rational coefficients.
    The implementation is deliberately exact and reference-oriented: its
    status decomposition is appropriate for low-dimensional reconstruction
    and validation, while the standalone C/D coefficient API remains the
    scalable building block for a later symmetry-compressed assembly.

    ``K`` follows the only dimensionally consistent reading of the printed
    definition: the undifferentiated ``t'_k`` condition is existential.  This
    convention is immaterial for Stadlmann's published certificate, where
    ``c2 = 0``.
    """
    support.validate()
    polynomials, k = _normalize_basis(basis)
    i_matrix = _symmetric_matrix(polynomials, lambda f, g: exact_i_entry(f, g, support))
    j_matrix = _symmetric_matrix(polynomials, lambda f, g: exact_j_entry(f, g, support, k))
    k_matrix = _symmetric_matrix(polynomials, lambda f, g: exact_k_entry(f, g, support, k))
    return IntegralMatrices(i_matrix, j_matrix, k_matrix)


def exact_i_entry(
    left: SparsePolynomial,
    right: SparsePolynomial,
    support: ExactSupportParameters,
) -> Fraction:
    """Exact bilinear ``I`` entry over ``T_k(delta,A,B,epsilon)``."""
    left_poly, k = _normalize_polynomial(left)
    right_poly, right_k = _normalize_polynomial(right)
    if k != right_k:
        raise ValueError("polynomial dimensions differ")
    result = Fraction(0)
    for a, ca in left_poly.items():
        for b, cb in right_poly.items():
            exponent = tuple(x + y for x, y in zip(a, b))
            result += ca * cb * _integrate_support_monomial(exponent, support)
    return result


def exact_j_entry(
    left: SparsePolynomial,
    right: SparsePolynomial,
    support: ExactSupportParameters,
    k: int | None = None,
) -> Fraction:
    """Exact bilinear ``J`` entry, with the last coordinate distinguished."""
    left_poly, left_k = _normalize_polynomial(left)
    right_poly, right_k = _normalize_polynomial(right)
    if left_k != right_k or (k is not None and left_k != k):
        raise ValueError("polynomial dimensions differ")
    result = Fraction(0)
    for a, ca in left_poly.items():
        for b, cb in right_poly.items():
            common = tuple(x + y for x, y in zip(a[:-1], b[:-1]))
            result += ca * cb * _j_monomial(common, a[-1], b[-1], support)
    return result


def exact_k_entry(
    left: SparsePolynomial,
    right: SparsePolynomial,
    support: ExactSupportParameters,
    k: int | None = None,
) -> Fraction:
    """Exact bilinear ``K`` entry under the documented existential convention."""
    left_poly, left_k = _normalize_polynomial(left)
    right_poly, right_k = _normalize_polynomial(right)
    if left_k != right_k or (k is not None and left_k != k):
        raise ValueError("polynomial dimensions differ")
    result = Fraction(0)
    for a, ca in left_poly.items():
        for b, cb in right_poly.items():
            common = tuple(x + y for x, y in zip(a[:-1], b[:-1]))
            result += ca * cb * _k_monomial(common, a[-1] + b[-1], support)
    return result


def _normalize_basis(
    basis: Sequence[SparsePolynomial],
) -> tuple[tuple[dict[Exponent, Fraction], ...], int]:
    if not basis:
        raise ValueError("basis must be non-empty")
    normalized = []
    dimension: int | None = None
    for polynomial in basis:
        current, k = _normalize_polynomial(polynomial)
        if dimension is None:
            dimension = k
        elif dimension != k:
            raise ValueError("all basis polynomials must have the same dimension")
        normalized.append(current)
    assert dimension is not None
    return tuple(normalized), dimension


def _normalize_polynomial(polynomial: SparsePolynomial) -> tuple[dict[Exponent, Fraction], int]:
    if not polynomial:
        raise ValueError("zero polynomials do not carry a dimension")
    result: dict[Exponent, Fraction] = {}
    dimension: int | None = None
    for raw_exponent, raw_coefficient in polynomial.items():
        exponent = tuple(raw_exponent)
        if dimension is None:
            dimension = len(exponent)
        elif len(exponent) != dimension:
            raise ValueError("inconsistent exponent dimensions")
        if not exponent:
            raise ValueError("polynomials must have at least one variable")
        if any(value < 0 for value in exponent):
            raise ValueError("exponents must be non-negative")
        coefficient = _fraction(raw_coefficient)
        if coefficient:
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    result = {exponent: value for exponent, value in result.items() if value}
    if not result or dimension is None:
        raise ValueError("zero polynomials are not valid basis functions")
    return result, dimension


def _symmetric_matrix(polynomials, entry):
    size = len(polynomials)
    result = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(i, size):
            value = entry(polynomials[i], polynomials[j])
            result[i][j] = value
            result[j][i] = value
    return tuple(tuple(row) for row in result)


def _integrate_support_monomial(
    exponents: Exponent, support: ExactSupportParameters
) -> Fraction:
    k = len(exponents)
    indices = tuple(range(k))
    answer = Fraction(0)
    for row_index, row in enumerate(support.B):
        lower = support.A[row_index] + support.epsilon
        upper = support.A[row_index + 1] + support.epsilon
        for large_count in range(k + 1):
            if large_count and large_count > len(row):
                continue
            cap = None if large_count == 0 else row[large_count - 1]
            for large in combinations(indices, large_count):
                answer += _mixed_status_cap(exponents, large, upper, cap, support.delta)
                answer -= _mixed_status_cap(exponents, large, lower, cap, support.delta)
    return answer


def _mixed_status_cap(
    exponents: Exponent,
    large: tuple[int, ...],
    total_cap: Fraction,
    large_cap: Fraction | None,
    delta: Fraction,
) -> Fraction:
    """Monomial integral for one large/small status and a total-sum cap."""
    if total_cap <= 0:
        return Fraction(0)
    k = len(exponents)
    large_set = frozenset(large)
    small = tuple(index for index in range(k) if index not in large_set)
    result = Fraction(0)
    # Inclusion--exclusion enforces the upper bound on every small coordinate.
    for shifted_small_count in range(len(small) + 1):
        for shifted_small in combinations(small, shifted_small_count):
            shifted = large_set | frozenset(shifted_small)
            sign = -1 if shifted_small_count % 2 else 1
            offset = (len(large) + shifted_small_count) * delta
            height = total_cap - offset
            if height <= 0:
                continue
            residual_large_cap = None if large_cap is None else large_cap - len(large) * delta
            if residual_large_cap is not None and residual_large_cap <= 0:
                continue
            for powers, coefficient in _shift_expansions(exponents, shifted, delta):
                result += sign * coefficient * _group_cap_integral(
                    powers, large_set, height, residual_large_cap
                )
    return result


def _shift_expansions(
    exponents: Exponent, shifted: frozenset[int], delta: Fraction
):
    choices = []
    for index, exponent in enumerate(exponents):
        if index in shifted:
            choices.append(
                tuple(
                    (power, Fraction(comb(exponent, power)) * delta ** (exponent - power))
                    for power in range(exponent + 1)
                )
            )
        else:
            choices.append(((exponent, Fraction(1)),))
    for selected in product(*choices):
        coefficient = Fraction(1)
        powers = []
        for power, multiplier in selected:
            powers.append(power)
            coefficient *= multiplier
        yield tuple(powers), coefficient


def _group_cap_integral(
    powers: Exponent,
    large: frozenset[int],
    height: Fraction,
    large_cap: Fraction | None,
) -> Fraction:
    large_powers = tuple(powers[index] for index in range(len(powers)) if index in large)
    small_powers = tuple(powers[index] for index in range(len(powers)) if index not in large)
    if not large_powers:
        return _simplex_monomial(small_powers, height)
    maximum = height if large_cap is None else min(height, large_cap)
    if maximum <= 0:
        return Fraction(0)
    if not small_powers:
        return _simplex_monomial(large_powers, maximum)

    x_power = sum(large_powers) + len(large_powers) - 1
    y_power = sum(small_powers) + len(small_powers) - 1
    density = Fraction(
        _factorial_product(large_powers) * _factorial_product(small_powers),
        factorial(x_power) * factorial(y_power),
    )
    integral = Fraction(0)
    for j in range(y_power + 2):
        integral += Fraction(
            (-1) ** j * comb(y_power + 1, j),
            (y_power + 1) * (x_power + j + 1),
        ) * height ** (y_power + 1 - j) * maximum ** (x_power + j + 1)
    return density * integral


def _simplex_monomial(powers: Exponent, height: Fraction) -> Fraction:
    if height < 0:
        return Fraction(0)
    if not powers:
        return Fraction(1)
    degree = len(powers) + sum(powers)
    return Fraction(_factorial_product(powers), factorial(degree)) * height**degree


def _factorial_product(values: Sequence[int]) -> int:
    result = 1
    for value in values:
        result *= factorial(value)
    return result


@dataclass(frozen=True)
class _SliceSpec:
    power: int
    large: bool
    lower: Fraction
    upper: Fraction
    support_limit: Fraction | None


def _j_monomial(
    common_exponents: Exponent,
    left_power: int,
    right_power: int,
    support: ExactSupportParameters,
) -> Fraction:
    answer = Fraction(0)
    common_count = len(common_exponents)
    indices = tuple(range(common_count))
    for left_band, left_row in enumerate(support.B):
        left_interval = (
            support.A[left_band] + support.epsilon,
            support.A[left_band + 1] + support.epsilon,
        )
        for right_band, right_row in enumerate(support.B):
            right_interval = (
                support.A[right_band] + support.epsilon,
                support.A[right_band + 1] + support.epsilon,
            )
            common_upper = max(
                support.A[left_band + 1] - support.epsilon,
                support.A[right_band + 1] - support.epsilon,
            )
            for large_count in range(common_count + 1):
                for large in combinations(indices, large_count):
                    for left_large in (False, True):
                        left_limit, left_valid = _status_limit(left_row, large_count, left_large)
                        if not left_valid:
                            continue
                        for right_large in (False, True):
                            right_limit, right_valid = _status_limit(
                                right_row, large_count, right_large
                            )
                            if not right_valid:
                                continue
                            specs = (
                                _SliceSpec(
                                    left_power,
                                    left_large,
                                    left_interval[0],
                                    left_interval[1],
                                    left_limit,
                                ),
                                _SliceSpec(
                                    right_power,
                                    right_large,
                                    right_interval[0],
                                    right_interval[1],
                                    right_limit,
                                ),
                            )
                            answer += _common_status_slices(
                                common_exponents,
                                large,
                                Fraction(0),
                                common_upper,
                                specs,
                                support.delta,
                            )
    return answer


def _k_monomial(
    common_exponents: Exponent, special_power: int, support: ExactSupportParameters
) -> Fraction:
    answer = Fraction(0)
    common_count = len(common_exponents)
    indices = tuple(range(common_count))
    for band, row in enumerate(support.B):
        interval = (
            support.A[band] + support.epsilon,
            support.A[band + 1] + support.epsilon,
        )
        for witness_band in range(len(support.B)):
            common_lower = max(
                support.A[band + 1] - support.epsilon,
                support.A[witness_band + 1] - support.epsilon,
            )
            # Existence of t'_k in its printed band is equivalent (up to
            # measure-zero boundaries) to common_sum < that band's upper end.
            common_upper = support.A[witness_band + 1] + support.epsilon
            if common_upper <= common_lower:
                continue
            for large_count in range(common_count + 1):
                for large in combinations(indices, large_count):
                    for special_large in (False, True):
                        limit, valid = _status_limit(row, large_count, special_large)
                        if not valid:
                            continue
                        spec = _SliceSpec(
                            special_power,
                            special_large,
                            interval[0],
                            interval[1],
                            limit,
                        )
                        answer += _common_status_slices(
                            common_exponents,
                            large,
                            common_lower,
                            common_upper,
                            (spec,),
                            support.delta,
                        )
    return answer


def _status_limit(
    row: tuple[Fraction, ...], common_large_count: int, special_large: bool
) -> tuple[Fraction | None, bool]:
    total_large = common_large_count + int(special_large)
    if total_large == 0:
        return None, True
    if total_large > len(row):
        return None, False
    return row[total_large - 1], True


def _common_status_slices(
    exponents: Exponent,
    large: tuple[int, ...],
    common_lower: Fraction,
    common_upper: Fraction,
    specs: tuple[_SliceSpec, ...],
    delta: Fraction,
) -> Fraction:
    large_set = frozenset(large)
    small = tuple(index for index in range(len(exponents)) if index not in large_set)
    result = Fraction(0)
    for shifted_small_count in range(len(small) + 1):
        for shifted_small in combinations(small, shifted_small_count):
            shifted = large_set | frozenset(shifted_small)
            sign = -1 if shifted_small_count % 2 else 1
            total_offset = (len(large) + shifted_small_count) * delta
            large_offset = len(large) * delta
            for powers, coefficient in _shift_expansions(exponents, shifted, delta):
                large_powers = tuple(
                    powers[index] for index in range(len(powers)) if index in large_set
                )
                small_powers = tuple(
                    powers[index] for index in range(len(powers)) if index not in large_set
                )
                large_sum_power, large_density = _sum_density(large_powers)
                small_sum_power, small_density = _sum_density(small_powers)
                result += (
                    sign
                    * coefficient
                    * large_density
                    * small_density
                    * _integrate_slice_cells(
                        large_sum_power,
                        small_sum_power,
                        common_lower,
                        common_upper,
                        total_offset,
                        large_offset,
                        specs,
                        delta,
                    )
                )
    return result


def _sum_density(powers: Exponent) -> tuple[int | None, Fraction]:
    if not powers:
        return None, Fraction(1)
    power = sum(powers) + len(powers) - 1
    return power, Fraction(_factorial_product(powers), factorial(power))


def _integrate_slice_cells(
    x_power: int | None,
    z_power: int | None,
    common_lower: Fraction,
    common_upper: Fraction,
    total_offset: Fraction,
    large_offset: Fraction,
    specs: tuple[_SliceSpec, ...],
    delta: Fraction,
) -> Fraction:
    residual_lower = max(Fraction(0), common_lower - total_offset)
    residual_upper = common_upper - total_offset
    if residual_upper <= residual_lower:
        return Fraction(0)
    lines = _unique_lines(
        line
        for spec in specs
        for line in _slice_lines(spec, total_offset, large_offset, delta)
    )
    # Add the lower common-sum face; the upper face is the initial simplex.
    if residual_lower > 0:
        lines = _unique_lines((*lines, (Fraction(1), Fraction(1), -residual_lower)))

    if x_power is None and z_power is None:
        if residual_lower > 0:
            return Fraction(0)
        sample = (Fraction(0), Fraction(0))
        value = Fraction(1)
        for spec in specs:
            poly = _slice_polynomial(spec, total_offset, large_offset, delta, sample)
            value *= poly.get((0, 0), Fraction(0))
        return value

    if x_power is None or z_power is None:
        variable_is_x = x_power is not None
        power = x_power if variable_is_x else z_power
        assert power is not None
        return _integrate_slice_intervals(
            variable_is_x,
            power,
            residual_lower,
            residual_upper,
            lines,
            total_offset,
            large_offset,
            specs,
            delta,
        )

    polygons: list[list[_Point]] = [
        [
            (Fraction(0), Fraction(0)),
            (residual_upper, Fraction(0)),
            (Fraction(0), residual_upper),
        ]
    ]
    for line in lines:
        next_polygons = []
        for polygon in polygons:
            positive, negative = _split_polygon(polygon, line)
            if positive:
                next_polygons.append(positive)
            if negative:
                next_polygons.append(negative)
        polygons = next_polygons

    result = Fraction(0)
    for polygon in polygons:
        sample = _polygon_centroid(polygon)
        if sample[0] + sample[1] < residual_lower:
            continue
        integrand: _Poly2 = {(x_power, z_power): Fraction(1)}
        for spec in specs:
            integrand = _poly_mul(
                integrand,
                _slice_polynomial(spec, total_offset, large_offset, delta, sample),
            )
            if not integrand:
                break
        result += _integrate_polygon(integrand, polygon)
    return result


def _integrate_slice_intervals(
    variable_is_x: bool,
    power: int,
    lower: Fraction,
    upper: Fraction,
    lines: tuple[_Affine, ...],
    total_offset: Fraction,
    large_offset: Fraction,
    specs: tuple[_SliceSpec, ...],
    delta: Fraction,
) -> Fraction:
    cuts = {lower, upper}
    for a, b, c in lines:
        slope = a if variable_is_x else b
        if slope:
            root = -c / slope
            if lower < root < upper:
                cuts.add(root)
    ordered = sorted(cuts)
    result = Fraction(0)
    for start, end in zip(ordered, ordered[1:]):
        sample_value = (start + end) / 2
        sample = (sample_value, Fraction(0)) if variable_is_x else (
            Fraction(0),
            sample_value,
        )
        polynomial: _Poly2 = {(power, 0) if variable_is_x else (0, power): Fraction(1)}
        for spec in specs:
            polynomial = _poly_mul(
                polynomial,
                _slice_polynomial(spec, total_offset, large_offset, delta, sample),
            )
        one_variable: dict[int, Fraction] = {}
        for (x_degree, z_degree), coefficient in polynomial.items():
            if variable_is_x and z_degree == 0:
                one_variable[x_degree] = one_variable.get(x_degree, Fraction(0)) + coefficient
            elif not variable_is_x and x_degree == 0:
                one_variable[z_degree] = one_variable.get(z_degree, Fraction(0)) + coefficient
        for degree, coefficient in one_variable.items():
            result += coefficient * (end ** (degree + 1) - start ** (degree + 1)) / (
                degree + 1
            )
    return result


def _slice_lines(
    spec: _SliceSpec, total_offset: Fraction, large_offset: Fraction, delta: Fraction
):
    common: _Affine = (Fraction(1), Fraction(1), total_offset)
    large_sum: _Affine = (Fraction(1), Fraction(0), large_offset)
    lower_band = _affine_sub((Fraction(0), Fraction(0), spec.lower), common)
    upper_band = _affine_sub((Fraction(0), Fraction(0), spec.upper), common)
    if spec.large:
        lowers = ((Fraction(0), Fraction(0), delta), lower_band)
        uppers = [(Fraction(0), Fraction(0), Fraction(1)), upper_band]
        if spec.support_limit is not None:
            uppers.append(
                _affine_sub(
                    (Fraction(0), Fraction(0), spec.support_limit), large_sum
                )
            )
    else:
        lowers = ((Fraction(0), Fraction(0), Fraction(0)), lower_band)
        uppers = [(Fraction(0), Fraction(0), delta), upper_band]
        if spec.support_limit is not None:
            yield _affine_sub(
                (Fraction(0), Fraction(0), spec.support_limit), large_sum
            )
    for first, second in combinations(lowers, 2):
        yield _affine_sub(first, second)
    for first, second in combinations(uppers, 2):
        yield _affine_sub(first, second)
    for lower in lowers:
        for upper in uppers:
            yield _affine_sub(upper, lower)


def _slice_polynomial(
    spec: _SliceSpec,
    total_offset: Fraction,
    large_offset: Fraction,
    delta: Fraction,
    sample: _Point,
) -> _Poly2:
    common: _Affine = (Fraction(1), Fraction(1), total_offset)
    large_sum: _Affine = (Fraction(1), Fraction(0), large_offset)
    lower_band = _affine_sub((Fraction(0), Fraction(0), spec.lower), common)
    upper_band = _affine_sub((Fraction(0), Fraction(0), spec.upper), common)
    if spec.large:
        lowers = ((Fraction(0), Fraction(0), delta), lower_band)
        uppers = [(Fraction(0), Fraction(0), Fraction(1)), upper_band]
        if spec.support_limit is not None:
            uppers.append(
                _affine_sub(
                    (Fraction(0), Fraction(0), spec.support_limit), large_sum
                )
            )
    else:
        if spec.support_limit is not None and _affine_value(large_sum, sample) > spec.support_limit:
            return {}
        lowers = ((Fraction(0), Fraction(0), Fraction(0)), lower_band)
        uppers = [(Fraction(0), Fraction(0), delta), upper_band]
    lower = max(lowers, key=lambda form: _affine_value(form, sample))
    upper = min(uppers, key=lambda form: _affine_value(form, sample))
    if _affine_value(upper, sample) <= _affine_value(lower, sample):
        return {}
    exponent = spec.power + 1
    return _poly_scale(
        _poly_add(_affine_power(upper, exponent), _poly_scale(_affine_power(lower, exponent), -1)),
        Fraction(1, exponent),
    )


def _affine_sub(left: _Affine, right: _Affine) -> _Affine:
    return tuple(x - y for x, y in zip(left, right))  # type: ignore[return-value]


def _affine_value(form: _Affine, point: _Point) -> Fraction:
    return form[0] * point[0] + form[1] * point[1] + form[2]


def _affine_power(form: _Affine, exponent: int) -> _Poly2:
    result: _Poly2 = {(0, 0): Fraction(1)}
    base: _Poly2 = {}
    if form[0]:
        base[(1, 0)] = form[0]
    if form[1]:
        base[(0, 1)] = form[1]
    if form[2]:
        base[(0, 0)] = form[2]
    for _ in range(exponent):
        result = _poly_mul(result, base)
    return result


def _poly_add(left: _Poly2, right: _Poly2) -> _Poly2:
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, Fraction(0)) + coefficient
        if not result[exponent]:
            del result[exponent]
    return result


def _poly_scale(polynomial: _Poly2, scalar: Rational) -> _Poly2:
    value = _fraction(scalar)
    return {exponent: coefficient * value for exponent, coefficient in polynomial.items() if coefficient * value}


def _poly_mul(left: _Poly2, right: _Poly2) -> _Poly2:
    result: _Poly2 = {}
    for (ax, az), ca in left.items():
        for (bx, bz), cb in right.items():
            exponent = (ax + bx, az + bz)
            result[exponent] = result.get(exponent, Fraction(0)) + ca * cb
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def _unique_lines(lines) -> tuple[_Affine, ...]:
    result = []
    seen = set()
    for line in lines:
        if line[0] == line[1] == 0:
            continue
        key = _line_key(line)
        if key not in seen:
            seen.add(key)
            result.append(line)
    return tuple(result)


def _line_key(line: _Affine) -> tuple[int, int, int]:
    denominator = 1
    for value in line:
        denominator = lcm(denominator, value.denominator)
    integers = [value.numerator * (denominator // value.denominator) for value in line]
    divisor = 0
    for value in integers:
        divisor = gcd(divisor, abs(value))
    integers = [value // divisor for value in integers]
    for value in integers:
        if value:
            if value < 0:
                integers = [-item for item in integers]
            break
    return tuple(integers)  # type: ignore[return-value]


def _split_polygon(
    polygon: list[_Point], line: _Affine
) -> tuple[list[_Point], list[_Point]]:
    values = [_affine_value(line, point) for point in polygon]
    if all(value >= 0 for value in values):
        return polygon, []
    if all(value <= 0 for value in values):
        return [], polygon
    return _clip_polygon(polygon, line, True), _clip_polygon(polygon, line, False)


def _clip_polygon(polygon: list[_Point], line: _Affine, keep_positive: bool) -> list[_Point]:
    output: list[_Point] = []
    for current, following in zip(polygon, polygon[1:] + polygon[:1]):
        current_value = _affine_value(line, current)
        following_value = _affine_value(line, following)
        current_inside = current_value >= 0 if keep_positive else current_value <= 0
        following_inside = following_value >= 0 if keep_positive else following_value <= 0
        if current_inside:
            output.append(current)
        if current_inside != following_inside:
            ratio = current_value / (current_value - following_value)
            output.append(
                (
                    current[0] + ratio * (following[0] - current[0]),
                    current[1] + ratio * (following[1] - current[1]),
                )
            )
    return _deduplicate_polygon(output)


def _deduplicate_polygon(polygon: list[_Point]) -> list[_Point]:
    result = []
    for point in polygon:
        if not result or point != result[-1]:
            result.append(point)
    if len(result) > 1 and result[0] == result[-1]:
        result.pop()
    return result if len(result) >= 3 else []


def _polygon_centroid(polygon: list[_Point]) -> _Point:
    count = len(polygon)
    return (
        sum(point[0] for point in polygon) / count,
        sum(point[1] for point in polygon) / count,
    )


def _integrate_polygon(polynomial: _Poly2, polygon: list[_Point]) -> Fraction:
    result = Fraction(0)
    for (x_power, z_power), coefficient in polynomial.items():
        boundary = Fraction(0)
        for start, end in zip(polygon, polygon[1:] + polygon[:1]):
            dx = end[0] - start[0]
            dz = end[1] - start[1]
            edge = Fraction(0)
            for i in range(x_power + 2):
                for j in range(z_power + 1):
                    edge += (
                        comb(x_power + 1, i)
                        * comb(z_power, j)
                        * start[0] ** (x_power + 1 - i)
                        * dx**i
                        * start[1] ** (z_power - j)
                        * dz**j
                        / (i + j + 1)
                    )
            boundary += dz * edge / (x_power + 1)
        result += coefficient * boundary
    return result
