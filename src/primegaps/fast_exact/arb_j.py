"""Boundary-only, cell-first rigorous interval contraction for exact J.

The discovery operator is floating point, while the final certificate only
needs a rigorous lower bound for one scalar.  This module keeps all polynomial
and geometry inputs exact, but performs the final dot products in Arb balls.
Geometry moments are constructed once per affected cell and shared by every
target density/candidate correlation on that cell.
"""

from __future__ import annotations

from dataclasses import dataclass
import time

from . import fast_j


def _ball_from_rational(field, value):
    numerator = value.numerator
    denominator = value.denominator
    if callable(numerator):
        numerator = numerator()
    if callable(denominator):
        denominator = denominator()
    return field(int(numerator)) / field(int(denominator))


@dataclass(frozen=True)
class BoundaryCell:
    large: int
    shifted: int
    left_large: bool
    right_large: bool
    left_legal: bool
    right_legal: bool
    left_limit: object | None
    right_limit: object | None
    total_offset: object
    large_offset: object
    kind: str
    cell: object
    sample: object

    @property
    def density_status(self):
        return (self.large, self.shifted)


def iter_boundary_cells(*, common_dimension: int, verifier):
    """Yield exactly the legal/unrestricted cells changed by a B cutoff."""
    common_dimension = int(common_dimension)
    maximum_offset = int(verifier.R // verifier.DELTA)
    for large in range(min(common_dimension, maximum_offset) + 1):
        for shifted in range(maximum_offset - large + 1):
            total_offset = (large + shifted) * verifier.DELTA
            large_offset = large * verifier.DELTA
            for left_large in (False, True):
                left_legal = large + int(left_large) <= len(verifier.B)
                left_limit = (
                    verifier._support_limit(large, left_large)
                    if left_legal else None
                )
                for right_large in (False, True):
                    right_legal = large + int(right_large) <= len(verifier.B)
                    right_limit = (
                        verifier._support_limit(large, right_large)
                        if right_legal else None
                    )
                    specs = tuple(
                        spec
                        for allowed, spec in (
                            (
                                left_legal,
                                verifier.RadialSlice(
                                    0, 0, left_large, support_limit=left_limit
                                ),
                            ),
                            (
                                right_legal,
                                verifier.RadialSlice(
                                    0, 0, right_large, support_limit=right_limit
                                ),
                            ),
                        )
                        if allowed
                    ) + (
                        verifier.RadialSlice(0, 0, left_large),
                        verifier.RadialSlice(0, 0, right_large),
                    )
                    kind, cells = verifier._slice_geometry(
                        large > 0,
                        common_dimension > large,
                        total_offset,
                        large_offset,
                        specs,
                    )
                    iterable = cells if kind != "point" else (cells[0],)
                    for cell in iterable:
                        if kind == "polygons":
                            _polygon, sample = cell
                        elif kind in ("xintervals", "zintervals"):
                            _start, _end, sample = cell
                        elif kind == "point":
                            sample = cell
                        else:
                            continue
                        unit_coefficients = {(0, 0): verifier.U / verifier.U}
                        legal_left = (
                            verifier._linear_slice_polynomial(
                                unit_coefficients,
                                left_large,
                                left_limit,
                                total_offset,
                                large_offset,
                                sample,
                            )
                            if left_legal else {}
                        )
                        legal_right = (
                            verifier._linear_slice_polynomial(
                                unit_coefficients,
                                right_large,
                                right_limit,
                                total_offset,
                                large_offset,
                                sample,
                            )
                            if right_legal else {}
                        )
                        full_left = verifier._linear_slice_polynomial(
                            unit_coefficients,
                            left_large,
                            None,
                            total_offset,
                            large_offset,
                            sample,
                        )
                        full_right = verifier._linear_slice_polynomial(
                            unit_coefficients,
                            right_large,
                            None,
                            total_offset,
                            large_offset,
                            sample,
                        )
                        if legal_left == full_left and legal_right == full_right:
                            continue
                        yield BoundaryCell(
                            large,
                            shifted,
                            left_large,
                            right_large,
                            left_legal,
                            right_legal,
                            left_limit,
                            right_limit,
                            total_offset,
                            large_offset,
                            kind,
                            cell,
                            sample,
                        )


def _slice_families(feature_groups, boundary_cell, *, verifier):
    """Construct each candidate marginal once on one geometry cell."""
    result = {}
    for signature, coefficients in feature_groups.items():
        legal_left = (
            verifier._linear_slice_polynomial(
                coefficients,
                boundary_cell.left_large,
                boundary_cell.left_limit,
                boundary_cell.total_offset,
                boundary_cell.large_offset,
                boundary_cell.sample,
            )
            if boundary_cell.left_legal else {}
        )
        legal_right = (
            verifier._linear_slice_polynomial(
                coefficients,
                boundary_cell.right_large,
                boundary_cell.right_limit,
                boundary_cell.total_offset,
                boundary_cell.large_offset,
                boundary_cell.sample,
            )
            if boundary_cell.right_legal else {}
        )
        full_left = verifier._linear_slice_polynomial(
            coefficients,
            boundary_cell.left_large,
            None,
            boundary_cell.total_offset,
            boundary_cell.large_offset,
            boundary_cell.sample,
        )
        full_right = verifier._linear_slice_polynomial(
            coefficients,
            boundary_cell.right_large,
            None,
            boundary_cell.total_offset,
            boundary_cell.large_offset,
            boundary_cell.sample,
        )
        result[tuple(signature)] = (
            legal_left,
            legal_right,
            full_left,
            full_right,
        )
    return result


def slice_regime_key(boundary_cell, *, verifier):
    """Identify cells having the same four marginal-slice operators.

    Within a fixed density status, each operator is determined by its active
    affine endpoints. Applying it to the constant polynomial records their
    affine difference, giving a small exact fingerprint for safe reuse of the
    much larger routed candidate polynomial.
    """
    unit = {(0, 0): verifier.U / verifier.U}

    def apply(large, support_limit, legal):
        if not legal:
            return ()
        polynomial = verifier._linear_slice_polynomial(
            unit,
            large,
            support_limit,
            boundary_cell.total_offset,
            boundary_cell.large_offset,
            boundary_cell.sample,
        )
        return tuple(sorted(polynomial.items()))

    return (
        boundary_cell.density_status,
        boundary_cell.left_large,
        boundary_cell.right_large,
        apply(
            boundary_cell.left_large,
            boundary_cell.left_limit,
            boundary_cell.left_legal,
        ),
        apply(
            boundary_cell.right_large,
            boundary_cell.right_limit,
            boundary_cell.right_legal,
        ),
        apply(boundary_cell.left_large, None, True),
        apply(boundary_cell.right_large, None, True),
    )


def target_correction_polynomials(
    pair_routes,
    feature_groups,
    boundary_cell,
    *,
    verifier,
    polynomial_backend,
    statistics=None,
):
    """Route each pair product once into exact target correction polynomials."""
    started = time.perf_counter()
    families = _slice_families(feature_groups, boundary_cell, verifier=verifier)
    family_seconds = time.perf_counter() - started
    encoded = {}

    def encode(polynomial):
        key = id(polynomial)
        cached = encoded.get(key)
        if cached is None or cached[0] is not polynomial:
            cached = (polynomial, polynomial_backend.from_dict(polynomial))
            encoded[key] = cached
        return cached[1]

    targets = {}
    for (left, right), outputs in pair_routes.items():
        left_legal, _left_other, left_full, _left_full_other = families[left]
        _right_other, right_legal, _right_full_other, right_full = families[right]
        if left_legal == left_full and right_legal == right_full:
            continue
        # Factor the legal-minus-unrestricted difference before multiplying.
        # Boundary cells commonly change only one distinguished coordinate;
        # those cells now need one product per signature pair instead of two.
        # When both sides change, the two-term identity remains exact:
        #   L_l R_l - L_u R_u
        #     = (L_l-L_u) R_l + L_u (R_l-R_u).
        left_changed = left_legal != left_full
        right_changed = right_legal != right_full
        if left_changed:
            left_difference = polynomial_backend.add_scaled(
                encode(left_legal), encode(left_full), -1
            )
            correction = polynomial_backend.multiply(
                left_difference, encode(right_legal)
            )
        else:
            correction = polynomial_backend.zero()
        if right_changed:
            right_difference = polynomial_backend.add_scaled(
                encode(right_legal), encode(right_full), -1
            )
            correction = polynomial_backend.add_scaled(
                correction,
                polynomial_backend.multiply(encode(left_full), right_difference),
                1,
            )
        if not correction:
            continue
        for target, structure in outputs:
            current = targets.get(target)
            if current is None:
                current = polynomial_backend.zero()
            targets[target] = polynomial_backend.add_scaled(
                current, correction, structure
            )
    answer = {target: polynomial for target, polynomial in targets.items() if polynomial}
    if statistics is not None:
        statistics.update({
            "slice_family_seconds": family_seconds,
            "pair_route_seconds": time.perf_counter() - started - family_seconds,
        })
    return answer


def contract_cell_real_ball(
    target_polynomials,
    target_densities,
    boundary_cell,
    *,
    verifier,
    polynomial_backend,
    precision: int,
):
    """Rigorous Arb contraction of all target rows on one boundary cell."""
    from sage.all import RealBallField

    field = (
        polynomial_backend.base_ring
        if getattr(polynomial_backend, "is_ball_backend", False)
        else RealBallField(int(precision))
    )
    answer = field.zero()
    raw_moments = {}
    ball_moments = {}
    converted_coefficients = {}
    density_status = boundary_cell.density_status
    active_targets = 0
    integrand_terms = 0
    for target, candidate in target_polynomials.items():
        density = target_densities.get(target, {}).get(density_status)
        if not density:
            continue
        active_targets += 1
        integrand = polynomial_backend.multiply(
            polynomial_backend.from_dict(density), candidate
        )
        for encoded_power, coefficient in integrand.dict().items():
            integrand_terms += 1
            x_power, z_power = divmod(
                int(encoded_power), polynomial_backend.stride
            )
            exponent = (x_power, z_power)
            try:
                moment = ball_moments[exponent]
            except KeyError:
                raw = fast_j.geometry_monomial_moment(
                    exponent,
                    kind=boundary_cell.kind,
                    cell=boundary_cell.cell,
                    verifier=verifier,
                    rational=polynomial_backend.rational,
                )
                raw_moments[exponent] = raw
                moment = _ball_from_rational(field, raw)
                ball_moments[exponent] = moment
            if getattr(polynomial_backend, "is_ball_backend", False):
                ball_coefficient = coefficient
            else:
                coefficient_key = (
                    int(coefficient.numerator()), int(coefficient.denominator())
                )
                ball_coefficient = converted_coefficients.get(coefficient_key)
                if ball_coefficient is None:
                    ball_coefficient = _ball_from_rational(field, coefficient)
                    converted_coefficients[coefficient_key] = ball_coefficient
            answer += ball_coefficient * moment
    return answer, {
        "active_targets": active_targets,
        "integrand_terms": integrand_terms,
        "raw_moments": len(raw_moments),
    }
