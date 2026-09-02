"""Pair-first exact J contraction for target-signature chunks.

The frozen verifier loops over target signatures and repeats each marginal
feature-pair polynomial product for every target in the monomial-symmetric
product. This module reverses those loops: one feature-pair product is formed
per support cell and distributed to all requested target signatures.
"""

from __future__ import annotations

from collections import defaultdict
import json
from typing import Iterable

from . import fast_i


def pair_routes(pair_groups, targets: Iterable[tuple[int, ...]]):
    """Invert target -> pair contributions for a target chunk."""
    routes = defaultdict(list)
    targets = tuple(targets)
    target_set = set(targets)
    for target in targets:
        for left, right, structure in pair_groups[target]:
            routes[(left, right)].append((target, structure))
    return targets, dict(routes)


def target_densities(
    targets,
    *,
    common_dimension: int,
    delta,
    max_large: int,
    max_offset_count: int,
    rational,
    orbit_size,
    positive_cache: dict | None = None,
):
    """Build normalized density polynomials for one target chunk."""
    answer = {}
    for target in targets:
        density = fast_i.orbit_status_densities(
            target,
            k=common_dimension,
            delta=delta,
            max_large=max_large,
            max_offset_count=max_offset_count,
            rational=rational,
            positive_cache=positive_cache,
        )
        orbit = orbit_size(target, common_dimension)
        grouped = defaultdict(dict)
        for (
            large,
            shifted,
            large_power,
            small_power,
        ), coefficient in fast_i.normalized_density_terms(
            density, common_dimension
        ):
            exponent = (
                0 if large_power is None else large_power,
                0 if small_power is None else small_power,
            )
            bucket = grouped[(large, shifted)]
            bucket[exponent] = bucket.get(exponent, rational(0)) + (
                orbit * coefficient
            )
        answer[target] = dict(grouped)
    return answer


def add_scaled_in_place(target, source, scalar, rational):
    """Accumulate scalar*source without temporary scale/add dictionaries."""
    for exponent, coefficient in source.items():
        value = target.get(exponent, rational(0)) + scalar * coefficient
        if value:
            target[exponent] = value
        else:
            target.pop(exponent, None)


def integrate_product_on_cell(
    density,
    candidate,
    *,
    kind,
    cell,
    verifier,
    rational,
):
    """Integrate one density/candidate polynomial product on one geometry cell."""
    if not density or not candidate:
        return rational(0)
    integrand = verifier.kernel._poly_mul(density, candidate)
    if kind == "polygons":
        polygon, _sample = cell
        return verifier._integrate_polynomial_polygon(integrand, polygon)
    if kind in ("xintervals", "zintervals"):
        start, end, _sample = cell
        variable_is_x = kind == "xintervals"
        answer = rational(0)
        for (x_power, z_power), coefficient in integrand.items():
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
    if kind == "point":
        return integrand.get((0, 0), rational(0))
    if kind == "empty":
        return rational(0)
    raise ValueError(f"unknown geometry kind {kind!r}")


def geometry_monomial_moment(
    exponent,
    *,
    kind,
    cell,
    verifier,
    rational,
):
    """Exact integral of one monomial on a single support cell."""
    x_power, z_power = exponent
    if kind == "polygons":
        polygon, _sample = cell
        return verifier._polygon_monomial_moment(x_power, z_power, polygon)
    if kind in ("xintervals", "zintervals"):
        start, end, _sample = cell
        variable_is_x = kind == "xintervals"
        if (variable_is_x and z_power) or (not variable_is_x and x_power):
            return rational(0)
        degree = x_power if variable_is_x else z_power
        return (end ** (degree + 1) - start ** (degree + 1)) / (degree + 1)
    if kind == "point":
        return rational(int(x_power == 0 and z_power == 0))
    if kind == "empty":
        return rational(0)
    raise ValueError(f"unknown geometry kind {kind!r}")


def density_weighted_moments(
    density,
    candidate_exponents,
    *,
    kind,
    cell,
    verifier,
    rational,
):
    """Build the reusable functional L(e)=integral density*x^e.

    The returned values depend on k, support, target signature and cell, but
    not on candidate coefficients.  They can therefore be persisted and
    extended when a higher degree asks for additional exponents.
    """
    answer = {}
    for exponent in candidate_exponents:
        x_power, z_power = exponent
        value = rational(0)
        for (density_x, density_z), coefficient in density.items():
            value += coefficient * geometry_monomial_moment(
                (x_power + density_x, z_power + density_z),
                kind=kind,
                cell=cell,
                verifier=verifier,
                rational=rational,
            )
        answer[tuple(exponent)] = value
    return answer


def evaluate_density_functional(candidate, moments, rational):
    """Contract a candidate polynomial against cached weighted moments."""
    return sum(
        (coefficient * moments[exponent]
         for exponent, coefficient in candidate.items()),
        rational(0),
    )


def evaluate_density_bilinear(
    left,
    right,
    density,
    *,
    kind,
    cell,
    verifier,
    rational,
    polynomial_backend=None,
):
    """Evaluate ``integral density*left*right`` without forming either product.

    The density-weighted monomial moments form a Hankel-like table indexed by
    exponent sums.  This is the exact scalar counterpart of the block-operator
    Hankel contraction.
    """
    if not left or not right or not density:
        return rational(0)
    required = tuple(sorted({
        (left_x + right_x, left_z + right_z)
        for left_x, left_z in left
        for right_x, right_z in right
    }))
    moments = density_weighted_moments(
        density,
        required,
        kind=kind,
        cell=cell,
        verifier=verifier,
        rational=rational,
    )
    if polynomial_backend is not None and hasattr(
        polynomial_backend, "contract_bilinear"
    ):
        return polynomial_backend.contract_bilinear(left, right, moments)
    return sum(
        (
            left_coefficient
            * right_coefficient
            * moments[(left_x + right_x, left_z + right_z)]
            for (left_x, left_z), left_coefficient in left.items()
            for (right_x, right_z), right_coefficient in right.items()
        ),
        rational(0),
    )


def _exact_payload(value):
    if isinstance(value, (tuple, list)):
        return [_exact_payload(item) for item in value]
    numerator = value.numerator
    denominator = value.denominator
    if callable(numerator):
        numerator = numerator()
    if callable(denominator):
        denominator = denominator()
    return [str(int(numerator)), str(int(denominator))]


def functional_id(target, status, kind, cell):
    """Return the stable cache key for one target/status/support cell."""
    return json.dumps(
        {
            "target": list(target),
            "status": list(status),
            "kind": kind,
            "cell": _exact_payload(cell),
        },
        sort_keys=True,
        separators=(",", ":"),
    )


def _candidate_dict(candidate, polynomial_backend):
    if polynomial_backend is None:
        return candidate
    return {
        (x_power, z_power): coefficient
        for x_power, z_power, coefficient in polynomial_backend.terms(candidate)
        if coefficient
    }


def evaluate_target_chunk_cached(
    targets,
    *,
    dimension: int,
    feature_groups,
    pair_groups,
    verifier,
    orbit_size,
    rational,
    functional_values,
    density_statuses,
    positive_cache: dict | None = None,
    polynomial_backend=None,
    slice_cache: dict | None = None,
    fresh_functional_sink=None,
    fresh_status_sink=None,
):
    """Evaluate a target chunk using persistent density/cell functionals.

    The returned fresh moments and density-status indexes are intended to be
    appended by the parent process.  A fully cached replay does not rebuild
    target densities or integrate any geometry moments; it only constructs the
    candidate polynomial and contracts its coefficients with cached values.
    """
    targets, routes = pair_routes(pair_groups, targets)
    common_dimension = dimension - 1
    maximum_offset = int(verifier.R // verifier.DELTA)
    known_statuses = {
        target: tuple(density_statuses[target])
        for target in targets
        if target in density_statuses
    }
    unknown = tuple(target for target in targets if target not in known_statuses)
    densities = {}
    fresh_statuses = {}
    if unknown:
        densities.update(target_densities(
            unknown,
            common_dimension=common_dimension,
            delta=verifier.DELTA,
            max_large=len(verifier.B),
            max_offset_count=maximum_offset,
            rational=rational,
            orbit_size=orbit_size,
            positive_cache=positive_cache,
        ))
        for target in unknown:
            statuses = tuple(sorted(densities[target]))
            known_statuses[target] = statuses
            if fresh_status_sink is None:
                fresh_statuses[target] = statuses
            else:
                fresh_status_sink(target, statuses)

    def density_for(target, density_status):
        if target not in densities:
            densities.update(target_densities(
                (target,),
                common_dimension=common_dimension,
                delta=verifier.DELTA,
                max_large=len(verifier.B),
                max_offset_count=maximum_offset,
                rational=rational,
                orbit_size=orbit_size,
                positive_cache=positive_cache,
            ))
        return densities[target][density_status]

    answers = {target: rational(0) for target in targets}
    fresh_functionals = {}
    for large in range(min(common_dimension, len(verifier.B)) + 1):
        for shifted in range(maximum_offset - large + 1):
            density_status = (large, shifted)
            active = tuple(
                target
                for target in targets
                if density_status in known_statuses[target]
            )
            if not active:
                continue
            total_offset = (large + shifted) * verifier.DELTA
            large_offset = large * verifier.DELTA
            for left_large in (False, True):
                if large + int(left_large) > len(verifier.B):
                    continue
                left_limit = verifier._support_limit(large, left_large)
                for right_large in (False, True):
                    if large + int(right_large) > len(verifier.B):
                        continue
                    right_limit = verifier._support_limit(large, right_large)
                    specs = (
                        verifier.RadialSlice(
                            0, 0, left_large, support_limit=left_limit
                        ),
                        verifier.RadialSlice(
                            0, 0, right_large, support_limit=right_limit
                        ),
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

                        left_polys = {}
                        right_polys = {}
                        candidate = {
                            target: (
                                {}
                                if polynomial_backend is None
                                else polynomial_backend.zero()
                            )
                            for target in active
                        }
                        active_set = set(active)
                        for (left, right), outputs in routes.items():
                            selected = tuple(
                                (target, structure)
                                for target, structure in outputs
                                if target in active_set
                            )
                            if not selected:
                                continue
                            left_poly = left_polys.get(left)
                            if left_poly is None:
                                slice_key = (
                                    left,
                                    left_large,
                                    left_limit,
                                    total_offset,
                                    large_offset,
                                    sample,
                                )
                                left_poly = (
                                    None
                                    if slice_cache is None
                                    else slice_cache.get(slice_key)
                                )
                                if left_poly is None:
                                    left_poly = verifier._linear_slice_polynomial(
                                        feature_groups[left],
                                        left_large,
                                        left_limit,
                                        total_offset,
                                        large_offset,
                                        sample,
                                    )
                                    if polynomial_backend is not None:
                                        left_poly = polynomial_backend.from_dict(
                                            left_poly
                                        )
                                    if slice_cache is not None:
                                        slice_cache[slice_key] = left_poly
                                left_polys[left] = left_poly
                            right_poly = right_polys.get(right)
                            if right_poly is None:
                                slice_key = (
                                    right,
                                    right_large,
                                    right_limit,
                                    total_offset,
                                    large_offset,
                                    sample,
                                )
                                right_poly = (
                                    None
                                    if slice_cache is None
                                    else slice_cache.get(slice_key)
                                )
                                if right_poly is None:
                                    right_poly = verifier._linear_slice_polynomial(
                                        feature_groups[right],
                                        right_large,
                                        right_limit,
                                        total_offset,
                                        large_offset,
                                        sample,
                                    )
                                    if polynomial_backend is not None:
                                        right_poly = polynomial_backend.from_dict(
                                            right_poly
                                        )
                                    if slice_cache is not None:
                                        slice_cache[slice_key] = right_poly
                                right_polys[right] = right_poly
                            product = (
                                verifier.kernel._poly_mul(left_poly, right_poly)
                                if polynomial_backend is None
                                else polynomial_backend.multiply(
                                    left_poly, right_poly
                                )
                            )
                            for target, structure in selected:
                                if polynomial_backend is None:
                                    add_scaled_in_place(
                                        candidate[target],
                                        product,
                                        structure,
                                        rational,
                                    )
                                else:
                                    candidate[target] = (
                                        polynomial_backend.add_scaled(
                                            candidate[target],
                                            product,
                                            structure,
                                        )
                                    )
                        status = (
                            large,
                            shifted,
                            int(left_large),
                            int(right_large),
                        )
                        for target in active:
                            polynomial = _candidate_dict(
                                candidate[target], polynomial_backend
                            )
                            if not polynomial:
                                continue
                            cache_id = functional_id(target, status, kind, cell)
                            moments = dict(functional_values.get(cache_id, {}))
                            missing = tuple(sorted(set(polynomial) - set(moments)))
                            if missing:
                                computed = density_weighted_moments(
                                    density_for(target, density_status),
                                    missing,
                                    kind=kind,
                                    cell=cell,
                                    verifier=verifier,
                                    rational=rational,
                                )
                                moments.update(computed)
                                if fresh_functional_sink is None:
                                    fresh_functionals[cache_id] = computed
                                else:
                                    fresh_functional_sink(cache_id, computed)
                            answers[target] += evaluate_density_functional(
                                polynomial, moments, rational
                            )
    return answers, fresh_functionals, fresh_statuses


def evaluate_target_chunk(
    targets,
    *,
    dimension: int,
    feature_groups,
    pair_groups,
    verifier,
    orbit_size,
    rational,
    positive_cache: dict | None = None,
    polynomial_backend=None,
    slice_cache: dict | None = None,
    control_variate: bool = False,
):
    """Return exact J contributions for a chunk of product signatures.

    With ``control_variate=True`` return only ``J_legal-J_unrestricted``.
    Geometry cells on which both marginal slice families agree are skipped.
    The unrestricted scalar can then be supplied by the closed simplex form.
    """
    targets, routes = pair_routes(pair_groups, targets)
    common_dimension = dimension - 1
    maximum_offset = int(verifier.R // verifier.DELTA)
    densities = target_densities(
        targets,
        common_dimension=common_dimension,
        delta=verifier.DELTA,
        max_large=(maximum_offset if control_variate else len(verifier.B)),
        max_offset_count=maximum_offset,
        rational=rational,
        orbit_size=orbit_size,
        positive_cache=positive_cache,
    )
    if polynomial_backend is not None:
        encoded_densities = {
            target: {
                context: polynomial_backend.from_dict(polynomial)
                for context, polynomial in grouped.items()
            }
            for target, grouped in densities.items()
        }
    else:
        encoded_densities = None
    answers = {
        target: (
            rational(0)
            if polynomial_backend is None
            else polynomial_backend.scalar_zero()
        )
        for target in targets
    }

    largest_common_count = min(
        common_dimension,
        maximum_offset if control_variate else len(verifier.B),
    )
    for large in range(largest_common_count + 1):
        for shifted in range(maximum_offset - large + 1):
            active = tuple(
                target
                for target in targets
                if (large, shifted) in densities[target]
            )
            if not active:
                continue
            total_offset = (large + shifted) * verifier.DELTA
            large_offset = large * verifier.DELTA
            for left_large in (False, True):
                left_legal = large + int(left_large) <= len(verifier.B)
                if not control_variate and not left_legal:
                    continue
                left_limit = (
                    verifier._support_limit(large, left_large)
                    if left_legal else None
                )
                for right_large in (False, True):
                    right_legal = large + int(right_large) <= len(verifier.B)
                    if not control_variate and not right_legal:
                        continue
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
                    )
                    geometry_specs = specs
                    if control_variate:
                        geometry_specs += (
                            verifier.RadialSlice(0, 0, left_large),
                            verifier.RadialSlice(0, 0, right_large),
                        )
                    kind, cells = verifier._slice_geometry(
                        large > 0,
                        common_dimension > large,
                        total_offset,
                        large_offset,
                        geometry_specs,
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

                        left_polys = {}
                        right_polys = {}
                        full_left_polys = {}
                        full_right_polys = {}
                        candidate = {
                            target: (
                                {}
                                if polynomial_backend is None
                                else polynomial_backend.zero()
                            )
                            for target in active
                        }
                        active_set = set(active)
                        for (left, right), outputs in routes.items():
                            selected = tuple(
                                (target, structure)
                                for target, structure in outputs
                                if target in active_set
                            )
                            if not selected:
                                continue
                            left_poly = left_polys.get(left)
                            if left_poly is None:
                                slice_key = (
                                    left,
                                    left_large,
                                    left_limit,
                                    total_offset,
                                    large_offset,
                                    sample,
                                )
                                left_poly = (
                                    None
                                    if slice_cache is None
                                    else slice_cache.get(slice_key)
                                )
                                if left_poly is None and not left_legal:
                                    left_poly = (
                                        {} if polynomial_backend is None
                                        else polynomial_backend.zero()
                                    )
                                elif left_poly is None:
                                    left_poly = verifier._linear_slice_polynomial(
                                        feature_groups[left],
                                        left_large,
                                        left_limit,
                                        total_offset,
                                        large_offset,
                                        sample,
                                    )
                                    if polynomial_backend is not None:
                                        left_poly = polynomial_backend.from_dict(
                                            left_poly
                                        )
                                    if slice_cache is not None:
                                        slice_cache[slice_key] = left_poly
                                left_polys[left] = left_poly
                            right_poly = right_polys.get(right)
                            if right_poly is None:
                                slice_key = (
                                    right,
                                    right_large,
                                    right_limit,
                                    total_offset,
                                    large_offset,
                                    sample,
                                )
                                right_poly = (
                                    None
                                    if slice_cache is None
                                    else slice_cache.get(slice_key)
                                )
                                if right_poly is None and not right_legal:
                                    right_poly = (
                                        {} if polynomial_backend is None
                                        else polynomial_backend.zero()
                                    )
                                elif right_poly is None:
                                    right_poly = verifier._linear_slice_polynomial(
                                        feature_groups[right],
                                        right_large,
                                        right_limit,
                                        total_offset,
                                        large_offset,
                                        sample,
                                    )
                                    if polynomial_backend is not None:
                                        right_poly = polynomial_backend.from_dict(
                                            right_poly
                                        )
                                    if slice_cache is not None:
                                        slice_cache[slice_key] = right_poly
                                right_polys[right] = right_poly
                            if control_variate:
                                full_left_poly = full_left_polys.get(left)
                                if full_left_poly is None:
                                    slice_key = (
                                        left,
                                        left_large,
                                        None,
                                        total_offset,
                                        large_offset,
                                        sample,
                                    )
                                    full_left_poly = (
                                        None
                                        if slice_cache is None
                                        else slice_cache.get(slice_key)
                                    )
                                    if full_left_poly is None:
                                        full_left_poly = verifier._linear_slice_polynomial(
                                            feature_groups[left],
                                            left_large,
                                            None,
                                            total_offset,
                                            large_offset,
                                            sample,
                                        )
                                        if polynomial_backend is not None:
                                            full_left_poly = polynomial_backend.from_dict(
                                                full_left_poly
                                            )
                                        if slice_cache is not None:
                                            slice_cache[slice_key] = full_left_poly
                                    full_left_polys[left] = full_left_poly
                                full_right_poly = full_right_polys.get(right)
                                if full_right_poly is None:
                                    slice_key = (
                                        right,
                                        right_large,
                                        None,
                                        total_offset,
                                        large_offset,
                                        sample,
                                    )
                                    full_right_poly = (
                                        None
                                        if slice_cache is None
                                        else slice_cache.get(slice_key)
                                    )
                                    if full_right_poly is None:
                                        full_right_poly = verifier._linear_slice_polynomial(
                                            feature_groups[right],
                                            right_large,
                                            None,
                                            total_offset,
                                            large_offset,
                                            sample,
                                        )
                                        if polynomial_backend is not None:
                                            full_right_poly = polynomial_backend.from_dict(
                                                full_right_poly
                                            )
                                        if slice_cache is not None:
                                            slice_cache[slice_key] = full_right_poly
                                    full_right_polys[right] = full_right_poly
                                if (
                                    left_poly == full_left_poly
                                    and right_poly == full_right_poly
                                ):
                                    continue
                            product = (
                                verifier.kernel._poly_mul(
                                    left_poly, right_poly
                                )
                                if polynomial_backend is None
                                else polynomial_backend.multiply(
                                    left_poly, right_poly
                                )
                            )
                            if control_variate:
                                full_product = (
                                    verifier.kernel._poly_mul(
                                        full_left_poly, full_right_poly
                                    )
                                    if polynomial_backend is None
                                    else polynomial_backend.multiply(
                                        full_left_poly, full_right_poly
                                    )
                                )
                                if polynomial_backend is None:
                                    add_scaled_in_place(
                                        product, full_product, -1, rational
                                    )
                                else:
                                    product = polynomial_backend.add_scaled(
                                        product, full_product, -1
                                    )
                            for target, structure in selected:
                                if polynomial_backend is None:
                                    add_scaled_in_place(
                                        candidate[target],
                                        product,
                                        structure,
                                        rational,
                                    )
                                else:
                                    candidate[target] = (
                                        polynomial_backend.add_scaled(
                                            candidate[target],
                                            product,
                                            structure,
                                        )
                                    )
                        for target in active:
                            if polynomial_backend is None:
                                answers[target] += integrate_product_on_cell(
                                    densities[target][(large, shifted)],
                                    candidate[target],
                                    kind=kind,
                                    cell=cell,
                                    verifier=verifier,
                                    rational=rational,
                                )
                            else:
                                contribution = (
                                    polynomial_backend.integrate_product(
                                        encoded_densities[target][
                                            (large, shifted)
                                        ],
                                        candidate[target],
                                        kind=kind,
                                        cell=cell,
                                        verifier=verifier,
                                    )
                                )
                                answers[target] = polynomial_backend.scalar_add(
                                    answers[target], contribution
                                )
    return answers


def evaluate_signature_pair_chunk_scalar(
    pairs,
    *,
    pair_route_map,
    dimension: int,
    feature_groups,
    verifier,
    orbit_size,
    rational,
    positive_cache: dict | None = None,
    target_density_cache=None,
    polynomial_backend=None,
    slice_cache: dict | None = None,
):
    """Contract a chunk of signature pairs directly to one exact J scalar.

    ``pair_route_map[(left, right)]`` contains ``(target, structure)`` rows.
    Target densities are combined *before* candidate slice multiplication and
    geometry integration.  The slice product is not formed: its coefficient
    vectors are contracted against density-weighted exponent-sum moments.

    This is the exact-certificate counterpart of :class:`JBlockOperator`: it
    contracts only the supplied candidate rather than materializing every
    entry of every candidate-independent block.
    """
    pairs = tuple((tuple(left), tuple(right)) for left, right in pairs)
    if not pairs:
        return rational(0)
    missing_pairs = set(pairs) - set(pair_route_map)
    if missing_pairs:
        raise KeyError(f"missing signature-pair routes: {sorted(missing_pairs)[:3]}")
    common_dimension = dimension - 1
    maximum_offset = int(verifier.R // verifier.DELTA)

    def density_groups(target):
        grouped = (
            None
            if target_density_cache is None
            else target_density_cache.get(target)
        )
        if grouped is None:
            grouped = target_densities(
                (target,),
                common_dimension=common_dimension,
                delta=verifier.DELTA,
                max_large=len(verifier.B),
                max_offset_count=maximum_offset,
                rational=rational,
                orbit_size=orbit_size,
                positive_cache=positive_cache,
            )[target]
            if target_density_cache is not None:
                target_density_cache[target] = grouped
        return grouped

    densities = {
        target: density_groups(target)
        for pair in pairs
        for target, _structure in pair_route_map[pair]
    }
    combined = {}
    for pair in pairs:
        grouped = {}
        for target, structure in pair_route_map[pair]:
            for status, density in densities[target].items():
                bucket = grouped.setdefault(status, {})
                add_scaled_in_place(bucket, density, structure, rational)
        combined[pair] = grouped
    answer = rational(0)
    for large in range(min(common_dimension, len(verifier.B)) + 1):
        for shifted in range(maximum_offset - large + 1):
            status = (large, shifted)
            active_pairs = tuple(
                pair for pair in pairs if status in combined[pair]
            )
            if not active_pairs:
                continue
            total_offset = (large + shifted) * verifier.DELTA
            large_offset = large * verifier.DELTA
            for left_large in (False, True):
                if large + int(left_large) > len(verifier.B):
                    continue
                left_limit = verifier._support_limit(large, left_large)
                for right_large in (False, True):
                    if large + int(right_large) > len(verifier.B):
                        continue
                    right_limit = verifier._support_limit(large, right_large)
                    specs = (
                        verifier.RadialSlice(
                            0, 0, left_large, support_limit=left_limit
                        ),
                        verifier.RadialSlice(
                            0, 0, right_large, support_limit=right_limit
                        ),
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
                        slice_polys = {}
                        for left, right in active_pairs:
                            polynomials = []
                            for signature, is_large, limit in (
                                (left, left_large, left_limit),
                                (right, right_large, right_limit),
                            ):
                                key = (
                                    signature,
                                    is_large,
                                    limit,
                                    total_offset,
                                    large_offset,
                                    sample,
                                )
                                polynomial = slice_polys.get(key)
                                if polynomial is None:
                                    polynomial = (
                                        None
                                        if slice_cache is None
                                        else slice_cache.get(key)
                                    )
                                    if polynomial is None:
                                        polynomial = verifier._linear_slice_polynomial(
                                            feature_groups[signature],
                                            is_large,
                                            limit,
                                            total_offset,
                                            large_offset,
                                            sample,
                                        )
                                        if slice_cache is not None:
                                            slice_cache[key] = polynomial
                                    slice_polys[key] = polynomial
                                polynomials.append(polynomial)
                            answer += evaluate_density_bilinear(
                                polynomials[0],
                                polynomials[1],
                                combined[(left, right)][status],
                                kind=kind,
                                cell=cell,
                                verifier=verifier,
                                rational=rational,
                                polynomial_backend=polynomial_backend,
                            )
    return answer
