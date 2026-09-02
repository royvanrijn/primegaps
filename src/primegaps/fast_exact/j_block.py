"""Candidate-independent block operator for repeated J evaluations.

The expensive support geometry is compiled into blocks indexed by marginal
signatures.  Repeated evaluation then consists only of a sparse marginal map,
block matrix-vector products, and the transpose sparse map.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np


Signature = tuple[int, ...]
FeatureKey = tuple[int, int]
BasisKey = tuple[Signature, int]


def canonical(values: Iterable[int]) -> Signature:
    return tuple(sorted((int(value) for value in values if value), reverse=True))


def marginal_images(signature: Signature, slack: int):
    """Yield the distinct marginal features produced by one basis term."""
    yield signature, (0, int(slack))
    for exponent in set(signature):
        erased = list(signature)
        erased.remove(exponent)
        yield canonical(erased), (int(exponent), int(slack))


@dataclass(frozen=True)
class MarginalMap:
    """Sparse map from candidate coefficients to marginal feature blocks."""

    basis: tuple[BasisKey, ...]
    feature_keys: Mapping[Signature, tuple[FeatureKey, ...]]
    routes: tuple[tuple[tuple[Signature, int], ...], ...]

    @classmethod
    def from_basis(cls, basis: Iterable[BasisKey]):
        ordered_basis = tuple((tuple(signature), int(slack)) for signature, slack in basis)
        feature_sets: dict[Signature, set[FeatureKey]] = {}
        images = []
        for signature, slack in ordered_basis:
            term_images = tuple(marginal_images(signature, slack))
            images.append(term_images)
            for marginal_signature, key in term_images:
                feature_sets.setdefault(marginal_signature, set()).add(key)
        feature_keys = {
            signature: tuple(sorted(keys))
            for signature, keys in sorted(feature_sets.items())
        }
        positions = {
            signature: {key: index for index, key in enumerate(keys)}
            for signature, keys in feature_keys.items()
        }
        routes = tuple(
            tuple(
                (signature, positions[signature][key])
                for signature, key in term_images
            )
            for term_images in images
        )
        return cls(ordered_basis, feature_keys, routes)

    def forward(self, coefficients, *, dtype=None):
        vector = np.asarray(coefficients, dtype=dtype)
        if vector.ndim != 1 or len(vector) != len(self.basis):
            raise ValueError("candidate vector has the wrong dimension")
        blocks = {
            signature: np.zeros(len(keys), dtype=vector.dtype)
            for signature, keys in self.feature_keys.items()
        }
        for coefficient, routes in zip(vector, self.routes):
            for signature, position in routes:
                blocks[signature][position] += coefficient
        return blocks

    def forward_matrix(self, coefficients, *, dtype=None):
        """Apply the sparse marginal map to several candidate vectors at once.

        The input has shape ``(candidate_dimension, projected_dimension)``.
        This is the useful orientation for projected J: each returned marginal
        block can be multiplied into its evaluated feature matrix, avoiding the
        construction of either signature-pair blocks or the full candidate
        feature matrix.
        """
        matrix = np.asarray(coefficients, dtype=dtype)
        if matrix.ndim != 2 or matrix.shape[0] != len(self.basis):
            raise ValueError("candidate matrix has the wrong dimension")
        blocks = {
            signature: np.zeros(
                (len(keys), matrix.shape[1]), dtype=matrix.dtype
            )
            for signature, keys in self.feature_keys.items()
        }
        for row, routes in zip(matrix, self.routes):
            for signature, position in routes:
                blocks[signature][position] += row
        return blocks

    def transpose(self, feature_blocks, *, dtype=None):
        if dtype is None:
            dtype = np.result_type(*(
                np.asarray(value).dtype for value in feature_blocks.values()
            ))
        result = np.zeros(len(self.basis), dtype=dtype)
        for basis_index, routes in enumerate(self.routes):
            for signature, position in routes:
                result[basis_index] += feature_blocks[signature][position]
        return result


class JBlockOperator:
    """Apply a symmetric J form without materializing its dense basis matrix.

    Blocks use the base symmetric convention: an off-diagonal signature pair
    is stored once without a factor of two.  Consequently ``quadratic(c)`` is
    exactly ``dot(c, matvec(c))``.
    """

    def __init__(self, marginal_map: MarginalMap, blocks):
        self.marginal_map = marginal_map
        normalized = {}
        for (raw_left, raw_right), raw_block in blocks.items():
            left, right = tuple(raw_left), tuple(raw_right)
            block = np.asarray(raw_block)
            if right < left:
                left, right = right, left
                block = block.T
            expected = (
                len(marginal_map.feature_keys[left]),
                len(marginal_map.feature_keys[right]),
            )
            if block.shape != expected:
                raise ValueError(
                    f"J block {(left, right)} has shape {block.shape}, expected {expected}"
                )
            if (left, right) in normalized:
                raise ValueError(f"duplicate J block {(left, right)}")
            if left == right and not np.array_equal(block, block.T):
                raise ValueError(f"diagonal J block {left} is not symmetric")
            normalized[(left, right)] = block
        self.blocks = normalized
        self.shape = (len(marginal_map.basis),) * 2
        self.dtype = np.result_type(*(block.dtype for block in normalized.values()))

    def feature_matvec(self, features):
        result = {
            signature: np.zeros(len(keys), dtype=self.dtype)
            for signature, keys in self.marginal_map.feature_keys.items()
        }
        for (left, right), block in self.blocks.items():
            if left == right:
                result[left] += block @ features[left]
            else:
                result[left] += block @ features[right]
                result[right] += block.T @ features[left]
        return result

    def matvec(self, coefficients):
        features = self.marginal_map.forward(coefficients, dtype=self.dtype)
        return self.marginal_map.transpose(
            self.feature_matvec(features), dtype=self.dtype
        )

    def quadratic(self, coefficients):
        vector = np.asarray(coefficients, dtype=self.dtype)
        return np.dot(vector, self.matvec(vector))

    def as_scipy_linear_operator(self):
        """Expose the block form directly to Davidson/Lanczos-style solvers."""
        from scipy.sparse.linalg import LinearOperator

        return LinearOperator(
            self.shape,
            matvec=self.matvec,
            dtype=self.dtype,
        )


def accumulate_feature_gram_blocks(
    marginal_map: MarginalMap,
    feature_values,
    weights,
    *,
    blocks=None,
):
    """Accumulate one numerical integration batch into signature blocks.

    ``feature_values[signature]`` has one row per integration point and one
    column per key in ``marginal_map.feature_keys[signature]``.  This is the
    stable numerical construction path: callers evaluate their preferred
    orthogonal basis directly, so no ill-conditioned monomial change of basis
    is required.  Passing the returned dictionary back as ``blocks`` streams
    arbitrarily many batches without ever constructing a dense J matrix.
    """
    weights = np.asarray(weights)
    if weights.ndim != 1:
        raise ValueError("integration weights must be one-dimensional")
    expected_signatures = set(marginal_map.feature_keys)
    if set(feature_values) != expected_signatures:
        missing = expected_signatures - set(feature_values)
        extra = set(feature_values) - expected_signatures
        raise ValueError(
            f"feature signature mismatch: missing={sorted(missing)}, "
            f"extra={sorted(extra)}"
        )
    values = {}
    for signature, keys in marginal_map.feature_keys.items():
        matrix = np.asarray(feature_values[signature])
        expected = (len(weights), len(keys))
        if matrix.shape != expected:
            raise ValueError(
                f"feature values for {signature} have shape {matrix.shape}, "
                f"expected {expected}"
            )
        values[signature] = matrix
    result = {} if blocks is None else blocks
    signatures = tuple(marginal_map.feature_keys)
    for left_index, left in enumerate(signatures):
        weighted_left_transpose = values[left].T * weights
        for right in signatures[left_index:]:
            key = (left, right)
            contribution = weighted_left_transpose @ values[right]
            if left == right:
                contribution = (contribution + contribution.T) / 2
            if key in result:
                if np.asarray(result[key]).shape != contribution.shape:
                    raise ValueError(f"existing J block {key} has the wrong shape")
                result[key] += contribution
            else:
                result[key] = contribution
    return result


def candidate_feature_values(marginal_map: MarginalMap, feature_values):
    """Assemble ``G = F M`` directly in candidate space.

    Rows are integration points and columns are candidate coefficients.  This
    discards the much larger intermediate signature-pair operator and lets one
    integration batch be accumulated with a single BLAS rank-k update.
    """
    expected_signatures = set(marginal_map.feature_keys)
    if set(feature_values) != expected_signatures:
        missing = expected_signatures - set(feature_values)
        extra = set(feature_values) - expected_signatures
        raise ValueError(
            f"feature signature mismatch: missing={sorted(missing)}, "
            f"extra={sorted(extra)}"
        )
    point_count = None
    dtype = None
    values = {}
    for signature, keys in marginal_map.feature_keys.items():
        matrix = np.asarray(feature_values[signature])
        if matrix.ndim != 2 or matrix.shape[1] != len(keys):
            raise ValueError(
                f"feature values for {signature} have shape {matrix.shape}, "
                f"expected (*, {len(keys)})"
            )
        if point_count is None:
            point_count = matrix.shape[0]
        elif matrix.shape[0] != point_count:
            raise ValueError("feature blocks have inconsistent row counts")
        dtype = matrix.dtype if dtype is None else np.result_type(dtype, matrix.dtype)
        values[signature] = matrix
    result = np.zeros((point_count or 0, len(marginal_map.basis)), dtype=dtype)
    for column, routes in enumerate(marginal_map.routes):
        for signature, position in routes:
            result[:, column] += values[signature][:, position]
    return result


def projected_feature_values(
    marginal_map: MarginalMap,
    feature_values,
    projection,
    *,
    mapped_projection=None,
):
    """Compute ``G Q`` without first materializing candidate-space ``G``.

    ``mapped_projection`` may be reused across all integration batches.  It is
    exactly ``marginal_map.forward_matrix(projection)``.
    """
    projection = np.asarray(projection)
    if projection.ndim != 2 or projection.shape[0] != len(marginal_map.basis):
        raise ValueError("projection has the wrong candidate dimension")
    mapped = (
        marginal_map.forward_matrix(projection)
        if mapped_projection is None
        else mapped_projection
    )
    if set(feature_values) != set(marginal_map.feature_keys):
        raise ValueError("feature signature mismatch")
    if set(mapped) != set(marginal_map.feature_keys):
        raise ValueError("mapped projection signature mismatch")
    point_count = None
    result_dtype = projection.dtype
    for signature, keys in marginal_map.feature_keys.items():
        values = np.asarray(feature_values[signature])
        coefficient_block = np.asarray(mapped[signature])
        if values.ndim != 2 or values.shape[1] != len(keys):
            raise ValueError(f"feature values for {signature} have the wrong shape")
        if coefficient_block.shape != (len(keys), projection.shape[1]):
            raise ValueError(f"mapped projection for {signature} has the wrong shape")
        if point_count is None:
            point_count = values.shape[0]
        elif values.shape[0] != point_count:
            raise ValueError("feature blocks have inconsistent row counts")
        result_dtype = np.result_type(
            result_dtype, values.dtype, coefficient_block.dtype
        )
    result = np.zeros(
        (point_count or 0, projection.shape[1]), dtype=result_dtype
    )
    for signature in marginal_map.feature_keys:
        result += np.asarray(feature_values[signature]) @ np.asarray(mapped[signature])
    return result


def accumulate_candidate_gram(values, weights, *, gram=None):
    """Accumulate ``values.T @ diag(weights) @ values`` with one rank-k update."""
    values = np.asarray(values)
    weights = np.asarray(weights)
    if values.ndim != 2 or weights.shape != (values.shape[0],):
        raise ValueError("values/weights have incompatible shapes")
    contribution = values.T @ (weights[:, None] * values)
    contribution = (contribution + contribution.T) / 2
    if gram is None:
        return contribution
    if np.asarray(gram).shape != contribution.shape:
        raise ValueError("existing Gram matrix has the wrong shape")
    gram += contribution
    return gram


def accumulate_gram_difference(legal, unrestricted, weights, *, gram=None):
    """Accumulate legal-minus-unrestricted using the symmetric cross form.

    Algebraically this is
    ``sym((legal + unrestricted).T W (legal - unrestricted))``.  Rows that
    agree exactly are omitted before the matrix multiplication; callers can
    therefore pass complete exact-m cells without paying for interior rows.
    """
    legal = np.asarray(legal)
    unrestricted = np.asarray(unrestricted)
    weights = np.asarray(weights)
    if legal.shape != unrestricted.shape or legal.ndim != 2:
        raise ValueError("legal and unrestricted values must have equal 2D shape")
    if weights.shape != (legal.shape[0],):
        raise ValueError("weights have the wrong shape")
    difference = legal - unrestricted
    active = np.any(difference != 0, axis=1)
    if np.any(active):
        left = legal[active] + unrestricted[active]
        right = difference[active]
        cross = left.T @ (weights[active, None] * right)
        contribution = (cross + cross.T) / 2
    else:
        contribution = np.zeros(
            (legal.shape[1], legal.shape[1]),
            dtype=np.result_type(legal.dtype, unrestricted.dtype, weights.dtype),
        )
    if gram is None:
        return contribution, int(active.sum())
    if np.asarray(gram).shape != contribution.shape:
        raise ValueError("existing Gram matrix has the wrong shape")
    gram += contribution
    return gram, int(active.sum())


def factorized_feature_values(
    marginal_map: MarginalMap,
    signature_values,
    radial_values,
    *,
    power_scale=1.0,
):
    """Assemble stable marginal feature columns from evaluated factors.

    ``signature_values[sigma]`` evaluates the common-variable symmetric
    factor, while ``radial_values[r, b]`` evaluates the already-integrated
    last-coordinate factor.  This matches the sparse marginal identity
    ``m_lambda -> m_lambda t^0 + sum_r m_(lambda-r) t^r`` and is suitable for
    direct Jacobi/Gauss evaluation.
    """
    missing = set(marginal_map.feature_keys) - set(signature_values)
    if missing:
        raise ValueError(f"missing signature values: {sorted(missing)}")
    radial_values = np.asarray(radial_values)
    if radial_values.ndim != 3:
        raise ValueError("radial values must have shape (power, slack, point)")
    point_count = radial_values.shape[2]
    result = {}
    for signature, keys in marginal_map.feature_keys.items():
        angular = np.asarray(signature_values[signature])
        if angular.shape != (point_count,):
            raise ValueError(
                f"signature values for {signature} have shape {angular.shape}, "
                f"expected {(point_count,)}"
            )
        columns = []
        for power, slack in keys:
            if power >= radial_values.shape[0] or slack >= radial_values.shape[1]:
                raise ValueError(f"radial values do not contain feature {(power, slack)}")
            columns.append(
                angular * (power_scale**power) * radial_values[power, slack]
            )
        result[signature] = np.column_stack(columns)
    return result


def polynomial_family_matrix(polynomials, *, dtype=np.float64):
    """Return a coefficient matrix and shared monomial column ordering."""
    polynomials = tuple(polynomials)
    exponents = tuple(sorted(set().union(*(set(poly) for poly in polynomials))))
    positions = {exponent: index for index, exponent in enumerate(exponents)}
    matrix = np.zeros((len(polynomials), len(exponents)), dtype=dtype)
    for row, polynomial in enumerate(polynomials):
        for exponent, coefficient in polynomial.items():
            matrix[row, positions[tuple(exponent)]] = coefficient
    return matrix, exponents


def hankel_moment_matrix(left_exponents, right_exponents, moments, *, dtype=np.float64):
    """Build H[i,j] = M[left_exponents[i] + right_exponents[j]]."""
    resolved_dtype = np.dtype(dtype)
    if not resolved_dtype.hasobject:
        left = np.asarray(left_exponents, dtype=np.intp)
        right = np.asarray(right_exponents, dtype=np.intp)
        maximum_x = int(left[:, 0].max() + right[:, 0].max())
        maximum_z = int(left[:, 1].max() + right[:, 1].max())
        table = np.empty((maximum_x + 1, maximum_z + 1), dtype=resolved_dtype)
        present = np.zeros(table.shape, dtype=bool)
        for (x_power, z_power), value in moments.items():
            if x_power <= maximum_x and z_power <= maximum_z:
                table[x_power, z_power] = value
                present[x_power, z_power] = True
        x_indices = left[:, 0, None] + right[None, :, 0]
        z_indices = left[:, 1, None] + right[None, :, 1]
        if not np.all(present[x_indices, z_indices]):
            missing_position = np.argwhere(~present[x_indices, z_indices])[0]
            exponent = (
                int(x_indices[tuple(missing_position)]),
                int(z_indices[tuple(missing_position)]),
            )
            raise KeyError(f"missing combined J moment {exponent}")
        return table[x_indices, z_indices]
    matrix = np.empty((len(left_exponents), len(right_exponents)), dtype=dtype)
    for row, (left_x, left_z) in enumerate(left_exponents):
        for column, (right_x, right_z) in enumerate(right_exponents):
            exponent = (left_x + right_x, left_z + right_z)
            try:
                matrix[row, column] = moments[exponent]
            except KeyError as error:
                raise KeyError(f"missing combined J moment {exponent}") from error
    return matrix


def contract_polynomial_families(
    left_polynomials,
    right_polynomials,
    moments,
    *,
    dtype=np.float64,
):
    """Evaluate every L(P_i Q_j) as one dense Hankel block contraction."""
    left, left_exponents = polynomial_family_matrix(left_polynomials, dtype=dtype)
    right, right_exponents = polynomial_family_matrix(right_polynomials, dtype=dtype)
    hankel = hankel_moment_matrix(
        left_exponents, right_exponents, moments, dtype=dtype
    )
    return left @ hankel @ right.T


def combine_target_moments(
    routes,
    functional_ids,
    functional_values,
    required_exponents,
    *,
    rational,
):
    """Compile sum_target route_weight * L_target on required monomials."""
    answer = {tuple(exponent): rational(0) for exponent in required_exponents}
    for target, structure in routes:
        functional_id = functional_ids[target]
        try:
            values = functional_values[functional_id]
        except KeyError as error:
            raise KeyError(f"missing target J functional {functional_id}") from error
        for exponent in answer:
            try:
                answer[exponent] += structure * values[exponent]
            except KeyError as error:
                raise KeyError(
                    f"missing target J moment {functional_id}, exponent={exponent}"
                ) from error
    return answer


def signature_pair_routes(pair_groups, left: Signature, right: Signature):
    """Extract base (non-signature-symmetrized) target routes for one block."""
    return signature_pair_route_index(pair_groups).get(
        tuple(sorted((tuple(left), tuple(right)))), ()
    )


def signature_pair_route_index(pair_groups):
    """Invert target routing once for all candidate-independent blocks."""
    index = defaultdict(list)
    for target, contributions in pair_groups.items():
        for raw_left, raw_right, structure in contributions:
            left, right = tuple(raw_left), tuple(raw_right)
            if right < left:
                left, right = right, left
            signature_symmetry = 1 if left == right else 2
            if structure % signature_symmetry:
                raise ArithmeticError("non-integral signature symmetry removal")
            index[(left, right)].append(
                (tuple(target), structure // signature_symmetry)
            )
    return {pair: tuple(routes) for pair, routes in index.items()}


def unrestricted_feature_pairing(
    common_dimension,
    outer_radius,
    marginal_radius,
    left_signature,
    left_key,
    right_signature,
    right_key,
):
    """Closed simplex value for one pair of already-erased marginal features."""
    from primegaps.symmetric import factorial_moment, radial_integral_factor
    from fractions import Fraction

    left_power, left_slack = left_key
    right_power, right_slack = right_key
    left_factor = Fraction(
        factorial(left_power) * factorial(left_slack),
        factorial(left_power + left_slack + 1),
    )
    right_factor = Fraction(
        factorial(right_power) * factorial(right_slack),
        factorial(right_power + right_slack + 1),
    )
    radial_degree = (
        left_power + left_slack + 1
        + right_power + right_slack + 1
    )
    q = common_dimension + sum(left_signature) + sum(right_signature)
    return (
        left_factor
        * right_factor
        * radial_integral_factor(
            q,
            marginal_radius,
            radial_degree,
            Fraction(outer_radius) - Fraction(marginal_radius),
        )
        * factorial_moment(
            common_dimension, left_signature, right_signature
        )
    )


def unrestricted_signature_pair_block(
    left_signature,
    right_signature,
    left_keys,
    right_keys,
    *,
    common_dimension,
    outer_radius,
    marginal_radius,
    dtype=np.float64,
):
    block = np.empty((len(left_keys), len(right_keys)), dtype=dtype)
    for row, left_key in enumerate(left_keys):
        for column, right_key in enumerate(right_keys):
            block[row, column] = unrestricted_feature_pairing(
                common_dimension,
                outer_radius,
                marginal_radius,
                left_signature,
                left_key,
                right_signature,
                right_key,
            )
    return block


def jacobi_monomial_coefficients(
    marginal_signature,
    feature_key,
    *,
    common_dimension,
    outer_radius,
    verifier,
    rational,
):
    """Expand one stable radial feature into exact slack-monomial coefficients."""
    power, radial_degree = feature_key
    beta = common_dimension + 2 * (sum(marginal_signature) + power)
    return {
        (power, slack): rational(coefficient) / outer_radius**slack
        for slack, coefficient in enumerate(
            verifier.jacobi_q_coefficients(radial_degree, beta)
        )
        if coefficient
    }


def unrestricted_jacobi_signature_pair_block(
    left_signature,
    right_signature,
    left_keys,
    right_keys,
    *,
    common_dimension,
    outer_radius,
    marginal_radius,
    verifier,
    rational,
    dtype=np.float64,
):
    """Closed unrestricted block in the stable source-vector radial basis."""
    block = np.empty((len(left_keys), len(right_keys)), dtype=dtype)
    left_expansions = tuple(
        jacobi_monomial_coefficients(
            left_signature,
            key,
            common_dimension=common_dimension,
            outer_radius=outer_radius,
            verifier=verifier,
            rational=rational,
        )
        for key in left_keys
    )
    right_expansions = tuple(
        jacobi_monomial_coefficients(
            right_signature,
            key,
            common_dimension=common_dimension,
            outer_radius=outer_radius,
            verifier=verifier,
            rational=rational,
        )
        for key in right_keys
    )
    for row, left_expansion in enumerate(left_expansions):
        for column, right_expansion in enumerate(right_expansions):
            value = rational(0)
            for left_key, left_coefficient in left_expansion.items():
                for right_key, right_coefficient in right_expansion.items():
                    value += (
                        left_coefficient
                        * right_coefficient
                        * unrestricted_feature_pairing(
                            common_dimension,
                            outer_radius,
                            marginal_radius,
                            left_signature,
                            left_key,
                            right_signature,
                            right_key,
                        )
                    )
            block[row, column] = value
    return block


def compile_signature_pair_block(
    left_signature,
    right_signature,
    *,
    left_keys,
    right_keys,
    pair_groups,
    functional_values,
    density_statuses,
    common_dimension,
    verifier,
    rational,
    dtype=np.float64,
    control_variate=False,
    radial_basis="monomial",
    route_index=None,
):
    """Compile one candidate-independent K block from target moment functionals.

    This is the offline target-indexed stage.  Slice products are never formed:
    every cell is contracted as ``A @ H @ B.T`` from combined exact moments.
    """
    from . import fast_j

    left_signature, right_signature = tuple(left_signature), tuple(right_signature)
    if route_index is None:
        routes = signature_pair_routes(
            pair_groups, left_signature, right_signature
        )
    else:
        routes = route_index.get(
            tuple(sorted((left_signature, right_signature))), ()
        )
    if not routes:
        return np.zeros((len(left_keys), len(right_keys)), dtype=dtype)
    if radial_basis not in ("monomial", "jacobi"):
        raise ValueError("radial_basis must be 'monomial' or 'jacobi'")
    if control_variate:
        unrestricted_builder = (
            unrestricted_signature_pair_block
            if radial_basis == "monomial"
            else unrestricted_jacobi_signature_pair_block
        )
        unrestricted_keywords = {
            "common_dimension": common_dimension,
            "outer_radius": verifier.U,
            "marginal_radius": verifier.R,
            "dtype": dtype,
        }
        if radial_basis == "jacobi":
            unrestricted_keywords.update({"verifier": verifier, "rational": rational})
        block = unrestricted_builder(
            left_signature,
            right_signature,
            left_keys,
            right_keys,
            **unrestricted_keywords,
        )
    else:
        block = np.zeros((len(left_keys), len(right_keys)), dtype=dtype)
    common_dimension = int(common_dimension)
    if common_dimension < max(len(left_signature), len(right_signature)):
        raise ValueError("common dimension is smaller than a marginal signature")
    maximum_offset = int(verifier.R // verifier.DELTA)
    largest_common_count = min(
        common_dimension,
        maximum_offset if control_variate else len(verifier.B),
    )
    for large in range(largest_common_count + 1):
        for shifted in range(maximum_offset - large + 1):
            density_status = (large, shifted)
            active_routes = tuple(
                (target, structure)
                for target, structure in routes
                if density_status in density_statuses.get(target, ())
            )
            if not active_routes:
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
                    geometry_specs = tuple(
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
                    status = (
                        large,
                        shifted,
                        int(left_large),
                        int(right_large),
                    )
                    for cell in iterable:
                        if kind == "polygons":
                            _polygon, sample = cell
                        elif kind in ("xintervals", "zintervals"):
                            _start, _end, sample = cell
                        elif kind == "point":
                            sample = cell
                        else:
                            continue
                        def feature_polynomial(signature, key, is_large, limit):
                            if radial_basis == "jacobi":
                                coefficients = jacobi_monomial_coefficients(
                                    signature,
                                    key,
                                    common_dimension=common_dimension,
                                    outer_radius=verifier.U,
                                    verifier=verifier,
                                    rational=rational,
                                )
                                return verifier._linear_slice_polynomial(
                                    coefficients,
                                    is_large,
                                    limit,
                                    total_offset,
                                    large_offset,
                                    sample,
                                )
                            power, slack = key
                            return verifier._slice_polynomial(
                                verifier.RadialSlice(
                                    power,
                                    slack,
                                    is_large,
                                    support_limit=limit,
                                ),
                                total_offset,
                                large_offset,
                                sample,
                            )

                        left_polynomials = (
                            tuple(
                                feature_polynomial(
                                    left_signature, key, left_large, left_limit
                                )
                                for key in left_keys
                            )
                            if left_legal else tuple({} for _key in left_keys)
                        )
                        right_polynomials = (
                            tuple(
                                feature_polynomial(
                                    right_signature, key, right_large, right_limit
                                )
                                for key in right_keys
                            )
                            if right_legal else tuple({} for _key in right_keys)
                        )
                        if control_variate:
                            full_left_polynomials = tuple(
                                feature_polynomial(
                                    left_signature, key, left_large, None
                                )
                                for key in left_keys
                            )
                            full_right_polynomials = tuple(
                                feature_polynomial(
                                    right_signature, key, right_large, None
                                )
                                for key in right_keys
                            )
                            if (
                                left_polynomials == full_left_polynomials
                                and right_polynomials == full_right_polynomials
                            ):
                                continue
                            left_family = left_polynomials + full_left_polynomials
                            right_family = right_polynomials + full_right_polynomials
                        else:
                            full_left_polynomials = full_right_polynomials = None
                            left_family = left_polynomials
                            right_family = right_polynomials
                        left_support = set().union(*(
                            set(polynomial) for polynomial in left_family
                        ))
                        right_support = set().union(*(
                            set(polynomial) for polynomial in right_family
                        ))
                        if not left_support or not right_support:
                            continue
                        required = tuple(sorted({
                            (left_x + right_x, left_z + right_z)
                            for left_x, left_z in left_support
                            for right_x, right_z in right_support
                        }))
                        functional_ids = {
                            target: fast_j.functional_id(target, status, kind, cell)
                            for target, _structure in active_routes
                        }
                        moments = combine_target_moments(
                            active_routes,
                            functional_ids,
                            functional_values,
                            required,
                            rational=rational,
                        )
                        legal_contribution = contract_polynomial_families(
                            left_polynomials, right_polynomials, moments, dtype=dtype
                        )
                        if control_variate:
                            full_contribution = contract_polynomial_families(
                                full_left_polynomials,
                                full_right_polynomials,
                                moments,
                                dtype=dtype,
                            )
                            block -= full_contribution - legal_contribution
                        else:
                            block += legal_contribution
    if left_signature == right_signature:
        block = (block + block.T) / 2
    return block


def file_sha256(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def save_block_operator(directory, operator: JBlockOperator, *, metadata=None):
    """Persist a block operator as one hashed, memory-mappable array."""
    directory = Path(directory)
    if directory.exists() and any(directory.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty directory: {directory}")
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    ordered_blocks = sorted(operator.blocks.items())
    offset = 0
    for (left, right), block in ordered_blocks:
        count = int(block.size)
        records.append({
            "left_signature": list(left),
            "right_signature": list(right),
            "offset": offset,
            "count": count,
            "shape": list(block.shape),
        })
        offset += count
    packed = np.concatenate([
        np.asarray(block).reshape(-1) for _pair, block in ordered_blocks
    ])
    packed_path = directory / "blocks.npy"
    temporary = packed_path.with_suffix(".npy.tmp")
    with temporary.open("wb") as stream:
        np.save(stream, packed, allow_pickle=False)
    temporary.replace(packed_path)
    manifest = {
        "schema": "primegaps-J-block-operator-v2",
        "basis": [
            {"signature": list(signature), "slack": slack}
            for signature, slack in operator.marginal_map.basis
        ],
        "dtype": str(operator.dtype),
        "expected_block_count": len(records),
        "metadata": {} if metadata is None else metadata,
        "packed_blocks": "blocks.npy",
        "packed_blocks_sha256": file_sha256(packed_path),
    }
    manifest_text = json.dumps(manifest, sort_keys=True) + "\n"
    index_text = "".join(json.dumps(record, sort_keys=True) + "\n" for record in records)
    for path, contents in (
        (directory / "manifest.json", manifest_text),
        (directory / "index.jsonl", index_text),
    ):
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(contents)
        temporary.replace(path)


def load_block_operator(directory, *, require_complete=True):
    """Load a checkpointed float block operator and verify every block hash."""
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    schema = manifest.get("schema")
    if schema not in {
        "primegaps-J-block-operator-v1",
        "primegaps-J-block-operator-v2",
    }:
        raise ValueError("unsupported J block operator schema")
    marginal_map = MarginalMap.from_basis(
        (tuple(item["signature"]), int(item["slack"]))
        for item in manifest["basis"]
    )
    records = {}
    index_path = directory / "index.jsonl"
    packed = None
    if schema == "primegaps-J-block-operator-v2":
        packed_path = directory / manifest["packed_blocks"]
        if file_sha256(packed_path) != manifest["packed_blocks_sha256"]:
            raise ValueError(f"J packed-block hash mismatch: {packed_path}")
        packed = np.load(packed_path, mmap_mode="r", allow_pickle=False)
    if index_path.exists():
        for line_number, line in enumerate(index_path.read_text().splitlines(), 1):
            record = json.loads(line)
            key = (tuple(record["left_signature"]), tuple(record["right_signature"]))
            if key in records:
                raise ValueError(f"duplicate J block index row {line_number}")
            if packed is None:
                path = directory / record["path"]
                if file_sha256(path) != record["sha256"]:
                    raise ValueError(f"J block hash mismatch: {path}")
                records[key] = np.load(path, allow_pickle=False)
            else:
                offset = int(record["offset"])
                count = int(record["count"])
                shape = tuple(int(value) for value in record["shape"])
                if count != int(np.prod(shape)) or offset < 0:
                    raise ValueError(f"invalid packed J block index row {line_number}")
                records[key] = packed[offset : offset + count].reshape(shape)
    expected = int(manifest["expected_block_count"])
    if require_complete and len(records) != expected:
        raise ValueError(f"incomplete J block operator: {len(records)}/{expected}")
    return JBlockOperator(marginal_map, records)
