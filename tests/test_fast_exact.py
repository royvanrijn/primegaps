from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from primegaps.fast_exact import fast_i as fast
from primegaps.fast_exact import fast_j, j_block, moment_cache
from primegaps.fast_exact import compiled_poly, modular_exact
from primegaps.symmetric import simplex_marginal_pairing


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


reference = load(
    "density_reference",
    ROOT / "reproduction" / "240" / "symmetry-assembler-design"
    / "orbit_status_densities.py",
)
verifier = load(
    "fast_j_frozen_verifier",
    ROOT / "reproduction" / "240" / "independent-reproducer"
    / "exact_symmetric_verifier.py",
)
reconstruct = load(
    "fast_exact_reconstruct_cli",
    ROOT / "scripts" / "reconstruct_modular_exact.py",
)


def test_closed_zero_dp_matches_coordinatewise_reference():
    cases = (
        ((), 3, 2, 2),
        ((2,), 4, 3, 3),
        ((4, 2), 6, 4, 5),
        ((6, 4, 2), 8, 5, 6),
        ((4, 4, 2, 2), 9, 6, 8),
    )
    for signature, k, max_large, max_offset in cases:
        expected = reference.orbit_status_densities(
            signature,
            k=k,
            delta=Fraction(7, 250),
            max_large=max_large,
            max_offset_count=max_offset,
        )
        actual = fast.orbit_status_densities(
            signature,
            k=k,
            delta=Fraction(7, 250),
            max_large=max_large,
            max_offset_count=max_offset,
            rational=Fraction,
        )
        assert actual == expected


def test_positive_core_is_reused_across_k():
    cache = {}
    for k in (8, 9, 10):
        actual = fast.orbit_status_densities(
            (6, 2, 2),
            k=k,
            delta=Fraction(7, 250),
            max_large=6,
            max_offset_count=9,
            rational=Fraction,
            positive_cache=cache,
        )
        expected = reference.orbit_status_densities(
            (6, 2, 2),
            k=k,
            delta=Fraction(7, 250),
            max_large=6,
            max_offset_count=9,
        )
        assert actual == expected
    assert len(cache) == 1


def test_validation():
    try:
        fast.orbit_status_densities(
            (2, 0), k=3, delta=Fraction(1, 4),
            max_large=2, max_offset_count=2, rational=Fraction,
        )
    except ValueError:
        pass
    else:
        raise AssertionError("zero exponent was accepted")

    try:
        fast.orbit_status_densities(
            (2, 2, 2), k=2, delta=Fraction(1, 4),
            max_large=2, max_offset_count=2, rational=Fraction,
        )
    except ValueError:
        pass
    else:
        raise AssertionError("overlong signature was accepted")


def test_pair_first_j_matches_target_first_reference():
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 3)
    actual = fast_j.evaluate_target_chunk(
        tuple(pair_groups),
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
    )
    expected = {
        signature: verifier.exact_j_signature_group(
            signature, pairs, feature_groups, 3
        )
        for signature, pairs in pair_groups.items()
    }
    assert actual == expected


def test_cached_pair_first_j_replays_without_density_reconstruction(monkeypatch):
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 3)
    targets = tuple(pair_groups)
    expected = fast_j.evaluate_target_chunk(
        targets,
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
    )
    actual, fresh, statuses = fast_j.evaluate_target_chunk_cached(
        targets,
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        functional_values={},
        density_statuses={},
    )
    assert actual == expected
    assert fresh
    assert set(statuses) == set(targets)

    def fail_density_reconstruction(*args, **kwargs):
        raise AssertionError("fully cached replay rebuilt target densities")

    monkeypatch.setattr(fast_j, "target_densities", fail_density_reconstruction)
    replay, replay_fresh, replay_statuses = fast_j.evaluate_target_chunk_cached(
        targets,
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        functional_values=fresh,
        density_statuses=statuses,
    )
    assert replay == expected
    assert not replay_fresh
    assert not replay_statuses


def test_flint_pair_first_j_matches_target_first_reference():
    pytest.importorskip("sage.all")
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 3)
    backend = compiled_poly.FlintEncodedPolynomialBackend(
        stride=64, rational=Fraction
    )
    slice_cache = {}
    actual = fast_j.evaluate_target_chunk(
        tuple(pair_groups),
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        polynomial_backend=backend,
        slice_cache=slice_cache,
    )
    expected = {
        signature: verifier.exact_j_signature_group(
            signature, pairs, feature_groups, 3
        )
        for signature, pairs in pair_groups.items()
    }
    assert actual == expected
    assert slice_cache
    assert fast_j.evaluate_target_chunk(
        tuple(pair_groups),
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        polynomial_backend=backend,
        slice_cache=slice_cache,
    ) == expected


def test_signature_pair_scalar_matches_target_group_sum():
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 3)
    targets, routes = fast_j.pair_routes(pair_groups, tuple(pair_groups))
    expected_groups = fast_j.evaluate_target_chunk(
        targets,
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
    )
    density_cache = {}
    actual = fast_j.evaluate_signature_pair_chunk_scalar(
        tuple(routes),
        pair_route_map=routes,
        dimension=3,
        feature_groups=feature_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        target_density_cache=density_cache,
        slice_cache={},
    )
    assert actual == sum(expected_groups.values(), Fraction())
    assert set(density_cache) == set(targets)


def test_flint_signature_pair_scalar_matches_target_group_sum():
    pytest.importorskip("sage.all")
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 3)
    targets, routes = fast_j.pair_routes(pair_groups, tuple(pair_groups))
    backend = compiled_poly.FlintEncodedPolynomialBackend(
        stride=64, rational=Fraction
    )
    expected = sum(
        fast_j.evaluate_target_chunk(
            targets,
            dimension=3,
            feature_groups=feature_groups,
            pair_groups=pair_groups,
            verifier=verifier,
            orbit_size=reference.monomial_symmetric_orbit_size,
            rational=Fraction,
            polynomial_backend=backend,
            slice_cache={},
        ).values(),
        Fraction(),
    )
    actual = fast_j.evaluate_signature_pair_chunk_scalar(
        tuple(routes),
        pair_route_map=routes,
        dimension=3,
        feature_groups=feature_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        polynomial_backend=backend,
        target_density_cache={},
        slice_cache={},
    )
    assert actual == expected


def test_exact_j_control_variate_matches_legal_minus_closed_unrestricted():
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    dimension = 3
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, dimension)
    legal = sum(
        fast_j.evaluate_target_chunk(
            tuple(pair_groups),
            dimension=dimension,
            feature_groups=feature_groups,
            pair_groups=pair_groups,
            verifier=verifier,
            orbit_size=reference.monomial_symmetric_orbit_size,
            rational=Fraction,
        ).values(),
        Fraction(),
    )
    correction = sum(
        fast_j.evaluate_target_chunk(
            tuple(pair_groups),
            dimension=dimension,
            feature_groups=feature_groups,
            pair_groups=pair_groups,
            verifier=verifier,
            orbit_size=reference.monomial_symmetric_orbit_size,
            rational=Fraction,
            control_variate=True,
        ).values(),
        Fraction(),
    )
    features = verifier.marginal_features(terms)
    unrestricted = Fraction()
    for left_index, left in enumerate(features):
        left_signature, left_power, left_slack, left_coefficient = left
        for right_index in range(left_index, len(features)):
            right_signature, right_power, right_slack, right_coefficient = features[right_index]
            feature_symmetry = 1 if left_index == right_index else 2
            unrestricted += (
                feature_symmetry
                * left_coefficient
                * right_coefficient
                * j_block.unrestricted_feature_pairing(
                    dimension - 1,
                    verifier.U,
                    verifier.R,
                    left_signature,
                    (left_power, left_slack),
                    right_signature,
                    (right_power, right_slack),
                )
            )
    assert legal == unrestricted + correction


def test_candidate_independent_moment_cache(tmp_path):
    signature = (2,)
    slacks = (0, 1, 3)
    moments = fast.signature_moments(
        signature,
        slacks,
        k=3,
        delta=verifier.DELTA,
        total_cap=verifier.U,
        large_caps=verifier.B,
        rational=Fraction,
        orbit_size=reference.monomial_symmetric_orbit_size,
        radial_moment=verifier._radial_group_integral,
    )
    coefficients = {
        (signature, slack): Fraction(slack + 2, slack + 1)
        for slack in slacks
    }
    expected = sum(
        coefficients[(signature, slack)] * moments[slack]
        for slack in slacks
    )
    context = {"kind": "I", "k": 3, "support": "test"}
    path = tmp_path / "moments.jsonl"
    cache = moment_cache.IMomentCache(
        path, context=context, rational=Fraction
    )
    assert cache.missing(signature, slacks) == slacks
    assert cache.append(signature, moments)
    assert cache.evaluate_atoms(coefficients) == expected

    reloaded = moment_cache.IMomentCache(
        path, context=context, rational=Fraction
    )
    assert not reloaded.missing(signature, slacks)
    assert reloaded.evaluate_atoms(coefficients) == expected
    assert not reloaded.append(signature, moments)


def test_crt_and_rational_reconstruction():
    primes = tuple(
        prime
        for _, prime in zip(
            range(3), modular_exact.descending_primes((1 << 31) - 1)
        )
    )
    assert len(set(primes)) == 3
    assert all(modular_exact.is_prime_64(prime) for prime in primes)
    numerator, denominator = -987654321, 1234567
    residues = tuple(
        modular_exact.rational_residue(numerator, denominator, prime)
        for prime in primes
    )
    residue, modulus = modular_exact.crt(residues, primes)
    assert modular_exact.rational_reconstruction(
        residue,
        modulus,
        numerator_bound=abs(numerator),
        denominator_bound=denominator,
    ) == (numerator, denominator)


def test_modular_pair_first_j_matches_exact_residue():
    pytest.importorskip("sage.all")
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 3)
    prime = 2_147_483_647
    backend = compiled_poly.FlintModularEncodedPolynomialBackend(
        prime, stride=64
    )
    actual = fast_j.evaluate_target_chunk(
        tuple(pair_groups),
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        polynomial_backend=backend,
    )
    expected = {
        signature: verifier.exact_j_signature_group(
            signature, pairs, feature_groups, 3
        )
        for signature, pairs in pair_groups.items()
    }
    assert actual == {
        signature: modular_exact.rational_residue(
            value.numerator, value.denominator, prime
        )
        for signature, value in expected.items()
    }

    second_prime = 2_147_483_629
    product_backend = compiled_poly.ProductPolynomialBackend((
        backend,
        compiled_poly.FlintModularEncodedPolynomialBackend(
            second_prime, stride=64
        ),
    ))
    batched = fast_j.evaluate_target_chunk(
        tuple(pair_groups),
        dimension=3,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        polynomial_backend=product_backend,
    )
    assert batched == {
        signature: tuple(
            modular_exact.rational_residue(
                value.numerator, value.denominator, modulus
            )
            for modulus in (prime, second_prime)
        )
        for signature, value in expected.items()
    }


def test_j_functional_cache_extends_across_degree(tmp_path):
    density = {(0, 0): Fraction(2), (1, 0): Fraction(-3, 5)}
    candidate = {(0, 0): Fraction(7, 3), (2, 0): Fraction(-4, 9)}
    cell = (Fraction(1, 7), Fraction(5, 11), (Fraction(1, 3), Fraction(0)))
    moments = fast_j.density_weighted_moments(
        density,
        candidate,
        kind="xintervals",
        cell=cell,
        verifier=verifier,
        rational=Fraction,
    )
    expected = fast_j.integrate_product_on_cell(
        density,
        candidate,
        kind="xintervals",
        cell=cell,
        verifier=verifier,
        rational=Fraction,
    )
    assert fast_j.evaluate_density_functional(
        candidate, moments, Fraction
    ) == expected

    right = {(0, 1): Fraction(5, 4), (1, 0): Fraction(-2, 7)}
    bilinear = fast_j.evaluate_density_bilinear(
        candidate,
        right,
        density,
        kind="xintervals",
        cell=cell,
        verifier=verifier,
        rational=Fraction,
    )
    explicit = fast_j.integrate_product_on_cell(
        density,
        verifier.kernel._poly_mul(candidate, right),
        kind="xintervals",
        cell=cell,
        verifier=verifier,
        rational=Fraction,
    )
    assert bilinear == explicit

    path = tmp_path / "j-moments.jsonl"
    context = {"kind": "J", "k": 3, "support": "test"}
    cache = moment_cache.JFunctionalCache(
        path, context=context, rational=Fraction
    )
    assert cache.append_density_statuses((2,), ((0, 0), (1, 2)))
    assert cache.append("target/cell", {(0, 0): moments[(0, 0)]})
    assert cache.missing("target/cell", candidate) == ((2, 0),)
    assert cache.append("target/cell", {(2, 0): moments[(2, 0)]})
    assert cache.evaluate("target/cell", candidate) == expected
    assert not cache.append("target/cell", moments)

    reloaded = moment_cache.JFunctionalCache(
        path, context=context, rational=Fraction
    )
    assert reloaded.density_statuses[(2,)] == ((0, 0), (1, 2))
    assert reloaded.evaluate("target/cell", candidate) == expected
    assert not reloaded.append_density_statuses((2,), ((1, 2), (0, 0)))

    shard = tmp_path / "worker-shard.jsonl"
    extra_id = "target/other-cell"
    extra_moments = {(1, 1): Fraction(5, 17)}
    shard.write_text(
        json.dumps(
            moment_cache.j_functional_record(
                reloaded.context_hash, extra_id, extra_moments
            ),
            sort_keys=True,
        )
        + "\n"
    )
    reloaded.ingest(shard, retain=False)
    assert extra_id not in reloaded.values
    ingested = moment_cache.JFunctionalCache(
        path, context=context, rational=Fraction
    )
    assert ingested.values[extra_id] == extra_moments


def test_batched_modular_checkpoint_loader(tmp_path):
    path = tmp_path / "batch.jsonl"
    rows = (
        {"signature": [], "primes": [101, 103], "residues": [7, 9]},
        {"signature": [2], "primes": [101, 103], "residues": [11, 13]},
    )
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    assert reconstruct.load_modular(path) == (
        (101, {(): 7, (2,): 11}),
        (103, {(): 9, (2,): 13}),
    )


def test_hankel_block_contraction_avoids_polynomial_product():
    left = (
        {(0, 0): Fraction(2), (1, 0): Fraction(-3)},
        {(0, 1): Fraction(5, 2)},
    )
    right = (
        {(0, 0): Fraction(-4), (0, 2): Fraction(7)},
        {(1, 1): Fraction(3, 5)},
    )
    moments = {
        (x_power, z_power): Fraction(2 * x_power - z_power + 7, 11)
        for x_power in range(3)
        for z_power in range(4)
    }
    actual = j_block.contract_polynomial_families(
        left, right, moments, dtype=object
    )
    expected = []
    for left_polynomial in left:
        row = []
        for right_polynomial in right:
            product = verifier.kernel._poly_mul(left_polynomial, right_polynomial)
            row.append(sum(
                coefficient * moments[exponent]
                for exponent, coefficient in product.items()
            ))
        expected.append(row)
    assert actual.tolist() == expected


def test_signature_pair_route_index_matches_reference_scan():
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 3)
    index = j_block.signature_pair_route_index(pair_groups)
    for left, right in index:
        signature_symmetry = 1 if left == right else 2
        expected = tuple(
            (tuple(target), structure // signature_symmetry)
            for target, contributions in pair_groups.items()
            for route_left, route_right, structure in contributions
            if route_left == left and route_right == right
        )
        assert index[(left, right)] == expected


def test_j_block_operator_matches_frozen_exact_j():
    dimension = 3
    basis = (((), 0), ((2,), 0), ((2, 2), 1))
    coefficients = np.asarray(
        [Fraction(2), Fraction(-3, 2), Fraction(1, 3)], dtype=object
    )
    marginal_map = j_block.MarginalMap.from_basis(basis)
    signatures = tuple(marginal_map.feature_keys)
    blocks = {}
    for left_index, left_signature in enumerate(signatures):
        left_keys = marginal_map.feature_keys[left_signature]
        for right_signature in signatures[left_index:]:
            right_keys = marginal_map.feature_keys[right_signature]
            block = np.empty((len(left_keys), len(right_keys)), dtype=object)
            routes = verifier.product_signatures(
                dimension - 1, left_signature, right_signature
            )
            for row, (left_power, left_slack) in enumerate(left_keys):
                for column, (right_power, right_slack) in enumerate(right_keys):
                    block[row, column] = sum(
                        structure * verifier.j_orbit(
                            dimension - 1,
                            target,
                            left_power,
                            left_slack,
                            right_power,
                            right_slack,
                        )
                        for target, structure in routes
                    )
            blocks[(left_signature, right_signature)] = block
    operator = j_block.JBlockOperator(marginal_map, blocks)
    terms = tuple(
        verifier.Term(signature, slack, coefficient)
        for (signature, slack), coefficient in zip(basis, coefficients)
    )
    expected = verifier.exact_j(terms, dimension)
    assert operator.quadratic(coefficients) == expected
    applied = operator.matvec(coefficients)
    assert sum(left * right for left, right in zip(coefficients, applied)) == expected


def test_streamed_feature_gram_blocks_match_dense_j_matrix():
    basis = (((), 0), ((2,), 0), ((2, 2), 1))
    marginal_map = j_block.MarginalMap.from_basis(basis)
    rng = np.random.default_rng(20260902)
    coefficients = rng.standard_normal(len(basis))
    all_weights = rng.random(11)
    all_values = {
        signature: rng.standard_normal((11, len(keys)))
        for signature, keys in marginal_map.feature_keys.items()
    }

    blocks = None
    for selection in (slice(0, 4), slice(4, 11)):
        blocks = j_block.accumulate_feature_gram_blocks(
            marginal_map,
            {
                signature: values[selection]
                for signature, values in all_values.items()
            },
            all_weights[selection],
            blocks=blocks,
        )
    operator = j_block.JBlockOperator(marginal_map, blocks)

    design = np.zeros((11, len(basis)))
    for basis_index, routes in enumerate(marginal_map.routes):
        for signature, position in routes:
            design[:, basis_index] += all_values[signature][:, position]
    dense_j = design.T @ (all_weights[:, None] * design)
    assert np.allclose(operator.matvec(coefficients), dense_j @ coefficients)
    assert np.allclose(
        operator.quadratic(coefficients),
        coefficients @ dense_j @ coefficients,
    )


def test_candidate_and_projected_batch_grams_match_signature_blocks():
    basis = (((), 0), ((2,), 0), ((2, 2), 1), ((4, 2), 2))
    marginal_map = j_block.MarginalMap.from_basis(basis)
    rng = np.random.default_rng(20260903)
    values = {
        signature: rng.standard_normal((17, len(keys)))
        for signature, keys in marginal_map.feature_keys.items()
    }
    weights = rng.random(17)
    projection = rng.standard_normal((len(basis), 3))

    candidate_values = j_block.candidate_feature_values(marginal_map, values)
    candidate_gram = j_block.accumulate_candidate_gram(candidate_values, weights)
    blocks = j_block.accumulate_feature_gram_blocks(
        marginal_map, values, weights
    )
    operator = j_block.JBlockOperator(marginal_map, blocks)
    expected = np.column_stack([
        operator.matvec(np.eye(len(basis))[:, index])
        for index in range(len(basis))
    ])
    assert np.allclose(candidate_gram, expected)

    mapped = marginal_map.forward_matrix(projection)
    projected = j_block.projected_feature_values(
        marginal_map, values, projection, mapped_projection=mapped
    )
    assert np.allclose(projected, candidate_values @ projection)
    projected_gram = j_block.accumulate_candidate_gram(projected, weights)
    assert np.allclose(projected_gram, projection.T @ expected @ projection)


def test_symmetric_cross_difference_matches_two_gram_updates_and_skips_equal_rows():
    rng = np.random.default_rng(20260904)
    unrestricted = rng.standard_normal((19, 5))
    legal = unrestricted.copy()
    legal[[2, 7, 13]] += rng.standard_normal((3, 5))
    weights = rng.random(19)
    correction, active = j_block.accumulate_gram_difference(
        legal, unrestricted, weights
    )
    expected = (
        legal.T @ (weights[:, None] * legal)
        - unrestricted.T @ (weights[:, None] * unrestricted)
    )
    assert active == 3
    assert np.allclose(correction, expected)


def test_j_block_scipy_linear_operator_matches_dense_eigenvalue():
    scipy = pytest.importorskip("scipy.sparse.linalg")
    marginal_map = j_block.MarginalMap.from_basis(
        (((), 0), ((2,), 0), ((2, 2), 1))
    )
    rng = np.random.default_rng(311)
    values = {
        signature: rng.standard_normal((12, len(keys)))
        for signature, keys in marginal_map.feature_keys.items()
    }
    weights = rng.random(12)
    operator = j_block.JBlockOperator(
        marginal_map,
        j_block.accumulate_feature_gram_blocks(marginal_map, values, weights),
    )
    dense = np.column_stack([
        operator.matvec(np.eye(operator.shape[0])[:, index])
        for index in range(operator.shape[0])
    ])
    eigenvalue = scipy.eigsh(
        operator.as_scipy_linear_operator(),
        k=1,
        which="LA",
        return_eigenvectors=False,
    )[0]
    assert np.isclose(eigenvalue, np.linalg.eigvalsh(dense)[-1])


def test_factorized_features_reconstruct_sparse_marginals():
    basis = (((3, 1), 0), ((2,), 1), ((), 2))
    marginal_map = j_block.MarginalMap.from_basis(basis)
    rng = np.random.default_rng(47)
    point_count = 7
    signature_values = {
        signature: rng.standard_normal(point_count)
        for signature in marginal_map.feature_keys
    }
    radial_values = rng.standard_normal((4, 3, point_count))
    power_scale = 1.75
    values = j_block.factorized_feature_values(
        marginal_map,
        signature_values,
        radial_values,
        power_scale=power_scale,
    )
    coefficients = rng.standard_normal(len(basis))
    reconstructed = np.zeros(point_count)
    for coefficient, routes in zip(coefficients, marginal_map.routes):
        for signature, position in routes:
            reconstructed += coefficient * values[signature][:, position]
    expected = np.zeros(point_count)
    for coefficient, (signature, slack) in zip(coefficients, basis):
        expected += coefficient * signature_values[signature] * radial_values[0, slack]
        for power in set(signature):
            erased = list(signature)
            erased.remove(power)
            erased = tuple(sorted(erased, reverse=True))
            expected += (
                coefficient
                * signature_values[erased]
                * power_scale**power
                * radial_values[power, slack]
            )
    assert np.allclose(reconstructed, expected)


def test_numerical_j_builder_features_match_dense_qmc_marginals():
    pytest.importorskip("scipy")
    builder = load(
        "test_numerical_j_block_builder",
        ROOT / "scripts" / "build_numerical_j_block_operator.py",
    )
    qmc = builder.qmc_verifier
    k = 7
    degree = 5
    basis = qmc.basis_indices(degree)
    marginal_map = j_block.MarginalMap.from_basis(basis)
    points = qmc.simplex_points(k - 1, qmc.R, 5, 73)
    coordinate_scale = k / (qmc.U * qmc.U)
    blocks = builder.evaluated_blocks(
        points,
        marginal_map,
        degree,
        k,
        coordinate_scale,
        legal=True,
    )
    reconstructed = np.zeros((len(points), len(basis)))
    for basis_index, routes in enumerate(marginal_map.routes):
        for signature, position in routes:
            reconstructed[:, basis_index] += blocks[signature][:, position]
    expected = qmc.marginal_features(
        points,
        basis,
        qmc.all_partitions(degree // 2),
        degree,
        k,
        coordinate_scale,
    )
    assert np.allclose(reconstructed, expected, rtol=2e-14, atol=2e-14)

    unrestricted = j_block.JBlockOperator(
        marginal_map,
        builder.unrestricted_blocks(marginal_map, degree, k),
    )
    assembled = np.column_stack([
        unrestricted.matvec(np.eye(len(basis))[:, index])
        for index in range(len(basis))
    ])
    assert np.allclose(
        assembled,
        qmc.unrestricted_marginal_gram(k, degree),
        rtol=2e-13,
        atol=2e-13,
    )

    total_room = np.asarray([0.24])
    lo = np.asarray([0.01])
    hi = np.asarray([0.20])
    high_degree = builder.integrated_jacobi_moments(
        total_room, lo, hi, 13, 27, 49
    )[13, 27, 0]
    nodes, weights = np.polynomial.legendre.leggauss(64)
    half = (hi[0] - lo[0]) / 2
    middle = (hi[0] + lo[0]) / 2
    coordinate = middle + half * nodes
    radial_coordinate = (total_room[0] - coordinate) / qmc.U
    reference = half * np.dot(
        weights,
        coordinate**26
        * qmc.eval_jacobi_basis(27, 49, radial_coordinate)[:, 27],
    )
    assert np.isclose(high_degree, reference, rtol=2e-13, atol=0.0)


def test_block_operator_persistence_roundtrip(tmp_path):
    marginal_map = j_block.MarginalMap.from_basis((((), 0), ((2,), 1)))
    rng = np.random.default_rng(11)
    values = {
        signature: rng.standard_normal((8, len(keys)))
        for signature, keys in marginal_map.feature_keys.items()
    }
    blocks = j_block.accumulate_feature_gram_blocks(
        marginal_map, values, rng.random(8)
    )
    original = j_block.JBlockOperator(marginal_map, blocks)
    directory = tmp_path / "operator"
    j_block.save_block_operator(directory, original, metadata={"seed": 11})
    loaded = j_block.load_block_operator(directory)
    candidate = np.asarray([0.75, -1.25])
    assert np.array_equal(loaded.matvec(candidate), original.matvec(candidate))
    assert json.loads((directory / "manifest.json").read_text())["metadata"] == {
        "seed": 11
    }


def test_unrestricted_feature_blocks_match_closed_simplex_pairing():
    dimension = 3
    basis = (((), 0), ((2,), 0), ((2, 2), 1))
    coefficients = np.asarray(
        [Fraction(2), Fraction(-3, 2), Fraction(1, 3)], dtype=object
    )
    marginal_map = j_block.MarginalMap.from_basis(basis)
    signatures = tuple(marginal_map.feature_keys)
    blocks = {
        (left, right): j_block.unrestricted_signature_pair_block(
            left,
            right,
            marginal_map.feature_keys[left],
            marginal_map.feature_keys[right],
            common_dimension=dimension - 1,
            outer_radius=verifier.U,
            marginal_radius=verifier.R,
            dtype=object,
        )
        for left_index, left in enumerate(signatures)
        for right in signatures[left_index:]
    }
    operator = j_block.JBlockOperator(marginal_map, blocks)
    expected = Fraction(0)
    for row, ((left_signature, left_slack), left_coefficient) in enumerate(
        zip(basis, coefficients)
    ):
        for column in range(row, len(basis)):
            right_signature, right_slack = basis[column]
            symmetry = 1 if row == column else 2
            expected += (
                symmetry
                * left_coefficient
                * coefficients[column]
                * simplex_marginal_pairing(
                    dimension,
                    verifier.U,
                    verifier.R,
                    left_slack,
                    left_signature,
                    right_slack,
                    right_signature,
                )
            )
    assert operator.quadratic(coefficients) == expected


def test_target_functionals_compile_to_float_j_blocks():
    dimension = 3
    terms = (
        verifier.Term((), 1, Fraction(2)),
        verifier.Term((2,), 0, Fraction(-3, 2)),
        verifier.Term((2, 2), 1, Fraction(1, 3)),
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, dimension)
    _values, functionals, statuses = fast_j.evaluate_target_chunk_cached(
        tuple(pair_groups),
        dimension=dimension,
        feature_groups=feature_groups,
        pair_groups=pair_groups,
        verifier=verifier,
        orbit_size=reference.monomial_symmetric_orbit_size,
        rational=Fraction,
        functional_values={},
        density_statuses={},
    )
    marginal_map = j_block.MarginalMap.from_basis(
        (term.signature, term.slack) for term in terms
    )
    signatures = tuple(marginal_map.feature_keys)
    blocks = {}
    for left_index, left_signature in enumerate(signatures):
        for right_signature in signatures[left_index:]:
            blocks[(left_signature, right_signature)] = (
                j_block.compile_signature_pair_block(
                    left_signature,
                    right_signature,
                    left_keys=marginal_map.feature_keys[left_signature],
                    right_keys=marginal_map.feature_keys[right_signature],
                    pair_groups=pair_groups,
                    functional_values=functionals,
                    density_statuses=statuses,
                    common_dimension=dimension - 1,
                    verifier=verifier,
                    rational=Fraction,
                )
            )
    operator = j_block.JBlockOperator(marginal_map, blocks)
    coefficients = np.asarray([float(term.coefficient) for term in terms])
    assert np.isclose(
        operator.quadratic(coefficients),
        float(verifier.exact_j(terms, dimension)),
        rtol=2e-14,
        atol=2e-14,
    )
