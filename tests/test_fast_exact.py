from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import pytest

from primegaps.fast_exact import fast_i as fast
from primegaps.fast_exact import fast_j, moment_cache
from primegaps.fast_exact import compiled_poly, modular_exact


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

    path = tmp_path / "j-moments.jsonl"
    context = {"kind": "J", "k": 3, "support": "test"}
    cache = moment_cache.JFunctionalCache(
        path, context=context, rational=Fraction
    )
    assert cache.append("target/cell", {(0, 0): moments[(0, 0)]})
    assert cache.missing("target/cell", candidate) == ((2, 0),)
    assert cache.append("target/cell", {(2, 0): moments[(2, 0)]})
    assert cache.evaluate("target/cell", candidate) == expected
    assert not cache.append("target/cell", moments)


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
