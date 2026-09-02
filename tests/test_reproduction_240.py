from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "reproduction" / "240"
VERIFIER = BUNDLE / "independent-reproducer" / "exact_symmetric_verifier.py"
DP = BUNDLE / "symmetry-assembler-design" / "orbit_status_densities.py"
SELFTEST = BUNDLE / "independent-reproducer" / "selftest_exact_bridge.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def test_frozen_evaluator_and_candidates_match_recorded_hashes():
    expected = {
        VERIFIER: "3739722c2e649daa8910e32071ff7dcad719f5dcc12106269a0a100cfabc5d8f",
        DP: "d025696eecb159cd8accb7bbeb133d345f08f073b1a0f7d2fad9d0373ac07740",
        SELFTEST: "5b641f530e8b0e93a6d8b5a04d5f3c5d8bffd3a68873bcc4416fbd1769703127",
        ROOT / "src" / "primegaps" / "integrals.py":
            "9d5e527f7479ec450a94cd1a82ba21ab1e04cc0cf10ec98d46926a9abb4fb531",
        BUNDLE / "independent-reproducer" / "exact_provenance.py":
            "e7cb7d8640225721b2966797cb3be5ce08118c487a638734bde5f23a6c84777d",
        BUNDLE / "independent-reproducer" / "run_parallel_exact.py":
            "fd55332bd6a74fd7b556ba1c527445a1c8707b06351b6431ab5fd8cfc33c1bb0",
        BUNDLE / "independent-reproducer" / "finalize_exact.py":
            "ef327e729f074ae7b8e5078e4423cc2436d119d38a224d9dd223efe765997b9f",
        BUNDLE / "independent-reproducer" / "candidate-k49-d21.json":
            "c840f99232b6c821b1f63aa81e496d1e850a4f5b482e5822fbf537c06be90815",
        BUNDLE / "independent-reproducer" / "candidate-k48-d21.json":
            "451301f90d6f5ded94f352a44cea935326fd9bb6dbac698fc45e829476f5479d",
        BUNDLE / "independent-reproducer" / "candidate-k49-binding.json":
            "cdd593ac275d90eb15e7ebb3f6472073a89b031faaa7a430891ca3d1911b2a61",
        BUNDLE / "independent-reproducer" / "candidate-k48-binding.json":
            "93b508adba400692b94cb47591e7719678de96dc2ec576496c4e8ac8cf234b42",
        BUNDLE / "exact-k49-d21-result.json":
            "2d6ff1239167c21ac79055deb7143d2bf52a24daafc82204d6529362ca0a81b5",
        BUNDLE / "exact-k48-d21-result.json":
            "fe186869bcbc51490bb7e2dbd499c80d7c5fa310fa2143fae0bb0e53548f02dc",
    }
    assert {path: _digest(path) for path in expected} == expected

    verifier = _load("frozen_exact_symmetric_verifier", VERIFIER)
    for k in (49, 48):
        candidate = BUNDLE / "independent-reproducer" / f"candidate-k{k}-d21.json"
        binding_path = BUNDLE / "independent-reproducer" / f"candidate-k{k}-binding.json"
        binding = json.loads(binding_path.read_text())
        assert binding["candidate_sha256"] == _digest(candidate)
        assert binding["exact_termwise_match"] is True
        assert len(verifier.rational_terms_from_candidate(candidate, k)) == 846


def test_recorded_k49_result_is_a_strict_exact_certificate():
    result_path = BUNDLE / "exact-k49-d21-result.json"
    assert _digest(result_path) == "2d6ff1239167c21ac79055deb7143d2bf52a24daafc82204d6529362ca0a81b5"
    result = json.loads(result_path.read_text())
    quotient = result["quotient_kJ_over_I"]
    difference = result["kJ_minus_I"]
    i_value = Fraction(int(result["I"]["numerator"]), int(result["I"]["denominator"]))
    j_value = Fraction(int(result["J"]["numerator"]), int(result["J"]["denominator"]))
    q_value = Fraction(int(quotient["numerator"]), int(quotient["denominator"]))
    diff_value = Fraction(int(difference["numerator"]), int(difference["denominator"]))
    assert q_value == 49 * j_value / i_value
    assert diff_value == 49 * j_value - i_value
    assert i_value > 0 and q_value > 1 and diff_value > 0


def test_recorded_k48_result_has_the_exact_reported_deficit():
    result_path = BUNDLE / "exact-k48-d21-result.json"
    assert _digest(result_path) == "fe186869bcbc51490bb7e2dbd499c80d7c5fa310fa2143fae0bb0e53548f02dc"
    result = json.loads(result_path.read_text())
    i_value = Fraction(int(result["I"]["numerator"]), int(result["I"]["denominator"]))
    j_value = Fraction(int(result["J"]["numerator"]), int(result["J"]["denominator"]))
    quotient = result["quotient_kJ_over_I"]
    deficit = result["deficit_one_minus_quotient"]
    q_value = Fraction(int(quotient["numerator"]), int(quotient["denominator"]))
    deficit_value = Fraction(int(deficit["numerator"]), int(deficit["denominator"]))
    assert q_value == 48 * j_value / i_value
    assert deficit_value == 1 - q_value
    assert i_value > 0 and q_value < 1 and deficit_value > 0
    assert result["certified_strictly_above_one"] is False


def test_status_density_dp_self_check():
    dp = _load("frozen_orbit_status_densities", DP)
    report = dp.self_check()
    assert report["explicit_density_cases"] == 4
    assert len(report["reference_integral_cases"]) == 5


def test_exact_i_j_bridge_against_reference_kernel():
    verifier = _load("exact_symmetric_verifier", VERIFIER)
    bridge = _load("frozen_selftest_exact_bridge", SELFTEST)
    terms = (verifier.Term((), 0, Fraction(2)),)
    polynomial = bridge.expand(terms, 2)
    parameters = bridge.support()
    assert verifier.exact_i_grouped(terms, 2) == bridge.exact_i_entry(
        polynomial, polynomial, parameters
    )
    ordinary_j = verifier.exact_j(terms, 2)
    assert ordinary_j == bridge.exact_j_entry(
        polynomial, polynomial, parameters, 2
    )
    feature_groups = verifier.grouped_marginal_coefficients(terms)
    pair_groups = verifier.grouped_signature_pairs(feature_groups, 2)
    grouped_j = sum(
        verifier.exact_j_signature_group(signature, pairs, feature_groups, 2)
        for signature, pairs in pair_groups.items()
    )
    assert grouped_j == ordinary_j
