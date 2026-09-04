import importlib.util
from fractions import Fraction
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_bf_dyadic_audit.py"
SPEC = importlib.util.spec_from_file_location("check_bf_dyadic_audit", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
audit = MODULE.audit
enumerate_blocks = MODULE.enumerate_blocks


def test_bf_audit_recovers_physical_exponents_and_sign_obstruction():
    result = audit()
    assert result["status"] == "refuted-proposed-bf"
    assert result["parameters"]["pair_modulus_exponent"] == "2742997/5000000"
    assert result["exact_checks"]["fi_v_upper_below_two_rough_primes"]
    assert result["exact_checks"]["physical_modulus_primes_below_rough_cutoff"]
    assert result["classification"]["potential_dispersion_or_trace"] == []


def test_declared_specialization_expands_every_v_q_pair():
    rows = enumerate_blocks(
        log2_x=64,
        d_exponent=Fraction(3, 4),
        capital_delta_exponent=Fraction(1, 8),
        log2_little_delta=4,
    )
    # Twelve V blocks times q=1 plus 36 positive dyadic q ranges.
    assert len(rows) == 12 * 37
    assert {row["parent_V_classification"] for row in rows} == {"impossible"}
    assert {row["q_slice_classification"] for row in rows} == {
        "undefined-without-physical-coefficient-norm"
    }
    assert rows[0]["V_log2_range"] == ["16", "17"]
    assert rows[-1]["q_log2_range"] == ["35", "2742997/78125"]
