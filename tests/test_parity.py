import json
from math import isclose
from pathlib import Path
import subprocess
import sys

import pytest

from primegaps.parity import (
    degree_two_prime_indicator,
    degree_three_prime_indicator,
    liouville_second_moment_prime_indicator,
    liouville_third_moment_prime_indicator,
    parity_contributions,
    parity_error_budget,
    rough_factor_constants,
)


def test_degree_three_factorial_and_liouville_identities_detect_primes():
    expected = (0, 1, 0, 0)
    assert tuple(degree_three_prime_indicator(n) for n in range(4)) == expected
    assert tuple(liouville_third_moment_prime_indicator(n) for n in range(4)) == expected


def test_degree_two_factorial_and_liouville_identities_detect_primes():
    expected = (0, 1, 0)
    assert tuple(degree_two_prime_indicator(n) for n in range(3)) == expected
    assert tuple(liouville_second_moment_prime_indicator(n) for n in range(3)) == expected


@pytest.mark.parametrize("omega", [-1, 4, 1.5, True])
def test_degree_three_identities_reject_values_outside_their_range(omega):
    with pytest.raises(ValueError):
        degree_three_prime_indicator(omega)
    with pytest.raises(ValueError):
        liouville_third_moment_prime_indicator(omega)


@pytest.mark.parametrize("omega", [-1, 3, 1.5, True])
def test_degree_two_identities_reject_values_outside_their_range(omega):
    with pytest.raises(ValueError):
        degree_two_prime_indicator(omega)
    with pytest.raises(ValueError):
        liouville_second_moment_prime_indicator(omega)


def test_beta_quarter_rough_factor_constants_reproduce_experiment():
    constants = rough_factor_constants(0.250001)
    assert isclose(constants.semiprime, 1.0986069553418876, abs_tol=2e-14)
    assert isclose(constants.triprime, 0.14721698020054072, abs_tol=2e-14)
    assert isclose(constants.omega_choose_1, 3.638864851285397, abs_tol=5e-14)
    assert isclose(constants.omega_choose_2, 1.5402578959435096, abs_tol=5e-14)
    assert isclose(constants.signed_identity, 1.0, abs_tol=2e-14)
    assert constants.detector_degree == 3
    assert isclose(constants.rough_carrier, 2.2458239355424283, abs_tol=5e-14)


def test_degree_two_rough_factor_constants():
    constants = rough_factor_constants(0.4)
    assert constants.detector_degree == 2
    assert constants.triprime == 0.0
    assert constants.omega_choose_3 == 0.0
    assert isclose(constants.semiprime, 0.4054651081081642, abs_tol=2e-14)
    assert isclose(constants.signed_identity, 1.0, abs_tol=2e-14)
    assert isclose(constants.gross_signed_condition, 2.621860432432657, abs_tol=2e-14)


def test_rough_constants_are_continuous_across_degree_transition():
    left = rough_factor_constants(1 / 3 - 1e-6)
    right = rough_factor_constants(1 / 3 + 1e-6)
    assert left.detector_degree == 3
    assert right.detector_degree == 2
    assert left.triprime < 2.1e-11
    assert abs(left.gross_signed_condition - right.gross_signed_condition) < 4e-5


@pytest.mark.parametrize("beta", [0.25, 0.5, 0.2, 0.6])
def test_rough_constants_reject_beta_outside_degree_two_three_range(beta):
    with pytest.raises(ValueError):
        rough_factor_constants(beta)


def test_k39_signed_contributions_and_error_budget():
    constants = rough_factor_constants(0.250001)
    terms = parity_contributions(1.0192774813942034, constants)
    assert isclose(terms.plus_omega, 3.709013000752072, abs_tol=2e-13)
    assert isclose(terms.minus_2_choose_2, -3.139900377749671, abs_tol=2e-13)
    assert isclose(terms.plus_3_choose_3, 0.45016485839180237, abs_tol=2e-13)
    eta, relative = parity_error_budget(terms)
    assert isclose(eta, 0.019277481394203377, abs_tol=2e-13)
    assert isclose(relative, 0.002641084362784935, abs_tol=2e-13)


def test_recorded_physical_parity_experiment_replays():
    root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [sys.executable, root / "scripts/check_physical_parity_viability.py"],
        cwd=root,
        text=True,
        capture_output=True,
        check=True,
    )
    replay = json.loads(completed.stdout)
    assert replay["status"] == "checked-exploratory-record"
    assert replay["dimensions"]["39"]["score"] > 1
    assert replay["dimensions"]["38"]["score"] > 1


def test_physical_parity_beta_sweep_crosses_degree_transition():
    root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [
            sys.executable,
            root / "scripts/sweep_physical_parity_beta.py",
            "--betas",
            "0.250001",
            "0.3333343333333333",
            "0.4",
        ],
        cwd=root,
        text=True,
        capture_output=True,
        check=True,
    )
    sweep = json.loads(completed.stdout)
    rows = sweep["rows"]
    assert [row["detector_degree"] for row in rows] == [3, 2, 2]
    assert rows[1]["bilinear_geometry"]["triprime_prime_times_semiprime"] is None
    assert rows[2]["constants"]["gross_signed_condition"] < rows[1]["constants"][
        "gross_signed_condition"
    ] < rows[0]["constants"]["gross_signed_condition"]
    assert rows[1]["dimension_budgets"]["39"]["ideal_score"] == rows[0][
        "dimension_budgets"
    ]["39"]["ideal_score"]
