import json
from math import isclose
from pathlib import Path
import subprocess
import sys

import pytest

from primegaps.parity import (
    degree_three_prime_indicator,
    liouville_third_moment_prime_indicator,
    parity_contributions,
    parity_error_budget,
    rough_factor_constants,
)


def test_degree_three_factorial_and_liouville_identities_detect_primes():
    expected = (0, 1, 0, 0)
    assert tuple(degree_three_prime_indicator(n) for n in range(4)) == expected
    assert tuple(liouville_third_moment_prime_indicator(n) for n in range(4)) == expected


@pytest.mark.parametrize("omega", [-1, 4, 1.5, True])
def test_degree_three_identities_reject_values_outside_their_range(omega):
    with pytest.raises(ValueError):
        degree_three_prime_indicator(omega)
    with pytest.raises(ValueError):
        liouville_third_moment_prime_indicator(omega)


def test_beta_quarter_rough_factor_constants_reproduce_experiment():
    constants = rough_factor_constants(0.250001)
    assert isclose(constants.semiprime, 1.0986069553418876, abs_tol=2e-14)
    assert isclose(constants.triprime, 0.14721698020054072, abs_tol=2e-14)
    assert isclose(constants.omega_choose_1, 3.638864851285397, abs_tol=5e-14)
    assert isclose(constants.omega_choose_2, 1.5402578959435096, abs_tol=5e-14)
    assert isclose(constants.signed_identity, 1.0, abs_tol=2e-14)


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
