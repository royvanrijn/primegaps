import pytest

from primegaps.minorants import (
    baker_irving_base_loss,
    baker_irving_parameters,
    discard_variants,
    regime_frontier,
    stadlmann_admissible,
    stadlmann_loss_components,
    stadlmann_xi1_interval,
)


def test_published_direct_parameters_are_admissible():
    assert stadlmann_admissible(0.38, 0.4, 0.4)
    assert stadlmann_loss_components(0.4) == (0.0, 0.0)


def test_xi1_interval_and_frontier_witness():
    xi2 = 0.40481
    lower, upper = stadlmann_xi1_interval(xi2)
    assert lower == pytest.approx(0.297595)
    assert upper == pytest.approx(0.35671)
    regimes = regime_frontier(xi2, 0.01)
    assert lower < regimes.xi1 < upper
    assert regimes.xi3 == xi2
    assert stadlmann_admissible(regimes.xi1, regimes.xi2, regimes.xi3)


def test_discard_variants_mass_and_pointwise_bounds():
    variants = discard_variants(0.01, 0.02)
    assert [(v.retained_mass, v.c2) for v in variants] == [
        (0.97, 24),
        (0.99, 4),
        (0.98, 20),
        (1.0, 0),
    ]
    assert variants[1].extra_required_regimes == ("exception-B",)
    assert variants[2].extra_required_regimes == ("exception-A",)


def test_baker_irving_mapping_and_loss_identity():
    eta = 0.00481
    regimes, beta, theta = baker_irving_parameters(eta)
    assert regimes.xi2 == pytest.approx(0.40481)
    assert beta == pytest.approx(0.19038)
    assert theta == pytest.approx(0.52401475)
    base = baker_irving_base_loss(eta, order=24)
    stadlmann_a, _ = stadlmann_loss_components(0.40481, order=24)
    assert base == pytest.approx(stadlmann_a, rel=1e-12)
    assert base > 0.0
