from fractions import Fraction

import primegaps.breakthrough as breakthrough
from primegaps.distribution import AnalyticSlack, Minorant, support_constraint_slacks
from primegaps.shadow_prices import ScoredSupport
from primegaps.support import SupportParameters


def _support(endpoint: float, b2: float = 0.21) -> SupportParameters:
    return SupportParameters(
        delta=0.2,
        epsilon=0.01,
        A=(-0.01, endpoint),
        B=((0.21, b2, 0.22, 0.22, 0.22),),
    )


def test_global_slacks_are_exact_and_simultaneous():
    support = SupportParameters(
        delta=0.028,
        epsilon=0.0075,
        A=(-0.0075, 0.254),
        B=((0.15, 0.15) + (0.17,) * 33,),
    )
    slacks = {
        item.constraint_id: item.slack
        for item in support_constraint_slacks(
            support, Minorant("0.38", "0.4", "0.4")
        )
    }
    assert slacks["P3.II.range"] == Fraction(8, 1000) - Fraction(1, 10**8)
    assert slacks["P3.II.delta"] == Fraction(25, 10000) + Fraction(1, 5 * 10**9)
    assert slacks["P3.I"] == 0
    assert slacks["P3.III"] == 0


def test_minimum_breakthrough_allows_joint_slack_and_score_gate(monkeypatch):
    candidates = (
        ScoredSupport("cheap", _support(0.2), 1.01, 0.002),
        ScoredSupport("expensive", _support(0.21), 1.03, 0.001),
        ScoredSupport("below", _support(0.22), 0.999, 0.0),
    )
    profiles = {
        0.2: (
            AnalyticSlack(
                "P3.II.delta", Fraction(1, 1000), "necessary-inequality"
            ),
            AnalyticSlack("P3.III", Fraction(2, 1000), "necessary-inequality"),
        ),
        0.21: (
            AnalyticSlack(
                "P3.II.delta", Fraction(6, 1000), "necessary-inequality"
            ),
        ),
        0.22: (),
    }
    monkeypatch.setattr(
        breakthrough,
        "support_constraint_slacks",
        lambda support, minorant: profiles[support.A[-1]],
    )
    result = breakthrough.minimum_breakthrough(
        candidates,
        Minorant("0.38", "0.4", "0.4"),
        {"P3.III": 2},
        score_standard_error_multiplier=2,
    )
    assert result.optimum_candidate_id == "cheap"
    assert result.optimum_weighted_cost == Fraction(5, 1000)
    by_id = {item.candidate_id: item for item in result.candidates}
    assert by_id["below"].reaches_target is False


def test_minimum_breakthrough_validates_weights():
    candidate = ScoredSupport("one", _support(0.2), 1.0)
    minorant = Minorant("0.38", "0.4", "0.4")
    try:
        breakthrough.minimum_breakthrough(
            (candidate,), minorant=minorant, weights={"unknown": 1}
        )
    except ValueError as exc:
        assert "unknown analytic constraints" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("unknown weight was accepted")
