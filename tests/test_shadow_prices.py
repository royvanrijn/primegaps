import primegaps.shadow_prices as shadows
from primegaps.distribution import ConstraintFailure, Minorant
from primegaps.shadow_prices import ScoredSupport, rank_constraint_relaxations
from primegaps.support import SupportParameters


def _support(endpoint: float) -> SupportParameters:
    return SupportParameters(
        delta=0.2,
        epsilon=0.01,
        A=(-0.01, endpoint),
        B=((0.21, 0.21, 0.21, 0.21, 0.21),),
    )


def test_shadow_price_reoptimizes_under_exactly_one_relaxation(monkeypatch):
    baseline = ScoredSupport("baseline", _support(0.2), 1.0, 0.01)
    type_ii = ScoredSupport("type-ii", _support(0.21), 1.2, 0.02)
    multiple = ScoredSupport("multiple", _support(0.22), 1.5, 0.03)

    failures = {
        0.2: (),
        0.21: (ConstraintFailure("P3.II.delta", "a", "b", "x", "necessary-inequality"),),
        0.22: (
            ConstraintFailure("P3.II.delta", "a", "b", "x", "necessary-inequality"),
            ConstraintFailure("P3.III", "a", "b", "x", "necessary-inequality"),
        ),
    }
    monkeypatch.setattr(
        shadows,
        "support_constraint_failures",
        lambda support, minorant: failures[support.A[-1]],
    )

    result = rank_constraint_relaxations(
        (baseline, type_ii, multiple),
        Minorant("0.38", "0.4", "0.4"),
        ("P3.II.delta", "P3.III"),
    )
    by_id = {item.constraint_id: item for item in result.constraints}
    priced = by_id["P3.II.delta"]
    assert priced.rank == 1
    assert priced.relaxed_candidate_id == "type-ii"
    assert abs(priced.delta_score - 0.2) < 1e-15
    assert priced.delta_standard_error_independent > 0.0
    assert by_id["P3.III"].relaxed_candidate_id == "baseline"
    assert by_id["P3.III"].delta_score == 0.0


def test_shadow_price_validates_ids_and_baseline(monkeypatch):
    candidate = ScoredSupport("candidate", _support(0.2), 1.0)
    monkeypatch.setattr(
        shadows,
        "support_constraint_failures",
        lambda support, minorant: (),
    )
    minorant = Minorant("0.38", "0.4", "0.4")
    try:
        rank_constraint_relaxations((candidate,), minorant, ("not-a-constraint",))
    except ValueError as exc:
        assert "unknown analytic constraints" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("unknown constraint was accepted")

    monkeypatch.setattr(
        shadows,
        "support_constraint_failures",
        lambda support, minorant: (
            ConstraintFailure("P3.I", "a", "b", "x", "necessary-inequality"),
        ),
    )
    try:
        rank_constraint_relaxations((candidate,), minorant, ("P3.I",))
    except ValueError as exc:
        assert "baseline" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("experiment without a baseline was accepted")
