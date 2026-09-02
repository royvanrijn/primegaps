import numpy as np

from primegaps.tradeoffs import (
    MinorantCandidate,
    TradeoffMeasurement,
    minorant_objective_matrix,
    optimize_minorant_score,
    pareto_frontier,
)


def test_complete_minorant_objective_includes_mass_and_k_penalty():
    candidate = MinorantCandidate("discard", 0.9, 2.0, "theorem-backed")
    J = np.diag([3.0, 4.0])
    K = np.diag([0.25, 0.5])
    objective = minorant_objective_matrix(J, K, k=5, candidate=candidate)
    assert np.allclose(objective, np.diag([11.0, 13.0]))


def test_minorant_score_reoptimizes_after_k_penalty():
    I = np.eye(2)
    J = np.diag([1.0, 0.8])
    K = np.diag([0.2, 0.0])
    direct = MinorantCandidate("direct", 1.0, 0.0, "theorem-backed")
    discard = MinorantCandidate("discard", 0.99, 2.0, "theorem-backed")
    direct_result = optimize_minorant_score(I, J, K, k=1, candidate=direct)
    discard_result = optimize_minorant_score(I, J, K, k=1, candidate=discard)
    assert direct_result.quotient == 1.0
    assert np.argmax(np.abs(direct_result.vector)) == 0
    assert abs(discard_result.quotient - 0.792) < 1e-12
    assert np.argmax(np.abs(discard_result.vector)) == 1


def test_pareto_frontier_keeps_real_tradeoffs_and_drops_dominated_points():
    def point(name, loss, unlocked, measure, delta):
        return TradeoffMeasurement(
            name, loss, unlocked, measure, 1.0, 1.0 + delta, delta, 0.0,
            "theorem-backed", name,
        )

    low_loss = point("low-loss", 0.0, 1, 0.0, 0.01)
    high_gain = point("high-gain", 0.1, 3, 0.2, 0.05)
    dominated = point("dominated", 0.2, 2, 0.1, 0.04)
    assert pareto_frontier((dominated, high_gain, low_loss)) == (low_loss, high_gain)
