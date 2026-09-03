import numpy as np
import pytest

from primegaps.sos import (
    best_rank_one,
    factor_psd,
    forbidden_component_pairs,
    maximal_cliques,
    solve_sparse_psd,
)


CYCLE4 = np.array(
    [
        [1, 1, 0, 1],
        [1, 1, 1, 0],
        [0, 1, 1, 1],
        [1, 0, 1, 1],
    ],
    dtype=bool,
)


def test_cycle_maximal_cliques_and_rank_one_support():
    assert maximal_cliques(CYCLE4) == ((0, 1), (0, 3), (1, 2), (2, 3))
    objective = np.zeros((4, 4))
    for left, right, sign in ((0, 1, 1), (1, 2, 1), (2, 3, -1), (3, 0, 1)):
        objective[left, right] = objective[right, left] = sign
    result = best_rank_one(np.eye(4), objective, CYCLE4, range(4))
    assert abs(result.value - 1.0) < 1e-12
    assert set(np.flatnonzero(result.vector)) <= set(result.clique)
    assert forbidden_component_pairs(CYCLE4, range(4)) == ((0, 2), (1, 3))


def test_psd_factor_round_trip():
    source = np.array([[1.0, 2.0, 0.0], [0.0, -1.0, 3.0]])
    matrix = source.T @ source
    factor, rank = factor_psd(matrix)
    assert rank == 2
    assert np.allclose(factor.T @ factor, matrix)


def test_sdp_can_strictly_beat_rank_one_on_a_nonchordal_mask():
    pytest.importorskip("cvxopt")
    objective = np.zeros((4, 4))
    for left, right, sign in ((0, 1, 1), (1, 2, 1), (2, 3, -1), (3, 0, 1)):
        objective[left, right] = objective[right, left] = sign
    rank_one = best_rank_one(np.eye(4), objective, CYCLE4, range(4))
    psd = solve_sparse_psd(np.eye(4), objective, CYCLE4, range(4))
    assert rank_one.value == pytest.approx(1.0)
    assert psd.value == pytest.approx(np.sqrt(2), rel=1e-7)
    assert psd.rank == 2
    assert psd.forbidden_max_abs < 1e-8
    assert psd.normalization == pytest.approx(1.0, abs=1e-7)
