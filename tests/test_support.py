import numpy as np

from primegaps.basis import symmetric_basis
from primegaps.support import contains, stadlmann_240_parameters


def test_published_parameters_validate():
    stadlmann_240_parameters().validate()


def test_support_boundaries_and_large_coordinate_constraint():
    p = stadlmann_240_parameters()
    points = np.array([
        [0.01, 0.01, 0.01],
        [0.2604, 0.0, 0.0],
        [0.14, 0.11, 0.0],
        [0.08, 0.07, 0.02],
        [0.10, 0.04, 0.03],
        [0.10, 0.04, 0.031],
        [0.01, 0.01, 0.241],
        [0.2605, 0.0, 0.0],
    ])
    assert contains(points, p).tolist() == [True, False, False, True, True, False, False, False]


def test_basis_sizes_for_k49():
    assert len(symmetric_basis(21, 49)) == 846
    assert len(symmetric_basis(27, 49)) == 2526
