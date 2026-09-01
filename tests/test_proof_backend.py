from dataclasses import replace
from fractions import Fraction
import json

import numpy as np
import pytest

from primegaps.certificate import (
    build_certificate,
    rationalize_candidate,
    verify_certificate,
)
from primegaps.eigen import solve_generalized_eigenproblem
from primegaps.exact_matrix import (
    ExactSymmetricMatrix,
    load_matrix_pair,
    load_matrix_pair_json,
    save_matrix_pair_json,
)
from primegaps.proof_backend import main, solve_and_certify


def test_exact_packed_and_sparse_storage_have_same_semantic_hash_and_form():
    dense = [
        [Fraction(1, 2), Fraction(1, 3), 0],
        [Fraction(1, 3), Fraction(2, 5), -1],
        [0, -1, 7],
    ]
    packed = ExactSymmetricMatrix.from_dense(dense)
    sparse = ExactSymmetricMatrix.from_sparse(
        3,
        [
            (0, 0, Fraction(1, 2)),
            (1, 0, Fraction(1, 3)),
            (1, 1, Fraction(2, 5)),
            (2, 1, -1),
            (2, 2, 7),
        ],
    )
    vector = (2, -3, 5)
    assert packed.semantic_sha256() == sparse.semantic_sha256()
    assert packed.quadratic_form(vector) == sparse.quadratic_form(vector)
    assert np.array_equal(packed.to_dense_float(), sparse.to_dense_float())


def test_dense_and_iterative_generalized_solvers_agree():
    rng = np.random.default_rng(260831126)
    raw = rng.standard_normal((12, 12))
    m1 = raw.T @ raw + np.diag(np.geomspace(0.5, 5.0, 12))
    raw2 = rng.standard_normal((12, 12))
    m2 = (raw2 + raw2.T) * 0.5
    dense = solve_generalized_eigenproblem(
        m1, m2, method="dense", exploit_blocks=False, diagnostic_iterations=4
    )
    iterative = solve_generalized_eigenproblem(
        m1,
        m2,
        method="iterative",
        exploit_blocks=False,
        max_iterations=12,
        tolerance=1e-12,
        diagnostic_iterations=4,
    )
    assert abs(dense.quotient - iterative.quotient) < 1e-9
    assert iterative.diagnostics.generalized_residual < 1e-9
    assert iterative.diagnostics.converged


def test_joint_block_structure_is_exploited():
    m1 = np.array([[2.0, 0.2, 0.0], [0.2, 1.0, 0.0], [0.0, 0.0, 3.0]])
    m2 = np.array([[1.0, 0.1, 0.0], [0.1, 1.0, 0.0], [0.0, 0.0, 6.0]])
    result = solve_generalized_eigenproblem(m1, m2, method="dense")
    assert result.quotient == pytest.approx(2.0)
    assert result.diagnostics.block_sizes == (2, 1)
    assert result.diagnostics.method == "block/dense"
    assert np.count_nonzero(result.vector[:2]) == 0


def test_exact_certificate_build_verify_and_tamper_rejection():
    m1 = ExactSymmetricMatrix.from_dense([[2, 0], [0, 1]])
    m2 = ExactSymmetricMatrix.from_dense([[5, 0], [0, 1]])
    candidate = rationalize_candidate(m1, m2, [1.0, 0.0], max_scale=1_000)
    certificate = build_certificate(m1, m2, candidate)
    verified = verify_certificate(m1, m2, certificate)
    assert verified.valid
    assert verified.quotient == Fraction(5, 2)
    assert verified.difference == 3

    tampered = replace(certificate, quotient=Fraction(999))
    with pytest.raises(ValueError, match="quotient"):
        verify_certificate(m1, m2, tampered)

    changed_m2 = ExactSymmetricMatrix.from_dense([[6, 0], [0, 1]])
    with pytest.raises(ValueError, match="hash"):
        verify_certificate(m1, changed_m2, certificate)


def test_end_to_end_backend_and_cli_exact_replay(tmp_path, capsys):
    m1 = ExactSymmetricMatrix.from_dense([[2, 1], [1, 2]])
    m2 = ExactSymmetricMatrix.from_dense([[4, 0], [0, 1]])
    result = solve_and_certify(m1, m2, method="dense", max_scale=10_000)
    assert result.numerical.quotient > 1
    assert result.certificate is not None
    assert verify_certificate(m1, m2, result.certificate).valid

    matrices_path = tmp_path / "matrices.json"
    certificate_path = tmp_path / "certificate.json"
    save_matrix_pair_json(matrices_path, m1, m2)
    loaded_m1, loaded_m2 = load_matrix_pair_json(matrices_path)
    assert loaded_m1.semantic_sha256() == m1.semantic_sha256()
    assert loaded_m2.semantic_sha256() == m2.semantic_sha256()

    assert main(
        [
            "solve",
            str(matrices_path),
            "--certificate",
            str(certificate_path),
            "--method",
            "dense",
        ]
    ) == 0
    solve_output = json.loads(capsys.readouterr().out)
    assert solve_output["certificate_produced"]
    assert certificate_path.exists()

    assert main(["verify", str(matrices_path), str(certificate_path)]) == 0
    verify_output = json.loads(capsys.readouterr().out)
    assert verify_output["valid"]


def test_backend_reports_when_no_strict_certificate_exists():
    m1 = ExactSymmetricMatrix.from_dense([[1, 0], [0, 1]])
    m2 = ExactSymmetricMatrix.from_dense([[1, 0], [0, Fraction(1, 2)]])
    result = solve_and_certify(m1, m2, method="dense")
    assert result.numerical.quotient == pytest.approx(1.0)
    assert result.certificate is None


def test_packed_npz_round_trip_is_exact(tmp_path):
    path = tmp_path / "matrices.npz"
    np.savez(
        path,
        m1_numerators=np.array([6, 2, 9], dtype=np.int64),
        m1_denominator=np.array(3, dtype=np.int64),
        m2_numerators=np.array([12, 1, 6], dtype=np.int64),
        m2_denominator=np.array(2, dtype=np.int64),
    )
    m1, m2 = load_matrix_pair(path)
    assert m1.quadratic_form((2, -1)) == Fraction(25, 3)
    assert m2.quadratic_form((2, -1)) == 25
    assert m1.storage == m2.storage == "packed-upper"


def test_solver_rejects_nonsymmetric_and_non_spd_inputs():
    with pytest.raises(ValueError, match="symmetric"):
        solve_generalized_eigenproblem([[1.0, 1.0], [0.0, 1.0]], np.eye(2))
    with pytest.raises(np.linalg.LinAlgError, match="positive"):
        solve_generalized_eigenproblem(
            [[1.0, 2.0], [2.0, 1.0]], np.eye(2), exploit_blocks=False
        )
