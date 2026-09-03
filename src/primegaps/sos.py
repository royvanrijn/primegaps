"""Sparse positive-semidefinite relaxations for sum-of-squares sieve weights.

The matrices supplied here are numerical screening matrices.  Exact sieve
certification remains a separate step.  Components may be grouped by analytic
support cell; a forbidden cell pair expands to a zero block in ``Q``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class RankOneResult:
    value: float
    clique: tuple[int, ...]
    vector: np.ndarray


@dataclass(frozen=True)
class SparsePSDResult:
    value: float
    matrix: np.ndarray
    factor: np.ndarray
    rank: int
    normalization: float
    forbidden_max_abs: float
    minimum_eigenvalue: float
    relative_gap: float
    status: str


def _validate_graph(allowed: np.ndarray) -> np.ndarray:
    graph = np.asarray(allowed, dtype=bool)
    if graph.ndim != 2 or graph.shape[0] != graph.shape[1]:
        raise ValueError("allowed must be a square matrix")
    if not np.array_equal(graph, graph.T):
        raise ValueError("allowed must be symmetric")
    if not np.all(np.diag(graph)):
        raise ValueError("every component group must have a certified diagonal")
    return graph


def maximal_cliques(allowed: np.ndarray) -> tuple[tuple[int, ...], ...]:
    """Enumerate maximal cliques with the Bron--Kerbosch pivot algorithm."""
    graph = _validate_graph(allowed)
    size = len(graph)

    def neighbors(vertex: int) -> set[int]:
        return {other for other in range(size) if other != vertex and graph[vertex, other]}

    found: list[tuple[int, ...]] = []

    def visit(chosen: set[int], candidates: set[int], excluded: set[int]) -> None:
        if not candidates and not excluded:
            found.append(tuple(sorted(chosen)))
            return
        pool = candidates | excluded
        pivot = max(pool, key=lambda item: len(candidates & neighbors(item))) if pool else None
        extension = candidates - (neighbors(pivot) if pivot is not None else set())
        for vertex in tuple(extension):
            linked = neighbors(vertex)
            visit(chosen | {vertex}, candidates & linked, excluded & linked)
            candidates.remove(vertex)
            excluded.add(vertex)

    visit(set(), set(range(size)), set())
    return tuple(sorted(found))


def _validate_forms(denominator: np.ndarray, objective: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    first = np.asarray(denominator, dtype=float)
    second = np.asarray(objective, dtype=float)
    if first.ndim != 2 or first.shape[0] != first.shape[1] or second.shape != first.shape:
        raise ValueError("denominator and objective must be equally sized square matrices")
    if not np.allclose(first, first.T, rtol=0.0, atol=1e-12):
        raise ValueError("denominator is not symmetric")
    if not np.allclose(second, second.T, rtol=0.0, atol=1e-12):
        raise ValueError("objective is not symmetric")
    return (first + first.T) / 2, (second + second.T) / 2


def _largest_generalized(
    denominator: np.ndarray,
    objective: np.ndarray,
    *,
    relative_cutoff: float,
) -> tuple[float, np.ndarray]:
    values, vectors = np.linalg.eigh(denominator)
    if values[-1] <= 0:
        raise np.linalg.LinAlgError("denominator has no positive direction")
    keep = values > values[-1] * relative_cutoff
    whitening = vectors[:, keep] / np.sqrt(values[keep])[None, :]
    reduced = whitening.T @ objective @ whitening
    eigenvalues, eigenvectors = np.linalg.eigh((reduced + reduced.T) / 2)
    vector = whitening @ eigenvectors[:, -1]
    vector /= np.sqrt(float(vector @ denominator @ vector))
    return float(eigenvalues[-1]), vector


def best_rank_one(
    denominator: np.ndarray,
    objective: np.ndarray,
    allowed: np.ndarray,
    group_of_component: Sequence[int],
    *,
    relative_cutoff: float = 1e-11,
) -> RankOneResult:
    """Best rank-one quotient satisfying every forbidden zero constraint.

    For ``Q=cc^T``, every two nonzero component groups must be adjacent, so the
    support of ``c`` lies in a clique.  It is therefore enough to solve one
    generalized eigenproblem per maximal clique.
    """
    denominator, objective = _validate_forms(denominator, objective)
    graph = _validate_graph(allowed)
    groups = np.asarray(group_of_component, dtype=int)
    if groups.shape != (len(denominator),):
        raise ValueError("group_of_component has the wrong length")
    if np.any(groups < 0) or np.any(groups >= len(graph)):
        raise ValueError("group_of_component contains an unknown group")

    best_value = -np.inf
    best_clique: tuple[int, ...] = ()
    best_vector = np.zeros(len(denominator))
    for clique in maximal_cliques(graph):
        indices = np.flatnonzero(np.isin(groups, clique))
        value, local = _largest_generalized(
            denominator[np.ix_(indices, indices)],
            objective[np.ix_(indices, indices)],
            relative_cutoff=relative_cutoff,
        )
        if value > best_value:
            best_value = value
            best_clique = clique
            best_vector = np.zeros(len(denominator))
            best_vector[indices] = local
    return RankOneResult(best_value, best_clique, best_vector)


def forbidden_component_pairs(
    allowed: np.ndarray, group_of_component: Sequence[int]
) -> tuple[tuple[int, int], ...]:
    """Expand forbidden group interactions to upper-triangular component pairs."""
    graph = _validate_graph(allowed)
    groups = tuple(int(value) for value in group_of_component)
    if any(value < 0 or value >= len(graph) for value in groups):
        raise ValueError("group_of_component contains an unknown group")
    return tuple(
        (left, right)
        for left in range(len(groups))
        for right in range(left + 1, len(groups))
        if not graph[groups[left], groups[right]]
    )


def factor_psd(matrix: np.ndarray, *, relative_cutoff: float = 1e-8) -> tuple[np.ndarray, int]:
    """Return ``V`` with ``Q approximately V.T @ V`` after spectral truncation."""
    symmetric = (np.asarray(matrix, dtype=float) + np.asarray(matrix, dtype=float).T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    scale = max(float(values[-1]), np.finfo(float).tiny)
    keep = values > scale * relative_cutoff
    factor = np.sqrt(values[keep])[:, None] * vectors[:, keep].T
    return factor, int(np.sum(keep))


def solve_sparse_psd(
    denominator: np.ndarray,
    objective: np.ndarray,
    allowed: np.ndarray,
    group_of_component: Sequence[int],
    *,
    tolerance: float = 1e-8,
    factor_cutoff: float = 1e-7,
    show_progress: bool = False,
) -> SparsePSDResult:
    """Solve the sparse SDP with CVXOPT's primal-dual interior-point solver.

    This uses the dual formulation ``min y`` subject to
    ``y*M_I + sum z_ab*E_ab - M_J >= 0``.  CVXOPT's PSD dual variable is the
    requested primal matrix ``Q``.  Install the optional ``sdp`` dependency or
    run under a Sage environment that provides CVXOPT.
    """
    try:
        from cvxopt import matrix as cvx_matrix
        from cvxopt import solvers
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("solve_sparse_psd requires CVXOPT") from exc

    denominator, objective = _validate_forms(denominator, objective)
    graph = _validate_graph(allowed)
    groups = tuple(int(value) for value in group_of_component)
    pairs = forbidden_component_pairs(graph, groups)
    size = len(denominator)
    constraints = [denominator]
    for left, right in pairs:
        selector = np.zeros((size, size))
        selector[left, right] = selector[right, left] = 1.0
        constraints.append(selector)
    cone_map = np.column_stack(
        [-item.reshape(-1, order="F") for item in constraints]
    )
    linear_objective = np.zeros(len(constraints))
    linear_objective[0] = 1.0
    solvers.options["show_progress"] = show_progress
    solvers.options["abstol"] = tolerance
    solvers.options["reltol"] = tolerance
    solvers.options["feastol"] = tolerance
    solution = solvers.sdp(
        cvx_matrix(linear_objective),
        Gs=[cvx_matrix(cone_map)],
        hs=[cvx_matrix(-objective)],
    )
    if solution["status"] != "optimal":
        raise RuntimeError(f"CVXOPT SDP status is {solution['status']!r}")
    q_matrix = np.asarray(solution["zs"][0], dtype=float)
    q_matrix = (q_matrix + q_matrix.T) / 2
    eigenvalues = np.linalg.eigvalsh(q_matrix)
    factor, rank = factor_psd(q_matrix, relative_cutoff=factor_cutoff)
    forbidden_max = max(
        (abs(float(q_matrix[left, right])) for left, right in pairs), default=0.0
    )
    normalization = float(np.sum(denominator * q_matrix))
    return SparsePSDResult(
        value=float(np.sum(objective * q_matrix)) / normalization,
        matrix=q_matrix,
        factor=factor,
        rank=rank,
        normalization=normalization,
        forbidden_max_abs=forbidden_max,
        minimum_eigenvalue=float(eigenvalues[0]),
        relative_gap=float(solution["relative gap"]),
        status=str(solution["status"]),
    )
