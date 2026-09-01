"""Numerical generalized-eigenvalue search and diagnostics.

Floating point locates a useful direction. Exact certificate verification is
implemented separately and does not trust any result from this module.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from typing import Callable, Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class ConditioningDiagnostics:
    dimension: int
    storage: str
    method: str
    block_sizes: tuple[int, ...]
    symmetry_error_m1: float
    symmetry_error_m2: float
    diagonal_min: float
    diagonal_max: float
    estimated_m1_eigenvalue_min: float
    estimated_m1_eigenvalue_max: float
    estimated_condition_m1: float
    estimated_condition_equilibrated_m1: float
    cholesky_diagonal_min: float | None
    cholesky_diagonal_max: float | None
    generalized_residual: float
    converged: bool
    iterations: int
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["block_sizes"] = list(self.block_sizes)
        result["warnings"] = list(self.warnings)
        return result


@dataclass(frozen=True)
class GeneralizedEigenResult:
    quotient: float
    vector: np.ndarray
    diagnostics: ConditioningDiagnostics


def _is_sparse(matrix: object) -> bool:
    try:
        from scipy.sparse import issparse
    except ImportError:
        return False
    return bool(issparse(matrix))


def _matrix_shape(matrix: object) -> tuple[int, ...]:
    return tuple(int(value) for value in getattr(matrix, "shape", ()))


def _max_abs(matrix) -> float:
    if _is_sparse(matrix):
        return float(np.max(np.abs(matrix.data))) if matrix.nnz else 0.0
    array = np.asarray(matrix)
    return float(np.max(np.abs(array))) if array.size else 0.0


def _symmetry_error(matrix) -> float:
    if _is_sparse(matrix):
        absolute = _max_abs(matrix - matrix.T)
    else:
        array = np.asarray(matrix, dtype=float)
        absolute = float(np.max(np.abs(array - array.T)))
    return absolute / max(1.0, _max_abs(matrix))


def _connected_blocks(m1, m2, tolerance: float) -> tuple[np.ndarray, ...]:
    """Connected components of the joint off-diagonal nonzero pattern."""
    dimension = _matrix_shape(m1)[0]
    if dimension == 1:
        return (np.array([0], dtype=int),)
    if _is_sparse(m1):
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components

        pattern = (abs(m1) > tolerance).astype(np.int8) + (abs(m2) > tolerance).astype(np.int8)
        pattern.setdiag(0)
        pattern.eliminate_zeros()
        count, labels = connected_components(csr_matrix(pattern), directed=False)
        return tuple(np.flatnonzero(labels == label) for label in range(count))

    parent = np.arange(dimension)

    def find(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = int(parent[item])
        return item

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    a = np.asarray(m1)
    b = np.asarray(m2)
    for i in range(dimension):
        linked = np.flatnonzero(
            (np.abs(a[i, i + 1 :]) > tolerance) | (np.abs(b[i, i + 1 :]) > tolerance)
        )
        for relative_j in linked:
            union(i, i + 1 + int(relative_j))
    groups: dict[int, list[int]] = {}
    for i in range(dimension):
        groups.setdefault(find(i), []).append(i)
    return tuple(np.asarray(indices, dtype=int) for indices in groups.values())


def _forward_substitution(lower: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    result = np.empty_like(rhs, dtype=float)
    for i in range(len(rhs)):
        result[i] = (rhs[i] - np.dot(lower[i, :i], result[:i])) / lower[i, i]
    return result


def _back_substitution_transpose(lower: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    result = np.empty_like(rhs, dtype=float)
    for i in range(len(rhs) - 1, -1, -1):
        result[i] = (rhs[i] - np.dot(lower[i + 1 :, i], result[i + 1 :])) / lower[i, i]
    return result


def _solve_cholesky(lower: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    return _back_substitution_transpose(lower, _forward_substitution(lower, rhs))


def _largest_spd_eigenvalue(
    operator: Callable[[np.ndarray], np.ndarray], dimension: int, iterations: int
) -> float:
    vector = np.arange(1, dimension + 1, dtype=float)
    vector /= np.linalg.norm(vector)
    value = 0.0
    for _ in range(max(2, iterations)):
        image = operator(vector)
        norm = float(np.linalg.norm(image))
        if not np.isfinite(norm) or norm == 0.0:
            return float("nan")
        vector = image / norm
        value = float(np.dot(vector, operator(vector)))
    return value


def _lanczos_largest(
    operator: Callable[[np.ndarray], np.ndarray],
    dimension: int,
    *,
    tolerance: float,
    max_iterations: int,
    seed: int,
) -> tuple[float, np.ndarray, bool, int]:
    if dimension == 1:
        unit = np.ones(1)
        return float(operator(unit)[0]), unit, True, 1
    iterations = min(max_iterations, dimension)
    rng = np.random.default_rng(seed)
    q = rng.standard_normal(dimension)
    q /= np.linalg.norm(q)
    previous = np.zeros(dimension)
    previous_beta = 0.0
    basis = np.empty((dimension, iterations), dtype=float)
    diagonal = np.empty(iterations, dtype=float)
    off_diagonal = np.empty(max(0, iterations - 1), dtype=float)
    best_value = float("nan")
    best_vector = q.copy()

    for j in range(iterations):
        basis[:, j] = q
        image = operator(q) - previous_beta * previous
        alpha = float(np.dot(q, image))
        diagonal[j] = alpha
        image -= alpha * q
        active_basis = basis[:, : j + 1]
        for _ in range(2):
            image -= active_basis @ (active_basis.T @ image)
        beta = float(np.linalg.norm(image))
        active = j + 1
        if active == 1:
            ritz_values = diagonal[:1].copy()
            ritz_vectors = np.ones((1, 1))
        else:
            tridiagonal = np.diag(diagonal[:active])
            links = off_diagonal[: active - 1]
            tridiagonal += np.diag(links, 1) + np.diag(links, -1)
            ritz_values, ritz_vectors = np.linalg.eigh(tridiagonal)
        best_value = float(ritz_values[-1])
        best_vector = active_basis @ ritz_vectors[:, -1]
        residual_estimate = abs(beta * float(ritz_vectors[-1, -1]))
        if beta <= np.finfo(float).eps or residual_estimate <= tolerance * max(1.0, abs(best_value)):
            best_vector /= np.linalg.norm(best_vector)
            return best_value, best_vector, True, active
        if j == iterations - 1:
            break
        off_diagonal[j] = beta
        previous, q = q, image / beta
        previous_beta = beta

    best_vector /= np.linalg.norm(best_vector)
    return best_value, best_vector, False, iterations


def _relative_generalized_residual(m1, m2, value: float, vector: np.ndarray) -> float:
    left = np.asarray(m2 @ vector).reshape(-1)
    right = np.asarray(m1 @ vector).reshape(-1)
    residual = np.linalg.norm(left - value * right)
    scale = np.linalg.norm(left) + abs(value) * np.linalg.norm(right)
    return float(residual / max(scale, np.finfo(float).tiny))


def _solve_dense_core(
    m1: np.ndarray,
    m2: np.ndarray,
    *,
    method: str,
    tolerance: float,
    max_iterations: int,
    seed: int,
    diagnostic_iterations: int,
) -> GeneralizedEigenResult:
    dimension = len(m1)
    diagonal = np.diag(m1)
    if np.any(diagonal <= 0.0):
        raise np.linalg.LinAlgError("M1 has a non-positive diagonal entry and cannot be SPD")
    scale = 1.0 / np.sqrt(diagonal)
    equilibrated_m1 = (scale[:, None] * m1) * scale[None, :]
    try:
        lower = np.linalg.cholesky(equilibrated_m1)
    except np.linalg.LinAlgError as exc:
        raise np.linalg.LinAlgError("M1 is not numerically positive definite") from exc
    del equilibrated_m1

    def solve_equilibrated(rhs: np.ndarray) -> np.ndarray:
        return _solve_cholesky(lower, rhs)

    raw_max = _largest_spd_eigenvalue(lambda x: m1 @ x, dimension, diagnostic_iterations)
    raw_inverse_max = _largest_spd_eigenvalue(
        lambda x: scale * solve_equilibrated(scale * x), dimension, diagnostic_iterations
    )
    equilibrated_max = _largest_spd_eigenvalue(
        lambda x: scale * (m1 @ (scale * x)), dimension, diagnostic_iterations
    )
    equilibrated_inverse_max = _largest_spd_eigenvalue(
        solve_equilibrated, dimension, diagnostic_iterations
    )

    if method == "dense":
        equilibrated_m2 = (scale[:, None] * m2) * scale[None, :]
        transformed_left = np.linalg.solve(lower, equilibrated_m2)
        del equilibrated_m2
        transformed = np.linalg.solve(lower, transformed_left.T).T
        del transformed_left
        transformed = (transformed + transformed.T) * 0.5
        values, vectors = np.linalg.eigh(transformed)
        y = vectors[:, -1]
        transformed_vector = np.linalg.solve(lower.T, y)
        iterations = dimension
        converged = True
    elif method == "iterative":
        def operator(y: np.ndarray) -> np.ndarray:
            transformed_vector = _back_substitution_transpose(lower, y)
            product = scale * (m2 @ (scale * transformed_vector))
            return _forward_substitution(lower, product)

        _, y, converged, iterations = _lanczos_largest(
            operator,
            dimension,
            tolerance=tolerance,
            max_iterations=max_iterations,
            seed=seed,
        )
        transformed_vector = _back_substitution_transpose(lower, y)
    else:
        raise ValueError(f"unsupported dense solve method {method!r}")

    vector = scale * transformed_vector
    norm = float(np.max(np.abs(vector)))
    if not np.isfinite(norm) or norm == 0.0:
        raise ArithmeticError("generalized eigensolver returned an invalid vector")
    vector /= norm
    numerator = float(vector @ (m2 @ vector))
    denominator = float(vector @ (m1 @ vector))
    quotient = numerator / denominator
    residual = _relative_generalized_residual(m1, m2, quotient, vector)
    condition = raw_max * raw_inverse_max
    equilibrated_condition = equilibrated_max * equilibrated_inverse_max
    warnings: list[str] = []
    if condition * np.finfo(float).eps > 1e-6:
        warnings.append("M1 is severely ill-conditioned at float64 precision")
    if not converged:
        warnings.append("iterative eigensolver reached its iteration limit")
    if residual > max(10.0 * tolerance, 1e-10):
        warnings.append("generalized eigenpair residual is larger than requested tolerance")
    diagnostics = ConditioningDiagnostics(
        dimension=dimension,
        storage="dense",
        method=method,
        block_sizes=(dimension,),
        symmetry_error_m1=_symmetry_error(m1),
        symmetry_error_m2=_symmetry_error(m2),
        diagonal_min=float(np.min(diagonal)),
        diagonal_max=float(np.max(diagonal)),
        estimated_m1_eigenvalue_min=float(1.0 / raw_inverse_max),
        estimated_m1_eigenvalue_max=float(raw_max),
        estimated_condition_m1=float(condition),
        estimated_condition_equilibrated_m1=float(equilibrated_condition),
        cholesky_diagonal_min=float(np.min(np.diag(lower))),
        cholesky_diagonal_max=float(np.max(np.diag(lower))),
        generalized_residual=residual,
        converged=converged,
        iterations=iterations,
        warnings=tuple(warnings),
    )
    return GeneralizedEigenResult(quotient, vector, diagnostics)


def _solve_sparse_core(m1, m2, *, tolerance: float, max_iterations: int) -> GeneralizedEigenResult:
    try:
        from scipy.sparse import diags
        from scipy.sparse.linalg import eigsh
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("install primegaps[sparse] for sparse numerical solving") from exc

    dimension = m1.shape[0]
    diagonal = np.asarray(m1.diagonal(), dtype=float)
    if np.any(diagonal <= 0.0):
        raise np.linalg.LinAlgError("M1 has a non-positive diagonal entry and cannot be SPD")
    scale = 1.0 / np.sqrt(diagonal)
    scaling = diags(scale)
    equilibrated_m1 = scaling @ m1 @ scaling
    equilibrated_m2 = scaling @ m2 @ scaling
    if dimension == 1:
        y = np.ones(1)
        iterations = 1
    else:
        _, vectors = eigsh(
            equilibrated_m2,
            k=1,
            M=equilibrated_m1,
            which="LA",
            tol=tolerance,
            maxiter=max_iterations,
        )
        y = vectors[:, 0]
        iterations = max_iterations
    vector = scale * y
    vector /= np.max(np.abs(vector))
    quotient = float((vector @ (m2 @ vector)) / (vector @ (m1 @ vector)))
    residual = _relative_generalized_residual(m1, m2, quotient, vector)
    try:
        largest = float(eigsh(m1, k=1, which="LA", return_eigenvectors=False)[0])
        smallest = float(eigsh(m1, k=1, which="SA", return_eigenvectors=False)[0])
        eq_largest = float(eigsh(equilibrated_m1, k=1, which="LA", return_eigenvectors=False)[0])
        eq_smallest = float(eigsh(equilibrated_m1, k=1, which="SA", return_eigenvectors=False)[0])
    except Exception:
        largest = smallest = eq_largest = eq_smallest = float("nan")
    warnings: list[str] = []
    if np.isfinite(smallest) and smallest <= 0.0:
        warnings.append("M1 is not numerically positive definite")
    if residual > max(10.0 * tolerance, 1e-10):
        warnings.append("generalized eigenpair residual is larger than requested tolerance")
    diagnostics = ConditioningDiagnostics(
        dimension=dimension,
        storage="sparse",
        method="scipy-eigsh",
        block_sizes=(dimension,),
        symmetry_error_m1=_symmetry_error(m1),
        symmetry_error_m2=_symmetry_error(m2),
        diagonal_min=float(np.min(diagonal)),
        diagonal_max=float(np.max(diagonal)),
        estimated_m1_eigenvalue_min=smallest,
        estimated_m1_eigenvalue_max=largest,
        estimated_condition_m1=largest / smallest,
        estimated_condition_equilibrated_m1=eq_largest / eq_smallest,
        cholesky_diagonal_min=None,
        cholesky_diagonal_max=None,
        generalized_residual=residual,
        converged=True,
        iterations=iterations,
        warnings=tuple(warnings),
    )
    return GeneralizedEigenResult(quotient, vector, diagnostics)


def solve_generalized_eigenproblem(
    m1,
    m2,
    *,
    method: str = "auto",
    dense_threshold: int = 1_200,
    tolerance: float = 1e-11,
    max_iterations: int = 180,
    seed: int = 0,
    diagnostic_iterations: int = 8,
    exploit_blocks: bool = True,
    block_tolerance: float = 0.0,
    symmetry_tolerance: float = 1e-12,
) -> GeneralizedEigenResult:
    """Return the largest generalized eigenpair of symmetric M2 and M1.

    M1 must be positive definite. Auto uses a complete dense solve up to
    dense_threshold and matrix-free Lanczos above it. SciPy sparse inputs use
    ARPACK. Exact simultaneous block structure is detected before solving.
    """
    if method not in {"auto", "dense", "iterative", "sparse"}:
        raise ValueError("method must be auto, dense, iterative, or sparse")
    if tolerance <= 0.0 or max_iterations < 1 or dense_threshold < 1:
        raise ValueError("solver tolerances and iteration limits must be positive")
    if not _matrix_shape(m1):
        m1 = np.asarray(m1, dtype=float)
    if not _matrix_shape(m2):
        m2 = np.asarray(m2, dtype=float)
    shape1, shape2 = _matrix_shape(m1), _matrix_shape(m2)
    if shape1 != shape2 or len(shape1) != 2 or shape1[0] != shape1[1] or shape1[0] < 1:
        raise ValueError("M1 and M2 must be non-empty square matrices of equal shape")
    sparse = _is_sparse(m1) or _is_sparse(m2)
    if _is_sparse(m1) != _is_sparse(m2):
        raise ValueError("M1 and M2 must use the same dense/sparse storage class")
    if not sparse:
        m1 = np.asarray(m1, dtype=float)
        m2 = np.asarray(m2, dtype=float)
        if not np.all(np.isfinite(m1)) or not np.all(np.isfinite(m2)):
            raise ValueError("matrices must contain only finite values")
    symmetry1, symmetry2 = _symmetry_error(m1), _symmetry_error(m2)
    if max(symmetry1, symmetry2) > symmetry_tolerance:
        raise ValueError(
            f"matrices are not symmetric within tolerance: errors {symmetry1:.3e}, {symmetry2:.3e}"
        )

    blocks = _connected_blocks(m1, m2, block_tolerance) if exploit_blocks else (
        np.arange(shape1[0]),
    )
    if len(blocks) > 1:
        results: list[GeneralizedEigenResult] = []
        for block_index, indices in enumerate(blocks):
            if sparse:
                block_m1 = m1[indices][:, indices]
                block_m2 = m2[indices][:, indices]
            else:
                block_m1 = m1[np.ix_(indices, indices)]
                block_m2 = m2[np.ix_(indices, indices)]
            results.append(
                solve_generalized_eigenproblem(
                    block_m1,
                    block_m2,
                    method=method,
                    dense_threshold=dense_threshold,
                    tolerance=tolerance,
                    max_iterations=max_iterations,
                    seed=seed + block_index,
                    diagnostic_iterations=diagnostic_iterations,
                    exploit_blocks=False,
                    symmetry_tolerance=symmetry_tolerance,
                )
            )
        best_index = max(range(len(results)), key=lambda index: results[index].quotient)
        best = results[best_index]
        vector = np.zeros(shape1[0])
        vector[blocks[best_index]] = best.vector
        minimum = min(result.diagnostics.estimated_m1_eigenvalue_min for result in results)
        maximum = max(result.diagnostics.estimated_m1_eigenvalue_max for result in results)
        residual = _relative_generalized_residual(m1, m2, best.quotient, vector)
        chol_mins = [
            value
            for result in results
            if (value := result.diagnostics.cholesky_diagonal_min) is not None
        ]
        chol_maxs = [
            value
            for result in results
            if (value := result.diagnostics.cholesky_diagonal_max) is not None
        ]
        diagnostics = ConditioningDiagnostics(
            dimension=shape1[0],
            storage="sparse" if sparse else "dense",
            method=f"block/{best.diagnostics.method}",
            block_sizes=tuple(sorted((len(block) for block in blocks), reverse=True)),
            symmetry_error_m1=symmetry1,
            symmetry_error_m2=symmetry2,
            diagonal_min=min(result.diagnostics.diagonal_min for result in results),
            diagonal_max=max(result.diagnostics.diagonal_max for result in results),
            estimated_m1_eigenvalue_min=minimum,
            estimated_m1_eigenvalue_max=maximum,
            estimated_condition_m1=maximum / minimum,
            estimated_condition_equilibrated_m1=max(
                result.diagnostics.estimated_condition_equilibrated_m1 for result in results
            ),
            cholesky_diagonal_min=min(chol_mins) if chol_mins else None,
            cholesky_diagonal_max=max(chol_maxs) if chol_maxs else None,
            generalized_residual=residual,
            converged=all(result.diagnostics.converged for result in results),
            iterations=sum(result.diagnostics.iterations for result in results),
            warnings=tuple(
                dict.fromkeys(warning for result in results for warning in result.diagnostics.warnings)
            ),
        )
        return GeneralizedEigenResult(best.quotient, vector, diagnostics)

    if sparse:
        if method in {"dense", "iterative"}:
            raise ValueError("use method auto or sparse with SciPy sparse matrices")
        return _solve_sparse_core(m1, m2, tolerance=tolerance, max_iterations=max_iterations)
    selected = method
    if method == "auto":
        selected = "dense" if shape1[0] <= dense_threshold else "iterative"
    if selected == "sparse":
        raise ValueError("method sparse requires SciPy sparse matrices")
    return _solve_dense_core(
        m1,
        m2,
        method=selected,
        tolerance=tolerance,
        max_iterations=max_iterations,
        seed=seed,
        diagnostic_iterations=diagnostic_iterations,
    )


def largest_generalized_eigenpair(m1: np.ndarray, m2: np.ndarray) -> tuple[float, np.ndarray]:
    """Backward-compatible complete dense generalized eigensolve."""
    result = solve_generalized_eigenproblem(m1, m2, method="dense", exploit_blocks=False)
    return result.quotient, result.vector


def rationalize_vector(
    vector: Iterable[float], max_denominator: int = 1_000_000
) -> tuple[Fraction, ...]:
    """Convert a floating eigenvector into a rational direction."""
    return tuple(Fraction(float(value)).limit_denominator(max_denominator) for value in vector)


def exact_quadratic_form(
    matrix: Sequence[Sequence[Fraction | int]], vector: Sequence[Fraction]
) -> Fraction:
    dimension = len(vector)
    if len(matrix) != dimension or any(len(row) != dimension for row in matrix):
        raise ValueError("matrix/vector dimensions differ")
    return sum(
        vector[i] * Fraction(matrix[i][j]) * vector[j]
        for i in range(dimension)
        for j in range(dimension)
    )


def exact_rayleigh_quotient(
    m1: Sequence[Sequence[Fraction | int]],
    m2: Sequence[Sequence[Fraction | int]],
    vector: Sequence[Fraction],
) -> Fraction:
    denominator = exact_quadratic_form(m1, vector)
    if denominator == 0:
        raise ZeroDivisionError("c M1 c^T is zero")
    return exact_quadratic_form(m2, vector) / denominator
