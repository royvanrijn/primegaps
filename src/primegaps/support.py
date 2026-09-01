from __future__ import annotations

from dataclasses import dataclass
from math import floor

import numpy as np


@dataclass(frozen=True)
class SupportParameters:
    """Parameters from Definition 1 of arXiv:2608.31126.

    A has n+1 entries with A[0] = -epsilon. B has n rows and floor(1/delta)
    entries; row j belongs to the interval [A[j]+epsilon, A[j+1]+epsilon).
    """

    delta: float
    epsilon: float
    A: tuple[float, ...]
    B: tuple[tuple[float, ...], ...]

    def validate(self) -> None:
        n = len(self.A) - 1
        if n <= 0 or len(self.B) != n:
            raise ValueError("B must have len(A)-1 rows")
        if abs(self.A[0] + self.epsilon) > 1e-12:
            raise ValueError("Definition 1 requires A[0] = -epsilon")
        if any(x >= y for x, y in zip(self.A, self.A[1:])):
            raise ValueError("A must be strictly increasing")
        width = floor(1.0 / self.delta)
        for row in self.B:
            if len(row) != width:
                raise ValueError(f"B rows must have {width} entries")
            for i, value in enumerate(row):
                if value <= self.delta:
                    raise ValueError("B[j,m] must exceed delta")
                if i and not (row[i - 1] <= value <= row[i - 1] + self.delta + 1e-12):
                    raise ValueError("B[j,m] must be monotone and grow by at most delta")


def stadlmann_240_parameters() -> SupportParameters:
    """Published parameters from the proof of Theorem 1 (k=49, H1 <= 240)."""
    epsilon = 0.0075
    delta = 0.028
    width = floor(1.0 / delta)
    row = (0.15, 0.15) + (0.17,) * (width - 2)
    p = SupportParameters(
        delta=delta,
        epsilon=epsilon,
        A=(-epsilon, 0.253),
        B=(row,),
    )
    p.validate()
    return p


def contains(points: np.ndarray, p: SupportParameters) -> np.ndarray:
    """Vectorized membership in T_k(delta,A,B,epsilon)."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2:
        raise ValueError("points must have shape (samples, k)")
    if np.any(points < 0.0) or np.any(points > 1.0):
        return np.zeros(points.shape[0], dtype=bool)

    total = points.sum(axis=1)
    large = points > p.delta
    m = large.sum(axis=1)
    large_sum = np.where(large, points, 0.0).sum(axis=1)

    result = np.zeros(points.shape[0], dtype=bool)
    for j, row in enumerate(p.B):
        lo = p.A[j] + p.epsilon
        hi = p.A[j + 1] + p.epsilon
        in_band = (total >= lo) & (total < hi)
        b_ok = m == 0
        nonempty = m > 0
        if np.any(nonempty):
            valid_m = nonempty & (m <= len(row))
            limits = np.zeros(points.shape[0])
            idx = np.flatnonzero(valid_m)
            limits[idx] = np.asarray(row)[m[idx] - 1]
            b_ok |= valid_m & (large_sum <= limits + 1e-15)
        result |= in_band & b_ok
    return result


def sample_uniform_simplex(k: int, total: float, samples: int, seed: int = 1) -> np.ndarray:
    """Uniform samples from {t_i >= 0, sum t_i <= total}."""
    rng = np.random.default_rng(seed)
    x = rng.exponential(size=(samples, k + 1))
    x /= x.sum(axis=1, keepdims=True)
    return x[:, :k] * total
