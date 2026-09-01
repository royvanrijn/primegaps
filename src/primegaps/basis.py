from __future__ import annotations

from collections.abc import Iterator


def partitions(n: int, max_part: int | None = None) -> Iterator[tuple[int, ...]]:
    """Integer partitions of n, represented in non-increasing order."""
    if n == 0:
        yield ()
        return
    max_part = n if max_part is None else min(max_part, n)
    for first in range(max_part, 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def symmetric_basis(degree: int, k: int) -> list[tuple[tuple[int, ...], int]]:
    """Indices for monomial-symmetric p_lambda times (1-sum(t))^b, 2a+b<=D."""
    basis: list[tuple[tuple[int, ...], int]] = []
    for a in range(degree // 2 + 1):
        for lam in partitions(a):
            if len(lam) > k:
                continue
            for b in range(degree - 2 * a + 1):
                basis.append((lam, b))
    return basis
