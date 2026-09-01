"""Compact exact symmetric matrices for proof certificates.

The numerical eigensolver is deliberately separate from this module.  Exact
verification needs only Python integer arithmetic and can use either a packed
upper triangle or a sparse list of upper-triangular entries.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt, lcm
from numbers import Integral
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np


SparseEntry = tuple[int, int, int]


def _packed_length(dimension: int) -> int:
    return dimension * (dimension + 1) // 2


def _infer_packed_dimension(length: int) -> int:
    if length < 1:
        raise ValueError("a packed symmetric matrix cannot be empty")
    dimension = (isqrt(8 * length + 1) - 1) // 2
    if _packed_length(dimension) != length:
        raise ValueError(f"packed length {length} is not triangular")
    return dimension


@dataclass(frozen=True)
class ExactSymmetricMatrix:
    """A common-denominator exact symmetric matrix.

    Exactly one of ``packed_upper`` and ``sparse_upper`` must be supplied.
    Packed values are ordered row by row: ``(0,0), (0,1), ..., (1,1), ...``.
    Sparse entries are sorted ``(row, column, numerator)`` triples with
    ``row <= column``; omitted entries are zero.
    """

    dimension: int
    denominator: int
    packed_upper: Sequence[int] | None = None
    sparse_upper: tuple[SparseEntry, ...] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.dimension, Integral) or not isinstance(
            self.denominator, Integral
        ):
            raise TypeError("matrix dimension and denominator must be integers")
        if self.dimension < 1:
            raise ValueError("matrix dimension must be positive")
        if self.denominator <= 0:
            raise ValueError("matrix denominator must be positive")
        if (self.packed_upper is None) == (self.sparse_upper is None):
            raise ValueError("supply exactly one matrix storage representation")
        if self.packed_upper is not None:
            expected = _packed_length(self.dimension)
            if len(self.packed_upper) != expected:
                raise ValueError(
                    f"packed upper triangle has {len(self.packed_upper)} values; expected {expected}"
                )
            if isinstance(self.packed_upper, np.ndarray):
                if self.packed_upper.dtype.kind not in "iu":
                    raise TypeError("packed numerators must be integers")
            elif any(not isinstance(value, Integral) for value in self.packed_upper):
                raise TypeError("packed numerators must be integers")
            return

        assert self.sparse_upper is not None
        previous: tuple[int, int] | None = None
        for raw_i, raw_j, raw_value in self.sparse_upper:
            if not all(isinstance(value, Integral) for value in (raw_i, raw_j, raw_value)):
                raise TypeError("sparse indices and numerators must be integers")
            i, j, value = int(raw_i), int(raw_j), int(raw_value)
            if not (0 <= i <= j < self.dimension):
                raise ValueError(f"invalid sparse upper entry ({i}, {j})")
            if value == 0:
                raise ValueError("sparse storage must omit zero entries")
            if previous is not None and (i, j) <= previous:
                raise ValueError("sparse upper entries must be unique and sorted")
            previous = (i, j)

    @classmethod
    def from_dense(
        cls, matrix: Sequence[Sequence[int | Fraction]], *, require_symmetric: bool = True
    ) -> "ExactSymmetricMatrix":
        """Pack an integer/Fraction matrix using one common denominator."""
        dimension = len(matrix)
        if dimension < 1 or any(len(row) != dimension for row in matrix):
            raise ValueError("matrix must be non-empty and square")
        upper: list[Fraction] = []
        denominator = 1
        for i in range(dimension):
            for j in range(i, dimension):
                value = Fraction(matrix[i][j])
                if require_symmetric and value != Fraction(matrix[j][i]):
                    raise ValueError(f"matrix is not exactly symmetric at ({i}, {j})")
                upper.append(value)
                denominator = lcm(denominator, value.denominator)
        numerators = tuple(value.numerator * (denominator // value.denominator) for value in upper)
        common = denominator
        for value in numerators:
            common = gcd(common, abs(value))
        if common > 1:
            denominator //= common
            numerators = tuple(value // common for value in numerators)
        return cls(dimension, denominator, packed_upper=numerators)

    @classmethod
    def from_sparse(
        cls,
        dimension: int,
        entries: Iterable[tuple[int, int, int | Fraction]],
    ) -> "ExactSymmetricMatrix":
        """Build sparse upper storage from exact entries, combining duplicates."""
        combined: dict[tuple[int, int], Fraction] = {}
        for raw_i, raw_j, raw_value in entries:
            i, j = sorted((int(raw_i), int(raw_j)))
            if not (0 <= i <= j < dimension):
                raise ValueError(f"invalid sparse entry ({raw_i}, {raw_j})")
            combined[(i, j)] = combined.get((i, j), Fraction(0)) + Fraction(raw_value)
        denominator = 1
        for value in combined.values():
            denominator = lcm(denominator, value.denominator)
        sparse = tuple(
            (i, j, value.numerator * (denominator // value.denominator))
            for (i, j), value in sorted(combined.items())
            if value
        )
        common = denominator
        for _, _, value in sparse:
            common = gcd(common, abs(value))
        if common > 1:
            denominator //= common
            sparse = tuple((i, j, value // common) for i, j, value in sparse)
        return cls(dimension, denominator, sparse_upper=sparse)

    @property
    def storage(self) -> str:
        return "packed-upper" if self.packed_upper is not None else "sparse-upper"

    @property
    def stored_entries(self) -> int:
        return len(self.packed_upper if self.packed_upper is not None else self.sparse_upper or ())

    @property
    def nonzero_entries(self) -> int:
        if self.packed_upper is not None:
            return sum(int(value) != 0 for value in self.packed_upper)
        return len(self.sparse_upper or ())

    def _iter_nonzero(self) -> Iterator[SparseEntry]:
        if self.sparse_upper is not None:
            for i, j, value in self.sparse_upper:
                yield int(i), int(j), int(value)
            return
        assert self.packed_upper is not None
        offset = 0
        for i in range(self.dimension):
            for j in range(i, self.dimension):
                value = int(self.packed_upper[offset])
                offset += 1
                if value:
                    yield i, j, value

    def semantic_sha256(self) -> str:
        """Hash exact entries, independent of packed versus sparse storage."""
        common = self.denominator
        for _, _, value in self._iter_nonzero():
            common = gcd(common, abs(value))
        denominator = self.denominator // common
        digest = sha256()
        digest.update(f"primegaps-exact-symmetric-v1\n{self.dimension}\n{denominator}\n".encode())
        for i, j, value in self._iter_nonzero():
            digest.update(f"{i},{j},{value // common}\n".encode())
        return digest.hexdigest()

    def quadratic_form(self, vector: Sequence[int]) -> Fraction:
        """Compute ``vector.T @ self @ vector`` with Python big integers."""
        if len(vector) != self.dimension:
            raise ValueError("matrix/vector dimensions differ")
        if any(not isinstance(value, Integral) for value in vector):
            raise TypeError("quadratic-form vector must contain integers")
        c = tuple(int(value) for value in vector)
        numerator = 0
        for i, j, value in self._iter_nonzero():
            term = value * c[i] * c[j]
            numerator += term if i == j else 2 * term
        return Fraction(numerator, self.denominator)

    def to_dense_float(self) -> np.ndarray:
        """Materialize a float64 dense matrix for numerical search."""
        result = np.zeros((self.dimension, self.dimension), dtype=float)
        reciprocal = 1.0 / self.denominator
        if self.packed_upper is not None:
            offset = 0
            for i in range(self.dimension):
                width = self.dimension - i
                row = np.asarray(self.packed_upper[offset : offset + width], dtype=float) * reciprocal
                result[i, i:] = row
                result[i:, i] = row
                offset += width
            return result
        assert self.sparse_upper is not None
        for i, j, value in self.sparse_upper:
            float_value = int(value) * reciprocal
            result[i, j] = float_value
            result[j, i] = float_value
        return result

    def to_scipy_sparse(self):
        """Materialize CSR storage; SciPy is an optional dependency."""
        try:
            from scipy.sparse import coo_matrix
        except ImportError as exc:  # pragma: no cover - depends on optional extra
            raise RuntimeError("SciPy is required for sparse numerical solving") from exc
        rows: list[int] = []
        columns: list[int] = []
        values: list[float] = []
        reciprocal = 1.0 / self.denominator
        for i, j, value in self._iter_nonzero():
            rows.append(i)
            columns.append(j)
            values.append(value * reciprocal)
            if i != j:
                rows.append(j)
                columns.append(i)
                values.append(value * reciprocal)
        return coo_matrix((values, (rows, columns)), shape=(self.dimension, self.dimension)).tocsr()

    def to_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "dimension": self.dimension,
            "denominator": self.denominator,
            "storage": self.storage,
        }
        if self.packed_upper is not None:
            payload["numerators"] = [int(value) for value in self.packed_upper]
        else:
            payload["entries"] = [list(entry) for entry in self.sparse_upper or ()]
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "ExactSymmetricMatrix":
        dimension = int(payload["dimension"])
        denominator = int(payload["denominator"])
        storage = payload.get("storage")
        if storage == "packed-upper":
            raw = payload.get("numerators")
            if not isinstance(raw, list):
                raise ValueError("packed matrix numerators must be a list")
            return cls(dimension, denominator, packed_upper=tuple(int(value) for value in raw))
        if storage == "sparse-upper":
            raw = payload.get("entries")
            if not isinstance(raw, list):
                raise ValueError("sparse matrix entries must be a list")
            entries = tuple((int(entry[0]), int(entry[1]), int(entry[2])) for entry in raw)
            return cls(dimension, denominator, sparse_upper=entries)
        raise ValueError(f"unsupported exact matrix storage {storage!r}")


def save_matrix_pair_json(
    path: str | Path, m1: ExactSymmetricMatrix, m2: ExactSymmetricMatrix
) -> None:
    if m1.dimension != m2.dimension:
        raise ValueError("matrix dimensions differ")
    payload = {"format": "primegaps-exact-matrix-pair-v1", "m1": m1.to_dict(), "m2": m2.to_dict()}
    Path(path).write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")


def load_matrix_pair_json(path: str | Path) -> tuple[ExactSymmetricMatrix, ExactSymmetricMatrix]:
    payload = json.loads(Path(path).read_text())
    if payload.get("format") != "primegaps-exact-matrix-pair-v1":
        raise ValueError("unsupported exact matrix-pair format")
    m1 = ExactSymmetricMatrix.from_dict(payload["m1"])
    m2 = ExactSymmetricMatrix.from_dict(payload["m2"])
    if m1.dimension != m2.dimension:
        raise ValueError("matrix dimensions differ")
    return m1, m2


def load_matrix_pair_npz(path: str | Path) -> tuple[ExactSymmetricMatrix, ExactSymmetricMatrix]:
    """Load compact int64 packed triangles from an NPZ container."""
    with np.load(path, allow_pickle=False) as payload:
        matrices: list[ExactSymmetricMatrix] = []
        for name in ("m1", "m2"):
            numerators = np.asarray(payload[f"{name}_numerators"])
            if numerators.ndim != 1 or numerators.dtype.kind not in "iu":
                raise ValueError(f"{name}_numerators must be a one-dimensional integer array")
            denominator = int(np.asarray(payload[f"{name}_denominator"]).item())
            dimension = _infer_packed_dimension(len(numerators))
            # Copy because NpzFile closes its backing ZipExtFile on return.
            matrices.append(
                ExactSymmetricMatrix(dimension, denominator, packed_upper=numerators.copy())
            )
    if matrices[0].dimension != matrices[1].dimension:
        raise ValueError("matrix dimensions differ")
    return matrices[0], matrices[1]


def load_matrix_pair(path: str | Path) -> tuple[ExactSymmetricMatrix, ExactSymmetricMatrix]:
    path = Path(path)
    if path.suffix.lower() == ".npz":
        return load_matrix_pair_npz(path)
    return load_matrix_pair_json(path)
