#!/usr/bin/env python3
"""Optimize admissible tuple cardinality up to a proposed diameter.

This is the finite endpoint check used in docs/shi-2025-audit.md.  It requires
NumPy and SciPy's HiGHS-backed ``milp`` implementation.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix


def primes_through(n: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, n + 1):
        if all(candidate % prime for prime in primes if prime * prime <= candidate):
            primes.append(candidate)
    return primes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diameter", type=int, required=True)
    parser.add_argument("--time-limit", type=float, default=300.0)
    args = parser.parse_args()
    if args.diameter < 0:
        raise SystemExit("diameter must be nonnegative")

    # Every admissible tuple omits one residue modulo 2, hence has one parity.
    # Translate its minimum to zero and search only the even representatives.
    points = list(range(0, args.diameter + 1, 2))
    primes = primes_through(len(points))[1:]
    x_index = {number: index for index, number in enumerate(points)}
    variable_count = len(points)
    y_index: dict[tuple[int, int], int] = {}
    for prime in primes:
        for residue in range(prime):
            y_index[prime, residue] = variable_count
            variable_count += 1

    # At least one y[p,r] is selected for every p.  If y[p,r]=1, no selected
    # integer may occupy residue r modulo p.
    rows: list[tuple[dict[int, float], float, float]] = []
    for prime in primes:
        rows.append(
            ({y_index[prime, residue]: 1.0 for residue in range(prime)}, 1.0, math.inf)
        )
        for number in points:
            rows.append(
                (
                    {
                        x_index[number]: 1.0,
                        y_index[prime, number % prime]: 1.0,
                    },
                    -math.inf,
                    1.0,
                )
            )

    matrix = lil_matrix((len(rows), variable_count), dtype=float)
    lower = np.empty(len(rows))
    upper = np.empty(len(rows))
    for row_number, (coefficients, lo, hi) in enumerate(rows):
        for column, value in coefficients.items():
            matrix[row_number, column] = value
        lower[row_number] = lo
        upper[row_number] = hi

    objective = np.zeros(variable_count)
    objective[: len(points)] = -1.0
    result = milp(
        objective,
        integrality=np.ones(variable_count),
        bounds=Bounds(np.zeros(variable_count), np.ones(variable_count)),
        constraints=LinearConstraint(matrix.tocsr(), lower, upper),
        options={"time_limit": args.time_limit, "mip_rel_gap": 0.0, "presolve": True},
    )
    if not result.success or result.mip_gap != 0.0 or result.x is None:
        raise SystemExit(f"solver did not certify optimality: {result.message}")

    selected = [number for number in points if result.x[x_index[number]] > 0.5]
    if not all(len({number % prime for number in selected}) < prime for prime in primes):
        raise SystemExit("returned set failed exact admissibility replay")
    if round(-result.mip_dual_bound) != len(selected):
        raise SystemExit("integral primal and dual cardinality bounds disagree")

    print(
        json.dumps(
            {
                "diameter_limit": args.diameter,
                "maximum_cardinality": len(selected),
                "maximizing_set": selected,
                "mip_dual_cardinality_bound": -float(result.mip_dual_bound),
                "mip_gap": float(result.mip_gap),
                "mip_node_count": int(result.mip_node_count),
                "solver_message": result.message,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
