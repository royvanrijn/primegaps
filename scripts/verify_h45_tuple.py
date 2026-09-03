#!/usr/bin/env python3
"""Exact finite verification of the H45 tuple displayed in BGP212."""

from __future__ import annotations

import json


H45 = (
    0, 2, 12, 14, 24, 26, 30, 36, 44, 50, 54, 56, 60, 66, 72,
    74, 80, 84, 92, 96, 102, 110, 114, 116, 122, 126, 134, 140,
    144, 150, 156, 162, 164, 170, 176, 180, 182, 186, 192, 194,
    200, 204, 206, 210, 212,
)


def primes_through(limit: int) -> list[int]:
    out: list[int] = []
    for candidate in range(2, limit + 1):
        if all(candidate % prime for prime in out if prime * prime <= candidate):
            out.append(candidate)
    return out


def main() -> None:
    assert tuple(sorted(set(H45))) == H45
    assert len(H45) == 45
    diameter = H45[-1] - H45[0]
    assert diameter == 212

    omissions = {}
    checked_primes = primes_through(len(H45))
    for prime in checked_primes:
        occupied = {value % prime for value in H45}
        missing = sorted(set(range(prime)) - occupied)
        assert missing, f"tuple covers every residue modulo {prime}"
        omissions[str(prime)] = missing

    print(json.dumps({
        "schema": "primegaps.bgp212-h45.v1",
        "cardinality": len(H45),
        "minimum": H45[0],
        "maximum": H45[-1],
        "diameter": diameter,
        "admissible": True,
        "checked_primes": checked_primes,
        "omitted_residues": omissions,
        "tuple": list(H45),
        "automatic_above": len(H45),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
