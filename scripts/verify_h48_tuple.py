#!/usr/bin/env python3
"""Exact finite replay for an admissible 48-tuple of diameter 236."""

from __future__ import annotations

import json


H48 = (
    0, 6, 8, 14, 18, 24, 26, 48, 50, 54, 56, 60,
    66, 68, 74, 78, 80, 84, 90, 96, 98, 104, 110, 116,
    120, 126, 134, 138, 144, 150, 158, 164, 168, 176, 180, 186,
    188, 194, 200, 204, 206, 210, 216, 224, 228, 230, 234, 236,
)


def primes_through(n: int) -> tuple[int, ...]:
    primes: list[int] = []
    for candidate in range(2, n + 1):
        if all(candidate % p for p in primes if p * p <= candidate):
            primes.append(candidate)
    return tuple(primes)


def main() -> None:
    if len(H48) != 48 or len(set(H48)) != 48:
        raise SystemExit("H48 does not contain 48 distinct elements")
    if tuple(sorted(H48)) != H48:
        raise SystemExit("H48 is not strictly increasing")
    if H48[0] != 0 or H48[-1] - H48[0] != 236:
        raise SystemExit("H48 does not have normalized diameter 236")

    rows = []
    for prime in primes_through(48):
        occupied = {value % prime for value in H48}
        omitted = sorted(set(range(prime)) - occupied)
        if not omitted:
            raise SystemExit(f"H48 occupies every residue modulo {prime}")
        rows.append(
            {
                "prime": prime,
                "occupied_residue_count": len(occupied),
                "omitted_residues": omitted,
            }
        )

    print(
        json.dumps(
            {
                "schema": "primegaps-h48-admissibility-v1",
                "tuple": list(H48),
                "cardinality": len(H48),
                "minimum": H48[0],
                "maximum": H48[-1],
                "diameter": H48[-1] - H48[0],
                "admissible": True,
                "checked_primes": rows,
                "large_prime_argument": (
                    "for p > 48, a 48-element set cannot occupy all p residue classes"
                ),
                "external_diameter_reference": "OEIS A008407: a(48)=236",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
