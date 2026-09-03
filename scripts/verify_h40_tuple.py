#!/usr/bin/env python3
"""Exact finite replay of the PrimeGaps186 H40 tuple and its H39 prefix."""

from __future__ import annotations

import json


H40 = (
    0, 2, 6, 12, 20, 26, 30, 32, 36, 42,
    48, 50, 56, 60, 68, 72, 78, 86, 90, 92,
    98, 102, 110, 116, 120, 126, 132, 138, 140, 146,
    152, 156, 158, 162, 168, 170, 176, 180, 182, 186,
)
H39 = H40[:-1]


def primes_through(limit: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, limit + 1):
        if all(candidate % prime for prime in primes if prime * prime <= candidate):
            primes.append(candidate)
    return primes


def verify(name: str, values: tuple[int, ...], expected_diameter: int) -> dict[str, object]:
    assert values == tuple(sorted(set(values))), f"{name} is not strictly increasing"
    assert values[0] == 0
    assert values[-1] - values[0] == expected_diameter

    omissions: dict[str, list[int]] = {}
    checked_primes = primes_through(len(values))
    for prime in checked_primes:
        occupied = {value % prime for value in values}
        missing = sorted(set(range(prime)) - occupied)
        assert missing, f"{name} covers all residues modulo {prime}"
        omissions[str(prime)] = missing

    return {
        "name": name,
        "cardinality": len(values),
        "minimum": values[0],
        "maximum": values[-1],
        "diameter": expected_diameter,
        "admissible": True,
        "checked_primes": checked_primes,
        "omitted_residues": omissions,
        "automatic_for_primes_above": len(values),
        "tuple": list(values),
    }


def main() -> None:
    assert len(H40) == 40
    assert len(H39) == 39
    # H39 is a subset of H40, hence admissibility also follows immediately;
    # replay it independently so the endpoint claim is self-contained.
    result = {
        "schema": "primegaps.primegaps186-tuples.v1",
        "source": "openai/PrimeGaps186@61340d0b74163003b32756bb16e91d9209a5e330",
        "H40": verify("H40", H40, 186),
        "H39_prefix": verify("H39_prefix", H39, 182),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
