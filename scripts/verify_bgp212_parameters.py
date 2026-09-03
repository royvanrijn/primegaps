#!/usr/bin/env python3
"""Replay the exact scalar parameter checks exposed by the BGP212 draft.

This verifies the published rational datum, support/Harman/global analytic
walls, and the explicit H45 tuple.  It deliberately does not claim to replay
the 455 continuum-packing certificates or the unpublished variational vector.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


def F(value: str | int) -> Fraction:
    return Fraction(value)


def primes_through(limit: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, limit + 1):
        if all(candidate % p for p in primes if p * p <= candidate):
            primes.append(candidate)
    return primes


def rational(value: Fraction) -> dict[str, str]:
    return {"fraction": str(value), "decimal": f"{float(value):.16g}"}


def check(name: str, left: Fraction, right: Fraction, *, strict: bool = False) -> dict:
    holds = left < right if strict else left <= right
    if not holds:
        relation = "<" if strict else "<="
        raise AssertionError(f"{name}: {left} not {relation} {right}")
    return {
        "name": name,
        "left": rational(left),
        "right": rational(right),
        "slack": rational(right - left),
        "strict": strict,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "manifest",
        nargs="?",
        type=Path,
        default=Path("reproduction/212/paper-parameters.json"),
    )
    args = parser.parse_args()
    data = json.loads(args.manifest.read_text())
    support = data["physical_support"]

    k = int(data["theorem"]["k"])
    omega = F(support["omega"])
    a0, a1 = map(F, support["A"])
    eps_s = F(support["support_epsilon"])
    delta = F(support["delta"])
    xi1, xi2, xi3 = map(F, support["xi"])
    eps_a = F(support["analytic_epsilon"])

    assert a0 == -eps_s
    assert omega == a1 - F(1) / 4
    assert F(support["total_cap"]) == a1 + eps_s
    assert F(support["marginal_cap"]) == a1 - eps_s
    assert F(data["rescaled_support"]["L"]) == 4 * (a1 + eps_s)
    assert F(data["rescaled_support"]["tau"]) == 4 * (a1 - eps_s)
    assert F(data["rescaled_support"]["g"]) == 4 * delta

    cap_data = support["rough_caps"]
    caps = [
        F(cap_data.get(str(index), cap_data["tail"]))
        for index in range(1, int(F(1) / delta) + 1)
    ]
    checks: list[dict] = []
    checks.append(check("support delta positive", F(0), delta, strict=True))
    checks.append(check("support A1 below half minus epsilon", a1, F(1) / 2 - eps_s, strict=True))
    for index, cap in enumerate(caps, 1):
        checks.append(check(f"B[{index}] above delta", delta, cap, strict=True))
        if index > 1:
            checks.append(check(f"B[{index - 1}] monotone", caps[index - 2], cap))
            checks.append(check(f"B[{index}] hereditary step", cap, caps[index - 2] + delta))

    # Harman/direct-prime conditions from Proposition 7.2.
    checks.extend([
        check("2 xi1 + 3 xi2 < 2", 2 * xi1 + 3 * xi2, F(2), strict=True),
        check("xi2 <= xi3", xi2, xi3),
        check("xi1 + 9 xi2 < 4", xi1 + 9 * xi2, F(4), strict=True),
        check("2 xi1 + xi2 > 1", F(1), 2 * xi1 + xi2, strict=True),
        check("17 xi2 < 7", 17 * xi2, F(7), strict=True),
        check("direct-prime endpoint xi2 <= 2/5", xi2, F(2) / 5),
    ])

    # Global walls from Proposition 7.8, with the fixed analytic retreat.
    type_i_1 = xi1 - 4 * a1 + F(2) / 3 - 2 * eps_a
    type_i_2 = F(9) / 7 - F(34) * a1 / 7 - 2 * eps_a
    type_ii_range = F(19) / 2 - 36 * a1 - 13 * delta - 9 * eps_a
    type_ii_1 = xi2 / 10 - F(32) * a1 / 10 + F(8) / 10 - 2 * eps_a
    type_ii_2 = xi2 / 4 + F(11) / 16 - 3 * a1 - 2 * eps_a
    type_iii = F(11) / 8 - F(7) * a1 / 2 - F(9) * xi3 / 8 - 2 * eps_a
    checks.extend([
        check("Type I wall 1", delta, type_i_1, strict=True),
        check("Type I wall 2", delta, type_i_2, strict=True),
        check("Type II range wall", F(0), type_ii_range),
        check("Type II cap wall 1", delta, type_ii_1),
        check("Type II cap wall 2", delta, type_ii_2),
        check("Type III wall", delta, type_iii, strict=True),
        check("dominant Type II ledger", 3 * a1 + delta + 2 * eps_a, F(63) / 80),
    ])

    h45 = tuple(int(value) for value in data["h45"])
    assert len(h45) == k
    assert tuple(sorted(set(h45))) == h45
    assert h45[-1] - h45[0] == int(data["theorem"]["tuple_diameter"])
    omissions = {}
    for prime in primes_through(k):
        missing = sorted(set(range(prime)) - {value % prime for value in h45})
        if not missing:
            raise AssertionError(f"H45 is not admissible modulo {prime}")
        omissions[str(prime)] = missing

    output = {
        "schema": "primegaps.bgp212.parameter-replay.v1",
        "source_pdf_sha256": data["source"]["pdf_sha256"],
        "checks_passed": len(checks),
        "checks": checks,
        "h45": {
            "cardinality": len(h45),
            "diameter": h45[-1] - h45[0],
            "admissible": True,
            "omitted_residues": omissions,
        },
        "not_replayed": data["unreleased_replay_inputs"],
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
