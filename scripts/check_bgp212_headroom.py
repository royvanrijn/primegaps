#!/usr/bin/env python3
"""Exact arithmetic comparison of BGP212's active Type-II walls.

This is a structural headroom calculation, not a proof that the experimental
incomplete-rectangle estimate applies to BGP212's parameter family.
"""

from __future__ import annotations

from fractions import Fraction
import json


def payload(value: Fraction) -> dict[str, str]:
    return {"fraction": str(value), "decimal": f"{float(value):.15g}"}


def main() -> None:
    xi2 = Fraction(2, 5)
    delta = Fraction(41, 2500)
    chosen_a = Fraction(257, 1000)
    epsilon_a = Fraction(1, 10_000_000_000)

    # Proposition 7.8(II), second branch / the terminal Type-IIc wall used by
    # BGP212. This is the wall with approximately 1e-4 slack in delta.
    old_type_iic_a_ceiling = (
        xi2 / 4 + Fraction(11, 16) - delta - 2 * epsilon_a
    ) / 3

    # Proposition 7.8(II), first branch. It becomes the next global Type-II
    # wall if the terminal wall is replaced by the experimental estimate below.
    first_type_ii_a_ceiling = (xi2 + 8 - 10 * delta) / 32

    # Counterfactual only: the structural, unpadded version of the
    # incomplete-rectangle estimate currently checked in this repository at
    # delta=7/250. It must be re-proved for BGP212's delta and exact divisor
    # intervals before this ceiling can be used.
    experimental_type_iic_a_ceiling = (
        Fraction(11, 36) * xi2 + Fraction(2, 3) - delta
    ) / 3

    counterfactual_a_ceiling = min(
        first_type_ii_a_ceiling, experimental_type_iic_a_ceiling
    )

    def level(a: Fraction) -> Fraction:
        # BGP212 writes omega=A-1/4, and generated moduli reach
        # x^(1/2+2*omega)=x^(2A).
        return 2 * a

    result = {
        "schema": "primegaps.bgp212-headroom.v1",
        "scope": (
            "exact structural comparison only; the experimental Type-IIc "
            "estimate is not yet licensed at delta=41/2500"
        ),
        "paper_parameters": {
            "k": 45,
            "A1": payload(chosen_a),
            "delta": payload(delta),
            "xi2": payload(xi2),
            "epsilon_analytic": payload(epsilon_a),
            "generated_modulus_level": payload(level(chosen_a)),
        },
        "ceilings": {
            "published_terminal_type_iic": payload(old_type_iic_a_ceiling),
            "first_type_ii_cap_wall": payload(first_type_ii_a_ceiling),
            "experimental_terminal_type_iic_unpadded": payload(
                experimental_type_iic_a_ceiling
            ),
            "counterfactual_combined": payload(counterfactual_a_ceiling),
        },
        "headroom_from_paper_A1": {
            "published_terminal_type_iic": payload(
                old_type_iic_a_ceiling - chosen_a
            ),
            "counterfactual_combined": payload(
                counterfactual_a_ceiling - chosen_a
            ),
            "generated_modulus_level_gain": payload(
                level(counterfactual_a_ceiling) - level(chosen_a)
            ),
        },
        "next_binding_wall_under_counterfactual": (
            "Proposition 7.8(II), first branch"
        ),
        "tuple_targets": {
            "45": 212,
            "44": 210,
            "43": 200,
        },
    }

    assert chosen_a == Fraction(257, 1000)
    assert level(chosen_a) == Fraction(257, 500)
    assert first_type_ii_a_ceiling == Fraction(2059, 8000)
    assert experimental_type_iic_a_ceiling == Fraction(17381, 67500)
    assert counterfactual_a_ceiling == Fraction(2059, 8000)
    assert counterfactual_a_ceiling - chosen_a == Fraction(3, 8000)
    assert level(counterfactual_a_ceiling) == Fraction(2059, 4000)

    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
