#!/usr/bin/env python3
"""Replay the exact parameter audit and optionally expand proposed B_F blocks."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARITY_RESULT = ROOT / "experiments/physical_parity_viability.json"
PHYSICAL_INPUT = ROOT / "reproduction/186/physical-parity-input.json"


def _fraction(value: str) -> Fraction:
    return Fraction(value)


def audit() -> dict[str, object]:
    parity = json.loads(PARITY_RESULT.read_text())
    physical = json.loads(PHYSICAL_INPUT.read_text())
    beta = Fraction(str(parity["beta"]))
    rho = Fraction(physical["rho_star"])
    outer_radius = Fraction(physical["outer_radius"])
    largest_fragment = Fraction(physical["cap_shell_data"]["full"][0][1])
    pair_exponent = 2 * rho * outer_radius
    modulus_prime_exponent = rho * largest_fragment
    return {
        "schema": "primegaps.bf-dyadic-audit.v1",
        "status": "refuted-proposed-bf",
        "parameters": {
            "beta": str(beta),
            "pair_modulus_exponent": str(pair_exponent),
            "pair_modulus_exponent_decimal": float(pair_exponent),
            "modulus_prime_exponent": str(modulus_prime_exponent),
            "modulus_prime_exponent_decimal": float(modulus_prime_exponent),
        },
        "exact_checks": {
            "beta_gt_one_quarter": beta > Fraction(1, 4),
            "fi_v_upper_below_two_rough_primes": Fraction(1, 2) < 2 * beta,
            "physical_modulus_primes_below_rough_cutoff": modulus_prime_exponent < beta,
            "gamma_at_C_equals_one": "gamma(v,1)=1",
            "active_inner_sign": "mu(uv)=-mu(u) when mu(uv)!=0",
        },
        "compressed_block_family": {
            "V": "V_j=2^j*sqrt(D)/Delta, 0<=j<J, V_j<sqrt(X)/delta",
            "U": "X/(2*V_j)<u<=2*X/V_j (split into at most three dyadic U intervals)",
            "q": "q=1 or 2^ell<q<=2^(ell+1), q<=X^(2742997/5000000)",
            "coefficient_norm": {
                "gamma_C_1_linf": 1,
                "gamma_C_1_l2_squared": "number of admissible primes v in the block",
                "physical_q_norm": "undefined until the discrete Lambda_d normalization and residue-coloured aggregation are specified",
            },
            "required_saving": "A_i(X)*(log X)^(-222) for every V and every 1<=C<=X/D",
        },
        "classification": {
            "already_controlled": "blocks with 2V<=X^beta are empty; blocks with C>=2V have gamma=0",
            "potential_dispersion_or_trace": [],
            "genuinely_new_parity": "not a B_F block: a replacement must retain cancellation between factor-count sectors before the outer absolute value",
            "impossible": "every nonempty prime-only V block at C=1 under the proposed rough sequence",
        },
    }


def enumerate_blocks(
    *,
    log2_x: int,
    d_exponent: Fraction,
    capital_delta_exponent: Fraction,
    log2_little_delta: int,
) -> list[dict[str, object]]:
    """Expand the Cartesian (V,q) family for one declared FI specialization."""

    base = audit()
    beta = Fraction(base["parameters"]["beta"])
    theta = Fraction(base["parameters"]["pair_modulus_exponent"])
    x_log = Fraction(log2_x)
    y_log = d_exponent * x_log / 2 - capital_delta_exponent * x_log
    upper_v_log = x_log / 2 - log2_little_delta
    z_log = beta * x_log
    q_top_log = theta * x_log
    v_logs: list[tuple[int, Fraction]] = []
    index = 0
    while y_log + index < upper_v_log:
        v_logs.append((index, y_log + index))
        index += 1
    q_logs: list[tuple[int, Fraction, Fraction]] = [(-1, Fraction(0), Fraction(0))]
    index = 0
    while Fraction(index) < q_top_log:
        q_logs.append((index, Fraction(index), min(Fraction(index + 1), q_top_log)))
        index += 1

    rows: list[dict[str, object]] = []
    for j, v_log in v_logs:
        if v_log + 1 <= z_log:
            support = "empty"
            classification = "already-controlled"
        elif v_log + 1 < 2 * z_log:
            support = "prime-only" if v_log >= z_log else "rough-boundary-prime-only"
            classification = "impossible"
        else:
            support = "composites-not-excluded"
            classification = "outside-audited-prime-only-range"
        for ell, q_low, q_high in q_logs:
            rows.append(
                {
                    "V_index": j,
                    "V_log2_range": [str(v_log), str(v_log + 1)],
                    "U_log2_range": [str(x_log - v_log - 1), str(x_log + 1 - v_log)],
                    "q_index": ell,
                    "q_log2_range": [str(q_low), str(q_high)],
                    "support": support,
                    "gamma_C_1_linf": 1,
                    "gamma_C_1_l2_squared": "number of admissible primes v",
                    "physical_q_norm": "undefined",
                    "required_saving": "A_i(X)*(log X)^(-222)",
                    "parent_V_classification": classification,
                    "q_slice_classification": "undefined-without-physical-coefficient-norm",
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log2-x", type=int)
    parser.add_argument("--d-exponent", type=_fraction)
    parser.add_argument("--capital-delta-exponent", type=_fraction)
    parser.add_argument("--log2-little-delta", type=int)
    parser.add_argument("--expand", action="store_true")
    args = parser.parse_args()
    result = audit()
    if args.expand:
        supplied = (
            args.log2_x,
            args.d_exponent,
            args.capital_delta_exponent,
            args.log2_little_delta,
        )
        if any(value is None for value in supplied):
            parser.error("--expand requires all four range parameters")
        result["declared_specialization"] = {
            "log2_x": args.log2_x,
            "D_exponent": str(args.d_exponent),
            "Delta_exponent": str(args.capital_delta_exponent),
            "log2_delta": args.log2_little_delta,
        }
        result["blocks"] = enumerate_blocks(
            log2_x=args.log2_x,
            d_exponent=args.d_exponent,
            capital_delta_exponent=args.capital_delta_exponent,
            log2_little_delta=args.log2_little_delta,
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
