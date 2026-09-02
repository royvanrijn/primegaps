#!/usr/bin/env python3
"""Run the rigorous cheap gate for surgical deletion of Type-IIc branches."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction

from primegaps.minorants import (
    optimistic_no_k_screen,
    type_iic_gamma_cutoff,
    type_iic_middle_high_loss_enclosure,
)


def fraction_json(value: Fraction) -> dict[str, str]:
    return {"fraction": str(value), "decimal": f"{float(value):.17g}"}


def bounded_fraction_json(
    value: Fraction, *, direction: str, places: int = 18
) -> dict[str, object]:
    """Serialize a huge exact fraction by an outward decimal and fingerprints."""
    if value < 0 or direction not in {"lower", "upper"}:
        raise ValueError("only nonnegative lower/upper bounds are supported")
    scale = 10**places
    quotient, remainder = divmod(value.numerator * scale, value.denominator)
    if direction == "upper" and remainder:
        quotient += 1
    whole, fractional = divmod(quotient, scale)

    def digest(integer: int) -> str:
        data = integer.to_bytes((integer.bit_length() + 7) // 8 or 1, "big")
        return hashlib.sha256(data).hexdigest()

    return {
        "decimal_outward": f"{whole}.{fractional:0{places}d}",
        "direction": direction,
        "exact_numerator_sha256": digest(value.numerator),
        "exact_denominator_sha256": digest(value.denominator),
        "numerator_bits": value.numerator.bit_length(),
        "denominator_bits": value.denominator.bit_length(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--support-max", default="913600001/3600000000")
    parser.add_argument("--delta", default="7/250")
    parser.add_argument("--epsilon", default="1/10000000000")
    parser.add_argument("--raw-score-cap", default="1.000670218")
    parser.add_argument("--parts", type=int, default=128)
    parser.add_argument("--log-terms", type=int, default=12)
    args = parser.parse_args()

    support_max = Fraction(args.support_max)
    delta = Fraction(args.delta)
    epsilon = Fraction(args.epsilon)
    raw_score_cap = Fraction(args.raw_score_cap)
    cutoff = type_iic_gamma_cutoff(support_max, delta, epsilon)
    loss = type_iic_middle_high_loss_enclosure(
        cutoff, parts=args.parts, log_terms=args.log_terms
    )
    screen = optimistic_no_k_screen(loss, raw_score_cap)

    payload = {
        "schema": "primegaps.surgical-type-iic-screen.v1",
        "scope": (
            "sign-preserving branch deletions of the current Buchstab identity; "
            "the raw score cap is caller-supplied"
        ),
        "inputs": {
            "support_max": fraction_json(support_max),
            "delta": fraction_json(delta),
            "theorem_epsilon": fraction_json(epsilon),
            "raw_no_k_score_cap": fraction_json(raw_score_cap),
            "parts": args.parts,
            "log_terms": args.log_terms,
        },
        "branch_trace": [
            {
                "id": "negative-one-prime-direct-type-II",
                "sign": -1,
                "source": "Stadlmann 2025 TeX line 1837, equation (buchstab1)",
                "selection_status": "cannot delete while preserving a pointwise minorant",
            },
            {
                "id": "positive-two-prime-middle-type-II-high-sum",
                "sign": 1,
                "source": "Stadlmann 2025 TeX line 1838, equation (buchstab2)",
                "selection_status": "deletable but mandatory subset already fails mass gate",
            },
        ],
        "variants": [
            {
                "name": "keep-current-core-branches",
                "pointwise_minorant": True,
                "bypasses_type_iic": False,
            },
            {
                "name": "delete-positive-high-sum-type-iic-slice",
                "pointwise_minorant": True,
                "bypasses_type_iic": False,
                "reason": "negative direct Type-II semiprime correction remains",
            },
            {
                "name": "delete-both-direct-type-iic-branches",
                "pointwise_minorant": False,
                "bypasses_type_iic": True,
                "reason": "deleting the negative branch is positive on rough semiprimes",
            },
            {
                "name": "hypothetical-free-negative-branch-reorganization",
                "pointwise_minorant": "conditional",
                "bypasses_type_iic": "conditional",
                "screen_survives": screen.survives,
                "reason": "optimistically charges only the mandatory positive slice",
            },
        ],
        "gamma_cutoff": fraction_json(cutoff),
        "mandatory_loss_enclosure": {
            "lower": bounded_fraction_json(loss.lower, direction="lower"),
            "upper": bounded_fraction_json(loss.upper, direction="upper"),
            "rigorous": True,
            "method": "rational atanh log bounds plus monotone interval rectangles",
        },
        "retained_mass_upper": bounded_fraction_json(
            screen.retained_mass_upper, direction="upper"
        ),
        "optimistic_no_k_score_upper": bounded_fraction_json(
            screen.optimistic_score_upper, direction="upper"
        ),
        "required_raw_score_lower": bounded_fraction_json(
            screen.required_raw_score_lower, direction="lower"
        ),
        "survives": screen.survives,
        "run_full_ijk": screen.survives,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
