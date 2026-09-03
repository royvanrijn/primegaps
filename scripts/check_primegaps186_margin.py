#!/usr/bin/env python3
"""Replay the exact printed margin and source-loss ledger of PrimeGaps186.

This checks rational bookkeeping copied from the numerical note. It does not
recompute the physical integrals and therefore does not discharge the Lean
`physical_integral_bounds` axiom.
"""

from __future__ import annotations

from fractions import Fraction
import json


F = Fraction

RHO_STAR = F(2_624_989, 10_000_000)
I_LOWER = F(23_685_317_816, 10**24)
I_UPPER = F(23_685_317_890, 10**24)
J_LOWER = F(90_248_755_123, 10**24)

GROUP_TOTALS = {
    "outer_order_2": 38_927_522,
    "outer_order_5_over_2": 622_829_241,
    "inner_old_order_2": 55_254,
    "inner_old_order_5_over_2": 435_544,
    "inner_new_order_2": 1_405_159,
    "inner_new_order_5_over_2": 32_422_390,
}
KIND_TOTALS = {
    "low": 79_845_205,
    "rank_two": 602_422_937,
    "high": 13_806_968,
}
TOP_COMPONENTS = [
    ("outer_order_5_over_2.P2", 72_504_766),
    ("outer_order_5_over_2.P3", 71_471_327),
    ("outer_order_5_over_2.P4", 64_233_828),
    ("outer_order_5_over_2.P1", 60_116_961),
    ("outer_order_5_over_2.P5", 54_915_393),
    ("outer_order_5_over_2.P9", 50_903_176),
    ("outer_order_5_over_2.P6", 45_639_918),
    ("outer_order_5_over_2.P7", 37_277_308),
]


def q(value: Fraction) -> dict[str, object]:
    return {
        "fraction": str(value),
        "decimal": f"{float(value):.18g}",
    }


def main() -> None:
    source_units = sum(GROUP_TOTALS.values())
    assert source_units == 696_075_110
    assert sum(KIND_TOTALS.values()) == source_units

    source_relative = F(source_units, 10**12)
    raw_cap_quotient = RHO_STAR * J_LOWER / I_UPPER
    final_quotient_using_printed_total = (
        RHO_STAR * (J_LOWER - I_LOWER * source_relative) / I_UPPER
    )
    final_margin_using_printed_total = final_quotient_using_printed_total - 1

    # Corollary 2.6 deliberately replaces the exact source total and raw cap
    # lower bound by simpler outward rational bounds.
    coarse_raw_cap = F(500_103, 500_000)
    coarse_source_relative = F(697, 10**6)
    certified_margin = coarse_raw_cap - 1 - RHO_STAR * coarse_source_relative
    assert certified_margin == F(230_382_667, 10**13)
    assert certified_margin > F(1, 50_000)

    # Tightening every source-loss estimate to zero cannot add more than this
    # to the final quotient while all other ingredients are frozen.
    maximum_gain_from_eliminating_all_printed_loss = RHO_STAR * source_relative

    result = {
        "schema": "primegaps.primegaps186-margin.v1",
        "source": "openai/PrimeGaps186@61340d0b74163003b32756bb16e91d9209a5e330",
        "scope": (
            "exact replay of printed rational endpoints and loss totals only; "
            "physical integrals are not recomputed"
        ),
        "rho_star": q(RHO_STAR),
        "cap_endpoints": {
            "I_lower": q(I_LOWER),
            "I_upper": q(I_UPPER),
            "J_lambda_lower": q(J_LOWER),
        },
        "raw_cap_quotient_lower_from_endpoints": q(raw_cap_quotient),
        "raw_cap_margin_lower_from_endpoints": q(raw_cap_quotient - 1),
        "source_loss": {
            "relative_units_1e12": source_units,
            "relative_bound": q(source_relative),
            "quotient_cost": q(RHO_STAR * source_relative),
            "group_totals": {
                name: {
                    "units": units,
                    "share": f"{units / source_units:.12%}",
                }
                for name, units in GROUP_TOTALS.items()
            },
            "kind_totals": {
                name: {
                    "units": units,
                    "share": f"{units / source_units:.12%}",
                }
                for name, units in KIND_TOTALS.items()
            },
            "largest_components": [
                {
                    "component": name,
                    "units": units,
                    "share": f"{units / source_units:.12%}",
                }
                for name, units in TOP_COMPONENTS
            ],
        },
        "final_from_printed_endpoints_and_exact_loss_total": {
            "quotient_lower": q(final_quotient_using_printed_total),
            "margin_lower": q(final_margin_using_printed_total),
        },
        "corollary_2_6_coarse_certificate": {
            "raw_cap_floor": q(coarse_raw_cap),
            "source_relative_ceiling": q(coarse_source_relative),
            "margin": q(certified_margin),
            "required_margin": q(F(1, 50_000)),
            "passes": certified_margin > F(1, 50_000),
        },
        "frozen-ingredient_headroom": {
            "maximum_quotient_gain_if_all_printed_source_loss_vanished": q(
                maximum_gain_from_eliminating_all_printed_loss
            ),
            "interpretation": (
                "source-bound tightening is valuable for robustness, but by "
                "itself can improve the quotient by less than 0.000183"
            ),
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
