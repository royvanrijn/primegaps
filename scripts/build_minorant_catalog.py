#!/usr/bin/env python3
"""Emit the optimizer-facing prime-minorant candidate catalog as JSON."""

from __future__ import annotations

import argparse
import json

from primegaps.minorants import (
    HarmanRegimes,
    baker_irving_base_loss,
    baker_irving_parameters,
    discard_variants,
    regime_frontier,
    stadlmann_loss_components,
    stadlmann_xi1_interval,
)


def regimes_json(regimes: HarmanRegimes) -> dict[str, object]:
    return {
        "xi1": regimes.xi1,
        "xi2": regimes.xi2,
        "xi3": regimes.xi3,
        "type_I": {
            "shape": "alpha * beta_smooth",
            "smooth_factor_gamma": list(regimes.type_i_smooth_gamma),
        },
        "type_II": {
            "shape": "alpha_SW * beta_SW",
            "factor_gamma": list(regimes.type_ii_gamma),
        },
        "type_III": {
            "shape": "alpha * psi1_smooth * psi2_smooth * psi3_smooth",
            "each_gamma": list(regimes.type_iii_each_gamma),
            "each_pair_sum_min": regimes.type_iii_pair_sum_min,
        },
    }


def variant_json(variant, delta: float) -> dict[str, object]:
    return {
        "decomposition": variant.name,
        "retained_mass": variant.retained_mass,
        "pointwise_negative_bound_c2": variant.c2,
        "extra_required_regimes": list(variant.extra_required_regimes),
        "quadrature_change_vs_comparison_order": delta,
        "mass_is_rigorous_enclosure": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=160)
    parser.add_argument("--comparison-order", type=int, default=96)
    args = parser.parse_args()

    catalog: dict[str, object] = {
        "schema": "primegaps.prime-minorant-catalog.v1",
        "scope": "prime minorants and required distribution regimes only",
        "quadrature": {
            "method": "analytic innermost integral plus nested Gauss-Legendre",
            "order": args.order,
            "comparison_order": args.comparison_order,
            "rigorous_enclosure": False,
        },
        "sources": [
            {
                "id": "stadlmann-2026v1",
                "url": "https://arxiv.org/abs/2608.31126v1",
                "locators": ["Definition 5", "Proposition 2", "eqs. (21),(22)"],
            },
            {
                "id": "baker-irving-2017",
                "url": "https://arxiv.org/abs/1505.01815v1",
                "locators": ["Lemma 2", "Lemmas 4--6", "Sections 3--5"],
            },
            {
                "id": "stadlmann-2025",
                "url": "https://arxiv.org/abs/2309.00425",
                "locators": ["Section 4", "Proposition 1"],
            },
        ],
        "modulus_families": {
            "Qstar-current": {
                "description": (
                    "Every squarefree q in the support-supplied Q*: q=e*e'*prod(f_i)*prod(f'_i), "
                    "e*e' x^delta-smooth, large factors >=x^delta, with the supplied A/B log-size bounds"
                ),
                "support_compatibility": "max_j B[j,1] < 1-2*xi2",
            },
            "BI-smooth": {
                "description": "squarefree x^delta-smooth q <= x^(theta-epsilon)",
            },
        },
        "exceptional_regimes": {
            "exception-A": {
                "shape": "five-prime ordered convolution",
                "polytope": [
                    "beta<a4<a3<a2<a1<xi2",
                    "a1+a2<xi2",
                    "a2+a3+a4>1-xi2",
                    "a5>a4",
                ],
                "pointwise_multiplicity_bound": 4,
            },
            "exception-B": {
                "shape": "five-prime reversal convolution",
                "polytope": [
                    "beta<a2,a3,a4,a5,a6<8*xi2-3",
                    "a2>a3 and a4>a3",
                    "a2+a4<xi2",
                    "a3+a4+a5>1-xi2",
                    "a6>a5",
                ],
                "pointwise_multiplicity_bound": 20,
            },
        },
        "families": [],
    }

    families = catalog["families"]
    assert isinstance(families, list)

    published = HarmanRegimes(0.38, 0.4, 0.4)
    families.append(
        {
            "family": "stadlmann-direct-primes",
            "status": "published anchor; exceptional regions empty",
            "modulus_family": "Qstar-current",
            "base_required_regimes": regimes_json(published),
            "sifted_support_exponent_beta": 0.2,
            "candidates": [
                {
                    "decomposition": "direct-prime-indicator",
                    "retained_mass": 1.0,
                    "pointwise_negative_bound_c2": 0,
                    "extra_required_regimes": [],
                    "mass_is_exact": True,
                }
            ],
        }
    )

    for xi2 in (0.4005, 0.4025, 0.40481, 0.4075, 0.409, 0.4105, 0.4115):
        lower, upper = stadlmann_xi1_interval(xi2)
        slack = 0.1 * (upper - lower)
        regimes = regime_frontier(xi2, slack)
        loss_a, loss_b = stadlmann_loss_components(xi2, args.order)
        old_a, old_b = stadlmann_loss_components(xi2, args.comparison_order)
        delta = abs((loss_a + loss_b) - (old_a + old_b))
        families.append(
            {
                "family": "stadlmann-xi",
                "status": "baseline published conditionally; partial-retention variants are derived conditional candidates",
                "modulus_family": "Qstar-current",
                "xi1_admissible_open_interval": [lower, upper],
                "xi1_frontier_note": "larger xi1 shrinks Type I; witness uses 10% of interval width as strict upper slack",
                "base_required_regimes": regimes_json(regimes),
                "sifted_support_exponent_beta": 1.0 - 2.0 * xi2,
                "loss_integrals": {"exception-A": loss_a, "exception-B": loss_b},
                "candidates": [
                    variant_json(variant, delta)
                    for variant in discard_variants(loss_a, loss_b)
                ],
            }
        )

    for eta in (0.001, 0.003, 0.00481, 0.006, 0.00667):
        regimes, beta, theta = baker_irving_parameters(eta)
        base = baker_irving_base_loss(eta, args.order)
        old_base = baker_irving_base_loss(eta, args.comparison_order)
        # Baker--Irving's two exceptional multiplicities have masses I and 5I.
        delta = 6.0 * abs(base - old_base)
        families.append(
            {
                "family": "baker-irving-eta",
                "status": "baseline published for smooth moduli; partial-retention variants are derived conditional candidates",
                "eta": eta,
                "modulus_family": "BI-smooth",
                "distribution_exponent_theta": theta,
                "base_required_regimes": regimes_json(regimes),
                "sifted_support_exponent_beta": beta,
                "loss_integrals": {"exception-A": base, "exception-B": 5.0 * base},
                "candidates": [
                    variant_json(variant, delta)
                    for variant in discard_variants(base, 5.0 * base)
                ],
            }
        )

    print(json.dumps(catalog, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
