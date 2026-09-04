#!/usr/bin/env python3
"""Map roughness, cancellation, and factorization geometry across beta."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import tempfile

from primegaps.parity import (
    parity_contributions,
    parity_error_budget,
    rough_factor_constants,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORE_RECORD = ROOT / "experiments/physical_parity_viability.json"
DEFAULT_PHYSICAL_INPUT = ROOT / "reproduction/186/physical-parity-input.json"
DEFAULT_BETAS = (
    0.250001,
    0.275,
    0.3,
    0.325,
    1.0 / 3.0 - 1e-6,
    1.0 / 3.0 + 1e-6,
    0.35,
    0.4,
    0.45,
    0.49,
)


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def portable_path(path: Path) -> str:
    """Prefer repository-relative provenance paths when possible."""

    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def _load_scores(path: Path) -> dict[int, float]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.physical-parity-production-extrapolation.v1":
        raise ValueError("unexpected physical parity score schema")
    return {
        int(row["dimension"]): float(
            row["quadratic_in_inverse_mesh"]["production_mesh_score"]
        )
        for row in payload["dimensions"]
    }


def _load_geometry(path: Path) -> dict[str, float | str]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.physical-parity-input.v1":
        raise ValueError("unexpected physical parity input schema")
    rho = Fraction(payload["rho_star"])
    outer_radius = Fraction(payload["outer_radius"])
    largest_target_fragment = max(
        Fraction(fragment) for _, fragment in payload["cap_shell_data"]["full"]
    )
    pair_support = 2 * rho * outer_radius
    target_fragment = rho * largest_target_fragment
    return {
        "rho_star": str(rho),
        "outer_radius": str(outer_radius),
        "crude_pair_support_exponent": float(pair_support),
        "largest_target_fragment_exponent": float(target_fragment),
    }


def _identity_record(degree: int) -> dict[str, object]:
    if degree == 3:
        return {
            "factorial": "A1 - 2*A2 + 3*A3",
            "liouville": "(A0 - L)/2 - A3",
            "coefficients": {"A1": 1, "A2": -2, "A3": 3},
        }
    return {
        "factorial": "A1 - 2*A2",
        "liouville": "(A0 - L)/2",
        "coefficients": {"A1": 1, "A2": -2},
    }


def sweep(
    score_record: Path,
    physical_input: Path,
    betas: tuple[float, ...],
) -> dict[str, object]:
    scores = _load_scores(score_record)
    geometry = _load_geometry(physical_input)
    if not betas:
        raise ValueError("at least one beta is required")
    if tuple(sorted(betas)) != betas or len(set(betas)) != len(betas):
        raise ValueError("betas must be strictly increasing")

    baseline = rough_factor_constants(betas[0]).gross_signed_condition
    rows = []
    for beta in betas:
        constants = rough_factor_constants(beta)
        if abs(constants.signed_identity - 1.0) > 5e-13:
            raise ArithmeticError(f"signed constants fail to reconstruct one at beta={beta}")
        degree = constants.detector_degree
        dimension_budgets = {}
        for dimension, score in sorted(scores.items(), reverse=True):
            contributions = parity_contributions(score, constants)
            reconstructed_eta, _ = parity_error_budget(contributions)
            eta = score - 1.0
            if abs(reconstructed_eta - eta) > 5e-13:
                raise ArithmeticError(
                    f"signed score reconstruction failed at beta={beta}, k={dimension}"
                )
            relative = eta / contributions.gross_absolute
            dimension_budgets[str(dimension)] = {
                "ideal_score": score,
                "absolute_error_budget_I": eta,
                "common_relative_error_budget": relative,
            }

        triprime_geometry = None
        if degree == 3:
            triprime_geometry = {
                "single_factor_exponent_range": [beta, 1.0 - 2.0 * beta],
                "complementary_semiprime_exponent_range": [
                    2.0 * beta,
                    1.0 - beta,
                ],
                "simplex_slack": 1.0 - 3.0 * beta,
            }
        pair_support = float(geometry["crude_pair_support_exponent"])
        target_fragment = float(geometry["largest_target_fragment_exponent"])
        rows.append(
            {
                "beta": beta,
                "detector_degree": degree,
                "identity": _identity_record(degree),
                "constants": {
                    **asdict(constants),
                    "rough_carrier": constants.rough_carrier,
                    "liouville_mean": (
                        -constants.prime
                        + constants.semiprime
                        - constants.triprime
                    ),
                },
                "gross_condition_relative_to_first_beta": (
                    constants.gross_signed_condition / baseline
                ),
                "bilinear_geometry": {
                    "semiprime_smaller_factor_exponent_range": [beta, 0.5],
                    "semiprime_strip_width": 0.5 - beta,
                    "triprime_prime_times_semiprime": triprime_geometry,
                },
                "roughness_geometry": {
                    "beta_in_rho_star_units": beta / float(Fraction(geometry["rho_star"])),
                    "headroom_above_largest_target_fragment_exponent": (
                        beta - target_fragment
                    ),
                    "crude_pair_support_exponent_in_beta_units": pair_support / beta,
                },
                "dimension_budgets": dimension_budgets,
            }
        )

    return {
        "schema": "primegaps.physical-parity-beta-sweep.v1",
        "status": "checked-elementary-tradeoff-on-exploratory-physical-scores",
        "score_independence": (
            "The full-face physical scores are held fixed; only the exact signed "
            "rough-almost-prime decomposition and its error amplification vary."
        ),
        "physical_geometry": {
            **geometry,
            "warning": (
                "The pair-support/beta ratio is a scale diagnostic, not a proved "
                "level of distribution for the rough sieve."
            ),
        },
        "inputs": {
            "score_record": portable_path(score_record),
            "score_record_sha256": file_hash(score_record),
            "physical_input": portable_path(physical_input),
            "physical_input_sha256": file_hash(physical_input),
        },
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-record", type=Path, default=DEFAULT_SCORE_RECORD)
    parser.add_argument("--physical-input", type=Path, default=DEFAULT_PHYSICAL_INPUT)
    parser.add_argument("--betas", type=float, nargs="+", default=DEFAULT_BETAS)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = sweep(args.score_record, args.physical_input, tuple(args.betas))
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        atomic_json(args.output, payload)
        print(rendered, end="")


if __name__ == "__main__":
    main()
