#!/usr/bin/env python3
"""Cheap replay of the exact H_{5/2} factorization-oracle comparison."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from primegaps.physical import (
    exact_group_failure,
    grouped_five_halves_failure,
    primegaps186_source_data,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECORD = ROOT / "experiments" / "physical_factorization_oracle.json"


def _file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "primegaps.physical-factorization-oracle.v1":
        raise ValueError("unexpected factorization-oracle schema")
    implementation = ROOT / "experiments" / "physical_factorization_oracle.py"
    if payload["implementation_sha256"] != _file_hash(implementation):
        raise ArithmeticError("factorization-oracle implementation hash changed")
    predicate = ROOT / payload["exact_inputs"]["predicate_module"]
    if payload["exact_inputs"]["predicate_module_sha256"] != _file_hash(predicate):
        raise ArithmeticError("factorization predicate module hash changed")
    old, new, groups, mesh = primegaps186_source_data()
    if (len(old), len(new)) != (
        payload["exact_inputs"]["old_source_rows"],
        payload["exact_inputs"]["new_source_rows"],
    ) or mesh != Fraction(payload["exact_inputs"]["mesh"]):
        raise ArithmeticError("exact source inputs changed")
    by_name = {group.name: group for group in groups}
    checked = {}
    for row in payload["groups"]:
        group = by_name[row["group"]]
        parameters = row["parameters"]
        for field in (
            "activation",
            "threshold",
            "radial_lower",
            "radial_upper",
            "hard_cap",
            "split",
        ):
            if getattr(group, field) != Fraction(parameters[field]):
                raise ArithmeticError(f"{group.name} parameter changed: {field}")
        witness = row["strict_overcover_witness"]
        fragments = tuple(Fraction(value) for value in witness["fragments"])
        total_mass = Fraction(witness["total_mass"])
        exact = exact_group_failure(group, fragments, total_mass)
        h25 = grouped_five_halves_failure(group, fragments, total_mass)
        if exact or not h25:
            raise ArithmeticError(f"{group.name} overcover witness did not replay")
        false_negatives = (
            row["critical_census"]["false_negatives"]
            + row["uniform_cell_stress_census"]["false_negatives"]
        )
        if false_negatives:
            raise ArithmeticError(f"{group.name} has a recorded H_5/2 false negative")
        checked[group.name] = {
            "tested_configurations": (
                row["critical_census"]["configurations"]
                + row["uniform_cell_stress_census"]["configurations"]
            ),
            "false_negatives": false_negatives,
            "strict_overcover_witness": True,
        }
    if set(checked) != {"outer_h25", "old_inner_h25", "new_inner_h25"}:
        raise ArithmeticError("order-5/2 group inventory changed")
    return {
        "schema": "primegaps.physical-factorization-oracle-replay.v1",
        "status": "checked-exact-predicate-record",
        "groups": checked,
        "claim_boundary": "finite census, not physical-law probabilities or a theorem",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("record", type=Path, nargs="?", default=DEFAULT_RECORD)
    args = parser.parse_args()
    print(json.dumps(check(args.record), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
