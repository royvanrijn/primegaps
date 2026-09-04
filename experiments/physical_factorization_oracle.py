#!/usr/bin/env python3
"""Compare exact order-three source predicates with grouped H_{5/2}.

All predicate decisions use ``fractions.Fraction``.  The finite census is a
structural experiment, not an integration against the physical fragment law.
It assumes the largest-fragment and opposite-root guards that precede the
H_{5/2} reduction in the PrimeGaps186 source proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import random
import tempfile

from primegaps.physical import (
    FiveHalvesGroup,
    exact_group_failure,
    grouped_five_halves_failure,
    primegaps186_source_data,
)


F = Fraction
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = Path(__file__).with_suffix(".json")
PREDICATE_SOURCE = ROOT / "src" / "primegaps" / "physical.py"
UPSTREAM_REVISION = "61340d0b74163003b32756bb16e91d9209a5e330"


OVERCOVER_WITNESSES = {
    "outer_h25": (
        "6745029623/32255864832",
        "101490889/2015991552",
        "737866193/32255864832",
        "19318927871/129023459328",
        "56744378939/129023459328",
        "4375080215/129023459328",
        "13053922723/258046918656",
        "161836823/2687988736",
    ),
    "old_inner_h25": (
        "1829578999/86015639552",
        "1105427791/86015639552",
        "2427552345/43007819776",
        "38871010487/258046918656",
        "29166287101/86015639552",
        "10398701627/258046918656",
        "7584386705/129023459328",
        "770782157/43007819776",
        "3157189547/64511729664",
        "48339836131/258046918656",
    ),
    "new_inner_h25": (
        "1160287731/43007819776",
        "194752787/8063966208",
        "30450009697/86015639552",
        "21705335261/258046918656",
        "2323318459/86015639552",
        "16636276805/86015639552",
        "419678541/2687988736",
        "1763747071/32255864832",
    ),
}


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def _critical_values(group: FiveHalvesGroup, mesh: F) -> tuple[F, ...]:
    values = {
        group.activation + mesh,
        2 * group.activation,
        group.split - mesh,
        group.split,
        group.split + mesh,
        group.threshold / F(7, 2) - mesh,
        group.threshold / F(7, 2) + mesh,
        group.hard_cap / 2,
        3 * group.hard_cap / 4,
        group.hard_cap - mesh,
        group.hard_cap,
    }
    stride = max(1, len(group.rows) // 4)
    values.update(row.activation + mesh for row in group.rows[::stride])
    return tuple(
        sorted(
            value
            for value in values
            if group.activation < value <= group.hard_cap
        )
    )


def _classify(group: FiveHalvesGroup, fragments: tuple[F, ...], mesh: F):
    active_mass = sum(fragments, F(0))
    # The remainder is realizable by fragments at or below the activation and
    # affects the radial condition without entering either obstruction.
    total_mass = max(active_mass, group.radial_lower + mesh)
    if total_mass > group.radial_upper:
        return None
    exact = exact_group_failure(group, fragments, total_mass)
    grouped = grouped_five_halves_failure(group, fragments, total_mass)
    return exact, grouped, total_mass


def _blank_counts() -> dict[str, int]:
    return {
        "configurations": 0,
        "exact_bad": 0,
        "h25_bad": 0,
        "both_bad": 0,
        "false_negatives": 0,
        "strict_overcovers": 0,
    }


def _add(counts: dict[str, int], exact: bool, grouped: bool) -> None:
    counts["configurations"] += 1
    counts["exact_bad"] += int(exact)
    counts["h25_bad"] += int(grouped)
    counts["both_bad"] += int(exact and grouped)
    counts["false_negatives"] += int(exact and not grouped)
    counts["strict_overcovers"] += int(grouped and not exact)


def critical_census(
    group: FiveHalvesGroup, mesh: F, max_fragments: int
) -> dict[str, object]:
    values = _critical_values(group, mesh)
    counts = _blank_counts()
    for size in range(2, max_fragments + 1):
        for fragments in combinations_with_replacement(values, size):
            classified = _classify(group, fragments, mesh)
            if classified is not None:
                _add(counts, classified[0], classified[1])
    return {
        **counts,
        "critical_value_count": len(values),
        "max_fragments": max_fragments,
    }


def random_cell_census(
    group: FiveHalvesGroup,
    mesh: F,
    *,
    samples: int,
    seed: int,
) -> dict[str, object]:
    generator = random.Random(seed)
    first = int(group.activation // mesh) + 1
    last = int(group.hard_cap // mesh)
    counts = _blank_counts()
    for _ in range(samples):
        size = generator.randint(0, 14)
        fragments = tuple(generator.randint(first, last) * mesh for _ in range(size))
        classified = _classify(group, fragments, mesh)
        if classified is not None:
            _add(counts, classified[0], classified[1])
    return {
        **counts,
        "draws": samples,
        "seed": seed,
        "fragment_count_range": [0, 14],
        "fragment_cell_range": [first, last],
    }


def exact_overcover_witness(group: FiveHalvesGroup, mesh: F) -> dict[str, object]:
    fragments = tuple(F(value) for value in OVERCOVER_WITNESSES[group.name])
    classified = _classify(group, fragments, mesh)
    if classified is None:
        raise ArithmeticError("recorded witness is outside its source window")
    exact, grouped, total_mass = classified
    if exact or not grouped:
        raise ArithmeticError("recorded witness is not a strict H_5/2 overcover")
    return {
        "fragments": [str(value) for value in fragments],
        "total_mass": str(total_mass),
        "exact_group_failure": exact,
        "grouped_h25_failure": grouped,
    }


def build_result(*, samples: int, seed: int, max_fragments: int) -> dict[str, object]:
    old, new, groups, mesh = primegaps186_source_data()
    rows = []
    for offset, group in enumerate(groups):
        critical = critical_census(group, mesh, max_fragments)
        random_census = random_cell_census(
            group,
            mesh,
            samples=samples,
            seed=seed + offset,
        )
        if critical["false_negatives"] or random_census["false_negatives"]:
            raise ArithmeticError("H_5/2 missed an exact source failure")
        rows.append(
            {
                "group": group.name,
                "role": group.role,
                "source_row_count": len(group.rows),
                "parameters": {
                    "activation": str(group.activation),
                    "threshold": str(group.threshold),
                    "radial_lower": str(group.radial_lower),
                    "radial_upper": str(group.radial_upper),
                    "hard_cap": str(group.hard_cap),
                    "split": str(group.split),
                },
                "critical_census": critical,
                "uniform_cell_stress_census": random_census,
                "strict_overcover_witness": exact_overcover_witness(group, mesh),
            }
        )
    return {
        "schema": "primegaps.physical-factorization-oracle.v1",
        "status": "checked-finite-exact-predicate-experiment",
        "upstream": {
            "repository": "https://github.com/openai/PrimeGaps186",
            "revision": UPSTREAM_REVISION,
            "source_sha256": "7f71bdefcfe3bb5ca76a143929b3cb3f4156c21dc483253cda3077420f1e5de4",
            "numerical_note": "equations (1.12), (1.21), (1.23)-(1.24), (1.34); Lemma 1.4",
        },
        "exact_inputs": {
            "old_source_rows": len(old),
            "new_source_rows": len(new),
            "mesh": str(mesh),
            "comparison": "union of exact order-three row failures versus grouped nonlargest H_{5/2}",
            "guards_assumed": [
                "group radial window",
                "largest-fragment cap",
                "opposite-root allocation bound from PrimeGaps186 Lemma 1.4",
            ],
            "predicate_module": str(PREDICATE_SOURCE.relative_to(ROOT)),
            "predicate_module_sha256": file_hash(PREDICATE_SOURCE),
        },
        "groups": rows,
        "conclusion": {
            "false_negatives": sum(
                row["critical_census"]["false_negatives"]
                + row["uniform_cell_stress_census"]["false_negatives"]
                for row in rows
            ),
            "strict_overcover_witnesses": len(rows),
            "statement": (
                "The exact rational oracle found no H_{5/2} false negative in either "
                "finite census and has an explicit strict-overcoverage witness in every "
                "group. This is consistent with the source majorization, while the exact "
                "witnesses prove that it is not an equivalence. The general implication "
                "comes from PrimeGaps186 Lemma 1.4, not this finite census, and census "
                "frequencies are not probabilities under the physical fragment law."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=186_250)
    parser.add_argument("--max-fragments", type=int, default=7)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build_result(
        samples=args.samples,
        seed=args.seed,
        max_fragments=args.max_fragments,
    )
    result["implementation_sha256"] = file_hash(Path(__file__))
    atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
