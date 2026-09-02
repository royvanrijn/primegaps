#!/usr/bin/env python3
"""Cheap exact finalization of complete accelerated I/J checkpoints."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


def file_hash(path):
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    rows = {}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        row = json.loads(line)
        signature = tuple(row["signature"])
        if signature in rows:
            raise ValueError(f"duplicate signature in {path}:{line_number}")
        rows[signature] = Fraction(
            int(row["numerator"]), int(row["denominator"])
        )
    return rows


def rational_payload(value):
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def decimal_string(value, digits):
    with localcontext() as context:
        context.prec = digits
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--i-groups", type=Path, required=True)
    parser.add_argument("--j-groups", type=Path, required=True)
    parser.add_argument("--expected-groups", type=int, default=2714)
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    i_rows, j_rows = load(args.i_groups), load(args.j_groups)
    if len(i_rows) != args.expected_groups or len(j_rows) != args.expected_groups:
        raise ValueError(
            f"expected {args.expected_groups} groups; got I={len(i_rows)}, "
            f"J={len(j_rows)}"
        )
    if set(i_rows) != set(j_rows):
        raise ValueError("I and J signature sets differ")
    i_value, j_value = sum(i_rows.values()), sum(j_rows.values())
    difference = args.k * j_value - i_value
    quotient = args.k * j_value / i_value
    deficit = 1 - quotient
    payload = {
        "schema": "primegaps-fast-exact-final-v1",
        "k": args.k,
        "groups": args.expected_groups,
        "i_groups_sha256": file_hash(args.i_groups),
        "j_groups_sha256": file_hash(args.j_groups),
        "I": rational_payload(i_value),
        "J": rational_payload(j_value),
        "kJ_minus_I": rational_payload(difference),
        "kJ_over_I": rational_payload(quotient),
        "kJ_over_I_decimal": decimal_string(quotient, args.digits),
        "one_minus_kJ_over_I": rational_payload(deficit),
        "one_minus_kJ_over_I_decimal": decimal_string(deficit, args.digits),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(encoded)
        temporary.replace(args.output)


if __name__ == "__main__":
    main()
