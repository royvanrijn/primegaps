#!/usr/bin/env python3
"""Compare accelerated exact rows with an independent oracle checkpoint."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


def load(path):
    rows = {}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        row = json.loads(line)
        signature = tuple(row["signature"])
        if signature in rows:
            raise ValueError(f"duplicate signature in {path}:{line_number}")
        rows[signature] = (str(row["numerator"]), str(row["denominator"]))
    if not rows:
        raise ValueError(f"empty checkpoint: {path}")
    return rows


def mapping_hash(rows):
    payload = [
        [list(signature), numerator, denominator]
        for signature, (numerator, denominator) in sorted(rows.items())
    ]
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()
    return sha256(encoded).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("actual", type=Path)
    parser.add_argument("oracle", type=Path)
    parser.add_argument("--allow-prefix", action="store_true")
    args = parser.parse_args()

    actual = load(args.actual)
    oracle = load(args.oracle)
    foreign = set(actual) - set(oracle)
    missing = set(oracle) - set(actual)
    mismatches = {
        signature: (actual[signature], oracle[signature])
        for signature in set(actual) & set(oracle)
        if actual[signature] != oracle[signature]
    }
    if foreign or mismatches or (missing and not args.allow_prefix):
        raise SystemExit(json.dumps({
            "status": "mismatch",
            "actual_rows": len(actual),
            "oracle_rows": len(oracle),
            "foreign": len(foreign),
            "missing": len(missing),
            "value_mismatches": len(mismatches),
        }, sort_keys=True))
    print(json.dumps({
        "status": "exact-prefix-match" if missing else "exact-complete-match",
        "actual_rows": len(actual),
        "oracle_rows": len(oracle),
        "mapping_sha256": mapping_hash(actual),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
