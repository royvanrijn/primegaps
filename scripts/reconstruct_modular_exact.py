#!/usr/bin/env python3
"""Reconstruct exact rational group rows from modular J checkpoints.

The bounds file is JSONL with ``signature``, ``numerator_bound`` and
``denominator_bound``.  Bounds are mandatory: successful reconstruction is a
proof only when ``2*N*D < product(primes)``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from primegaps.fast_exact.modular_exact import crt, rational_reconstruction


def load_unique(path, required):
    rows = {}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        row = json.loads(line)
        signature = tuple(row["signature"])
        if signature in rows:
            raise ValueError(f"duplicate signature in {path}:{line_number}")
        missing = required - set(row)
        if missing:
            raise ValueError(f"missing {sorted(missing)} in {path}:{line_number}")
        rows[signature] = row
    if not rows:
        raise ValueError(f"empty input: {path}")
    return rows


def load_modular(path):
    """Return one or more per-prime row mappings from a checkpoint."""
    raw = load_unique(path, set())
    first = next(iter(raw.values()))
    if "prime" in first and "residue" in first:
        primes = (int(first["prime"]),)
        mappings = ({},)
        scalar = True
    elif "primes" in first and "residues" in first:
        primes = tuple(int(prime) for prime in first["primes"])
        mappings = tuple({} for _ in primes)
        scalar = False
    else:
        raise ValueError(f"missing modular values in {path}")
    for signature, row in raw.items():
        if scalar:
            if int(row.get("prime", -1)) != primes[0]:
                raise ValueError(f"inconsistent prime in {path}")
            mappings[0][signature] = int(row["residue"])
        else:
            row_primes = tuple(int(prime) for prime in row.get("primes", ()))
            residues = tuple(int(value) for value in row.get("residues", ()))
            if row_primes != primes or len(residues) != len(primes):
                raise ValueError(f"inconsistent batched primes in {path}")
            for mapping, residue in zip(mappings, residues):
                mapping[signature] = residue
    return tuple(zip(primes, mappings))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--bounds", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    images = [image for path in args.inputs for image in load_modular(path)]
    expected = set(images[0][1])
    if any(set(rows) != expected for _prime, rows in images[1:]):
        raise ValueError("modular checkpoints have different signature sets")
    bounds = load_unique(
        args.bounds, {"numerator_bound", "denominator_bound"}
    )
    if set(bounds) != expected:
        raise ValueError("bounds and checkpoints have different signature sets")

    primes = [prime for prime, _rows in images]
    if len(set(primes)) != len(primes):
        raise ValueError("duplicate CRT prime")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    with temporary.open("w") as stream:
        for signature in sorted(expected):
            residue, modulus = crt(
                [rows[signature] for _prime, rows in images], primes
            )
            bound = bounds[signature]
            numerator, denominator = rational_reconstruction(
                residue,
                modulus,
                numerator_bound=int(bound["numerator_bound"]),
                denominator_bound=int(bound["denominator_bound"]),
            )
            stream.write(json.dumps({
                "signature": list(signature),
                "numerator": str(numerator),
                "denominator": str(denominator),
                "crt_modulus": str(modulus),
            }, sort_keys=True) + "\n")
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
