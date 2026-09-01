#!/usr/bin/env python3
"""Replay a PrimeGapsLib-style sparse simplex certificate in Python."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from primegaps.symmetric import (
    evaluate_sparse_symmetric_certificate,
    load_sparse_symmetric_terms,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--dimension", type=int, default=50)
    parser.add_argument("--epsilon-denominator", type=int, default=25)
    parser.add_argument("--degree-bound", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    source = args.certificate.read_bytes()
    started = perf_counter()
    terms = load_sparse_symmetric_terms(args.certificate)
    result = evaluate_sparse_symmetric_certificate(
        terms,
        args.dimension,
        args.epsilon_denominator,
        degree_bound=args.degree_bound,
    )
    payload = {
        "format": "primegaps-sparse-symmetric-replay-v1",
        "input_sha256": sha256(source).hexdigest(),
        "dimension": result.dimension,
        "epsilon_denominator": result.epsilon_denominator,
        "degree_bound": result.degree_bound,
        "term_count": result.term_count,
        "signature_count": result.signature_count,
        "mass_group_count": result.mass_group_count,
        "mass_transform_count": result.mass_transform_count,
        "marginal_feature_count": result.marginal_feature_count,
        "marginal_group_count": result.marginal_group_count,
        "marginal_transform_count": result.marginal_transform_count,
        "mass": result.mass,
        "marginal": result.marginal,
        "difference": result.difference,
        "quotient": {
            "numerator": result.quotient.numerator,
            "denominator": result.quotient.denominator,
        },
        "certified": result.certified,
        "elapsed_seconds": perf_counter() - started,
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output is not None:
        args.output.write_text(serialized)
    else:
        print(serialized, end="")
    return 0 if result.certified else 2


if __name__ == "__main__":
    raise SystemExit(main())
