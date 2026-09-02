#!/usr/bin/env python3
"""Convert one normalized-Jacobi vector into an exact-decimal candidate."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in (HERE, *HERE.parents) if (parent / "pyproject.toml").is_file())
SOURCE = ROOT / "scripts/build_numerical_j_block_operator.py"
spec = importlib.util.spec_from_file_location("candidate_conversion_basis", SOURCE)
builder = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = builder
spec.loader.exec_module(builder)
q = builder.qmc_verifier


def jacobi_q_coefficients(degree, beta):
    answer = [0] * (degree + 1)
    for middle in range(degree + 1):
        base = (
            math.comb(degree, middle)
            * math.comb(degree + beta, degree - middle)
            * (-1) ** (degree - middle)
        )
        for offset in range(middle + 1):
            answer[degree - middle + offset] += (
                base * math.comb(middle, offset) * (-1) ** offset
            )
    return answer


def atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name, delete=False
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vectors", type=Path, required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--support-config", type=Path, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--k", type=int, default=48)
    parser.add_argument("--digits", type=int, default=15)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    support_payload = json.loads(args.support_config.read_text())
    support = builder.load_support(args.support_config)
    with np.load(args.vectors, allow_pickle=False) as payload:
        vector = np.asarray(payload[args.key], dtype=float)
    basis = tuple(q.basis_indices(args.degree))
    if vector.shape != (len(basis),):
        raise ValueError("vector and basis dimensions differ")
    outer = Fraction(support_payload["U"])
    marginal = Fraction(support_payload["R"])
    epsilon = (outer - marginal) / 2
    coordinate_scale = Fraction(args.k) / outer**2
    combined = defaultdict(float)
    for coefficient, (partition, radial_degree) in zip(vector, basis):
        if coefficient == 0:
            continue
        signature = tuple(2 * value for value in partition)
        normalization = math.sqrt((2 * radial_degree + args.k) / args.k)
        angular_scale = float(coordinate_scale ** sum(partition))
        for slack, jacobi_coefficient in enumerate(
            jacobi_q_coefficients(radial_degree, args.k - 1)
        ):
            combined[(signature, slack)] += (
                coefficient * angular_scale * normalization
                * jacobi_coefficient / float(outer) ** slack
            )
    terms = []
    for (signature, slack), coefficient in sorted(combined.items()):
        rational = Fraction(format(coefficient, f".{args.digits}e"))
        if rational:
            terms.append({
                "signature": list(signature),
                "slack_power": slack,
                "numerator": rational.numerator,
                "denominator": rational.denominator,
            })
    bands = support_payload["bands"]
    payload = {
        "schema": "primegaps-stadlmann-rational-candidate-v1",
        "k": args.k,
        "degree": args.degree,
        "basis_dimension": len(terms),
        "radial_coordinate": f"(U-sum(t)) with U={outer}",
        "rationalization_significant_digits": args.digits,
        "source_coordinate_convention": (
            "m_lambda((k/U^2)t_i^2) times normalized "
            "P_b^(0,k-1)(1-2(U-sum(t))/U)"
        ),
        "support": {
            "delta": support_payload["delta"],
            "epsilon": str(epsilon),
            "A": [str(-epsilon)] + [
                str(Fraction(band["upper"]) - epsilon) for band in bands
            ],
            "B": [band["large_caps"] for band in bands],
        },
        "terms": terms,
    }
    atomic_json(args.output, payload)
    print(json.dumps({
        "output": str(args.output), "terms": len(terms),
        "degree": args.degree, "digits": args.digits,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
