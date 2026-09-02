#!/usr/bin/env python3
"""Measure repeated J/Jc cost for the candidate-independent block layout."""

from __future__ import annotations

import argparse
import json
import time

import numpy as np

from primegaps.basis import symmetric_basis
from primegaps.fast_exact.j_block import JBlockOperator, MarginalMap


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--k", type=int, default=49)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    basis = [
        (tuple(2 * value for value in signature), slack)
        for signature, slack in symmetric_basis(args.degree, args.k)
    ]
    marginal_map = MarginalMap.from_basis(basis)
    signatures = tuple(marginal_map.feature_keys)
    rng = np.random.default_rng(20260902)
    blocks = {}
    stored_entries = 0
    for left_index, left in enumerate(signatures):
        left_size = len(marginal_map.feature_keys[left])
        for right in signatures[left_index:]:
            right_size = len(marginal_map.feature_keys[right])
            block = rng.standard_normal((left_size, right_size))
            if left == right:
                block = (block + block.T) / 2
            blocks[(left, right)] = block
            stored_entries += block.size
    operator = JBlockOperator(marginal_map, blocks)
    vector = rng.standard_normal(len(basis))
    operator.matvec(vector)
    started = time.perf_counter()
    for _ in range(args.repeats):
        result = operator.matvec(vector)
    elapsed = time.perf_counter() - started
    quadratic = float(vector @ result)
    print(json.dumps({
        "degree": args.degree,
        "k": args.k,
        "basis_dimension": len(basis),
        "marginal_signature_blocks": len(signatures),
        "marginal_feature_dimension": sum(
            len(keys) for keys in marginal_map.feature_keys.values()
        ),
        "stored_blocks": len(blocks),
        "stored_float64_entries": stored_entries,
        "stored_mebibytes": stored_entries * 8 / 2**20,
        "repeats": args.repeats,
        "seconds_per_matvec": elapsed / args.repeats,
        "finite_quadratic_smoke": bool(np.isfinite(quadratic)),
        "random_seed": 20260902,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
